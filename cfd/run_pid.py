#!/usr/bin/env python3
"""
Python co-simulation: chtMultiRegionFoam + PI(D) controller on the coolant temperature.

Every `chunk_s` seconds of simulated time the script
  1. runs chtMultiRegionFoam for one chunk,
  2. reads the bed (fluid-region) temperature field and computes mean / max / min,
  3. updates the coolant temperature Ta with a PI(D) law that holds the mean bed
     temperature at the setpoint,
  4. recomputes the reaction heat release of every bed zone with the exact Arrhenius
     law at that zone's current temperature (explicit coupling, one update per chunk),
  5. at `disturbance_time_s` applies a feed-flow step (+10 % by default): the inlet
     velocity and the reaction heat release are both scaled by `feed_flow_step`
     (heat scaled 1:1 with flow = conservative upper bound).

Outputs (in ./results):
  pid_log.csv          time history
  pid_response.png     bed temperatures and coolant temperature vs time
  summary.txt          spec check: max (T_max - T_min) in the bed <= 15 K

Usage:
  python3 setup_case.py     # once
  python3 run_pid.py        # add --open-loop to keep the coolant fixed (no control)
"""
import argparse
import csv
import gzip
import json
import os
import re
import subprocess
import sys

import setup_case as sc

HERE = os.path.dirname(os.path.abspath(__file__))


# ------------------------------------------------------------------ file helpers
def _open(path):
    if os.path.exists(path + ".gz"):
        return path + ".gz", True
    if os.path.exists(path):
        return path, False
    raise FileNotFoundError(path)


def read_text(path):
    p, gz = _open(path)
    if gz:
        with gzip.open(p, "rt") as f:
            return f.read(), p, gz
    with open(p) as f:
        return f.read(), p, gz


def write_text(p, gz, txt):
    if gz:
        with gzip.open(p, "wt") as f:
            f.write(txt)
    else:
        with open(p, "w") as f:
            f.write(txt)


def internal_field(path):
    txt, _, _ = read_text(path)
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(([^)]*)\)", txt, re.S)
    if m:
        return [float(v) for v in m.group(2).split()]
    m = re.search(r"internalField\s+uniform\s+([0-9.eE+-]+)", txt)
    return [float(m.group(1))]


def time_dirs(case):
    out = []
    for d in os.listdir(case):
        try:
            out.append((float(d), d))
        except ValueError:
            pass
    return sorted(out)


def set_coolant(case, tdir, Ta):
    txt, p, gz = read_text(os.path.join(case, tdir, "solid", "T"))
    new, n = re.subn(r"(coolant\s*\{[^}]*?\bTa\s+(?:uniform|constant)\s+)[0-9.eE+-]+", r"\g<1>%.4f" % Ta, txt, flags=re.S)
    if n != 1:
        sys.exit("could not find coolant Ta entry in %s" % p)
    write_text(p, gz, new)


def set_inlet_velocity(case, tdir, u):
    txt, p, gz = read_text(os.path.join(case, tdir, "fluid", "U"))
    new, n = re.subn(r"(inlet\s*\{[^}]*?\bvalue\s+uniform\s*)\([^)]*\)", r"\g<1>(%.6g 0 0)" % u, txt, flags=re.S)
    if n != 1:
        sys.exit("could not find inlet velocity in %s" % p)
    write_text(p, gz, new)


def read_cell_set(path):
    txt, _, _ = read_text(path)
    body = txt[txt.index("FoamFile"):]
    body = body[body.index("}") + 1:]
    m = re.search(r"(\d+)\s*\(([^)]*)\)", body, re.S)
    return [int(v) for v in m.group(2).split()]


def patch_values(path, patch):
    txt, _, _ = read_text(path)
    i = txt.index("boundaryField")
    m = re.search(r"\n\s*%s\s*\{(.*?)\n\s*\}" % re.escape(patch), txt[i:], re.S)
    blk = m.group(1)
    mm = re.search(r"\bvalue\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(([^)]*)\)", blk, re.S)
    if mm:
        return [float(v) for v in mm.group(2).split()]
    return [float(re.search(r"\bvalue\s+uniform\s+([0-9.eE+-]+)", blk).group(1))]


def energy_balance(s, case, tdir, q_zone, zone_cells, Ta, feed):
    """Steady-state check: heat generated = heat to coolant + gas sensible heat."""
    g, m, op = s["geometry"], s["mesh"], s["operating"]
    L, hb, H = g["length_m"], g["bed_half_thickness_m"], g["slot_height_m"]
    v_cell = L * hb * H / (m["cells_x"] * m["cells_bed_y"] * m["cells_z"])
    Q_gen = sum(q_zone[n] * len(c) * v_cell for n, c in zone_cells.items())
    Tw = patch_values(os.path.join(case, tdir, "solid", "T"), "coolant")
    A_f = L * H / len(Tw)
    Q_cool = sum(op["coolant_h_W_m2K"] * (T - Ta) * A_f for T in Tw)
    Tout = patch_values(os.path.join(case, tdir, "fluid", "T"), "outlet")
    Tin = op["T_inlet_K"]
    rho_in = op["P_Pa"] * op["molWeight_kg_kmol"] / 1000.0 / (8.314 * Tin)
    mdot = rho_in * op["superficial_velocity_m_s"] * feed * hb * H
    Q_gas = mdot * s["gas"]["Cp_J_kgK"] * (sum(Tout) / len(Tout) - Tin)
    # axial conduction entering through the fixed-temperature inlet face
    # (cells are ordered x fastest, then y, then z, as written by blockMesh)
    T = internal_field(os.path.join(case, tdir, "fluid", "T"))
    nx, ny, nz = m["cells_x"], m["cells_bed_y"], m["cells_z"]
    Q_in = 0.0
    if len(T) == nx * ny * nz:
        A_in, dx = (hb / ny) * (H / nz), L / nx
        Q_in = sum(s["gas"]["k_eff_W_mK"] * A_in * (Tin - T[i]) / (0.5 * dx) for i in range(0, nx * ny * nz, nx))
    return Q_gen, Q_in, Q_cool, Q_gas


def set_end_time(case, t_end):
    p = os.path.join(case, "system", "controlDict")
    txt = open(p).read()
    txt = re.sub(r"\nendTime\s+[^;]+;", "\nendTime         %g;" % t_end, txt)
    open(p, "w").write(txt)


# ------------------------------------------------------------------ main loop
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default=os.path.join(HERE, "case"))
    ap.add_argument("--settings", default=os.path.join(HERE, "case_settings.json"))
    ap.add_argument("--open-loop", action="store_true", help="no control: coolant stays at its initial value")
    ap.add_argument("--no-disturbance", action="store_true")
    ap.add_argument("--t-end", type=float, default=None)
    ap.add_argument("--results", default=os.path.join(HERE, "results"))
    a = ap.parse_args()

    case = os.path.abspath(a.case)
    sc.SETTINGS = os.path.abspath(a.settings)
    sc.CASE = case
    s, hp = sc.load()
    c, op = s["control"], s["operating"]
    if not os.path.isdir(os.path.join(case, "constant", "fluid", "polyMesh")):
        sys.exit("mesh not found - run setup_case.py first")
    os.makedirs(a.results, exist_ok=True)

    dt = c["chunk_s"]
    t_end = a.t_end if a.t_end is not None else c["t_end_s"]
    Tsp = c["setpoint_mean_bed_T_K"]
    Ta = op["coolant_T_initial_K"]
    Ta0 = Ta
    integ, e_prev = 0.0, None
    feed = 1.0
    rows = []

    sets_dir = os.path.join(case, "constant", "fluid", "polyMesh", "sets")
    zone_cells = {n: read_cell_set(os.path.join(sets_dir, n)) for n, _, _ in sc.zone_names(s, hp)}
    n_cells = sum(len(v) for v in zone_cells.values())
    t = time_dirs(case)[-1][0]
    print("%d heat-release zones covering %d bed cells" % (len(zone_cells), n_cells))
    print("co-simulation from t = %g s to %g s, chunk %g s, %s"
          % (t, t_end, dt, "OPEN LOOP" if a.open_loop else "PID on coolant T"))
    with open(os.path.join(case, "log.chtMultiRegionFoam"), "a") as log:
        while t < t_end - 1e-9:
            tdir = time_dirs(case)[-1][1]
            # feed disturbance
            if (not a.no_disturbance) and feed == 1.0 and t >= c["disturbance_time_s"] - 1e-9:
                feed = c["feed_flow_step"]
                set_inlet_velocity(case, tdir, op["superficial_velocity_m_s"] * feed)
                print("  t = %6.1f s  feed-flow step x%.2f applied" % (t, feed))
            # reaction heat for this chunk: exact Arrhenius at the current zone temperatures
            T_now = internal_field(os.path.join(case, tdir, "fluid", "T"))
            zone_T = {n: (sum(T_now[k] for k in cells) / len(cells)) if len(T_now) > 1 else T_now[0]
                      for n, cells in zone_cells.items()}
            q_zone = sc.zone_heat(s, hp, feed, zone_T)
            Ta_used = Ta
            with open(os.path.join(case, "constant", "fluid", "fvOptions"), "w") as f:
                f.write(sc.fv_options(s, hp, feed_scale=feed, zone_T=zone_T))
            set_coolant(case, tdir, Ta)
            set_end_time(case, t + dt)
            r = subprocess.run(["chtMultiRegionFoam"], cwd=case, stdout=log, stderr=subprocess.STDOUT)
            if r.returncode != 0:
                sys.exit("chtMultiRegionFoam failed - see %s/log.chtMultiRegionFoam" % case)
            t, tdir = time_dirs(case)[-1]
            T = internal_field(os.path.join(case, tdir, "fluid", "T"))   # uniform cells -> plain mean = volume mean
            Tavg, Tmax, Tmin = sum(T) / len(T), max(T), min(T)
            # PI(D) with clamping anti-windup
            if not a.open_loop:
                e = Tsp - Tavg
                d = 0.0 if e_prev is None else (e - e_prev) / dt
                trial_int = integ + e * dt
                u = Ta0 + c["Kp"] * e + c["Ki"] * trial_int + c["Kd"] * d
                if c["coolant_T_min_K"] < u < c["coolant_T_max_K"]:
                    integ = trial_int
                Ta = min(max(Ta0 + c["Kp"] * e + c["Ki"] * integ + c["Kd"] * d, c["coolant_T_min_K"]), c["coolant_T_max_K"])
                e_prev = e
            rows.append({"time_s": t, "feed_scale": feed, "coolant_T_C": Ta - 273.15,
                         "bed_T_mean_C": Tavg - 273.15, "bed_T_max_C": Tmax - 273.15,
                         "bed_T_min_C": Tmin - 273.15, "bed_dT_K": Tmax - Tmin})
            print("  t = %6.1f s  bed mean %.2f C  max %.2f C  min %.2f C  dT %.2f K  coolant %.2f C"
                  % (t, Tavg - 273.15, Tmax - 273.15, Tmin - 273.15, Tmax - Tmin, Ta - 273.15))

    with open(os.path.join(a.results, "pid_log.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    dT_max = max(r["bed_dT_K"] for r in rows)
    after = [r for r in rows if r["time_s"] >= c["disturbance_time_s"]]
    lines = [
        "Reactor-slot CFD summary",
        "  bed half-thickness          : %.2f mm (bed thickness %.2f mm)" % (s["geometry"]["bed_half_thickness_m"] * 1e3, 2e3 * s["geometry"]["bed_half_thickness_m"]),
        "  mode                        : %s" % ("open loop" if a.open_loop else "PID on coolant temperature"),
        "  max bed dT over the run     : %.2f K   (spec: <= 15 K)  -> %s" % (dT_max, "PASS" if dT_max <= 15 else "FAIL"),
        "  final mean bed temperature  : %.2f C  (setpoint %.2f C)" % (rows[-1]["bed_T_mean_C"], Tsp - 273.15),
        "  final coolant temperature   : %.2f C" % rows[-1]["coolant_T_C"],
    ]
    if after:
        lines.append("  peak bed T after disturbance: %.2f C" % max(r["bed_T_max_C"] for r in after))
    Qg, Qi, Qc, Qs = energy_balance(s, case, tdir, q_zone, zone_cells, Ta_used, feed)
    lines += [
        "  energy balance, last chunk (half slot):  generated + inlet conduction = coolant + gas",
        "     heat generated by reaction       : %.4f W" % Qg,
        "     conduction in through inlet face : %.4f W" % Qi,
        "     heat removed by coolant          : %.4f W" % Qc,
        "     gas enthalpy change (out - in)   : %.4f W" % Qs,
        "     imbalance                        : %.2f %% of generation (should be small at steady state)" % (100 * (Qg + Qi - Qc - Qs) / Qg),
    ]
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(a.results, "summary.txt"), "w").write(txt + "\n")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        tt = [r["time_s"] for r in rows]
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
        ax1.plot(tt, [r["bed_T_max_C"] for r in rows], label="bed max")
        ax1.plot(tt, [r["bed_T_mean_C"] for r in rows], label="bed mean")
        ax1.plot(tt, [r["bed_T_min_C"] for r in rows], label="bed min")
        ax1.axhline(Tsp - 273.15, ls="--", lw=0.8, color="k", label="setpoint")
        ax1.set_ylabel("Temperature [C]")
        ax1.legend(fontsize=8)
        ax2.plot(tt, [r["coolant_T_C"] for r in rows], color="tab:purple")
        ax2.set_ylabel("Coolant T [C]")
        ax2.set_xlabel("Time [s]")
        for ax in (ax1, ax2):
            if not a.no_disturbance:
                ax.axvline(c["disturbance_time_s"], color="gray", lw=0.8, ls=":")
            ax.grid(alpha=0.3)
        ax1.set_title("Reactor slot: bed temperature control (+%d%% feed step at %g s)"
                      % (round((c["feed_flow_step"] - 1) * 100), c["disturbance_time_s"]))
        fig.tight_layout()
        fig.savefig(os.path.join(a.results, "pid_response.png"), dpi=150)
        print("plot written to", os.path.join(a.results, "pid_response.png"))
    except ImportError:
        print("matplotlib not installed - skipped plot")


if __name__ == "__main__":
    main()
