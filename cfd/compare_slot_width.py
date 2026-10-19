#!/usr/bin/env python3
"""
Design study: why the bed must be micro-structured.

Runs the same reactor slot (same catalyst, same GHSV, same coolant temperature,
no controller) for several bed thicknesses and reports the bed temperature
difference. A 25 mm bed corresponds to a conventional packed tube.

Usage:  python3 compare_slot_width.py            (about 5-15 min on a laptop)
Output: results/slot_width_study.csv and results/slot_width_study.png
"""
import copy
import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.join(HERE, "study")
RES = os.path.join(HERE, "results")

# bed thickness [mm] -> cells across the half bed (multiple of n_radial_bands = 4)
CASES = [(4.0, 16), (10.0, 20), (25.0, 48)]
T_END = 30.0     # s, enough for these cases to reach steady state


def main():
    base = json.load(open(os.path.join(HERE, "case_settings.json")))
    os.makedirs(STUDY, exist_ok=True)
    os.makedirs(RES, exist_ok=True)
    rows = []
    for thick_mm, ny in CASES:
        s = copy.deepcopy(base)
        s["geometry"]["bed_half_thickness_m"] = thick_mm / 2000.0
        s["mesh"]["cells_bed_y"] = ny
        s["control"]["chunk_s"] = 1.0
        tag = "bed_%gmm" % thick_mm
        sp = os.path.join(STUDY, tag + ".json")
        json.dump(s, open(sp, "w"), indent=2)
        case = os.path.join(STUDY, tag)
        print("=== %s" % tag)
        subprocess.run([sys.executable, os.path.join(HERE, "setup_case.py"), "--settings", sp, "--case", case], check=True)
        out = os.path.join(STUDY, tag + "_results")
        subprocess.run([sys.executable, os.path.join(HERE, "run_pid.py"), "--settings", sp, "--case", case,
                        "--open-loop", "--no-disturbance", "--t-end", str(T_END), "--results", out], check=True)
        last = list(csv.DictReader(open(os.path.join(out, "pid_log.csv"))))[-1]
        rows.append({"bed_thickness_mm": thick_mm, "bed_T_max_C": float(last["bed_T_max_C"]),
                     "bed_T_mean_C": float(last["bed_T_mean_C"]), "bed_T_min_C": float(last["bed_T_min_C"]),
                     "bed_dT_K": float(last["bed_dT_K"]), "coolant_T_C": float(last["coolant_T_C"])})

    with open(os.path.join(RES, "slot_width_study.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("\nbed thickness [mm]   max bed T [C]   bed dT [K]   spec dT <= 15 K")
    for r in rows:
        print("   %6.1f            %7.2f        %7.2f       %s" % (r["bed_thickness_mm"], r["bed_T_max_C"], r["bed_dT_K"],
                                                                  "PASS" if r["bed_dT_K"] <= 15 else "FAIL"))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 4))
        x = [r["bed_thickness_mm"] for r in rows]
        ax.plot(x, [r["bed_dT_K"] for r in rows], "o-")
        ax.axhline(15, ls="--", color="k", lw=0.8, label="spec limit 15 K")
        ax.set_xlabel("Bed thickness [mm]")
        ax.set_ylabel("Bed temperature difference [K]")
        ax.set_title("Hot-spot size vs. bed thickness (same coolant, no control)")
        ax.grid(alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(RES, "slot_width_study.png"), dpi=150)
    except ImportError:
        pass


if __name__ == "__main__":
    main()
