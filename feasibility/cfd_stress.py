#!/usr/bin/env python3
"""
Spec 3 stress tests (bed temperature difference <= 15 K) with the KINETIC heat profile.

Cases
  pid_plus10   PID on, +10 % feed step at 60 s          (spec case)
  pid_minus10  PID on, -10 % feed step at 60 s
  mesh_*       open loop, coarse / base / fine mesh     (mesh independence)
  keff_*       open loop, effective bed conductivity 0.3 and 1.5 W/(m K)
  h_*          open loop, coolant heat-transfer coefficient 1000 and 10000 W/(m2 K)
  worst        open loop, k_eff 0.3 and h 1000 together
  bed10mm      open loop, 10 mm bed (how much margin the 4 mm choice has)

Usage:  python3 cfd_stress.py [--jobs 2]
Writes results/cfd_stress.csv and a summary to stdout.
"""
import argparse
import copy
import csv
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
CFD = os.path.join(REPO, "cfd")
RUNS = os.path.join(HERE, "cfd_runs")


def cases(base):
    out = []

    def add(name, mode, t_end, **mods):
        s = copy.deepcopy(base)
        for path, v in mods.items():
            sec, key = path.split("__")
            s[sec][key] = v
        out.append((name, mode, t_end, s))

    open_loop_T = base["control"]["setpoint_mean_bed_T_K"] - 1.0
    ol = {"operating__coolant_T_initial_K": open_loop_T}
    add("pid_plus10", "pid", 160.0)
    add("pid_minus10", "pid", 160.0, control__feed_flow_step=0.90)
    add("mesh_coarse", "open", 40.0, mesh__cells_x=45, mesh__cells_bed_y=8, **ol)
    add("mesh_base", "open", 40.0, **ol)
    add("mesh_fine", "open", 40.0, mesh__cells_x=150, mesh__cells_bed_y=32, **ol)
    add("keff_0.3", "open", 40.0, gas__k_eff_W_mK=0.3, **ol)
    add("keff_1.5", "open", 40.0, gas__k_eff_W_mK=1.5, **ol)
    add("h_1000", "open", 40.0, operating__coolant_h_W_m2K=1000.0, **ol)
    add("h_10000", "open", 40.0, operating__coolant_h_W_m2K=10000.0, **ol)
    add("worst", "open", 40.0, gas__k_eff_W_mK=0.3, operating__coolant_h_W_m2K=1000.0, **ol)
    add("bed10mm", "open", 40.0, geometry__bed_half_thickness_m=0.005, mesh__cells_bed_y=20, **ol)
    return out


def run(case):
    name, mode, t_end, s = case
    d = os.path.join(RUNS, name)
    os.makedirs(d, exist_ok=True)
    sp = os.path.join(d, "settings.json")
    json.dump(s, open(sp, "w"), indent=2)
    log = open(os.path.join(d, "driver.log"), "w")
    subprocess.run([sys.executable, os.path.join(CFD, "setup_case.py"), "--settings", sp, "--case", os.path.join(d, "case")],
                   check=True, stdout=log, stderr=subprocess.STDOUT)
    args = [sys.executable, os.path.join(CFD, "run_pid.py"), "--settings", sp, "--case", os.path.join(d, "case"),
            "--results", os.path.join(d, "results"), "--t-end", str(t_end)]
    if mode == "open":
        args += ["--open-loop", "--no-disturbance"]
    subprocess.run(args, check=True, stdout=log, stderr=subprocess.STDOUT)
    rows = list(csv.DictReader(open(os.path.join(d, "results", "pid_log.csv"))))
    summ = open(os.path.join(d, "results", "summary.txt")).read()
    imb = [l for l in summ.splitlines() if "imbalance" in l]
    return {"case": name, "mode": mode,
            "max_dT_K": max(float(r["bed_dT_K"]) for r in rows),
            "final_dT_K": float(rows[-1]["bed_dT_K"]),
            "final_Tmax_C": float(rows[-1]["bed_T_max_C"]),
            "final_Tmean_C": float(rows[-1]["bed_T_mean_C"]),
            "final_coolant_C": float(rows[-1]["coolant_T_C"]),
            "energy_imbalance": imb[0].split(":")[1].strip() if imb else ""}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=2)
    a = ap.parse_args()
    base = json.load(open(os.path.join(HERE, "results", "settings_kinetic.json")))
    cs = cases(base)
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        res = list(ex.map(run, cs))
    with open(os.path.join(HERE, "results", "cfd_stress.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(res[0].keys()))
        w.writeheader()
        w.writerows(res)
    print("%-12s %-5s %9s %9s %9s %9s  %s" % ("case", "mode", "max dT K", "Tmax C", "Tmean C", "cool C", "energy imbalance"))
    for r in res:
        print("%-12s %-5s %9.2f %9.2f %9.2f %9.2f  %s %s" % (r["case"], r["mode"], r["max_dT_K"], r["final_Tmax_C"],
              r["final_Tmean_C"], r["final_coolant_C"], r["energy_imbalance"], "PASS" if r["max_dT_K"] <= 15 else "FAIL"))


if __name__ == "__main__":
    main()
