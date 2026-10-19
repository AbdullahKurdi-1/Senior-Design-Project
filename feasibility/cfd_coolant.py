#!/usr/bin/env python3
"""
Item 1 test: boiling Dowtherm A coolant (about 3.4 bar at 319 C) instead of boiling water (113 bar).
Organic fluids boil with a lower heat-transfer coefficient than water, so h = 1,500 and 2,500 W/(m2 K)
are tested. The plate is thickened to 1.5 mm for the 3.3 MPa worst-case pressure difference.
"""
import copy, csv, json, os, sys
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cfd_stress as cs

base = json.load(open(os.path.join(HERE, "results", "settings_kinetic.json")))
base["geometry"]["wall_thickness_m"] = 0.0015
base["mesh"]["cells_wall_y"] = 6
ol = base["control"]["setpoint_mean_bed_T_K"] - 1.0
cases = []
def add(name, mode, t_end, h, keff=None, step=None):
    s = copy.deepcopy(base)
    s["operating"]["coolant_h_W_m2K"] = h
    if keff: s["gas"]["k_eff_W_mK"] = keff
    if step: s["control"]["feed_flow_step"] = step
    if mode == "open": s["operating"]["coolant_T_initial_K"] = ol
    cases.append((name, mode, t_end, s))
add("dow_h1500_pid_plus10", "pid", 160.0, 1500.0)
add("dow_h1500_pid_minus10", "pid", 160.0, 1500.0, step=0.90)
add("dow_h2500_pid_plus10", "pid", 160.0, 2500.0)
add("dow_h1500_worst_keff03", "open", 40.0, 1500.0, keff=0.3)
with ThreadPoolExecutor(max_workers=2) as ex:
    res = list(ex.map(cs.run, cases))
with open(os.path.join(HERE, "results", "cfd_coolant.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(res[0].keys())); w.writeheader(); w.writerows(res)
for r in res:
    print("%-24s max dT %.2f K  Tmax %.2f C  Tmean %.2f C  coolant %.2f C  %s  %s" % (r["case"], r["max_dT_K"], r["final_Tmax_C"],
          r["final_Tmean_C"], r["final_coolant_C"], r["energy_imbalance"].split(" of")[0], "PASS" if r["max_dT_K"] <= 15 else "FAIL"))
