#!/usr/bin/env python3
"""
Heat-release profile along one 150 mm slot from the calibrated kinetics (instead of the
first-order placeholder), at the GHSV that gives 35 % single-pass conversion.
Writes results/heat_profile_kinetic.json (same format as reactor_model/heat_profile.json)
and results/settings_kinetic.json (CFD settings with the matching gas velocity).
"""
import copy
import json
import math
import os
import sys

import numpy as np

trapz = getattr(np, "trapezoid", None) or np.trapz   # numpy >= 2.0 / older

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import kinetics as kin

A_CASE, SCO_CASE = 1.0, 0.05      # base kinetic case (see kinetics.py)


def main():
    s = json.load(open(os.path.join(HERE, "..", "cfd", "case_settings.json")))
    L, rho_b = s["geometry"]["length_m"], s["bed"]["bulk_density_kg_m3"]
    T, P = s["reaction"]["T_ref_K"], s["operating"]["P_Pa"] / 1e5
    dist = kin.distribution()
    d = dist["d_extra_H2"]
    k1, k2 = kin.calibrate(S_co_cal=SCO_CASE, a=A_CASE, d=d)
    g = kin.ghsv_for_conversion(0.35, k1, k2, A_CASE, d)
    F0 = kin.ghsv_feed(g)                       # mol/s per kg catalyst
    sol = kin.pfr(F0, 1.0, P, T, k1, k2, A_CASE, d, n_out=1501)
    W = sol.t
    q = []
    for i in range(len(W)):
        r1, r2 = kin.rates(sol.y[:, i], P, T, k1, k2, A_CASE, d)
        q.append(-(r1 * 41.2e3 + r2 * (-(165e3 + 41e3 * d))) * rho_b)   # W/m3 bed
    q = np.array(q)
    z = W / 1.0 * L                              # 1 kg spans the full slot length
    V = 1.0 / rho_b
    A = V / L
    Fg = F0[0] + F0[1]
    u = Fg * 22.414e-3 * (T / 273.15) * (1.01325 / P) / A
    n = s["reaction"]["n_zones"]
    edges = np.linspace(0, L, n + 1)
    zones = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (z >= a) & (z <= b)
        zones.append({"x0": float(a), "x1": float(b), "q_W_m3": float(trapz(q[m], z[m]) / (z[m][-1] - z[m][0]))})
    X = 1 - sol.y[0, -1] / F0[0]
    out = {"note": "Heat-release profile from calibrated RWGS+FT kinetics (feasibility/heat_profile_kinetic.py)",
           "GHSV_mL_g_h": g, "outlet_conversion": float(X), "q_inlet_W_m3": float(q[0]),
           "q_max_W_m3": float(q.max()), "z_at_q_max_m": float(z[q.argmax()]),
           "q_average_W_m3": float(trapz(q, z) / L), "superficial_velocity_m_s": float(u), "zones": zones}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    json.dump(out, open(os.path.join(HERE, "results", "heat_profile_kinetic.json"), "w"), indent=2)
    s2 = copy.deepcopy(s)
    s2["operating"]["superficial_velocity_m_s"] = float(u)
    s2["heat_profile_file"] = "feasibility/results/heat_profile_kinetic.json"
    json.dump(s2, open(os.path.join(HERE, "results", "settings_kinetic.json"), "w"), indent=2)
    print("GHSV %.0f mL/(g h), X %.3f, u %.4f m/s" % (g, X, u))
    print("q inlet %.0f, max %.0f at z = %.0f mm, average %.0f W/m3"
          % (q[0], q.max(), 1000 * z[q.argmax()], out["q_average_W_m3"]))
    print("first-order placeholder for comparison: inlet 397,515 / average 322,971 W/m3 at GHSV 6,000")


if __name__ == "__main__":
    main()
