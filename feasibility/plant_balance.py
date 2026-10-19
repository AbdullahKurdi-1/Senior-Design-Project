#!/usr/bin/env python3
"""
Plant mass balance with recycle: tests Spec 1, Spec 2, Integrated Spec 1 and the
+-0.5 % balance-closure requirement, and gives H2 use for the LCA.

Flowsheet (lumped):
  fresh CO2 + fresh H2 + recycle -> micro-structured reactor (320 C, 3.0 MPa, inlet H2/CO2 = 3)
  -> cooling: water and hydrocarbons removed (CH4 to fuel gas, C2+ to product)
  -> unreacted CO2 / CO / H2 recycled; a small purge removes impurities.
The catalyst mass is sized so that the single-pass CO2 conversion equals X_design.
"""
import json
import math
import os
import sys

from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kinetics as kin

HOURS = 8000.0


def solve_loop(F0_co2, k1, k2, a, dist, X_design=0.35, purge=0.02, T=593.15, P=30.0, tol=1e-10):
    d = dist["d_extra_H2"]

    def loop(W):
        R = [0.0, 0.0, 0.0]                       # recycle CO2, H2, CO
        for _ in range(500):
            co2_in = F0_co2 + R[0]
            h2_fresh = max(3.0 * co2_in - R[1], 0.0)
            Fin = [co2_in, R[1] + h2_fresh, R[2], 0.0, 0.0]
            s = kin.pfr(Fin, W, P, T, k1, k2, a, d, n_out=2)
            Fo = s.y[:, -1]
            Rn = [(1 - purge) * Fo[0], (1 - purge) * Fo[1], (1 - purge) * Fo[2]]
            if max(abs(x - y) for x, y in zip(Rn, R)) < tol * F0_co2:
                R = Rn
                break
            R = Rn
        X = 1 - Fo[0] / Fin[0]
        return X, Fin, Fo, h2_fresh

    W = brentq(lambda w: loop(w)[0] - X_design, 1e-3 * F0_co2, 1e6 * F0_co2, xtol=1e-9 * F0_co2)
    X, Fin, Fo, h2_fresh = loop(W)
    pf = purge
    hc = Fo[4]
    out = {
        "catalyst_kg": W, "single_pass_X": X, "h2_fresh_mol_s": h2_fresh,
        "inlet": Fin, "outlet": Fo, "hc_molC_s": hc,
        "purge_mol_s": [pf * Fo[0], pf * Fo[1], pf * Fo[2]],
    }
    # atom balances (C, H, O)
    C_in, C_out = F0_co2, hc + pf * (Fo[0] + Fo[2])
    H_in, H_out = 2 * h2_fresh, 2 * Fo[3] + hc * (2 + 2 * d) + 2 * pf * Fo[1]
    O_in, O_out = 2 * F0_co2, Fo[3] + pf * (2 * Fo[0] + Fo[2])
    out["closure_%"] = {k: 100 * (o - i) / i for k, i, o in (("C", C_in, C_out), ("H", H_in, H_out), ("O", O_in, O_out))}
    # total gas feed to the reactor -> GHSV
    Ftot = Fin[0] + Fin[1] + Fin[2]
    out["GHSV_mL_g_h"] = Ftot * 22414.0 * 3600.0 / (W * 1000.0)
    return out


def annual(res, dist, F0_co2):
    ts = HOURS * 3600.0 / 1e6                          # mol/s -> (t per yr) / (g/mol)
    hc = res["hc_molC_s"]
    eth = hc * dist["C2H4"] * 14.0 * ts
    ole = hc * dist["C2-C4 olefins"] * 14.0 * ts
    hc_mass = hc * (1 - dist["CH4"]) * 14.0 * ts + hc * dist["CH4"] * 16.0 * ts
    return {
        "CO2_feed_t_yr": F0_co2 * 44.0 * ts,
        "H2_fresh_t_yr": res["h2_fresh_mol_s"] * 2.016 * ts,
        "ethylene_t_yr": eth, "C2C4_olefins_t_yr": ole, "hydrocarbons_t_yr": hc_mass,
        "carbon_to_HC_%": 100 * hc / F0_co2,
        "H2_kg_per_kg_HC": res["h2_fresh_mol_s"] * 2.016 / (hc_mass / ts),
    }


def main():
    co2_t_yr = 150000.0
    F0 = co2_t_yr / 44.0 / (HOURS * 3600.0) * 1e6       # mol/s
    dist = kin.distribution()
    out = {"distribution": dist, "cases": []}
    for a in (0.3, 1.0, 3.0):
        for sco in (0.05, 0.10):
            try:
                k1, k2 = kin.calibrate(S_co_cal=sco, a=a, d=dist["d_extra_H2"])
            except RuntimeError:
                continue
            res = solve_loop(F0, k1, k2, a, dist)
            ann = annual(res, dist, F0)
            case = {"a": a, "S_CO_cal": sco, "k1": k1, "k2": k2,
                    "catalyst_t": res["catalyst_kg"] / 1000, "GHSV_mL_g_h": res["GHSV_mL_g_h"],
                    "single_pass_X": res["single_pass_X"], "closure_%": res["closure_%"], **ann}
            # CO2 needed for exactly 5,000 t/yr ethylene (linear in scale)
            case["CO2_needed_for_5000t_ethylene_t_yr"] = co2_t_yr * 5000.0 / ann["ethylene_t_yr"]
            out["cases"].append(case)
            print("a=%.1f S_CO=%.2f | catalyst %.0f t, GHSV %.0f | X %.3f | ethylene %.0f t/yr, "
                  "C2-C4= %.0f t/yr, HC %.0f t/yr | H2 %.0f t/yr (%.3f kg/kg HC) | C->HC %.1f %% | "
                  "closure C %.1e H %.1e O %.1e %%"
                  % (a, sco, case["catalyst_t"], case["GHSV_mL_g_h"], case["single_pass_X"], ann["ethylene_t_yr"],
                     ann["C2C4_olefins_t_yr"], ann["hydrocarbons_t_yr"], ann["H2_fresh_t_yr"], ann["H2_kg_per_kg_HC"],
                     ann["carbon_to_HC_%"], *res["closure_%"].values()))
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "plant_balance.json"), "w"), indent=2)


if __name__ == "__main__":
    os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), "results"), exist_ok=True)
    main()
