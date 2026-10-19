#!/usr/bin/env python3
"""
1-D plug-flow model of one catalyst slot of the micro-structured fixed-bed reactor.

Purpose
-------
* Gives a quick, hand-checkable estimate of CO2 conversion and of the volumetric
  heat-release profile q(z) along the bed.
* Writes that heat-release profile to heat_profile.json, which the OpenFOAM case
  (../cfd) uses as its reaction heat source.

Model (deliberately simple, every assumption is listed)
------------------------------------------------------
1. Pseudo-homogeneous plug flow at constant T_ref and constant pressure.
2. Lumped first-order CO2 consumption:  r = k(T) * C_CO2   [mol / (m3_bed s)].
3. The rate constant k_ref is CALIBRATED so that the outlet conversion equals
   X_target at T_ref. Replace this with the fitted literature kinetics
   (Student 1 task) as soon as they are available.
4. Fixed selectivities: S_CO to CO (RWGS, endothermic, +41 kJ/mol) and
   1 - S_CO to hydrocarbons, taken as -CH2- units (CO2 + 3H2 -> CH2 + 2H2O,
   about -124 kJ/mol CH2, i.e. RWGS +41 combined with FT -165 kJ/mol).
5. Gas-volume change on reaction is neglected (conservative for heat release).

Run:  python3 plug_flow_model.py            (uses settings in ../cfd/case_settings.json)
"""
import json
import math
import os

import numpy as np

R = 8.314  # J/(mol K)
HERE = os.path.dirname(os.path.abspath(__file__))


def load_settings():
    with open(os.path.join(HERE, "..", "cfd", "case_settings.json")) as f:
        return json.load(f)


def run(s, verbose=True):
    op, geo, rx = s["operating"], s["geometry"], s["reaction"]
    L = geo["length_m"]
    T = rx["T_ref_K"]
    P = op["P_Pa"]
    u = op["superficial_velocity_m_s"]          # at T_ref, P
    y_co2 = 1.0 / (1.0 + op["H2_to_CO2_ratio"])
    X_target = rx["X_target_at_Tref"]
    S_co = rx["CO_selectivity"]

    C0 = y_co2 * P / (R * T)                       # mol CO2 / m3 gas
    k = -u * math.log(1.0 - X_target) / L          # 1/s, calibrated (assumption 3)
    dH = (1.0 - S_co) * 124e3 - S_co * 41e3        # J released per mol CO2 converted

    z = np.linspace(0.0, L, 601)
    X = 1.0 - np.exp(-k * z / u)
    q = dH * k * C0 * np.exp(-k * z / u)           # W/m3 of bed at T_ref

    # zone averages for the CFD heat source
    n = rx["n_zones"]
    edges = np.linspace(0.0, L, n + 1)
    zones = []
    for a, b in zip(edges[:-1], edges[1:]):
        # exact integral of the exponential profile over the zone
        qa = dH * k * C0 * (u / k) * (math.exp(-k * a / u) - math.exp(-k * b / u)) / (b - a)
        zones.append({"x0": a, "x1": b, "q_W_m3": qa})

    # GHSV check (mL STP per g catalyst per h)
    rho_bed = s["bed"]["bulk_density_kg_m3"]
    u_stp = u * (273.15 / T) * (P / 101325.0)
    ghsv = u_stp / L * 3600.0 * 1e6 / (rho_bed * 1e3)

    trapz = getattr(np, "trapezoid", None) or np.trapz   # numpy >= 2.0 / older
    q_avg = float(trapz(q, z) / L)
    out = {
        "note": "Heat-release profile at T_ref from reactor_model/plug_flow_model.py",
        "T_ref_K": T,
        "Ea_J_mol": rx["Ea_J_mol"],
        "k_ref_1_s": k,
        "C_CO2_inlet_mol_m3": C0,
        "heat_per_mol_CO2_J": dH,
        "outlet_conversion": float(X[-1]),
        "q_inlet_W_m3": float(q[0]),
        "q_average_W_m3": q_avg,
        "GHSV_mL_gcat_h": ghsv,
        "zones": zones,
    }
    if verbose:
        print("1-D plug-flow model (isothermal at T_ref = %.1f C)" % (T - 273.15))
        print("  GHSV                     : %8.0f mL/(g_cat h)" % ghsv)
        print("  CO2 inlet concentration  : %8.1f mol/m3" % C0)
        print("  calibrated k_ref         : %8.4f 1/s" % k)
        print("  outlet CO2 conversion    : %8.3f" % X[-1])
        print("  heat per mol CO2 conv.   : %8.1f kJ/mol" % (dH / 1e3))
        print("  q at inlet / average     : %8.0f / %8.0f W/m3" % (q[0], q_avg))
        # adiabatic temperature rise: why cooling is essential
        M = op["molWeight_kg_kmol"] / 1000.0
        rho_g = P * M / (R * T)
        dT_ad = q_avg * L / (rho_g * u * s["gas"]["Cp_J_kgK"])
        print("  adiabatic temperature rise (no cooling): %6.0f K" % dT_ad)
    return out


if __name__ == "__main__":
    s = load_settings()
    res = run(s)
    path = os.path.join(HERE, "heat_profile.json")
    with open(path, "w") as f:
        json.dump(res, f, indent=2)
    print("wrote", path)
