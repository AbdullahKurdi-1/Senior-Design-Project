#!/usr/bin/env python3
"""
Semi-mechanistic CO2-FTS kinetics for the feasibility study.

Reactions (lumped, Riedel-type rate forms for K-promoted Fe):
  (1) RWGS   CO2 + H2  <-> CO + H2O                     dH = +41.2 kJ/mol
  (2) FT     CO  + (2+d) H2 -> "CH2+d" + H2O            dH = -165 - 41*d kJ/mol (approx.)
      d = extra H2 per C from methane and paraffins (product-distribution dependent)

Rates [mol/(kg_cat s)], partial pressures in bar:
  r1 = k1 * (pCO2 - pCO*pH2O/(K1*pH2))              reversible RWGS, equilibrium K1(T)
  r2 = k2 * pH2*pCO / (pCO + a*pH2O + b*pCO2)       FT with water/CO2 competition (Riedel form), b = 0.1

k1, k2 are CALIBRATED to one published data point measured at the project design point
(320 C, 3.0 MPa, H2/CO2 = 3): CO2 conversion X_cal at GHSV_cal with CO selectivity S_CO_cal.
The water-inhibition constant `a` is not known for our catalyst, so every result is
repeated for a = 0.3, 1 and 3 (sensitivity).

Product distribution (carbon basis, CO-free hydrocarbons):
  CH4 selectivity s1 fixed; C2+ follows Anderson-Schulz-Flory with growth probability alpha;
  olefin fraction of C2-C4 from the olefin/paraffin ratio; ethylene from its own O/P ratio.
"""
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, fsolve

R = 8.314
SP = ["CO2", "H2", "CO", "H2O", "HC"]
B_CO2 = 0.1   # CO2 term in the FT denominator (Riedel form); keeps the rate finite at the inlet   # HC = moles of carbon in hydrocarbons


def K_rwgs(T):
    """RWGS equilibrium constant (inverse of Moe's WGS correlation)."""
    return 1.0 / math.exp(4577.8 / T - 4.33)


# ------------------------------------------------------------------ product distribution
def distribution(alpha=0.5, s_ch4=0.12, op_c3c4=4.0, op_c2=1.5, nmax=30):
    """Carbon-based selectivities within hydrocarbons (sum = 1)."""
    n = np.arange(2, nmax + 1)
    w = n * (1 - alpha) ** 2 * alpha ** (n - 1)
    w = w / w.sum() * (1 - s_ch4)            # C2+ share renormalised after methane
    f_ole_c2 = op_c2 / (1 + op_c2)
    f_ole = op_c3c4 / (1 + op_c3c4)
    c2, c3, c4 = w[0], w[1], w[2]
    out = {
        "CH4": s_ch4,
        "C2H4": c2 * f_ole_c2, "C2H6": c2 * (1 - f_ole_c2),
        "C3-C4 olefins": (c3 + c4) * f_ole, "C3-C4 paraffins": (c3 + c4) * (1 - f_ole),
        "C5+": float(w[3:].sum()),
    }
    out["C2-C4 olefins"] = out["C2H4"] + out["C3-C4 olefins"]
    # extra H2 per mol C relative to -CH2- : CH4 +1, paraffin CnH2n+2 +1/n, olefins 0, C5+ ~ +1/8
    d = s_ch4 * 1.0 + c2 * (1 - f_ole_c2) / 2 + c3 * (1 - f_ole) / 3 + c4 * (1 - f_ole) / 4 + out["C5+"] / 8
    out["d_extra_H2"] = d
    # mass of hydrocarbon product per mol C (kg): CH2 = 14 g; add 2 g per extra H2
    out["kg_HC_per_molC"] = (14.0 + 2.0 * d) / 1000.0
    return out


# ------------------------------------------------------------------ reactor
def rates(F, P_bar, T, k1, k2, a, d):
    Ft = F[0] + F[1] + F[2] + F[3]               # gas moles (HC removed from gas count: small)
    Ft = max(Ft, 1e-30)
    p = [max(F[i], 0.0) / Ft * P_bar for i in range(4)]
    pCO2, pH2, pCO, pH2O = p
    r1 = k1 * (pCO2 - pCO * pH2O / (K_rwgs(T) * max(pH2, 1e-12)))
    r2 = k2 * pH2 * pCO / (pCO + a * pH2O + B_CO2 * pCO2 + 1e-12)
    return r1, r2


def pfr(F0, W, P_bar, T, k1, k2, a, d, n_out=200):
    """Isothermal plug flow; F0 molar flows [mol/s]; W catalyst [kg]."""
    def rhs(w, F):
        r1, r2 = rates(F, P_bar, T, k1, k2, a, d)
        return [-r1, -r1 - (2 + d) * r2, r1 - r2, r1 + r2, r2]
    sol = solve_ivp(rhs, (0, W), F0, method="Radau", rtol=1e-8, atol=1e-12,
                    t_eval=np.linspace(0, W, n_out))
    return sol


def ghsv_feed(ghsv_ml_g_h, h2_co2=3.0):
    """Feed per kg catalyst [mol/s] for a given GHSV (mL STP per g per h)."""
    Ftot = ghsv_ml_g_h * 1000.0 / 22414.0 / 3600.0   # mL/(g h) -> mol/(kg s)
    co2 = Ftot / (1 + h2_co2)
    return [co2, co2 * h2_co2, 0.0, 0.0, 0.0]


def outlet(ghsv, k1, k2, a, d, T=593.15, P=30.0):
    F0 = ghsv_feed(ghsv)
    s = pfr(F0, 1.0, P, T, k1, k2, a, d, n_out=2)
    F = s.y[:, -1]
    X = 1 - F[0] / F0[0]
    conv = F0[0] - F[0]
    S_co = F[2] / conv if conv > 0 else 0.0
    return X, S_co, F


def calibrate(X_cal=0.472, S_co_cal=0.10, ghsv_cal=6000.0, a=1.0, d=0.2, T=593.15, P=30.0):
    """Find k1, k2 that reproduce X_cal and S_co_cal at ghsv_cal (nested 1-D root finding)."""
    def k1_for_X(lk2):
        f = lambda lk1: outlet(ghsv_cal, math.exp(lk1), math.exp(lk2), a, d, T, P)[0] - X_cal
        lo, hi = math.log(1e-6), math.log(1.0)
        if f(hi) < 0:            # X_cal not reachable with this k2
            return None
        return brentq(f, lo, hi, xtol=1e-6)

    def g(lk2):
        lk1 = k1_for_X(lk2)
        if lk1 is None:
            return 1.0           # treat as "CO selectivity too high"
        return outlet(ghsv_cal, math.exp(lk1), math.exp(lk2), a, d, T, P)[1] - S_co_cal

    lo, hi = math.log(1e-5), math.log(1.0)
    if g(lo) < 0 or g(hi) > 0:
        raise RuntimeError("target not reachable for a = %g, S_CO = %g" % (a, S_co_cal))
    lk2 = brentq(g, lo, hi, xtol=1e-6)
    lk1 = k1_for_X(lk2)
    return math.exp(lk1), math.exp(lk2)


def ghsv_for_conversion(X_target, k1, k2, a, d, T=593.15, P=30.0):
    f = lambda g: outlet(g, k1, k2, a, d, T, P)[0] - X_target
    return brentq(f, 200.0, 2.0e5)


if __name__ == "__main__":
    dist = distribution()
    print("Product distribution (alpha 0.5, CH4 12 %, O/P 4, C2 O/P 1.5):")
    for k, v in dist.items():
        print("  %-16s %.3f" % (k, v))
    for a in (0.3, 1.0, 3.0):
        for sco in (0.05, 0.10):
            try:
                k1, k2 = calibrate(S_co_cal=sco, a=a, d=dist["d_extra_H2"])
            except RuntimeError:
                print("a = %.1f  S_CO = %.2f : not reachable together with X = 47 %%" % (a, sco))
                continue
            g35 = ghsv_for_conversion(0.35, k1, k2, a, dist["d_extra_H2"])
            print("a = %.1f  S_CO = %.2f  k1 = %.3e  k2 = %.3e  GHSV for X = 35 %%: %.0f mL/(g h)" % (a, sco, k1, k2, g35))
