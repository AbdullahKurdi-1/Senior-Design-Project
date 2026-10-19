#!/usr/bin/env python3
"""
Feasibility checks that do not need CFD:
  Spec 4  catalyst BET >= 300 m2/g        - dilution / pore-blockage estimate vs MCF quality and loading
  Int 3   >= 8,000 h at X >= 35 %         - catalyst over-design needed for different end-of-run activities
  C1      ASME VIII shell thickness and API 520 relief area
  C1/S3   coolant: saturation pressure of water at the coolant temperature
  C6      CO2 / H2 demand vs Saudi supply
  C7      lab budget estimate
All inputs are stated in the printout; prices are ESTIMATES to be replaced with supplier quotes.
"""
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = []


def p(s=""):
    print(s)
    OUT.append(s)


# ------------------------------------------------------------------ Spec 4: surface area
def composition(fe_wt, fe_zn_molar=2.0, k_wt=1.5):
    fe2o3 = fe_wt * 159.69 / (2 * 55.845)
    zn = fe_wt / 55.845 / fe_zn_molar * 65.38
    zno = zn * 81.38 / 65.38
    k2o = k_wt * 94.2 / (2 * 39.098)
    oxides = fe2o3 + zno + k2o
    return oxides, 100 - oxides


p("SPEC 4 - estimated BET of the finished catalyst [m2/g] (silica dilution x pore blockage)")
p("  Fe wt% | oxide wt% | MCF 500 m2/g      | MCF 650 m2/g      | MCF 800 m2/g")
p("         |           | block 10% / 25%   | block 10% / 25%   | block 10% / 25%")
for fe in (10.0, 15.0):
    ox, si = composition(fe)
    cells = []
    for s_mcf in (500, 650, 800):
        vals = [si / 100 * s_mcf * (1 - b) for b in (0.10, 0.25)]
        cells.append("%5.0f / %5.0f     " % tuple(vals))
    p("  %5.1f  |   %5.1f   | %s| %s| %s" % (fe, ox, *cells))
for fe in (10.0, 15.0):
    ox, si = composition(fe)
    p("  minimum MCF support area for >= 300 m2/g at %.0f wt%% Fe (20 %% blockage): %.0f m2/g" % (fe, 300 / (si / 100 * 0.8)))

# ------------------------------------------------------------------ Int 3: run length
p("\nINTEGRATED SPEC 3 - catalyst over-design to keep X >= 35 % for 8,000 h")
p("  conversion depends on (activity x catalyst mass / feed), so the extra catalyst needed is 1 / (end-of-run activity)")
for a_end in (0.8, 0.6, 0.5, 0.3):
    p("  end-of-run activity %3.0f %% of fresh -> install %.2f x the fresh-catalyst mass (+%.0f %%)" % (100 * a_end, 1 / a_end, 100 * (1 / a_end - 1)))
p("  literature: Zn-promoted Fe stable > 200 h; Ba-Na-Fe > 500 h; no published 8,000 h CO2-FTS test was found")

# ------------------------------------------------------------------ C1: pressure design
p("\nCONSTRAINT 1 - pressure design")
P_set = 3.3e6
S_allow, E, CA = 100e6, 0.85, 3.0
for D in (0.5, 1.0, 1.5):
    t = P_set * (D / 2) / (S_allow * E - 0.6 * P_set) * 1000 + CA
    p("  module shell ID %.1f m: ASME VIII UG-27 t = PR/(SE-0.6P) + CA = %.1f mm (S = 100 MPa SS316 at ~350 C assumed, E = 0.85, CA = 3 mm)" % (D, t))
# API 520 critical gas flow, SI form: A[mm2] = W / (C Kd P1 Kb Kc) * sqrt(T Z / M)
W = 18750 + 2560                 # kg/h fresh feed (blocked-outlet case)
n_co2, n_h2 = 18750 / 44.01, 2560 / 2.016
M = W / (n_co2 + n_h2)
k = 1.38
C = 0.03948 * math.sqrt(k * (2 / (k + 1)) ** ((k + 1) / (k - 1)))
P1 = P_set * 1.10 / 1000 + 101.3      # kPa abs, 10 % overpressure
A = W / (C * 0.975 * P1 * 1.0 * 1.0) * math.sqrt(593.15 * 1.0 / M)
orifices = [("J", 830), ("K", 1186), ("L", 1841), ("M", 2323), ("N", 2800)]
letter = next(n for n, a in orifices if a >= A)
p("  API 520 relief, blocked outlet with full fresh feed %.0f kg/h (M = %.1f, k = %.2f): A = %.0f mm2 -> API orifice '%s'" % (W, M, k, A, letter))

# ------------------------------------------------------------------ coolant pressure
p("\nCOOLANT CHOICE (affects Spec 3 and Constraint 1)")
psat = {280: 6.42, 300: 8.59, 310: 9.87, 320: 11.29}     # MPa, IAPWS steam tables
for T, ps in psat.items():
    p("  boiling water at %d C needs %.1f MPa steam-drum pressure (process side 3.0 MPa)" % (T, ps))
p("  -> boiling water at ~319 C means coolant-side pressure ~3.7x the process pressure. Alternatives: an organic heat-transfer")
p("     fluid (Dowtherm A / Therminol type; low-pressure, check vendor vapour-pressure data) with h ~ 1,000-3,000 W/m2K, which the")
p("     CFD stress tests cover through the h = 1,000 case.")

# ------------------------------------------------------------------ C6: supply
p("\nCONSTRAINT 6 - feedstock vs Saudi supply")
pb = os.path.join(HERE, "results", "plant_balance.json")
if os.path.exists(pb):
    cases = json.load(open(pb))["cases"]
    co2 = [c["CO2_needed_for_5000t_ethylene_t_yr"] for c in cases]
    h2 = [c["H2_fresh_t_yr"] * c["CO2_needed_for_5000t_ethylene_t_yr"] / c["CO2_feed_t_yr"] for c in cases]
    p("  CO2 needed for 5,000 t/yr ethylene: %.0f - %.0f kt/yr (all kinetic cases)" % (min(co2) / 1e3, max(co2) / 1e3))
    p("  H2 needed                         : %.1f - %.1f kt/yr" % (min(h2) / 1e3, max(h2) / 1e3))
p("  SABIC/Linde Jubail purification ~1,500 t/d (~550 kt/yr, already supplied to methanol/urea); Aramco-Linde-SLB CCS hub up to 9 Mt/yr by 2027")
p("  NEOM Green Hydrogen: up to 219 kt/yr (contracted as ammonia export); first production expected around 2027")

# ------------------------------------------------------------------ C7: budget
p("\nCONSTRAINT 7 - laboratory budget ESTIMATE (SAR, catalogue-type prices, replace with quotes)")
items = [("Pluronic P123, 250 g", 650), ("TEOS 98 %, 1 L", 450), ("1,3,5-trimethylbenzene, 500 mL", 300),
         ("Fe(NO3)3.9H2O, 500 g", 250), ("Zn(NO3)2.6H2O, 500 g", 250), ("KNO3 or K2CO3, 500 g", 150),
         ("HCl 37 %, 1 L", 150), ("crucibles, filter paper, glassware, gas for aging (N2)", 700)]
tot = 0
for n, c in items:
    p("  %-58s %6d" % (n, c))
    tot += c
p("  %-58s %6d" % ("TOTAL chemicals and consumables", tot))
p("  characterization (XRD, BET, SEM/EDS, TGA) for 3 batches + 2 aged samples = about 20-25 measurements:")
p("  internal KFUPM rates decide whether SAR 5,000 is enough; ask the lab managers for their student rates")

open(os.path.join(HERE, "results", "catalyst_and_constraints.txt"), "w").write("\n".join(OUT) + "\n")
