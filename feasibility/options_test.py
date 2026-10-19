#!/usr/bin/env python3
"""
Tests of the three open decisions, and the most favourable option for each.

  Item 1  coolant          boiling water vs boiling Dowtherm A (pressure, plate thickness, flow, hazards)
  Item 2  Spec 2 target    40 % vs 35 % C2-C4 olefins against published 2-3 MPa data; effect on ethylene output
  Item 3  catalyst loading 8 / 10 / 12 / 15 wt % Fe: surface area (Spec 4), activity, catalyst inventory

Every number that is an assumption is labelled ASSUMED in the printout.
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import kinetics as kin

OUT = []


def p(s=""):
    print(s)
    OUT.append(s)


# ====================================================================== ITEM 1
p("ITEM 1 - COOLANT AT 319 C")
# Dowtherm A: Clausius-Clapeyron through two published points (NOAA data sheet):
# 3.008 psia at 380 F (466.5 K) and the normal boiling point 257 C (530.15 K, 1.01325 bar)
P1, T1 = 3.008 * 0.0689476, (380 - 32) / 1.8 + 273.15
P2, T2 = 1.01325, 530.15
B = math.log(P2 / P1) / (1 / T1 - 1 / T2)          # dH_vap / R  [K]
p_dow = lambda T: P2 * math.exp(B * (1 / T2 - 1 / T))
Tc = 319 + 273.15
dHvap_dow = B * 8.314 / 0.166 / 1000                # kJ/kg, M = 0.166 kg/mol (eutectic)
p_water = 11.10                                     # MPa at 319 C (IAPWS: 11.29 at 320, 9.87 at 310)
p("  saturation pressure at 319 C : water %.1f MPa | Dowtherm A %.2f bar (Clausius-Clapeyron fit, dHvap %.0f kJ/kg)"
  % (p_water, p_dow(Tc), dHvap_dow))
p("  check of the fit at 400 C    : %.1f bar (Dow literature value is about 10 bar)" % p_dow(673.15))

# plate between process slot and coolant: clamped strip, sigma = dP b^2 / (2 t^2)
S, b = 100e6, 0.010                                  # ASSUMED allowable stress SS316 ~350 C, rib span = slot height 10 mm
for name, dp in (("water (11.1 MPa coolant, process depressurised)", p_water * 1e6),
                 ("Dowtherm A (process 3.3 MPa at relief, coolant 0.34 MPa)", 3.3e6 - p_dow(Tc) * 1e5)):
    t = math.sqrt(dp * b ** 2 / (2 * S))
    p("  plate thickness needed, %-58s: %.2f mm (worst dP %.1f MPa)" % (name, 1000 * t, dp / 1e6))

# reactor duty from the plant balance
pb = json.load(open(os.path.join(HERE, "results", "plant_balance.json")))
c0 = next(c for c in pb["cases"] if c["a"] == 1.0 and c["S_CO_cal"] == 0.05)
d = pb["distribution"]["d_extra_H2"]
molC = c0["hydrocarbons_t_yr"] * 1e6 / (14.0 + 2 * d) / (8000 * 3600)       # mol C/s to hydrocarbons
duty = molC * (165e3 + 41e3 * d - 41.2e3)                                    # W (FT exotherm minus RWGS)
p("  reactor heat duty            : %.1f MW (from the plant balance)" % (duty / 1e6))
p("  coolant evaporated           : water %.1f kg/s (dHvap ~1,240 kJ/kg) | Dowtherm A %.0f kg/s" % (duty / 1.24e6, duty / (dHvap_dow * 1e3)))
p("  heat use                     : water -> 111-bar steam (high value) | Dowtherm vapour -> condensed in a steam")
p("                                 generator making ~ 10-40 bar steam")
p("  hazards                      : water - none (but 111-bar drum) | Dowtherm A - flash point 124 C, autoignition 621 C,")
p("                                 runs above its flash point: leaks are a fire hazard -> HAZOP node, leak detection")
p("  CFD check (bed dT with boiling Dowtherm, h 1,500-2,500 W/m2K): see results/cfd_coolant.csv")

# ====================================================================== ITEM 2
p("\nITEM 2 - SPEC 2 TARGET (C2-C4 olefins, CO-free hydrocarbons) vs PUBLISHED DATA AT >= 1 MPa")
data = [
    ("Na/Fe@Co core-shell", 330, 3.0, 50.3, "33.8 (CO 4.6 %, so ~35 CO-free)", 35.4),
    ("FeCoK/ZrO2", 350, 3.0, 36.0, "28", 28.0),
    ("Catalyst in review (320 C, 3 MPa, 12,000 mL/g h)", 320, 3.0, 48.0, "29.0 mmol C/(g h) = ~45 of converted C", 45.0),
    ("(Fe1-xCox)5C2 Na-promoted / CNT (Kim et al.)", 340, 1.0, 34.0, "~39", 39.0),
    ("ZnCoxFe2-xO4 spinel (Xu et al.)", None, None, 49.6, "36.1", 36.1),
    ("Fe-Zn-K (Zhang et al.) - low pressure", 320, 0.5, 51.0, "53.6-57.4", 55.0),
    ("Fe7C3-KMg (2025)", 340, 2.0, 41.5, "67.1 high-value olefins incl. C5+ alpha-olefins", None),
]
p("  %-48s %5s %5s %6s  %s" % ("catalyst", "T C", "MPa", "X %", "light-olefin selectivity %"))
vals = []
for n, T, P, X, s, v in data:
    p("  %-48s %5s %5s %6.1f  %s" % (n, T or "n/r", P or "n/r", X, s))
    if v is not None and (P is None or P >= 1.0):
        vals.append(v)
vals.sort()
med = vals[len(vals) // 2]
p("  data at >= 1 MPa (or unstated): %s  -> median %.0f %%" % (", ".join("%.0f" % v for v in vals), med))
p("  share of these points meeting 40 %%: %d of %d   meeting 35 %%: %d of %d"
  % (sum(v >= 40 for v in vals), len(vals), sum(v >= 35 for v in vals), len(vals)))
for tgt in (0.40, 0.35):
    amax = max(a / 100 for a in range(40, 75) if kin.distribution(alpha=a / 100, s_ch4=0.12, op_c3c4=4.0)["C2-C4 olefins"] >= tgt)
    p("  target %.0f %%: met up to chain-growth alpha = %.2f (CH4 12 %%, O/P 4)" % (100 * tgt, amax))
# ethylene output if the catalyst sits at the 35 % edge
base = c0
for a in (0.50, 0.60, 0.62):
    dist = kin.distribution(alpha=a, s_ch4=0.12, op_c3c4=4.0)
    hc = base["hydrocarbons_t_yr"]
    eth = hc * dist["C2H4"] * 14 / (14 + 2 * d)
    co2_5000 = 150000 * 5000 / eth
    p("  alpha %.2f: C2-C4 olefins %.1f %%, ethylene %.1f %% -> %.0f t/yr from 150 kt CO2 (CO2 for 5,000 t: %.0f kt/yr)"
      % (a, 100 * dist["C2-C4 olefins"], 100 * dist["C2H4"], eth, co2_5000 / 1e3))

# ====================================================================== ITEM 3
p("\nITEM 3 - CATALYST LOADING (Fe:Zn = 2:1 molar, 1.5 wt % K)")
p("  ASSUMED: 20 % pore blockage; Fe-oxide crystallite 6/7/8/9 nm at 8/10/12/15 wt % Fe (confined MCF literature 6-10 nm);")
p("           activity per gram ~ Fe loading / crystallite size (surface Fe); reference = 15 wt % (6-9 t catalyst at that activity)")


def oxides(fe):
    fe2o3 = fe * 159.69 / (2 * 55.845)
    zno = fe / 55.845 / 2 * 81.38
    k2o = 1.5 * 94.2 / (2 * 39.098)
    return fe2o3 + zno + k2o


ref_act = 15 / 9.0
rows = []
p("  Fe wt% | BET @ MCF 500 / 600 / 700 m2/g | min MCF for 300 | rel. activity | catalyst t | crystallite <= 15 nm")
for fe, dnm in ((8, 6.0), (10, 7.0), (12, 8.0), (15, 9.0)):
    si = (100 - oxides(fe)) / 100
    bets = [si * s * 0.8 for s in (500, 600, 700)]
    act = (fe / dnm) / ref_act
    cat = 9.0 / act
    rows.append((fe, bets, act, cat))
    p("  %5d  |      %4.0f / %4.0f / %4.0f           |      %4.0f       |     %4.2f      |   %4.1f     |  yes (%d nm)"
      % (fe, *bets, 300 / (si * 0.8), act, cat, dnm))
best = max((r for r in rows if r[1][0] >= 300), key=lambda r: r[2])
p("  most favourable: %d wt %% Fe - highest activity that still passes Spec 4 even on a low-end 500 m2/g support"
  % best[0])
p("  (12 wt %% Fe gives only %.0f %% more activity but fails on a 500 m2/g support; use it only if the bare MCF measures >= 600 m2/g)"
  % (100 * (rows[2][2] / best[2] - 1)))

open(os.path.join(HERE, "results", "options_test.txt"), "w").write("\n".join(OUT) + "\n")
