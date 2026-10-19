#!/usr/bin/env python3
"""
Spec 2 feasibility: C2-C4 olefin selectivity >= 40 % (carbon basis, CO-free hydrocarbons).

The selectivity is set by the catalyst's product distribution, not by reactor design.
We scan the literature range for K-promoted Fe catalysts:
  chain-growth probability alpha 0.40 - 0.70
  methane selectivity        0.08 - 0.25
  C3-C4 olefin/paraffin      2 - 7
and report where the spec is met, plus the ethylene share of the hydrocarbons.
"""
import itertools
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kinetics as kin

alphas = np.round(np.arange(0.40, 0.701, 0.05), 2)
ch4s = [0.08, 0.12, 0.16, 0.20, 0.25]
ops = [2.0, 3.0, 4.0, 5.0, 7.0]

rows = []
for a, m, op in itertools.product(alphas, ch4s, ops):
    d = kin.distribution(alpha=a, s_ch4=m, op_c3c4=op, op_c2=max(1.0, op / 2.5))
    rows.append({"alpha": float(a), "CH4": m, "OP": op, "C2C4_olefins": d["C2-C4 olefins"], "ethylene": d["C2H4"]})

ok = [r for r in rows if r["C2C4_olefins"] >= 0.40]
print("grid points: %d, meeting >= 40 %%: %d (%.0f %%)" % (len(rows), len(ok), 100 * len(ok) / len(rows)))
print("\nC2-C4 olefin selectivity [%] at CH4 = 12 %, rows = alpha, cols = O/P", ops)
for a in alphas:
    vals = [next(r for r in rows if r["alpha"] == a and r["CH4"] == 0.12 and r["OP"] == op)["C2C4_olefins"] for op in ops]
    print("  alpha %.2f : %s" % (a, "  ".join("%5.1f" % (100 * v) for v in vals)))
print("\nmaximum alpha that still meets 40 % (CH4 = 12 %, O/P = 4):",
      max(r["alpha"] for r in rows if r["CH4"] == 0.12 and r["OP"] == 4.0 and r["C2C4_olefins"] >= 0.40))
eth = [r["ethylene"] for r in ok]
print("ethylene share of hydrocarbons where the spec is met: %.1f - %.1f %%" % (100 * min(eth), 100 * max(eth)))
base = kin.distribution()
print("base case (alpha 0.5, CH4 12 %%, O/P 4): C2-C4 olefins %.1f %%, ethylene %.1f %%"
      % (100 * base["C2-C4 olefins"], 100 * base["C2H4"]))
print("literature anchor: Fe-Zn-K at 320 C reported 53.6-57.4 % C2-C4 olefins in hydrocarbons")
os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), "results"), exist_ok=True)
json.dump(rows, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "selectivity_grid.json"), "w"), indent=1)
