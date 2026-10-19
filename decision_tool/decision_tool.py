#!/usr/bin/env python3
"""
ISE decision-support tool (starter) for the CO2-to-ethylene plant.

Modules
  1. Supply-chain LP      : cheapest mix of CO2 and green-H2 sources that meets plant demand
  2. Supply reliability   : Monte Carlo year-long simulation of source outages -> buffer storage
                            needed for >= 98 % feed availability (Spec 7)
  3. Techno-economics     : Monte Carlo NPV / IRR / unit production cost + tornado sensitivity
  4. Life-cycle GHG       : cradle-to-gate kg CO2e per kg product vs naphtha cracker (Int. Spec 2)
  5. SPC                  : I-MR control limits and RSD for the MSE catalyst batches (Spec 8)

All inputs live in inputs.json. Values marked PLACEHOLDER must be replaced with sourced data.
Only numpy, scipy and matplotlib are needed (the LP uses SciPy's built-in HiGHS solver).

Usage:  python3 decision_tool.py            -> prints a report, writes results/*.png and results/report.txt
"""
import json
import os

import numpy as np
from scipy.optimize import brentq, linprog

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
rng = np.random.default_rng(42)


# ------------------------------------------------------------------ 1. supply-chain LP
def supply_lp(sources, demand, unit_key, transport_key):
    """min sum c_i x_i  s.t. sum x_i = demand, 0 <= x_i <= capacity_i (t/h)."""
    c = np.array([s[unit_key] + s[transport_key] for s in sources])
    res = linprog(c, A_eq=[np.ones(len(sources))], b_eq=[demand],
                  bounds=[(0, s["capacity_t_per_h"]) for s in sources], method="highs")
    if not res.success:
        raise RuntimeError("LP infeasible: total capacity below demand?")
    return res.x, c


# ------------------------------------------------------------------ 2. reliability simulation
def outage_profile(avail, hours, mean_outage_h=12.0):
    """Hourly up/down sequence for one source with the given long-run availability."""
    mtbf = mean_outage_h * avail / (1 - avail)
    up = np.ones(hours, dtype=bool)
    t = rng.exponential(mtbf)
    while t < hours:
        d = rng.exponential(mean_outage_h)
        up[int(t):int(min(hours, t + d))] = False
        t += d + rng.exponential(mtbf)
    return up


def availability_with_buffer(sources, alloc, demand, buffer_t, hours=8760, years=20):
    """Fraction of hours the plant feed is fully met, with a storage buffer (t)."""
    caps = np.array([s["capacity_t_per_h"] for s in sources])
    met = 0
    for _ in range(years):
        up = np.array([outage_profile(s["availability"], hours) for s in sources])
        store = buffer_t
        for h in range(hours):
            supply = np.sum(caps * up[:, h])          # spare capacity of running sources can cover gaps
            surplus = supply - demand
            if surplus >= 0:
                store = min(buffer_t, store + surplus)
                met += 1
            elif store + surplus >= 0:
                store += surplus
                met += 1
            else:
                store = 0.0
    return met / (hours * years)


# ------------------------------------------------------------------ 3. techno-economics
def tri(rng_, v, n):
    lo, mode, hi = v
    return rng_.triangular(lo, mode, hi, n)


def npv(rate, capex, cash, life):
    yrs = np.arange(1, life + 1)
    return -capex + np.sum(cash / (1 + rate) ** yrs)


def irr(capex, cash, life):
    f = lambda r: npv(r, capex, cash, life)
    try:
        return brentq(f, -0.99, 1.0)
    except ValueError:
        return np.nan


def economics(inp, h2_transport_usd_per_kg, co2_transport_usd_per_t):
    e, p = inp["economics"], inp["plant"]
    n, life, r = e["iterations"], e["plant_life_years"], e["discount_rate"]
    hrs = p["operating_hours_per_year"]
    h2_t = p["H2_demand_t_per_h"] * hrs
    co2_t = p["CO2_demand_t_per_h"] * hrs
    olef, hc = p["olefin_product_t_per_year"], p["hydrocarbon_product_t_per_year"]
    s = {k: tri(rng, e[k], n) for k in ("capex_usd", "fixed_opex_fraction_of_capex", "h2_price_usd_per_kg",
                                        "co2_price_usd_per_t", "olefin_price_usd_per_t",
                                        "byproduct_price_usd_per_t", "utilities_usd_per_t_product",
                                        "catalyst_cost_usd_per_year")}

    def one(i, v=None):
        g = (lambda k: v[k]) if v else (lambda k: s[k][i])
        opex = (h2_t * 1000 * (g("h2_price_usd_per_kg") + h2_transport_usd_per_kg)
                + co2_t * (g("co2_price_usd_per_t") + co2_transport_usd_per_t)
                + g("fixed_opex_fraction_of_capex") * g("capex_usd")
                + g("utilities_usd_per_t_product") * hc + g("catalyst_cost_usd_per_year"))
        rev = olef * g("olefin_price_usd_per_t") + (hc - olef) * g("byproduct_price_usd_per_t")
        cash = rev - opex
        crf = r * (1 + r) ** life / ((1 + r) ** life - 1)
        lcop = (g("capex_usd") * crf + opex - (hc - olef) * g("byproduct_price_usd_per_t")) / olef
        return npv(r, g("capex_usd"), cash, life), irr(g("capex_usd"), cash, life), lcop

    res = np.array([one(i) for i in range(n)])
    # tornado: each input low / high, others at mode
    mode = {k: e[k][1] for k in s}
    base = one(0, mode)[0]
    torn = []
    for k in s:
        lo, hi = dict(mode), dict(mode)
        lo[k], hi[k] = e[k][0], e[k][2]
        torn.append((k, one(0, lo)[0] - base, one(0, hi)[0] - base))
    # break-even H2 price (NPV = 0, other inputs at mode)
    def f(h2):
        m = dict(mode)
        m["h2_price_usd_per_kg"] = h2
        return one(0, m)[0]
    try:
        be = brentq(f, -10, 20)
    except ValueError:
        be = np.nan
    return res, torn, base, be


# ------------------------------------------------------------------ 4. LCA
def lca(inp):
    l, p = inp["lca"], inp["plant"]
    hc_kg = p["hydrocarbon_product_t_per_year"] * 1000
    hrs = p["operating_hours_per_year"]
    h2 = p["H2_demand_t_per_h"] * hrs * 1000 / hc_kg            # kg H2 / kg product (mass allocation)
    co2 = p["CO2_demand_t_per_h"] * hrs * 1000 / hc_kg          # kg CO2 fed / kg product
    uptake = 44.0 / 14.0 if l["credit_co2_uptake"] else 0.0     # carbon embodied in -CH2- product
    out = {}
    for name, ef in (("renewable", l["renewable_kgCO2e_per_kWh"]), ("grid", l["grid_kgCO2e_per_kWh"])):
        burden = (h2 * l["electricity_kWh_per_kg_H2"] + l["other_electricity_kWh_per_kg_product"]) * ef \
                 + co2 * l["co2_capture_kgCO2e_per_kg_co2"]
        out[name] = {"gross": burden, "net": burden - uptake}
    base = l["naphtha_cracker_baseline_kgCO2e_per_kg"]
    # break-even electricity intensity for a 45 % reduction
    kwh = h2 * l["electricity_kWh_per_kg_H2"] + l["other_electricity_kWh_per_kg_product"]
    be = (0.55 * base + uptake - co2 * l["co2_capture_kgCO2e_per_kg_co2"]) / kwh
    return out, base, h2, co2, be


# ------------------------------------------------------------------ 5. SPC
def spc(values, limit):
    x = np.array(values, float)
    mr = np.abs(np.diff(x))
    mrbar = mr.mean() if len(mr) else 0.0
    return {"mean": x.mean(), "rsd": x.std(ddof=1) / x.mean(), "UCL": x.mean() + 2.66 * mrbar,
            "LCL": x.mean() - 2.66 * mrbar, "MR_UCL": 3.267 * mrbar, "pass": x.std(ddof=1) / x.mean() <= limit}


def main():
    inp = json.load(open(os.path.join(HERE, "inputs.json")))
    os.makedirs(OUT, exist_ok=True)
    p = inp["plant"]
    rep = []

    # 1. LP
    x_c, c_c = supply_lp(inp["co2_sources"], p["CO2_demand_t_per_h"], "cost_usd_per_t", "transport_usd_per_t")
    x_h, c_h = supply_lp(inp["h2_sources"], p["H2_demand_t_per_h"], "cost_usd_per_kg", "transport_usd_per_kg")
    rep.append("1. SUPPLY-CHAIN LP (minimum delivered cost)")
    for s, x in zip(inp["co2_sources"], x_c):
        rep.append("   CO2  %-45s %6.2f t/h" % (s["name"], x))
    for s, x in zip(inp["h2_sources"], x_h):
        rep.append("   H2   %-45s %6.2f t/h" % (s["name"], x))
    co2_tr = float(np.dot(x_c, [s["transport_usd_per_t"] for s in inp["co2_sources"]]) / x_c.sum())
    h2_tr = float(np.dot(x_h, [s["transport_usd_per_kg"] for s in inp["h2_sources"]]) / x_h.sum())

    # 2. reliability -> buffer size
    rep.append("\n2. FEED AVAILABILITY vs H2 BUFFER STORAGE (target >= 98 %)")
    chosen = None
    for buf in (0, 5, 10, 20, 40, 80):
        a = availability_with_buffer(inp["h2_sources"], x_h, p["H2_demand_t_per_h"], buf, years=5)
        rep.append("   buffer %5.0f t H2  -> availability %.2f %%" % (buf, 100 * a))
        if chosen is None and a >= 0.98:
            chosen = buf
    rep.append("   smallest buffer meeting 98 %%: %s" % ("%g t" % chosen if chosen is not None else "none in range"))

    # 3. economics
    inp["economics"].setdefault("byproduct_price_usd_per_t", [400.0, 600.0, 800.0])
    res, torn, base, be = economics(inp, h2_tr, co2_tr)
    q = lambda a, pct: np.nanpercentile(a, pct)
    rep.append("\n3. TECHNO-ECONOMICS (Monte Carlo, %d runs)" % len(res))
    rep.append("   NPV  P10 / P50 / P90 : %8.1f / %8.1f / %8.1f  MUSD" % tuple(q(res[:, 0], k) / 1e6 for k in (10, 50, 90)))
    rep.append("   P(NPV > 0)           : %.1f %%" % (100 * np.mean(res[:, 0] > 0)))
    rep.append("   unit production cost P50 : %.0f USD/t olefin (net of by-product credit)" % q(res[:, 2], 50))
    rep.append("   break-even green-H2 price (NPV = 0): %.2f USD/kg" % be)

    # 4. LCA
    out, base_l, h2, co2, be_ef = lca(inp)
    rep.append("\n4. CRADLE-TO-GATE GHG (kg CO2e / kg product, mass allocation)")
    rep.append("   inventory: %.3f kg H2 and %.2f kg CO2 per kg product" % (h2, co2))
    for k, v in out.items():
        red = 100 * (1 - v["net"] / base_l)
        rep.append("   %-9s net %6.2f  gross %5.2f   reduction vs naphtha (%.2f): %6.1f %%" % (k, v["net"], v["gross"], base_l, red))
    rep.append("   break-even electricity intensity for a 45 %% reduction: %.3f kg CO2e/kWh" % be_ef)
    if inp["lca"]["credit_co2_uptake"]:
        rep.append("   NOTE: CO2 feedstock uptake credited (state this in the LCA goal & scope); set credit_co2_uptake false to see the result without it.")

    # 5. SPC
    sp = spc(inp["spc"]["bet_m2_per_g"], inp["spc"]["spec_limit_rsd"])
    rep.append("\n5. SPC - catalyst BET surface area")
    rep.append("   mean %.1f m2/g, RSD %.1f %% (limit %.0f %%) -> %s" % (sp["mean"], 100 * sp["rsd"], 100 * inp["spc"]["spec_limit_rsd"], "PASS" if sp["pass"] else "FAIL"))
    rep.append("   I-chart limits: LCL %.1f  UCL %.1f m2/g   MR-chart UCL %.1f" % (sp["LCL"], sp["UCL"], sp["MR_UCL"]))

    txt = "\n".join(rep)
    print(txt)
    open(os.path.join(OUT, "report.txt"), "w").write(txt + "\n")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].hist(res[:, 0] / 1e6, bins=60, color="tab:blue", alpha=0.8)
        ax[0].axvline(0, color="k", lw=0.8)
        ax[0].set_xlabel("NPV [MUSD]"); ax[0].set_ylabel("count"); ax[0].set_title("Monte Carlo NPV")
        torn.sort(key=lambda t: abs(t[2] - t[1]))
        names = [t[0] for t in torn]
        ax[1].barh(names, [t[1] / 1e6 for t in torn], color="tab:red", alpha=0.7, label="input low")
        ax[1].barh(names, [t[2] / 1e6 for t in torn], color="tab:green", alpha=0.7, label="input high")
        ax[1].set_xlabel("change in NPV [MUSD]"); ax[1].set_title("Tornado"); ax[1].legend(fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(OUT, "economics.png"), dpi=150)

        fig, ax = plt.subplots(figsize=(6, 3.5))
        v = inp["spc"]["bet_m2_per_g"]
        ax.plot(range(1, len(v) + 1), v, "o-")
        for y, ls in ((sp["UCL"], "--"), (sp["mean"], "-"), (sp["LCL"], "--")):
            ax.axhline(y, ls=ls, color="gray", lw=0.8)
        ax.set_xlabel("batch"); ax.set_ylabel("BET [m2/g]"); ax.set_title("I-chart: catalyst surface area")
        fig.tight_layout(); fig.savefig(os.path.join(OUT, "spc_bet.png"), dpi=150)
        print("plots written to", OUT)
    except ImportError:
        pass


if __name__ == "__main__":
    main()
