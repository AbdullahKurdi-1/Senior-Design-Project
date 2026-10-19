# `decision_tool/` – ISE decision-support tool (Student 5)

One script, one input file:

```bash
python3 decision_tool.py        # reads inputs.json, prints the report, writes results/*.png
```

| Module | Method | Spec it verifies |
|---|---|---|
| 1. Supply-chain optimization | Linear program (SciPy HiGHS): cheapest CO₂ and H₂ source mix that meets plant demand | Constraint 6 |
| 2. Supply reliability | Monte Carlo year-long outage simulation → smallest H₂ buffer giving ≥ 98 % feed availability | Spec 7 |
| 3. Techno-economics | 10 000-run Monte Carlo: NPV, IRR, unit production cost, tornado chart, break-even H₂ price | — (reported, not a pass/fail spec) |
| 4. Life-cycle GHG | Cradle-to-gate, mass allocation, renewable vs grid electricity, vs naphtha cracker (1.61 kg CO₂e/kg) | Integrated Spec 2, Constraint 5 |
| 5. SPC | I-MR control limits and RSD of the catalyst batch surface areas | Spec 8 |

## Important

* **All numbers in `inputs.json` marked PLACEHOLDER are illustrative.** Replace them with sourced values
  and cite each one. The plant block must come from the CHE mass balance; the batch data must come from
  the MSE BET measurements.
* With the placeholder numbers the NPV is negative. The literature shows the same for CO₂-derived olefins,
  so the report should present break-even H₂ and carbon prices, not assume profitability.
* The ≥ 45 % GHG reduction holds only with renewable electricity **and** with the CO₂ feedstock uptake
  credited (`credit_co2_uptake`). State this choice in the LCA goal-and-scope section, and also report
  the result without the credit.

## To do (Student 5)

1. Source every input (H₂ and CO₂ suppliers, distances, prices, emission factors).
2. Validate the LP on a hand-solvable 2-source case and show the Monte Carlo has converged (P50 stable within ±1 %).
3. Run the gauge R&R on repeated BET/XRD measurements and add it to the SPC section.
4. Optional: wrap the script in a simple Streamlit page for the final demo.
