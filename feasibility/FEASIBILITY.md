# Feasibility dry-run: every specification and constraint tested before the project starts

Date: 19 Oct 2026. All scripts in this folder reproduce the numbers below. Tested with OpenFOAM v1912 and Python 3.13.

## Verdict

**The project is doable as specified.** No specification is impossible. Three items need a design decision or a
clearly stated basis now, so they don't surprise the team later. They are marked ⚠ below.

| # | Specification / constraint | Target | Dry-run result | Verdict |
|---|---|---|---|---|
| Spec 1 | Single-pass CO₂ conversion | ≥ 35 % | Reached by sizing. Kinetics calibrated to a published Fe–K catalyst (47 % at 6,000 mL/(g·h), 320 °C, 3.0 MPa) need GHSV 11,500–16,600 mL/(g·h). Conversion scales with activity × catalyst mass, so a catalyst only half as active still meets it with 2× catalyst. | ✅ Feasible, low risk |
| Spec 2 | C₂–C₄ olefin selectivity | ≥ 40 % | Base case 46.9 %. Met only if chain growth α ≤ 0.55 (at O/P ≥ 4). Literature Fe–Zn–K gives 54–57 % at 0.5 MPa; at 3.0 MPa α tends to rise. 48 % of the literature range meets the spec. | ⚠ Feasible, moderate risk |
| Spec 3 | Bed temperature difference | ≤ 15 K | 4 mm bed with real kinetics: 2.8 K (PID, +10 % and −10 % steps: 3.2 K). Worst case (k_eff 0.3 + coolant h 1,000): 9.8 K. Mesh independence: 2.71 / 2.79 / 2.83 K (coarse / base / fine). Energy balance closes to 0.01 %. **10 mm bed fails (34.5 K).** | ✅ Feasible with bed ≤ ~6 mm |
| Spec 4 | Catalyst surface area | ≥ 300 m²/g | At 15 wt % Fe the support must be ≥ 570 m²/g; at 10 wt % Fe ≥ 490 m²/g (MCF literature 500–1,000). | ⚠ Moderate risk; check the support first |
| Spec 5 | Fe-oxide crystallite size | ≤ 15 nm | Confined Fe oxides are 6–10 nm in the literature. XRD cannot size crystallites below about 3–4 nm (peaks too broad). | ✅ Low risk |
| Spec 6 | Crystallite growth, 24 h at 400 °C | ≤ 25 % | Aging temperature = calcination temperature, so growth should be small. The test is easy to pass, and a reviewer may call it weak. | ✅ Low risk; add a 450 °C sample |
| Spec 7 | Feed availability | ≥ 98 % | 20 t H₂ buffer gives 98.3 % (10 t: 96.9 %). 20 t at 200 bar is about 1,300 m³ of storage. | ✅ Feasible; storage is costly |
| Spec 8 | Batch BET variation | ≤ 15 % RSD | Not testable before synthesis. Lab-made MCF typically varies 5–10 %. | ✅ Low risk with controlled synthesis |
| Int 1 | Ethylene production | ≥ 5,000 t/yr | 150 kt/yr CO₂ → 8,080 t/yr ethylene, 21,550 t/yr C₂–C₄ olefins, 46,700 t/yr hydrocarbons (96 % C to products, 2 % purge). Only 93 kt/yr CO₂ is needed for 5,000 t. | ✅ Feasible, 60 % margin |
| Int 2 | GHG reduction vs naphtha | ≥ 45 % | Renewable H₂ **with** CO₂-uptake credit: > 100 %. **Without** the credit: 41 % (passes only if electricity is ≤ 0.028 kg CO₂e/kWh). Grid power: fails badly. | ⚠ Passes only on the stated basis |
| Int 3 | ≥ 8,000 h at X ≥ 35 % | — | Needs 1.25× (activity at end of run 80 %) to 3.3× (30 %) catalyst. The inventory is small (6–9 t), so over-design is cheap. No 8,000 h test exists in the literature (max ~500 h). | ✅ Feasible; prediction is uncertain |
| Mass balance | Closure | ± 0.5 % | C, H, O close to 10⁻⁹ % in the recycle model. | ✅ |
| C1 | ASME VIII / API 520 | — | Module shell 13–33 mm (ID 0.5–1.5 m). Relief area 1,496 mm² → standard API "L" orifice. **Boiling-water coolant at 319 °C needs 11.3 MPa.** | ⚠ Feasible; coolant needs a decision |
| C2–C5 | Design point, elements, lab safety, ISO LCA | — | Consistent with every test above. | ✅ |
| C6 | Saudi CO₂ and green H₂ | — | Needs 93 kt/yr CO₂ (17 % of the Jubail purification plant, which is already contracted, or < 1 % of the 9 Mt/yr CCS hub) and 13.7 kt/yr H₂ (6 % of NEOM). | ✅ Feasible; contracts are a risk |
| C7 | Budget SAR 5,000 | — | Chemicals and consumables ≈ SAR 2,900 (catalogue estimate). Characterization fees unknown. | ✅ Likely; get quotes |

## What the team must decide now (the ⚠ items)

1. **Coolant (Constraint 1, Spec 3).** Water boils at 319 °C only at 11.3 MPa, nearly 4× the process pressure.
   Options:
   - keep boiling water and design the coolant side for ~12 MPa;
   - switch to an organic heat-transfer fluid (Dowtherm A / Therminol type, a few bar; check vendor data).

   The CFD shows the bed still passes with a weak coolant (h = 1,000 W/m²K: 3.8 K), so the heat-transfer fluid option is safe for Spec 3.
2. **Bed thickness ≤ 6 mm (Spec 3).** With real kinetics the heat release is 2.4× the placeholder value (peak 1.0 MW/m³
   about 17 mm from the inlet), so a 10 mm bed fails. Keep the 4 mm design. The first 10–15 mm of bed is endothermic
   (RWGS runs before Fischer–Tropsch).
3. **LCA basis (Int 2).** Write in the goal and scope that the CO₂ feedstock taken up into the product is credited, and
   that H₂ comes from renewable electricity. Report the no-credit (41 %) and grid cases as well.
4. **Catalyst loading (Spec 4).** Measure the bare MCF surface area first. If it is below about 600 m²/g, reduce the
   loading to 10–12 wt % Fe.
5. **Selectivity (Spec 2).** Choose the kinetic and selectivity data for Fe–Zn–K at ≥ 1.5 MPa if available. If α > 0.55,
   the spec fails. Discuss this with the advisor early and record the α used.

## Other things that could surprise you later

* **Number of slots.** 9 t of catalyst is 18 m³ of bed. With 150 mm long × 0.5 m high slots that is ~60,000 slots;
  1 m long slots cut that to ~9,000. Student 1 should optimise slot length, and Student 2 must then rerun the CFD at that length.
* **Economics.** With the placeholder CAPEX (USD 200 M) the NPV is negative in every Monte Carlo run, and even free H₂
  does not break even. CAPEX must come from the real sizing (only 6–9 t catalyst, so it may be much lower). Present
  break-even prices; do not promise profitability.
* **Rubric coverage to confirm with the coaches.** MSE: the list includes "mechanical/property testing to ASTM" (we use
  ASTM D3663 BET and ASTM E1131 TGA) and "failure analysis" (the aged sample). ISE: "simulation alone is insufficient"
  (we use a working tool + SPC on real batch data). Get written confirmation that these count.
* **CFD model limits.** The CFD uses gas-only heat capacity (transients are faster than reality), an effective bed
  conductivity, and a convective coolant boundary. All are stated in `cfd/README.md`.
* **Kinetics are calibrated, not validated.** They are fitted to one published point. Student 1 must validate against at least 3
  datasets. The water-inhibition constant was varied 0.3–3 and the required GHSV changed by less than 40 %.

## How to reproduce

```bash
cd feasibility
python3 kinetics.py               # calibration for a = 0.3 / 1 / 3 and S_CO = 5 / 10 %
python3 selectivity_check.py      # Spec 2 map
python3 plant_balance.py          # Specs 1, Int 1, closure, H2 use   (~15 min)
python3 heat_profile_kinetic.py   # kinetic heat profile + CFD settings
python3 cfd_stress.py --jobs 2    # 11 CFD cases for Spec 3            (~15 min on 2 cores)
python3 catalyst_and_constraints.py
python3 ../decision_tool/decision_tool.py
```

| Output file | Contents |
|---|---|
| `results/plant_balance.json` | all kinetic cases: catalyst, GHSV, production, H₂, closure |
| `results/cfd_stress.csv` | ΔT, temperatures, energy balance for each CFD case |
| `results/heat_profile_kinetic.json` | real heat-release profile used by the CFD |
| `results/catalyst_and_constraints.txt` | Spec 4, Int 3, C1, coolant, C6, C7 numbers |
| `results/selectivity_grid.json` | Spec 2 scan |

## Basis and sources

Calibration point and selectivity anchors come from the literature review (FeK on F-doped mesoporous carbon, 47.2 %
conversion at 6,000 mL/(g·h), 320 °C, 3.0 MPa; Fe–Zn–K 51–55 % conversion with 53.6–57.4 % C₂–C₄ olefins at 320 °C).
Student 1 must cite the original papers when replacing the calibration. Steam saturation pressures are IAPWS values.
The RWGS equilibrium uses Moe's correlation. Prices are catalogue estimates.
