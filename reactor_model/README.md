# `reactor_model/` – 1-D plug-flow model of one catalyst slot

`plug_flow_model.py` reads the operating conditions from `../cfd/case_settings.json` and prints:

* GHSV, CO₂ inlet concentration, outlet conversion
* heat released per mole of CO₂ converted (85 % to hydrocarbons, –124 kJ/mol CH₂; 15 % to CO, +41 kJ/mol)
* the inlet and average volumetric heat release, and the **adiabatic temperature rise without cooling (≈ 260 K)**.
  This number shows why the reactor must be cooled.

It writes `heat_profile.json`, the axial heat-release profile that the CFD case uses.

```bash
python3 plug_flow_model.py
```

## Hand check (base case)

| Quantity | Formula | Value |
|---|---|---|
| CO₂ concentration | y·P/(R·T) = 0.25 × 3.0 MPa / (8.314 × 593 K) | 152 mol/m³ |
| Rate constant | k = –u·ln(1–X)/L = 0.00917 × 0.431 / 0.15 | 0.0263 s⁻¹ |
| Inlet heat release | ΔH·k·C₀ = 99.2 kJ/mol × 0.0263 × 152 | ≈ 0.40 MW/m³ |
| GHSV | u_STP/L per kg of catalyst (bed 500 kg/m³) | ≈ 6000 mL/(g·h) |

## To do (Student 1)

1. Replace the calibrated first-order rate with **literature Fe–Zn–K kinetics** (RWGS + FT with olefin
   selectivity) and validate against at least 3 published datasets.
2. Use the **same kinetics** in Aspen Plus (RPlug) and check that the Python and Aspen conversions agree.
3. Add the selectivity outputs (C₂–C₄ olefins ≥ 40 %, Spec 2) once the kinetics are in.
