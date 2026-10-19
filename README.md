# CO₂-to-Ethylene Micro-Structured Reactor: Digital Twin

KFUPM TEAM Design Project (Term 261)
**Design and Simulation of a Micro-Structured Fixed-Bed Reactor for Carbon-Circularity Ethylene Production via Catalytic CO₂ Hydrogenation**

This repository is the **computer-based half of the team's hybrid prototype**. The physical half is the
Fe–Zn–K/MCF catalyst that the MSE students synthesize and test in the lab. Everything here is open-source
code that the team writes and owns. Commercial simulators (Aspen Plus) are used only to *validate* it.

| Folder | Owner | Program | What it does | What you get |
|---|---|---|---|---|
| [`reactor_model/`](reactor_model) | CHE Student 1 (with Student 2) | Python | 1-D plug-flow model of one catalyst slot: CO₂ conversion and heat-release profile | conversion and heat-release numbers, `heat_profile.json` for the CFD |
| [`cfd/`](cfd) | CHE Student 2 | OpenFOAM + Python | 3-D conjugate heat transfer of one slot (porous bed + steel wall + coolant) with a PID controller on the coolant temperature | temperature maps (ParaView), bed-temperature control plot, **ΔT ≤ 15 K check** |
| [`decision_tool/`](decision_tool) | ISE Student 5 | Python | supply-chain LP, supply-reliability simulation, Monte Carlo NPV/IRR, LCA, SPC of catalyst batches | supply plan, buffer size, cost and NPV distributions, GHG comparison, control charts |

Measured catalyst data from the MSE students (surface area, crystallite size, bed density, aging) enter through
`cfd/case_settings.json` (bed) and `decision_tool/inputs.json` (SPC, catalyst life and cost).

## Quick start (Ubuntu 24.04 or Windows + WSL)

```bash
sudo apt update && sudo apt install -y openfoam python3-numpy python3-scipy python3-matplotlib paraview git
echo 'source /usr/share/openfoam/etc/bashrc 2>/dev/null' >> ~/.bashrc && source ~/.bashrc

git clone https://github.com/AbdullahKurdi-1/Senior-Design-Project.git && cd Senior-Design-Project
python3 reactor_model/plug_flow_model.py     # heat-release profile
cd cfd && python3 setup_case.py && python3 run_pid.py && cd ..   # CFD + PID (≈ 5–10 min)
python3 decision_tool/decision_tool.py       # ISE analyses (seconds)
```

Reference outputs from a test run are in `cfd/results/` and `decision_tool/results/`.

## Feasibility dry-run

[`feasibility/FEASIBILITY.md`](feasibility/FEASIBILITY.md) tests every specification and constraint (calibrated kinetics, plant mass balance with recycle, 11 CFD stress cases, catalyst/pressure/budget checks) and lists the decisions the team must make early.

## Status: starter, not final

* The code runs end to end. It was tested with OpenFOAM v1912 and Python 3.12.
* **Values marked PLACEHOLDER** (kinetics calibration, costs, emission factors, batch data) must be
  replaced with sourced or measured data. Each owner's README lists the remaining tasks that make up
  their senior-design contribution (validation, mesh study, sensitivity, tuning).
