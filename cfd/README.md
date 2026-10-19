# `cfd/` – 3-D reactor-slot CFD with PID coolant control (Student 2, CHE)

This folder models **one catalyst slot** of the micro-structured fixed-bed reactor in OpenFOAM:

```
 coolant (boiling Dowtherm A ≈ 3.4 bar, h = 1500 W/m²K, temperature Ta set by the PID controller)
 ┌──────────────────────────────────────────────────────────┐
 │ steel plate, 1.5 mm          (region "solid")             │
 ├──────────────────────────────────────────────────────────┤ ← conjugate heat transfer
 │ porous catalyst bed, half of 4 mm   (region "fluid")      │
 │ H2/CO2 = 3 gas, 320 °C, 3.0 MPa → reaction heat in 60 zones│ → outlet
 └──────────────────────────────────────────────────────────┘   y = 0 symmetry plane
   inlet            length 150 mm, slot height 10 mm (z, 3-D)
```

* **Solver:** `chtMultiRegionFoam` (conjugate heat transfer, gas + steel).
* **Porous bed:** Darcy–Forchheimer resistance from the Ergun equation (particle size, porosity).
* **Reaction heat:** the axial profile comes from the calibrated RWGS + Fischer–Tropsch kinetics (`../reactor_model/heat_profile_kinetic.json`). Before every
  2-second chunk, each of the 60 bed zones (15 axial × 4 across the bed) gets the **exact Arrhenius**
  heat release for its current temperature. This is how a hot spot can feed on itself.
* **Control:** `run_pid.py` is a Python co-simulation. It runs OpenFOAM in 2 s chunks, reads the bed
  temperature, and moves the coolant temperature with a PI(D) law that holds the mean bed temperature
  at 320 °C. At t = 60 s it applies a **+10 % feed-flow step**.
* **No compiled code:** only standard OpenFOAM utilities plus Python, so it runs on any recent
  OpenFOAM.com (ESI) version. It was **tested with OpenFOAM v1912** (the Ubuntu 24.04 package).

## 1. Install (Windows → WSL Ubuntu, or native Ubuntu)

```bash
# Windows only, in PowerShell (admin):  wsl --install -d Ubuntu-24.04   then open Ubuntu
sudo apt update
sudo apt install -y openfoam python3-numpy python3-scipy python3-matplotlib paraview
echo 'source /usr/share/openfoam/etc/bashrc 2>/dev/null' >> ~/.bashrc   # load OpenFOAM in every terminal
source ~/.bashrc
blockMesh -help | head -3        # should print usage, not "command not found"
```

## 2. Run

```bash
git clone https://github.com/AbdullahKurdi-1/Senior-Design-Project.git && cd Senior-Design-Project
python3 reactor_model/plug_flow_model.py   # 1) heat-release profile -> reactor_model/heat_profile.json
cd cfd
python3 setup_case.py                      # 2) writes ./case and meshes it (≈ 15 000 cells, seconds)
python3 run_pid.py                         # 3) 160 s with PID + feed step (≈ 5–10 min on a laptop)
python3 compare_slot_width.py              # 4) optional design study: 4 / 10 / 25 mm beds (≈ 10–20 min)
```

Results land in `cfd/results/`:

| File | What it shows |
|---|---|
| `pid_response.png` | bed max / mean / min temperature and coolant temperature vs time, with the feed step |
| `pid_log.csv` | the same data, one row per 2 s |
| `summary.txt` | **spec check: max bed ΔT ≤ 15 K**, plus an energy balance (generation = coolant + gas heat) |
| `slot_width_study.png/.csv` | hot-spot size vs bed thickness: why the bed must be micro-structured |

**3-D pictures:** open `cfd/case/case.foam` in ParaView. Under *Mesh Regions* tick `fluid/internalMesh`
and `solid/internalMesh`, colour by `T`, and use a *Slice* normal to z.

To change anything (bed thickness, coolant h, conductivity, setpoint, PID gains, disturbance), edit
**`case_settings.json`** and run steps 2–3 again. Never edit files inside `case/`, because they are regenerated.

## 3. What Student 2 still has to do (do not skip – this is your contribution)

1. **Replace the placeholder kinetics.** In `reactor_model`, `k_ref` is calibrated so conversion = 35 %.
   Swap in the literature Fe–Zn–K kinetics that Student 1 validates.
2. **Mesh independence:** run with `cells_bed_y` = 8 / 16 / 32 and `cells_x` = 45 / 75 / 150, and show that ΔT changes by less than about 2 %.
3. **Validation by hand:** the radial temperature rise must match `q·b²/(2·k_eff)` (b = half-thickness).
   The energy balance in `summary.txt` must close within a few percent at steady state.
4. **Sensitivity:** vary `k_eff_W_mK` from 0.3 to 1.5 and `coolant_h_W_m2K` from 1000 to 10000. Report the worst case.
5. **Controller tuning:** justify `Kp`, `Ki` (e.g., from a step test with `--open-loop`), and show the
   response to +10 % and –10 % feed steps.

## 4. Model assumptions (state these in the report)

* The heat-release profile comes from a 1-D model with first-order lumped kinetics and fixed CO selectivity (15 %).
* The bed's heat capacity is the **gas only**. The real bed (particles) heats up about 20× more slowly, so
  transients here are faster than reality. Steady-state temperatures are unaffected.
* The coolant is a convective boundary (boiling Dowtherm A, h = 1,500 W/m²K); the coolant channel itself is not meshed.
  Boiling water was rejected: it needs ≈ 111 bar at 319 °C.
* The effective bed conductivity (0.8 W/m·K) is an assumption; see the sensitivity task above.
* The Arrhenius update is explicit, once per chunk (2 s). Reduce `chunk_s` if temperatures oscillate.
* The feed step scales the heat release 1:1 with flow, a conservative upper bound.
* The Arrhenius multiplier is capped at 1/X_target (35 % → 2.86×), so heat release never exceeds full
  conversion of the CO₂ fed. Without this cap a thick (runaway) bed heats without limit.
