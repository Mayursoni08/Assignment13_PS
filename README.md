# Stepped Distance Protection — EE 315 Assignment 13

Design and verification of a three-zone mho distance relay for a 220 kV, 50 Hz
transmission line, built in Python with an interactive web demo.

**Project page (live demo):** https://Mayursoni08.github.io/Assignment13_PS/

## Group

| Member | Roll Number | Name |
|---|---|---|
| Member 1 | 240002007 | Akarsh J |
| Member 2 | 240002073 | Sreehari J. |
| Member 3 | 240002043 | Naman V Shetty |
| Member 4 | 240002037 | **Mayur Soni** (Project Lead) |
| Member 5 | 240002017 | Bhasuru Nikhil |
| Member 6 | 240002069 | Ridhi Shakkar |

## Problem

Implement the stepped distance protection scheme for a power system and verify
it using simulation data. A distance relay measures the apparent impedance
Z = V/I at its location; this is proportional to the distance to the fault, so
the relay can decide whether the fault lies inside its reach. "Stepped" means
several zones of increasing reach and increasing time delay: Zone 1 clears
faults on the protected line instantly, Zones 2 and 3 give time-delayed backup
for the line and the adjacent lines.

**Test system (assumed standard example system):**
- Relay R at Bus A protects line AB (80 km); adjacent lines BC (60 km) and BD (100 km)
- 220 kV (line-line), 50 Hz, source level 5000 MVA (Zs1 = 9.68∠85° Ω, Zs0 = 12.58∠80° Ω)
- Line impedance: Z1 = 0.08 + j0.40 Ω/km (78.7°), Z0 = 0.30 + j1.20 Ω/km
- CT 1000/1 A, PT 220 kV/110 V (secondary impedance = 0.5 × primary)
- 32 samples/cycle (1.6 kHz)

## Zone settings

| Zone | Rule | Reach (primary) | Reach (secondary) | Reach (km) | Delay | Trip time* |
|---|---|---|---|---|---|---|
| Z1 | 80 % of line AB | 26.11 Ω | 13.05 Ω | 64 | 0 s | 0.04 s |
| Z2 | AB + 50 % of shortest adjacent line (BC) | 44.87 Ω | 22.44 Ω | 110 | 0.4 s | 0.44 s |
| Z3 | AB + 120 % of longest adjacent line (BD) | 81.58 Ω | 40.79 Ω | 200 | 1.0 s | 1.04 s |

\* includes 2 cycles (40 ms) for the data window and decision. The relay uses a
mho characteristic (circle through the origin) with its maximum-torque angle
equal to the line angle, 78.7°.

## Files

| File | Description |
|---|---|
| `ee315_assignment13.py` | Main program: fault calculation, waveform synthesis, DFT phasor estimation, mho zone decision, plots |
| `results_sweep.csv` | Results for all 204 simulated faults (location, R, X, zone, ideal zone, trip time) |
| `results_meta.json` | Zone settings and summary statistics |
| `figs/` | Generated figures (system, waveforms, R-X diagram, staircase, accuracy, fault-resistance map) |
| `docs/index.html` | Interactive web demo (served by GitHub Pages) |
| `report/` | LaTeX report (`main.tex`, figures, tables) — upload this folder to Overleaf |
| `requirements.txt` | Python dependencies (numpy, matplotlib) |

## How to run

```
pip install -r requirements.txt
python ee315_assignment13.py
```

This regenerates every figure in `figs/` and the files `results_sweep.csv` and
`results_meta.json`. The run is reproducible (fixed random seed).

## How it works

1. **Fault calculation** — sequence networks give relay voltages and currents for
   A-G, B-C and A-B-C faults at any distance and fault resistance.
2. **Waveform synthesis** — sampled three-phase signals at the CT/PT secondary:
   two pre-fault cycles, then the fault with a decaying DC offset, 2 % third
   harmonic, measurement noise (0.1 % of rated peak) and a random inception angle.
3. **Phasor estimation** — full-cycle DFT (Lec 33) in the window 20–40 ms after
   inception. A digital mimic filter removes the DC offset.
4. **Apparent impedance** — line-line loop for phase faults; phase-ground loop
   with zero-sequence compensation k0 = (Z0 − Z1)/(3 Z1) for A-G faults.
5. **Zone decision and trip time** — mho test against Z1, Z2, Z3; trip time =
   40 ms + zone delay.

## Web demo

`docs/index.html` reproduces the sequence-network model in JavaScript. Choose the
fault type, distance and fault resistance and see the apparent impedance on the
R-X diagram, the zone that picks up, the trip time and the time-distance
staircase. It shows the steady-state result, so it has no noise or DC offset.

## Results

204 faults (A-G, B-C, A-B-C) were applied at 5 % steps on lines AB and BC and up
to 140 % on BD (the part beyond Bus D is a hypothetical extension used to test
non-operation).

| Fault location | Distance | Zone picked | Trip time |
|---|---|---|---|
| Line AB, 10 % | 8 km | Z1 | 0.04 s |
| Line AB, 50 % | 40 km | Z1 | 0.04 s |
| Line AB, 90 % | 72 km | Z2 | 0.44 s |
| Line BC, 100 % | 140 km | Z3 | 1.04 s |
| Line BD, 100 % | 180 km | Z3 | 1.04 s |
| Beyond Zone 3 (BD, 125 %) | 205 km | none | no trip |

| Check | Result |
|---|---|
| Zone matches ideal zone | 94.1 % of 204 faults (all mismatches lie exactly on a reach boundary) |
| Impedance error, DFT only | mean ≈ 4.6 %, max ≈ 8.2 % |
| Impedance error, DFT + mimic filter | mean ≈ 0.1 %, max ≈ 0.8 % |
| Max Rf for instantaneous Z1 trip (A-G) | 20 Ω at 16 and 32 km, 15 Ω at 48 km |

The staircase 0.04 s / 0.44 s / 1.04 s appears at the Z1, Z2 and Z3 reaches, and
faults beyond Zone 3 are not tripped. Fault resistance is the main cause of
under-reach: it moves the measured point away from the line-angle ray, so a
slower zone operates or the relay does not operate.

## Assumptions and limitations

Single source behind the relay, no load current, homogeneous lines, ideal CT/PT,
and ideal phase selection (only the faulted loop is evaluated). Remote infeed,
shunt capacitance, CT saturation, load encroachment and power-swing blocking are
not modelled.

## Reference

EE 315 Power Systems, Module Lecture Material (Lec 33 phasor estimation, Lec 36
distance relaying), Department of Electrical Engineering.
