# EE 315 Assignment 13 - Stepped Distance Protection

Design and verification of a 3-zone mho distance relay for a 220 kV line (relay at Bus A, line AB = 80 km, adjacent lines BC = 60 km, BD = 100 km).

## Contents
| Path | Purpose |
|---|---|
| `ee315_assignment13.py` | Fault calculation, waveform synthesis, DFT phasor estimation, mho zone decision, plots |
| `results_sweep.csv`, `results_meta.json` | Numerical results (204 faults) |
| `figs/` | Generated figures |
| `docs/index.html` | Interactive web demo (GitHub Pages) |
| `report/` | LaTeX report (upload this folder to Overleaf, open `main.tex`) |

## Run
```
pip install -r requirements.txt
python ee315_assignment13.py
```

## Web demo (GitHub Pages)
Repository Settings -> Pages -> Source: `main` branch, folder `/docs`. The demo will be at `https://mayursoni08.github.io/Assignment13_PS/`.
It shows the steady-state result (no noise or DC offset).

## Report
Zip the `report/` folder and upload it to Overleaf (New Project -> Upload Project). Edit names, roll number and links in the "EDIT HERE" block at the top of `main.tex`.
