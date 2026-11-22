#!/usr/bin/env python3
"""
EE315 Power Systems - Assignment 13
Stepped distance protection (3-zone mho) for a 220 kV transmission line,
verified with simulated fault data.

Pipeline
  1. Sequence-network fault calculation (AG, BC, ABC faults, any location, any Rf)
  2. Time-domain synthesis of relay voltages/currents (decaying DC offset,
     3rd harmonic, measurement noise, random fault-inception angle)
  3. Phasor estimation: full-cycle DFT (+ optional digital mimic filter)
  4. Apparent impedance (loop impedance incl. zero-sequence compensation)
  5. Mho zone decision + stepped time delays -> trip time
  6. Plots and tables

Run:  python3 ee315_assignment13.py
"""
import os, csv, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figs")
os.makedirs(FIG, exist_ok=True)
rng = np.random.default_rng(315)

# --------------------------------------------------------------------------
# 1. SYSTEM DATA  (assumed standard example system, 50 Hz, 220 kV)
# --------------------------------------------------------------------------
f0 = 50.0
w = 2 * np.pi * f0
N = 32                       # samples per cycle (1.6 kHz sampling)
fs = N * f0
Ts = 1.0 / fs
a = np.exp(2j * np.pi / 3)

VLL = 220e3
E = VLL / np.sqrt(3)         # source phase voltage (rms), 1.0 pu
CTR = 1000.0                 # CT 1000/1 A
VTR = 220e3 / 110.0          # PT 220 kV / 110 V (line-line)  -> 2000
KSEC = CTR / VTR             # Z_secondary = Z_primary * CTR / VTR

z1 = complex(0.08, 0.40)     # ohm/km  positive = negative sequence
z0 = complex(0.30, 1.20)     # ohm/km  zero sequence
LEN = {"AB": 80.0, "BC": 60.0, "BD": 100.0}      # km

S_SC = 5000e6                # source short-circuit level (MVA) at bus A
Zs1 = (VLL ** 2 / S_SC) * np.exp(1j * np.deg2rad(85))   # ohm
Zs0 = 1.3 * Zs1 * np.exp(-1j * np.deg2rad(5))

k0 = (z0 - z1) / (3 * z1)    # residual compensation factor
theta = np.angle(z1)         # line angle = mho maximum-torque angle

# --------------------------------------------------------------------------
# 2. ZONE SETTINGS
# --------------------------------------------------------------------------
Z_AB = abs(z1) * LEN["AB"]
Z_BC = abs(z1) * LEN["BC"]   # shortest adjacent line
Z_BD = abs(z1) * LEN["BD"]   # longest adjacent line

ZR_PRI = {1: 0.80 * Z_AB,
          2: Z_AB + 0.50 * min(Z_BC, Z_BD),
          3: Z_AB + 1.20 * max(Z_BC, Z_BD)}
ZR_SEC = {z: ZR_PRI[z] * KSEC for z in ZR_PRI}
T_DELAY = {1: 0.0, 2: 0.4, 3: 1.0}       # s, intentional delays
T_MEAS = 2 * (1.0 / f0)                  # 2 cycles: data window + decision
TAU_MIMIC = 0.020                        # s, nominal L/R for mimic filter


def zone_of(Z):
    """Mho test (circles through origin at the line angle). 0 = no pickup."""
    for z in (1, 2, 3):
        r = ZR_SEC[z]
        centre = 0.5 * r * np.exp(1j * theta)
        if abs(Z - centre) <= 0.5 * r:
            return z
    return 0


def trip_time(z):
    return np.nan if z == 0 else T_MEAS + T_DELAY[z]


def expected_zone(d_km):
    """Ideal zone from true electrical distance (bolted fault)."""
    zt = abs(z1) * d_km
    for z in (1, 2, 3):
        if zt <= ZR_PRI[z]:
            return z
    return 0


# --------------------------------------------------------------------------
# 3. FAULT CALCULATION (sequence networks, single source behind the relay)
# --------------------------------------------------------------------------
def fault_phasors(ftype, d_km, Rf):
    """Return relay-point phase voltages and currents (primary, rms phasors)."""
    Z1 = Zs1 + z1 * d_km
    Z2 = Z1
    Z0 = Zs0 + z0 * d_km
    if ftype == "AG":
        I1 = I2 = I0 = E / (Z1 + Z2 + Z0 + 3 * Rf)
    elif ftype == "BC":
        I1 = E / (Z1 + Z2 + Rf)
        I2, I0 = -I1, 0
    elif ftype == "ABC":
        I1 = E / (Z1 + Rf)
        I2 = I0 = 0
    else:
        raise ValueError(ftype)
    V1 = E - I1 * Zs1
    V2 = -I2 * Zs1
    V0 = -I0 * Zs0
    to_ph = lambda x0, x1, x2: np.array([x0 + x1 + x2,
                                         x0 + a ** 2 * x1 + a * x2,
                                         x0 + a * x1 + a ** 2 * x2])
    tau = np.imag(Z1) / (w * np.real(Z1))
    return to_ph(V0, V1, V2), to_ph(I0, I1, I2), tau


# --------------------------------------------------------------------------
# 4. TIME-DOMAIN DATA GENERATION
# --------------------------------------------------------------------------
PRE, POST = 2, 6          # cycles before / after fault inception


def make_waves(ftype, d_km, Rf, alpha=None, noise=True):
    """Secondary-side sampled waveforms (3 x n) for v and i. Fault at index PRE*N."""
    Vf, If, tau = fault_phasors(ftype, d_km, Rf)
    if alpha is None:
        alpha = rng.uniform(0, 2 * np.pi)
    n = np.arange(-PRE * N, POST * N)
    t = n * Ts
    v = np.zeros((3, n.size))
    i = np.zeros((3, n.size))
    Vpre = E * np.array([1, a ** 2, a])
    for k in range(3):
        v_pre = np.sqrt(2) * abs(Vpre[k]) * np.cos(w * t + alpha + np.angle(Vpre[k]))
        v_flt = np.sqrt(2) * abs(Vf[k]) * np.cos(w * t + alpha + np.angle(Vf[k]))
        ph = alpha + np.angle(If[k])
        i_flt = np.sqrt(2) * abs(If[k]) * (np.cos(w * t + ph) - np.cos(ph) * np.exp(-t / tau)
                                          + 0.02 * np.cos(3 * w * t + 0.7))   # DC offset + 3rd harmonic
        v[k] = np.where(t < 0, v_pre, v_flt)
        i[k] = np.where(t < 0, 0.0, i_flt)
    v, i = v / VTR, i / CTR                                # to secondary (phase quantities)
    if noise:
        v += rng.normal(0, 1e-3 * np.sqrt(2) * (110 / np.sqrt(3)), v.shape)
        i += rng.normal(0, 1e-3 * np.sqrt(2) * 1.0, i.shape)
    return v, i


# --------------------------------------------------------------------------
# 5. PHASOR ESTIMATION  (Lec 33: full-cycle DFT) + digital mimic filter
# --------------------------------------------------------------------------
_W = np.exp(-2j * np.pi * np.arange(N) / N)


def mimic(x, tau_m=TAU_MIMIC):
    """y[n] = x[n] + (tau/Ts)(x[n]-x[n-1]) : removes a decaying DC term.
    Applied to BOTH v and i so the V/I ratio of the fundamental is unchanged."""
    y = x.copy()
    y[:, 1:] = x[:, 1:] + (tau_m / Ts) * (x[:, 1:] - x[:, :-1])
    return y


def dft(x, start):
    """Full-cycle DFT: rms phasor of the fundamental from x[:, start:start+N]."""
    return np.sqrt(2) / N * (x[:, start:start + N] @ _W)


def loop_impedance(ftype, V, I):
    if ftype == "AG":
        return V[0] / (I[0] + k0 * (I[0] + I[1] + I[2]))
    if ftype == "BC":
        return (V[1] - V[2]) / (I[1] - I[2])
    if ftype == "ABC":
        return (V[0] - V[1]) / (I[0] - I[1])


def estimate_Z(ftype, v, i, use_mimic=True, window=1):
    """Apparent impedance from the window starting `window` cycles after inception."""
    if use_mimic:
        v, i = mimic(v), mimic(i)
    s = PRE * N + window * N
    return loop_impedance(ftype, dft(v, s), dft(i, s))


def run_case(ftype, d_km, Rf, use_mimic=True, noise=True, alpha=None):
    v, i = make_waves(ftype, d_km, Rf, alpha, noise)
    Z = estimate_Z(ftype, v, i, use_mimic)
    z = zone_of(Z)
    return Z, z, trip_time(z)


# --------------------------------------------------------------------------
# 6. EXPERIMENTS
# --------------------------------------------------------------------------
RF_PRI = {"AG": 5.0, "BC": 2.0, "ABC": 1.0}      # ohm (primary) for main sweep
FTYPES = ["AG", "BC", "ABC"]


def location(line, m):
    """Distance from relay (km) of a fault at fraction m of `line`."""
    base = {"AB": 0.0, "BC": LEN["AB"], "BD": LEN["AB"]}[line]
    return base + m * LEN[line]


def sweep():
    rows = []
    grid = ([("AB", m) for m in np.arange(0.05, 1.0001, 0.05)] +
            [("BC", m) for m in np.arange(0.05, 1.0001, 0.05)] +
            [("BD", m) for m in np.arange(0.05, 1.4001, 0.05)])
    for ft in FTYPES:
        for line, m in grid:
            d = location(line, m)
            Z, z, tt = run_case(ft, d, RF_PRI[ft])
            rows.append(dict(ftype=ft, line=line, m=round(m, 2), d_km=d,
                             Rf=RF_PRI[ft], R=Z.real, X=Z.imag, zone=z,
                             expected=expected_zone(d), t_trip=tt))
    return rows


def accuracy_study(ntrial=300):
    out = {}
    for ft in FTYPES:
        err = {"plain": [], "mimic": []}
        for _ in range(ntrial):
            d = rng.uniform(5, 150)
            alpha = rng.uniform(0, 2 * np.pi)
            v, i = make_waves(ft, d, 0.0, alpha, True)
            zt = z1 * d * KSEC          # true line impedance, secondary ohm
            for key, um in (("plain", False), ("mimic", True)):
                Z = estimate_Z(ft, v, i, um)
                err[key].append(abs(Z - zt) / abs(zt) * 100)
        out[ft] = {k: np.array(vv) for k, vv in err.items()}
    return out


def rf_map():
    d_axis = np.arange(5, 221, 5.0)
    rf_axis = np.arange(0, 101, 5.0)
    M = np.zeros((rf_axis.size, d_axis.size), int)
    for r, Rf in enumerate(rf_axis):
        for c, d in enumerate(d_axis):
            Z, z, _ = run_case("AG", d, Rf)
            M[r, c] = z
    return d_axis, rf_axis, M


# --------------------------------------------------------------------------
# 7. PLOTS
# --------------------------------------------------------------------------
COL = {1: "#c0392b", 2: "#e67e22", 3: "#2471a3"}
FCOL = {"AG": "#1b7f3b", "BC": "#8e44ad", "ABC": "#d35400"}


def draw_mho(ax, lw=1.6):
    for z in (3, 2, 1):
        r = ZR_SEC[z]
        c = 0.5 * r * np.exp(1j * theta)
        ax.add_patch(Circle((c.real, c.imag), r / 2, fill=False, ec=COL[z], lw=lw))
        tip = r * np.exp(1j * theta)
        ax.annotate(f"Z{z}", (tip.real, tip.imag), textcoords="offset points",
                    xytext=(6, -3), color=COL[z], fontsize=10, fontweight="bold")
    ax.axhline(0, color="k", lw=.6)
    ax.axvline(0, color="k", lw=.6)
    ax.set_aspect("equal")
    ax.set_xlabel("R  (ohm, secondary)")
    ax.set_ylabel("X  (ohm, secondary)")
    ax.grid(alpha=.3)


def fig_system():
    fig, ax = plt.subplots(figsize=(9, 3.6))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 5)
    ax.axis("off")
    buses = {"A": (2.0, 2.5), "B": (6.0, 2.5), "C": (10.0, 3.8), "D": (10.0, 1.2)}
    ax.add_patch(Circle((0.7, 2.5), .45, fc="white", ec="k"))
    ax.text(0.7, 2.5, "~", ha="center", va="center", fontsize=18)
    ax.plot([1.15, 2.0], [2.5, 2.5], "k")
    ax.text(1.55, 2.75, "Zs", ha="center", fontsize=9)
    for nm, (x, y) in buses.items():
        ax.plot([x, x], [y - .5, y + .5], "k", lw=5)
        ax.text(x, y + .75, f"Bus {nm}", ha="center", fontweight="bold")
    ax.plot([2.0, 6.0], [2.5, 2.5], "k")
    ax.plot([6.0, 10.0], [2.5, 3.8], "k")
    ax.plot([6.0, 10.0], [2.5, 1.2], "k")
    ax.add_patch(FancyBboxPatch((2.15, 2.25), .5, .5, boxstyle="round", fc="#f6d5d5", ec="#c0392b"))
    ax.text(2.4, 2.5, "R", ha="center", va="center", color="#c0392b", fontweight="bold")
    ax.text(4.0, 2.05, "Line AB, 80 km", ha="center", fontsize=9)
    ax.text(8.0, 3.55, "Line BC, 60 km", ha="center", fontsize=9, rotation=18)
    ax.text(8.0, 1.25, "Line BD, 100 km", ha="center", fontsize=9, rotation=-18)
    ax.text(4.2, 4.25, "Relay R protects line AB\nCT 1000/1 A, PT 220 kV/110 V\n220 kV, 50 Hz",
            fontsize=9, ha="center", va="center",
            bbox=dict(boxstyle="round", fc="#f4f6f7", ec="#7f8c8d"))
    fig.savefig(os.path.join(FIG, "fig1_system.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_rx(rows):
    fig, ax = plt.subplots(figsize=(8, 7.6))
    draw_mho(ax)
    mk = {"AG": "o", "BC": "s", "ABC": "^"}
    for ft in FTYPES:
        pts = [(r["R"], r["X"]) for r in rows if r["ftype"] == ft]
        ax.scatter(*zip(*pts), s=26, marker=mk[ft], color=FCOL[ft], alpha=.85,
                   label=f"{ft} fault apparent impedance", zorder=3)
    ax.plot([0, 12 * np.cos(theta)], [0, 12 * np.sin(theta)], "k--", lw=.7)
    ax.set_xlim(-18, 28)
    ax.set_ylim(-6, 48)
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title("Mho zones and estimated fault impedances (R-X plane)")
    fig.savefig(os.path.join(FIG, "fig2_rx_diagram.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_staircase(rows):
    fig, ax = plt.subplots(figsize=(9, 4.6))
    off = {"AG": -0.8, "BC": 0.0, "ABC": 0.8}
    for ft in FTYPES:
        # path A->B->D only (the extended line beyond D included)
        pts = sorted([(r["d_km"], r["t_trip"]) for r in rows
                      if r["ftype"] == ft and r["line"] in ("AB", "BD")])
        d = np.array([p[0] for p in pts]) + off[ft]
        t = np.array([p[1] for p in pts])
        ax.plot(d, np.where(np.isnan(t), np.nan, t), "-o", ms=4, color=FCOL[ft], label=ft)
        for dd, tt in zip(d, t):
            if np.isnan(tt):
                ax.plot(dd, 1.2, "x", color=FCOL[ft], ms=5)
    reach_km = {z: ZR_PRI[z] / abs(z1) for z in (1, 2, 3)}
    for z in (1, 2, 3):
        ax.axvline(reach_km[z], color=COL[z], ls="--", lw=1)
        ax.text(reach_km[z] - 2, 1.27, f"Z{z} reach\n{reach_km[z]:.0f} km", color=COL[z],
                ha="right", fontsize=8)
    ax.axvline(80, color="gray", ls=":", lw=1)
    ax.axvline(180, color="gray", ls=":", lw=1)
    ax.text(80, -0.09, "Bus B", ha="center", fontsize=8)
    ax.text(180, -0.09, "Bus D", ha="center", fontsize=8)
    ax.set_ylim(-0.05, 1.4)
    ax.set_xlabel("Fault distance from relay R along A-B-D  (km)")
    ax.set_ylabel("Relay trip time (s)")
    ax.set_title("Stepped time-distance characteristic ( x = no trip, fault beyond Z3 )")
    ax.grid(alpha=.3)
    ax.legend(title="Fault type", loc="upper left")
    fig.savefig(os.path.join(FIG, "fig3_staircase.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_waveforms():
    d = 0.60 * LEN["AB"]
    # inception angle giving the maximum DC offset in phase A (worst case)
    alpha = -np.angle(fault_phasors("AG", d, 5.0)[1][0])
    v, i = make_waves("AG", d, 5.0, alpha, True)
    n = np.arange(-PRE * N, POST * N)
    t = n * Ts * 1e3
    fig, axs = plt.subplots(2, 1, figsize=(9, 5.2), sharex=True)
    for k, (nm, c) in enumerate(zip("ABC", ["#c0392b", "#27ae60", "#2471a3"])):
        axs[0].plot(t, v[k], color=c, lw=1, label=f"v{nm}")
        axs[1].plot(t, i[k], color=c, lw=1, label=f"i{nm}")
    for ax in axs:
        ax.axvline(0, color="k", ls="--", lw=.8)
        ax.grid(alpha=.3)
        ax.legend(loc="upper right", ncol=3, fontsize=8)
    axs[0].set_ylabel("Voltage (V, sec.)")
    axs[1].set_ylabel("Current (A, sec.)")
    axs[1].set_xlabel("Time (ms)   [fault inception at 0 ms]")
    axs[0].set_title("Simulated relay waveforms: A-G fault at 60 % of line AB, Rf = 5 ohm")
    fig.savefig(os.path.join(FIG, "fig4_waveforms.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)

    # impedance trajectory (sliding window)
    fig, axs = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw={"width_ratios": [1, 1.1]})
    ax = axs[0]
    draw_mho(ax)
    for um, col, lab in ((False, "#7f8c8d", "DFT only"), (True, "#1b7f3b", "DFT + mimic filter")):
        vv, ii = (mimic(v), mimic(i)) if um else (v, i)
        Zs_ = []
        for s in range(PRE * N, PRE * N + 4 * N):
            Zs_.append(loop_impedance("AG", dft(vv, s), dft(ii, s)))
        Zs_ = np.array(Zs_)
        ax.plot(Zs_.real, Zs_.imag, "-", color=col, lw=1, alpha=.8)
        ax.scatter(Zs_.real, Zs_.imag, s=8, color=col, label=lab, zorder=3)
    zt = z1 * d * KSEC
    ax.plot(zt.real, zt.imag, "k*", ms=11, label="true line impedance", zorder=4)
    ax.set_xlim(-6, 14)
    ax.set_ylim(-4, 28)
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title("Impedance trajectory after fault")
    ax = axs[1]
    ts = np.arange(PRE * N, PRE * N + 4 * N)
    for um, col, lab in ((False, "#7f8c8d", "DFT only"), (True, "#1b7f3b", "DFT + mimic filter")):
        vv, ii = (mimic(v), mimic(i)) if um else (v, i)
        mag = [abs(loop_impedance("AG", dft(vv, s), dft(ii, s))) for s in ts]
        ax.plot((ts - PRE * N) * Ts * 1e3, mag, color=col, label=lab)
    ax.axhline(abs(zt), color="k", ls=":", label="true |Z|")
    ax.axhline(ZR_SEC[1], color=COL[1], ls="--", lw=1)
    ax.text(0.5, ZR_SEC[1] + .3, "Z1 reach", color=COL[1], fontsize=8)
    ax.set_xlabel("Window start after inception (ms)")
    ax.set_ylabel("|Z apparent| (ohm, sec.)")
    ax.set_title("Convergence of the impedance estimate")
    ax.set_ylim(0, 30)
    ax.grid(alpha=.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig5_trajectory.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_accuracy(acc):
    fig, axs = plt.subplots(1, 3, figsize=(10, 3.8), sharey=True)
    for ax, ft in zip(axs, FTYPES):
        ax.boxplot([acc[ft]["plain"], acc[ft]["mimic"]], tick_labels=["DFT", "DFT+mimic"],
                   showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor="#d6eaf8"))
        ax.set_title(f"{ft} faults")
        ax.grid(alpha=.3)
    axs[0].set_ylabel("Impedance error (%)")
    fig.suptitle("Phasor-estimation accuracy (bolted faults, random location and inception angle)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig6_accuracy.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_rfmap(d_axis, rf_axis, M):
    fig, ax = plt.subplots(figsize=(9, 4.6))
    cmap = ListedColormap(["#ecf0f1", COL[1], COL[2], COL[3]])
    ax.imshow(M, origin="lower", aspect="auto", cmap=cmap, vmin=0, vmax=3,
              extent=[d_axis[0] - 2.5, d_axis[-1] + 2.5, rf_axis[0] - 2.5, rf_axis[-1] + 2.5])
    ax.axvline(80, color="k", ls=":", lw=1)
    ax.axvline(180, color="k", ls=":", lw=1)
    ax.text(81, 96, "Bus B", fontsize=8)
    ax.text(181, 96, "Bus D", fontsize=8)
    ax.set_xlabel("Fault distance from relay along A-B-D (km); > 180 km = hypothetical extension")
    ax.set_ylabel("Fault resistance Rf (ohm, primary)")
    ax.set_title("A-G fault: zone that picks up vs. distance and fault resistance")
    ax.legend(handles=[Patch(fc=COL[1], label="Zone 1 (instantaneous)"),
                       Patch(fc=COL[2], label="Zone 2 (0.4 s)"),
                       Patch(fc=COL[3], label="Zone 3 (1.0 s)"),
                       Patch(fc="#ecf0f1", ec="gray", label="No pickup")],
              loc="upper right", fontsize=8, framealpha=.95)
    fig.savefig(os.path.join(FIG, "fig7_rf_map.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------
# 8. MAIN
# --------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("SYSTEM / SETTINGS")
    print(f"  |Z1 line| = {abs(z1):.4f} ohm/km, angle = {np.degrees(theta):.1f} deg")
    print(f"  Z_AB={Z_AB:.2f}  Z_BC={Z_BC:.2f}  Z_BD={Z_BD:.2f} ohm (primary)")
    print(f"  CT {CTR:.0f}/1, PT {VTR:.0f}:1  -> Z_sec = {KSEC:.3f} x Z_pri")
    for z in (1, 2, 3):
        print(f"  Zone {z}: {ZR_PRI[z]:.2f} ohm pri = {ZR_SEC[z]:.2f} ohm sec "
              f"(= {ZR_PRI[z] / abs(z1):.1f} km),  t = {T_DELAY[z]:.1f} s")

    rows = sweep()
    with open(os.path.join(HERE, "results_sweep.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)

    agree = np.mean([r["zone"] == r["expected"] for r in rows]) * 100
    print(f"\nSweep: {len(rows)} faults, zone agrees with ideal in {agree:.1f} %")
    for ft in FTYPES:
        sub = [r for r in rows if r["ftype"] == ft]
        bad = [r for r in sub if r["zone"] != r["expected"]]
        print(f"  {ft}: {len(sub) - len(bad)}/{len(sub)} match;"
              f" mismatches at d(km)={[round(b['d_km'], 1) for b in bad]}"
              f" got {[b['zone'] for b in bad]} exp {[b['expected'] for b in bad]}")

    acc = accuracy_study()
    stats = {}
    for ft in FTYPES:
        stats[ft] = {k: dict(mean=float(np.mean(v)), p95=float(np.percentile(v, 95)),
                             max=float(np.max(v))) for k, v in acc[ft].items()}
        print(f"  Accuracy {ft}: plain mean {stats[ft]['plain']['mean']:.2f} % / max "
              f"{stats[ft]['plain']['max']:.1f} %; mimic mean {stats[ft]['mimic']['mean']:.2f} % / max "
              f"{stats[ft]['mimic']['max']:.1f} %")

    d_axis, rf_axis, M = rf_map()

    fig_system()
    fig_rx(rows)
    fig_staircase(rows)
    fig_waveforms()
    fig_accuracy(acc)
    fig_rfmap(d_axis, rf_axis, M)

    # max Rf for Zone 1 pickup at a few distances (AG)
    zmax = {}
    for d in (16, 32, 48, 64):
        c = int(np.argmin(abs(d_axis - d)))
        ok = [rf_axis[r] for r in range(rf_axis.size) if M[r, c] == 1]
        zmax[d] = max(ok) if ok else None
    print("  Max Rf (ohm primary) for Z1 pickup, AG:", zmax)

    meta = dict(ZR_PRI=ZR_PRI, ZR_SEC=ZR_SEC, T_DELAY=T_DELAY, KSEC=KSEC,
                Z_AB=Z_AB, Z_BC=Z_BC, Z_BD=Z_BD, agree=agree, stats=stats,
                theta_deg=float(np.degrees(theta)), zmax_rf={str(k): v for k, v in zmax.items()},
                Zs1=[abs(Zs1), float(np.degrees(np.angle(Zs1)))],
                Zs0=[abs(Zs0), float(np.degrees(np.angle(Zs0)))])
    with open(os.path.join(HERE, "results_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=1, default=float)
    print("\nDone. Figures in", FIG)


if __name__ == "__main__":
    main()
