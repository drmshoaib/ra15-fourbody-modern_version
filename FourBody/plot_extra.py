"""
Additional orbit figures for the 4-body RA15 integrator output.

Produces:
  angular_momentum.png  – total Lz vs time (conservation check)
  phase_portrait.png    – x vs vx phase-space for each body
  chaos_transition.png  – zoomed orbit paths across the stability break
  centre_of_mass.png    – CoM position vs time (should be ~constant)

Run:
  python plot_extra.py [result_file]  (defaults to result1.txt)
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

RESULT_FILE = sys.argv[1] if len(sys.argv) > 1 else "result1.txt"
OUT_DIR     = os.path.dirname(os.path.abspath(RESULT_FILE))

BODY_COLORS  = ["#E63946", "#2196F3", "#4CAF50", "#FF9800"]
BODY_LABELS  = ["Body 0", "Body 1", "Body 2", "Body 3"]
BODY_MARKERS = ["o", "s", "^", "D"]

STYLE = {
    "figure.facecolor": "#0d1117",
    "axes.facecolor":   "#0d1117",
    "axes.edgecolor":   "#30363d",
    "axes.labelcolor":  "#c9d1d9",
    "axes.titlecolor":  "#c9d1d9",
    "xtick.color":      "#8b949e",
    "ytick.color":      "#8b949e",
    "grid.color":       "#21262d",
    "grid.linewidth":   0.5,
    "text.color":       "#c9d1d9",
    "legend.facecolor": "#161b22",
    "legend.edgecolor": "#30363d",
    "font.size":        11,
}

# ── load ─────────────────────────────────────────────────────────────────────
print(f"Reading {RESULT_FILE} ...")
data = np.loadtxt(RESULT_FILE, comments="#")

t  = data[:, 0]
x  = [data[:, 1 + 2*i] for i in range(4)]
y  = [data[:, 2 + 2*i] for i in range(4)]
vx = [data[:, 9 + 2*i] for i in range(4)]
vy = [data[:,10 + 2*i] for i in range(4)]
M  = np.ones(4)   # masses (equal, from initialcond.txt)
N  = len(t)

# ── helpers ───────────────────────────────────────────────────────────────────
def fading_trail(ax, xs, ys, color, lw=1.2):
    pts  = np.array([xs, ys]).T.reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    alphas = np.linspace(0.05, 0.9, len(segs))
    rgba   = np.array([matplotlib.colors.to_rgba(color, a) for a in alphas])
    ax.add_collection(LineCollection(segs, colors=rgba, linewidths=lw))

# ─────────────────────────────────────────────────────────────────────────────
#  Figure 1 – Angular momentum  Lz = sum_i m_i (x_i * vy_i - y_i * vx_i)
# ─────────────────────────────────────────────────────────────────────────────
Lz = sum(M[i] * (x[i] * vy[i] - y[i] * vx[i]) for i in range(4))
L0 = Lz[0]
rel_Ldrift = (Lz - L0) / abs(L0)

with plt.rc_context(STYLE):
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True,
                             gridspec_kw={"height_ratios": [2, 1], "hspace": 0.08})

    ax = axes[0]
    ax.plot(t, Lz, color="#58a6ff", lw=1.3)
    ax.axhline(L0, color="#30363d", lw=0.8, ls="--")
    ax.set_ylabel("Total Lz")
    ax.set_title("Angular Momentum Conservation", fontsize=13)
    ax.grid(True)

    ax2 = axes[1]
    ax2.plot(t, rel_Ldrift, color="#f0883e", lw=1.0)
    ax2.axhline(0, color="#30363d", lw=0.8, ls="--")
    ax2.fill_between(t, rel_Ldrift, 0, alpha=0.2, color="#f0883e")
    ax2.set_ylabel("(Lz - Lz0) / |Lz0|")
    ax2.set_xlabel("Time")
    ax2.grid(True)

    out = os.path.join(OUT_DIR, "angular_momentum.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")

# ─────────────────────────────────────────────────────────────────────────────
#  Figure 2 – Phase portrait  (x vs vx, y vs vy) — coloured by time
# ─────────────────────────────────────────────────────────────────────────────
cmap = plt.get_cmap("plasma")
norm = matplotlib.colors.Normalize(vmin=t[0], vmax=t[-1])

with plt.rc_context(STYLE):
    fig, axes = plt.subplots(2, 4, figsize=(16, 7))
    fig.suptitle("Phase Portraits  (colour = time, dark=early, bright=late)",
                 fontsize=13, y=1.01)

    for i in range(4):
        # x vs vx
        ax = axes[0, i]
        sc = ax.scatter(x[i], vx[i], c=t, cmap="plasma", s=4, lw=0, norm=norm)
        ax.set_xlabel("x")
        ax.set_ylabel("vx")
        ax.set_title(f"{BODY_LABELS[i]}", color=BODY_COLORS[i])
        ax.grid(True)

        # y vs vy
        ax = axes[1, i]
        ax.scatter(y[i], vy[i], c=t, cmap="plasma", s=4, lw=0, norm=norm)
        ax.set_xlabel("y")
        ax.set_ylabel("vy")
        ax.grid(True)

    fig.colorbar(sc, ax=axes, orientation="vertical",
                 fraction=0.015, pad=0.02, label="Time")
    fig.tight_layout()

    out = os.path.join(OUT_DIR, "phase_portrait.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")

# ─────────────────────────────────────────────────────────────────────────────
#  Figure 3 – Chaos transition zoom  (t = 35 to 65)
# ─────────────────────────────────────────────────────────────────────────────
mask = (t >= 35) & (t <= 65)

with plt.rc_context(STYLE):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("Stability Breakdown  (t = 35 to 65)", fontsize=13)

    # Left: just the transition window
    ax = axes[0]
    for i in range(4):
        fading_trail(ax, x[i][mask], y[i][mask], BODY_COLORS[i])
        ax.plot(x[i][mask][0],  y[i][mask][0],  BODY_MARKERS[i],
                color=BODY_COLORS[i], ms=9, mfc="none", mew=1.8)
        ax.plot(x[i][mask][-1], y[i][mask][-1], BODY_MARKERS[i],
                color=BODY_COLORS[i], ms=9, label=BODY_LABELS[i])
    ax.set_aspect("equal"); ax.autoscale()
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_title("Orbit paths during transition")
    ax.legend(fontsize=9); ax.grid(True)
    ax.text(0.02, 0.97, "Hollow = t=35   Filled = t=65",
            transform=ax.transAxes, fontsize=8, color="#8b949e", va="top")

    # Right: separations in the same window
    pairs  = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    p_cols = ["#E63946","#2196F3","#4CAF50","#FF9800","#9C27B0","#00BCD4"]
    p_labs = ["r01","r02","r03","r12","r13","r23"]
    ax2 = axes[1]
    for (pi, pj), col, lab in zip(pairs, p_cols, p_labs):
        dist = np.sqrt((x[pi]-x[pj])**2 + (y[pi]-y[pj])**2)
        ax2.plot(t[mask], dist[mask], color=col, lw=1.2, label=lab)
    ax2.set_xlabel("Time"); ax2.set_ylabel("Separation")
    ax2.set_title("Pairwise separations during transition")
    ax2.legend(ncol=2, fontsize=9); ax2.grid(True)
    # Mark the onset
    ax2.axvline(44, color="#8b949e", lw=0.8, ls=":")
    ax2.text(44.5, ax2.get_ylim()[1]*0.95, "onset ~t=44",
             color="#8b949e", fontsize=8, va="top")

    fig.tight_layout()
    out = os.path.join(OUT_DIR, "chaos_transition.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")

# ─────────────────────────────────────────────────────────────────────────────
#  Figure 4 – Centre of mass  (should stay near zero — barycentric coords)
# ─────────────────────────────────────────────────────────────────────────────
Mtot = M.sum()
cx = sum(M[i] * x[i] for i in range(4)) / Mtot
cy = sum(M[i] * y[i] for i in range(4)) / Mtot
cvx = sum(M[i] * vx[i] for i in range(4)) / Mtot
cvy = sum(M[i] * vy[i] for i in range(4)) / Mtot

with plt.rc_context(STYLE):
    fig, axes = plt.subplots(2, 2, figsize=(12, 6))
    fig.suptitle("Centre of Mass (should be ~0 in barycentric coordinates)", fontsize=13)

    for ax, data_y, label, color in [
        (axes[0,0], cx,  "CoM x",  "#58a6ff"),
        (axes[0,1], cy,  "CoM y",  "#4CAF50"),
        (axes[1,0], cvx, "CoM vx", "#E63946"),
        (axes[1,1], cvy, "CoM vy", "#FF9800"),
    ]:
        ax.plot(t, data_y, color=color, lw=1.0)
        ax.axhline(0, color="#30363d", lw=0.8, ls="--")
        ax.set_ylabel(label)
        ax.set_xlabel("Time")
        ax.grid(True)
        # Note the scale — if near zero it confirms barycentric coords
        rng = np.ptp(data_y)
        ax.text(0.98, 0.95, f"range = {rng:.2e}",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=8, color="#8b949e")

    fig.tight_layout()
    out = os.path.join(OUT_DIR, "centre_of_mass.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")

print("\nAll extra figures saved to:", OUT_DIR)
