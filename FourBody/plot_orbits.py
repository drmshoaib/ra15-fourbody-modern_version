"""
Orbit visualisation for the 4-body RA15 integrator output.

Reads result1.txt and produces four publication-quality figures:
  1. orbit_paths.png       – full trajectories of all 4 bodies
  2. orbit_animated.gif    – animated orbit (optional, can be slow)
  3. energy_drift.png      – fractional energy drift vs time
  4. speed_profile.png     – speed of each body vs time
  5. orbit_phases.png      – side-by-side early / mid / late snapshots

Run:
  python plot_orbits.py [result_file]  (defaults to result1.txt)
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")           # no display needed — saves directly to PNG
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle

# ── configuration ────────────────────────────────────────────────────────────
RESULT_FILE = sys.argv[1] if len(sys.argv) > 1 else "result1.txt"
OUT_DIR     = os.path.dirname(os.path.abspath(RESULT_FILE))

BODY_COLORS  = ["#E63946", "#2196F3", "#4CAF50", "#FF9800"]   # red, blue, green, orange
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

# ── load data ─────────────────────────────────────────────────────────────────
print(f"Reading {RESULT_FILE} ...")
data = np.loadtxt(RESULT_FILE, comments="#")

t      = data[:, 0]
x      = [data[:, 1 + 2*i] for i in range(4)]   # x-coords for each body
y      = [data[:, 2 + 2*i] for i in range(4)]   # y-coords
vx     = [data[:, 9 + 2*i] for i in range(4)]
vy     = [data[:,10 + 2*i] for i in range(4)]
energy = data[:, 17]

N = len(t)
print(f"  {N} snapshots,  t = {t[0]:.2f} to {t[-1]:.2f}")

E0       = energy[0]
rel_drift = (energy - E0) / abs(E0)
speed    = [np.sqrt(vx[i]**2 + vy[i]**2) for i in range(4)]

# ─────────────────────────────────────────────────────────────────────────────
#  Figure 1 — Orbit paths  (fading colour trail + current-position dot)
# ─────────────────────────────────────────────────────────────────────────────
def make_fading_trail(ax, xs, ys, color, alpha_min=0.05, alpha_max=0.9, lw=1.2):
    """Draw a trajectory whose opacity fades from old (dim) to new (bright)."""
    points  = np.array([xs, ys]).T.reshape(-1, 1, 2)
    segs    = np.concatenate([points[:-1], points[1:]], axis=1)
    alphas  = np.linspace(alpha_min, alpha_max, len(segs))
    rgba    = np.array([matplotlib.colors.to_rgba(color, a) for a in alphas])
    lc      = LineCollection(segs, colors=rgba, linewidths=lw)
    ax.add_collection(lc)

with plt.rc_context(STYLE):
    fig, ax = plt.subplots(figsize=(9, 9))

    for i in range(4):
        make_fading_trail(ax, x[i], y[i], BODY_COLORS[i])
        # Start position (hollow marker)
        ax.plot(x[i][0], y[i][0], BODY_MARKERS[i], color=BODY_COLORS[i],
                ms=10, mfc="none", mew=1.8, zorder=5)
        # End position (filled marker)
        ax.plot(x[i][-1], y[i][-1], BODY_MARKERS[i], color=BODY_COLORS[i],
                ms=10, zorder=6, label=f"{BODY_LABELS[i]}")

    ax.set_aspect("equal")
    ax.autoscale()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"4-Body Orbits  (RA15 integrator,  T = {t[-1]:.0f})",
                 fontsize=13, pad=12)
    ax.legend(loc="upper right")
    ax.grid(True)

    # Annotation
    ax.text(0.02, 0.02,
            "Hollow = start   •   Filled = end",
            transform=ax.transAxes, fontsize=9, color="#8b949e")

    out = os.path.join(OUT_DIR, "orbit_paths.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")

# ─────────────────────────────────────────────────────────────────────────────
#  Figure 2 — Energy drift
# ─────────────────────────────────────────────────────────────────────────────
with plt.rc_context(STYLE):
    fig, ax = plt.subplots(figsize=(10, 4))

    ax.plot(t, rel_drift, color="#58a6ff", lw=1.2)
    ax.axhline(0, color="#30363d", lw=0.8, ls="--")
    ax.fill_between(t, rel_drift, 0, alpha=0.15, color="#58a6ff")

    ax.set_xlabel("Time")
    ax.set_ylabel("(E - E0) / |E0|")
    ax.set_title("Fractional Energy Drift", fontsize=13)
    ax.grid(True)

    # Annotate max drift
    idx_max = np.argmax(np.abs(rel_drift))
    ax.annotate(f"max |dE/E0| = {np.abs(rel_drift[idx_max]):.2e}",
                xy=(t[idx_max], rel_drift[idx_max]),
                xytext=(t[idx_max] + (t[-1]-t[0])*0.05, rel_drift[idx_max] * 1.3),
                color="#f0883e",
                arrowprops=dict(arrowstyle="->", color="#f0883e", lw=1),
                fontsize=9)

    out = os.path.join(OUT_DIR, "energy_drift.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")

# ─────────────────────────────────────────────────────────────────────────────
#  Figure 3 — Speed profiles
# ─────────────────────────────────────────────────────────────────────────────
with plt.rc_context(STYLE):
    fig, ax = plt.subplots(figsize=(10, 4))

    for i in range(4):
        ax.plot(t, speed[i], color=BODY_COLORS[i], lw=1.1, label=BODY_LABELS[i])

    ax.set_xlabel("Time")
    ax.set_ylabel("|v|")
    ax.set_title("Speed of Each Body vs Time", fontsize=13)
    ax.legend()
    ax.grid(True)

    out = os.path.join(OUT_DIR, "speed_profile.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")

# ─────────────────────────────────────────────────────────────────────────────
#  Figure 4 — Phase snapshots (early / middle / late)
# ─────────────────────────────────────────────────────────────────────────────
slices = [
    (0,           N // 4,  "Early"),
    (N // 4,      3*N//4,  "Middle"),
    (3*N // 4,    N,       "Late"),
]

with plt.rc_context(STYLE):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("4-Body Orbit — Phase Snapshots", fontsize=14, y=1.01)

    for ax, (s, e, label) in zip(axes, slices):
        for i in range(4):
            ax.plot(x[i][s:e], y[i][s:e], color=BODY_COLORS[i],
                    lw=1.1, alpha=0.85, label=BODY_LABELS[i])
            ax.plot(x[i][s], y[i][s], BODY_MARKERS[i],
                    color=BODY_COLORS[i], ms=7, mfc="none", mew=1.5)
            ax.plot(x[i][e-1], y[i][e-1], BODY_MARKERS[i],
                    color=BODY_COLORS[i], ms=7)
        ax.set_aspect("equal")
        ax.autoscale()
        ax.set_title(f"{label}  (t = {t[s]:.1f} – {t[e-1]:.1f})", fontsize=11)
        ax.grid(True)
        if ax is axes[0]:
            ax.legend(fontsize=8, loc="upper right")

    fig.tight_layout()
    out = os.path.join(OUT_DIR, "orbit_phases.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")

# ─────────────────────────────────────────────────────────────────────────────
#  Figure 5 — Separation distances between all pairs over time
# ─────────────────────────────────────────────────────────────────────────────
pairs  = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
p_cols = ["#E63946","#2196F3","#4CAF50","#FF9800","#9C27B0","#00BCD4"]
p_labs = ["r01", "r02", "r03", "r12", "r13", "r23"]

with plt.rc_context(STYLE):
    fig, ax = plt.subplots(figsize=(10, 4))

    for (i, j), col, lab in zip(pairs, p_cols, p_labs):
        dist = np.sqrt((x[i]-x[j])**2 + (y[i]-y[j])**2)
        ax.plot(t, dist, color=col, lw=1.0, label=lab)

    ax.set_xlabel("Time")
    ax.set_ylabel("Separation")
    ax.set_title("Pairwise Separations vs Time", fontsize=13)
    ax.legend(ncol=3, fontsize=9)
    ax.grid(True)

    out = os.path.join(OUT_DIR, "separations.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")

print("\nAll figures saved to:", OUT_DIR)
