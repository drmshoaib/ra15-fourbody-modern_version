"""
Generate orbit_phases and phase_portrait figures for three different initial
condition sets, providing a comparative survey of 4-body dynamics.

Sets:
  1. square_large  – square R=2, long-period orbit  (T=200)
  2. square_tight  – square R=0.5, fast orbit        (T=30)
  3. hierarchical  – two close binaries + mutual orbit (T=300)

Run:
  python plot_examples.py
(exe path and output directory are hard-coded relative to this script's location)
"""

import os
import sys
import subprocess
import tempfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection

# ── paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.dirname(SCRIPT_DIR)
EXE        = os.path.join(REPO_ROOT, "x64", "Release", "FourBody.exe")
OUT_DIR    = SCRIPT_DIR

if not os.path.isfile(EXE):
    sys.exit(f"Executable not found: {EXE}\nBuild the Release|x64 configuration first.")

# ── style ────────────────────────────────────────────────────────────────────
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
    "font.size":        10,
}

BODY_COLORS  = ["#E63946", "#2196F3", "#4CAF50", "#FF9800"]
BODY_LABELS  = ["Body 0", "Body 1", "Body 2", "Body 3"]
BODY_MARKERS = ["o", "s", "^", "D"]

# ── initial condition sets ────────────────────────────────────────────────────
# Each entry: (tag, title_label, ic_file)
# IC files live in data/ alongside initialcond.txt.

DATA_DIR = os.path.join(SCRIPT_DIR, "data")

EXAMPLES = [
    ("square_large",  "Larger square  R = 2  (T = 200)",
     os.path.join(DATA_DIR, "square_large.txt")),
    ("square_tight",  "Tight square  R = 0.5  (T = 30)",
     os.path.join(DATA_DIR, "square_tight.txt")),
    ("hierarchical",  "Hierarchical binary-binary  (T = 300)",
     os.path.join(DATA_DIR, "hierarchical.txt")),
]

# ── helpers ───────────────────────────────────────────────────────────────────
def fading_trail(ax, xs, ys, color, lw=1.2):
    pts  = np.array([xs, ys]).T.reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    alphas = np.linspace(0.05, 0.90, len(segs))
    rgba   = np.array([mcolors.to_rgba(color, a) for a in alphas])
    ax.add_collection(LineCollection(segs, colors=rgba, linewidths=lw))


def run_integrator(ic_path, tag):
    """Run the exe with an IC file, return path to result file."""
    if not os.path.isfile(ic_path):
        print(f"  IC file not found: {ic_path}")
        return None
    tmp_dir  = tempfile.mkdtemp()
    out_path = os.path.join(tmp_dir, f"result_{tag}.txt")

    result = subprocess.run(
        [EXE, ic_path, out_path],
        capture_output=True, text=True
    )
    print(f"  [{tag}] {result.stdout.strip()}")
    if result.returncode != 0:
        print(f"  STDERR: {result.stderr.strip()}")
        return None
    return out_path


def load_data(path):
    data = np.loadtxt(path, comments="#")
    t  = data[:, 0]
    x  = [data[:, 1 + 2*i] for i in range(4)]
    y  = [data[:, 2 + 2*i] for i in range(4)]
    vx = [data[:, 9 + 2*i] for i in range(4)]
    vy = [data[:,10 + 2*i] for i in range(4)]
    return t, x, y, vx, vy


# ── figure generators ─────────────────────────────────────────────────────────
def make_orbit_phases(t, x, y, title, tag):
    N = len(t)
    slices = [
        (0,       N // 4,  "Early"),
        (N // 4,  3*N//4,  "Middle"),
        (3*N//4,  N,       "Late"),
    ]
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        fig.suptitle(f"4-Body Orbit Phases — {title}", fontsize=13, y=1.01)

        for ax, (s, e, label) in zip(axes, slices):
            for i in range(4):
                ax.plot(x[i][s:e], y[i][s:e],
                        color=BODY_COLORS[i], lw=1.1, alpha=0.85,
                        label=BODY_LABELS[i])
                ax.plot(x[i][s], y[i][s], BODY_MARKERS[i],
                        color=BODY_COLORS[i], ms=7, mfc="none", mew=1.5)
                ax.plot(x[i][e-1], y[i][e-1], BODY_MARKERS[i],
                        color=BODY_COLORS[i], ms=7)
            ax.set_aspect("equal")
            ax.autoscale()
            ax.set_title(f"{label}  (t = {t[s]:.1f} to {t[e-1]:.1f})", fontsize=10)
            ax.grid(True)
            if ax is axes[0]:
                ax.legend(fontsize=8, loc="upper right")

        fig.tight_layout()
        out = os.path.join(OUT_DIR, f"orbit_phases_{tag}.png")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved {out}")
    return out


def make_phase_portrait(t, x, y, vx, vy, title, tag):
    norm = mcolors.Normalize(vmin=t[0], vmax=t[-1])
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 4, figsize=(16, 7))
        fig.suptitle(
            f"Phase Portraits — {title}\n(colour = time: dark=early, bright=late)",
            fontsize=12, y=1.02
        )
        sc_ref = None
        for i in range(4):
            # x vs vx
            ax = axes[0, i]
            sc = ax.scatter(x[i], vx[i], c=t, cmap="plasma", s=3, lw=0, norm=norm)
            if sc_ref is None:
                sc_ref = sc
            ax.set_xlabel("x");  ax.set_ylabel("vx")
            ax.set_title(BODY_LABELS[i], color=BODY_COLORS[i])
            ax.grid(True)

            # y vs vy
            ax = axes[1, i]
            ax.scatter(y[i], vy[i], c=t, cmap="plasma", s=3, lw=0, norm=norm)
            ax.set_xlabel("y");  ax.set_ylabel("vy")
            ax.grid(True)

        fig.colorbar(sc_ref, ax=axes, orientation="vertical",
                     fraction=0.012, pad=0.02, label="Time")
        fig.tight_layout()
        out = os.path.join(OUT_DIR, f"phase_portrait_{tag}.png")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved {out}")
    return out


# ── main loop ─────────────────────────────────────────────────────────────────
generated = {}   # tag -> {orbit_phases, phase_portrait}

for tag, title, ic_path in EXAMPLES:
    print(f"\nRunning integrator: {title}")
    result_path = run_integrator(ic_path, tag)
    if result_path is None:
        print(f"  SKIPPED (integrator failed)")
        continue

    t, x, y, vx, vy = load_data(result_path)
    print(f"  Loaded {len(t)} snapshots  t = {t[0]:.1f} to {t[-1]:.1f}")

    op  = make_orbit_phases(t, x, y, title, tag)
    pp  = make_phase_portrait(t, x, y, vx, vy, title, tag)
    generated[tag] = {"orbit_phases": op, "phase_portrait": pp}

print("\nAll example figures saved to:", OUT_DIR)
