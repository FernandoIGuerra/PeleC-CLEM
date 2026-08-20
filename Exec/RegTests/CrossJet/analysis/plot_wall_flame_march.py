"""Wall-adjacent temperature along x (y = 0 wall row, z = 0.25 cm) at a
sequence of times from 0.15 to 0.45 ms — the upstream march of the flame."""
import os
import glob
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from plot_inlet_pressure import read_plotfile_component

CDIR = "/home/pc08/Fernando/CLEMNew/PeleC-CLEM/Exec/RegTests/CrossJet"
RDIR = "/home/pc08/PeleC/PeleC/Exec/RegTests/CrossJet"
TIMES = np.arange(0.15, 0.4501, 0.05) * 1e-3     # s
XJET = 4.2


def time_of(d):
    f = open(os.path.join(d, "Header")).read().split("\n")
    return float(f[2 + int(f[1]) + 1])


def wall_line(pltdir):
    """T in the wall-adjacent cell row (j = 0) on the z = 0.25 plane."""
    T, t, lo, dx = read_plotfile_component(pltdir, "Temp")
    zc = lo[2] + (np.arange(T.shape[2]) + 0.5) * dx[2]
    k = np.argsort(np.abs(zc - 0.25))[:2]
    xc = lo[0] + (np.arange(T.shape[0]) + 0.5) * dx[0]
    return xc, T[:, 0, k].mean(axis=1), t


CASES = (("CLEM  (SGS on)", CDIR), ("no-model reference", RDIR))
norm = Normalize(TIMES[0] * 1e3, TIMES[-1] * 1e3)
cmap = plt.get_cmap("viridis")
INK, MUTED = "#1a1a1a", "#6b6b6b"

fig, axes = plt.subplots(2, 1, figsize=(10.0, 6.4), sharex=True, sharey=True,
                         constrained_layout=True)
for ax, (label, d) in zip(axes, CASES):
    dirs = glob.glob(f"{d}/plt_React*")
    for tt in TIMES:
        best = min(dirs, key=lambda p: abs(time_of(p) - tt))
        xc, T, t = wall_line(best)
        ax.plot(xc, T, color=cmap(norm(t * 1e3)), lw=1.6)
    ax.axvline(XJET, color=MUTED, lw=0.9, ls=":")
    ax.set_title(label, loc="left", fontsize=10, color=INK, pad=4)
    ax.set_ylabel("T  [K]", color=INK, fontsize=9)
    ax.grid(color="#eeeeee", lw=0.7)
    ax.tick_params(colors=MUTED, labelsize=8, length=3)
    for s in ax.spines.values():
        s.set_color("#cccccc")
    print(f"{label}: T range over the window "
          f"[{min(T):.0f}, {max(T):.0f}] K at the last time")
axes[1].set_xlabel("x  [cm]      (dotted: H$_2$ injector at x = 4.2 cm)",
                   color=INK, fontsize=9)
axes[0].set_xlim(0, 12.5)

cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=axes,
                  location="right", shrink=0.85, pad=0.015)
cb.set_label("t  [ms]", color=INK, fontsize=9)
cb.ax.tick_params(colors=MUTED, labelsize=8)
cb.outline.set_edgecolor("#cccccc")

fig.suptitle("Wall-row temperature along x  (y = 0 wall, z = 0.25 cm), "
             "t = 0.15 → 0.45 ms", fontsize=11.5, color=INK, x=0.01, ha="left")
out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "wall_flame_march.png")
fig.savefig(out, dpi=170, facecolor="white")
print("wrote", out)
