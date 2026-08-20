"""Three-way comparison: no-model reference, old CLEM (Smagorinsky),
new CLEM (WALE + wall model). Pressure at a pre-flame time, temperature once
the flame exists, and the near-wall Mach diagnostics that decide whether a
flashback channel exists."""
import os
import glob
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_inlet_pressure import read_plotfile_component as rd

C = "/home/pc08/Fernando/CLEMNew/PeleC-CLEM/Exec/RegTests/CrossJet"
R = "/home/pc08/PeleC/PeleC/Exec/RegTests/CrossJet"
P0, XMAX = 1.01325e6, 5.5
CASES = [("no-model reference", R + "/plt_React*", "#2166ac"),
         ("CLEM + Smagorinsky (old)", C + "/plt_React*", "#b2182b"),
         ("CLEM + WALE + wall model", C + "/plt_WM*", "#1b7837")]
INK, MUTED = "#1a1a1a", "#6b6b6b"


def tof(d):
    f = open(os.path.join(d, "Header")).read().split("\n")
    return float(f[2 + int(f[1]) + 1])


def near(pat, t):
    c = [p for p in glob.glob(pat) if ".old." not in p]
    return min(c, key=lambda p: abs(tof(p) - t))


def mid(d, var, xmax=XMAX):
    a, t, lo, dx = rd(d, var)
    zc = lo[2] + (np.arange(a.shape[2]) + 0.5) * dx[2]
    k = np.argsort(np.abs(zc - 0.25))[:2]
    nx = int(round(xmax / dx[0]))
    return a[:nx, :, k].mean(axis=2), t, lo, dx


def field_fig(var, tt, cmap, vmin, vmax, label, title, fname, pct=False):
    fig, axes = plt.subplots(3, 1, figsize=(11.0, 6.2), sharex=True,
                             sharey=True, constrained_layout=True)
    for ax, (lab, pat, _c) in zip(axes, CASES):
        d = near(pat, tt)
        f, t, lo, dx = mid(d, var)
        if pct:
            f = 100.0 * (f - P0) / P0
        xe = lo[0] + np.arange(f.shape[0] + 1) * dx[0]
        ye = lo[1] + np.arange(f.shape[1] + 1) * dx[1]
        m = ax.pcolormesh(xe, ye, f.T, cmap=cmap, vmin=vmin, vmax=vmax,
                          shading="flat", rasterized=True)
        ax.set_aspect("equal")
        ax.axvline(4.2, color="w", lw=0.9, ls=":")
        ax.set_ylabel("y [cm]", color=INK, fontsize=9)
        ax.set_title(f"{lab}   —   t = {t*1e3:.3f} ms", loc="left",
                     fontsize=10, color=INK, pad=3)
        ax.tick_params(colors=MUTED, labelsize=8, length=3)
        for s in ax.spines.values():
            s.set_color("#cccccc")
    axes[-1].set_xlabel("x [cm]   (dotted: injector)", color=INK, fontsize=9)
    cb = fig.colorbar(m, ax=axes, location="right", shrink=0.85, pad=0.015)
    cb.set_label(label, color=INK, fontsize=9)
    cb.ax.tick_params(colors=MUTED, labelsize=8)
    cb.outline.set_edgecolor("#cccccc")
    fig.suptitle(title, fontsize=11.5, color=INK, x=0.01, ha="left")
    fig.savefig(fname, dpi=165, facecolor="white")
    print("wrote", fname)


field_fig("pressure", 1.8e-4, "coolwarm", -25, 25,
          r"$(p-p_\infty)/p_\infty$  [%]",
          "Inlet pressure, pre-flame — the wave lattice is gone",
          "wm_pressure.png", pct=True)
field_fig("Temp", 5.6e-4, "inferno", 300, 2800, "T [K]",
          "Temperature once the flame exists", "wm_temperature.png")

# --- near-wall Mach diagnostics -------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(11.5, 3.9), constrained_layout=True)
for lab, pat, c in CASES:
    d = near(pat, 1.8e-4)
    M, t, lo, dx = mid(d, "MachNumber", 12.5)
    xc = lo[0] + (np.arange(M.shape[0]) + 0.5) * dx[0]
    ax[0].plot(xc, M[:, 0], color=c, lw=1.8, label=lab)
    ts, ms = [], []
    for T0 in np.arange(0.02, 0.57, 0.04) * 1e-3:
        dd = near(pat, T0)
        MM, tt2, lo2, dx2 = mid(dd, "MachNumber", 12.5)
        i3 = int(3.0 / dx2[0])
        ts.append(tt2 * 1e3)
        ms.append(MM[i3, 0])
    ax[1].plot(ts, ms, color=c, lw=1.8, label=lab)
ax[0].axhline(1.0, color=MUTED, lw=0.9, ls="--")
ax[1].axhline(1.0, color=MUTED, lw=0.9, ls="--")
ax[0].set_xlabel("x [cm]", color=INK, fontsize=9)
ax[0].set_ylabel("Mach, wall cell", color=INK, fontsize=9)
ax[0].set_title("wall-row Mach along x,  t = 0.18 ms (pre-flame)", loc="left",
                fontsize=10, color=INK)
ax[0].axvline(4.2, color=MUTED, lw=0.8, ls=":")
ax[1].set_xlabel("t [ms]", color=INK, fontsize=9)
ax[1].set_ylabel("Mach, wall cell at x = 3 cm", color=INK, fontsize=9)
ax[1].set_title("history at x = 3 cm — a subsonic wall cell is a "
                "flashback path", loc="left", fontsize=10, color=INK)
for a in ax:
    a.grid(color="#eeeeee", lw=0.7)
    a.legend(fontsize=8, frameon=False)
    a.tick_params(colors=MUTED, labelsize=8, length=3)
    for s in a.spines.values():
        s.set_color("#cccccc")
fig.suptitle("Does a subsonic near-wall channel exist? (dashed: M = 1)",
             fontsize=11.5, color=INK, x=0.01, ha="left")
fig.savefig("wm_mach.png", dpi=165, facecolor="white")
print("wrote wm_mach.png")
