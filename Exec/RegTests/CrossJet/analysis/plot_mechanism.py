"""Upstream-only pressure comparison plus the boundary-layer profile that
explains it. Same mid-plane, x restricted to the region upstream of the jet."""
import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_inlet_pressure import read_plotfile_component

CLEM = "/home/pc08/Fernando/CLEMNew/PeleC-CLEM/Exec/RegTests/CrossJet/plt_React15631"
REF = "/home/pc08/PeleC/PeleC/Exec/RegTests/CrossJet/plt_React15200"
P0, U0, XMAX = 1.01325e6, 9.30e4, 4.0
CASES = (("CLEM  (pelec.do_les = 1)", CLEM, "#b2182b"),
         ("no-model reference  (pelec.do_les = 0)", REF, "#2166ac"))


def midplane(pltdir, var):
    a, t, lo, dx = read_plotfile_component(pltdir, var)
    zc = lo[2] + (np.arange(a.shape[2]) + 0.5) * dx[2]
    k = np.argsort(np.abs(zc - 0.25))[:2]
    nx = int(round(XMAX / dx[0]))
    return a[:nx, :, k].mean(axis=2), t, lo, dx


data = {}
for label, d, _ in CASES:
    p, t, lo, dx = midplane(d, "pressure")
    u, _, _, _ = midplane(d, "x_velocity")
    data[label] = (p, u, t, lo, dx)
    f = 100.0 * (p - P0) / P0
    print(f"{label}\n  upstream x<{XMAX} cm:  mean {f.mean():+.2f} %   "
          f"rms {f.std():.2f} %   min {f.min():+.2f} %   max {f.max():+.2f} %")

lim = 25.0
INK, MUTED = "#1a1a1a", "#6b6b6b"
fig = plt.figure(figsize=(12.2, 6.0), constrained_layout=True)
gs = fig.add_gridspec(2, 2, width_ratios=[3.05, 1.0])
axp = [fig.add_subplot(gs[0, 0]), None]
axp[1] = fig.add_subplot(gs[1, 0], sharex=axp[0], sharey=axp[0])
axu = fig.add_subplot(gs[:, 1])

for ax, (label, d, _c) in zip(axp, CASES):
    p, u, t, lo, dx = data[label]
    xe = lo[0] + np.arange(p.shape[0] + 1) * dx[0]
    ye = lo[1] + np.arange(p.shape[1] + 1) * dx[1]
    m = ax.pcolormesh(xe, ye, (100 * (p - P0) / P0).T, cmap="coolwarm",
                      vmin=-lim, vmax=lim, shading="flat", rasterized=True)
    ax.set_aspect("equal")
    ax.set_ylabel("y  [cm]", color=INK, fontsize=9)
    ax.set_title(f"{label}   —   t = {t*1e3:.4f} ms", loc="left",
                 fontsize=10, color=INK, pad=4)
    ax.tick_params(colors=MUTED, labelsize=8, length=3)
    for s in ax.spines.values():
        s.set_color("#cccccc")
axp[1].set_xlabel("x  [cm]", color=INK, fontsize=9)
cb = fig.colorbar(m, ax=axp, location="bottom", shrink=0.6, pad=0.02,
                  aspect=45)
cb.set_label(r"$(p-p_\infty)/p_\infty$   [%]", color=INK, fontsize=9)
cb.ax.tick_params(colors=MUTED, labelsize=8)
cb.outline.set_edgecolor("#cccccc")

# Boundary layer at x = 2 cm, well upstream of the injector
XPROBE = 2.0
for label, d, c in CASES:
    p, u, t, lo, dx = data[label]
    i = int(XPROBE / dx[0])
    yc = lo[1] + (np.arange(u.shape[1]) + 0.5) * dx[1]
    prof = u[i, :] / U0
    axu.plot(prof, yc, color=c, lw=2.0, label=label.split("  (")[0])
    # The CLEM profile has no uniform core left, so a delta_99 threshold is
    # ill-posed. Report delta_90 (a robust thickness) and the DISPLACEMENT
    # thickness, which is what actually turns the flow and launches the waves.
    half = yc < 0.75
    d90 = np.interp(0.90, prof[half], yc[half])
    rho, _, _, _ = midplane(d, "density")
    rhou = rho[i, :] * u[i, :]
    dstar = np.trapezoid(1.0 - rhou[half] / rhou.max(), yc[half])
    axu.plot([0, 0.90], [d90, d90], color=c, lw=0.9, ls=":", alpha=0.8)
    print(f"{label}:  u_core = {prof.max():.3f} u_inf,  "
          f"delta90 = {d90*10:.3f} mm ({d90/dx[1]:.1f} cells),  "
          f"delta* = {dstar*10:.3f} mm")
axu.set_xlabel(r"$u/u_\infty$   at  x = 2.0 cm", color=INK, fontsize=9)
axu.set_ylabel("y  [cm]", color=INK, fontsize=9)
axu.set_title("boundary layer (dotted: $\\delta_{90}$)", loc="left",
              fontsize=10, color=INK, pad=4)
axu.legend(fontsize=8, frameon=False, loc="lower left")
axu.grid(color="#eeeeee", lw=0.7)
axu.tick_params(colors=MUTED, labelsize=8, length=3)
for s in axu.spines.values():
    s.set_color("#cccccc")

fig.suptitle("Upstream of the injector: the CLEM run's wave train, and the "
             "boundary layer that drives it", fontsize=11.5, color=INK,
             x=0.01, ha="left")
out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "inlet_mechanism.png")
fig.savefig(out, dpi=170, facecolor="white")
print("wrote", out)
