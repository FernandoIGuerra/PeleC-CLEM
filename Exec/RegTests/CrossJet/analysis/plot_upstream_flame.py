"""Upstream ignition / flashback diagnostic.

Question: the CLEM run ignites upstream of the injector and the flame walks
toward the inlet; the no-model reference does not. Is that the SGS model?

The discriminator is whether a SUBSONIC channel exists along the wall. A flame
cannot travel upstream through a M = 1.5 core; it can only crawl up a subsonic
near-wall layer. The artificially thick Smagorinsky boundary layer makes that
layer much deeper, so we measure it directly, alongside T and Y(H2).
"""
import os
import glob
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_inlet_pressure import read_plotfile_component

CDIR = "/home/pc08/Fernando/CLEMNew/PeleC-CLEM/Exec/RegTests/CrossJet"
RDIR = "/home/pc08/PeleC/PeleC/Exec/RegTests/CrossJet"
XMAX, XJET = 5.5, 4.2


def time_of(d):
    f = open(os.path.join(d, "Header")).read().split("\n")
    return float(f[2 + int(f[1]) + 1])


def match(dirpat, t):
    c = [(abs(time_of(d) - t), d) for d in glob.glob(dirpat)]
    return min(c)[1]


clem = max(glob.glob(f"{CDIR}/plt_React*"), key=time_of)
tc = time_of(clem)
ref = match(f"{RDIR}/plt_React*", tc)
print(f"CLEM {os.path.basename(clem)}  t = {tc:.6e}")
print(f"REF  {os.path.basename(ref)}  t = {time_of(ref):.6e}\n")

CASES = (("CLEM  (SGS on)", clem, "#b2182b"),
         ("no-model reference", ref, "#2166ac"))


def mid(pltdir, var):
    a, t, lo, dx = read_plotfile_component(pltdir, var)
    zc = lo[2] + (np.arange(a.shape[2]) + 0.5) * dx[2]
    k = np.argsort(np.abs(zc - 0.25))[:2]
    nx = int(round(XMAX / dx[0]))
    return a[:nx, :, k].mean(axis=2), lo, dx


D = {}
for label, d, _ in CASES:
    T, lo, dx = mid(d, "Temp")
    YH2, _, _ = mid(d, "Y(H2)")
    M, _, _ = mid(d, "MachNumber")
    D[label] = (T, YH2, M, lo, dx)
    iu = slice(0, int(4.0 / dx[0]))          # upstream of the injector
    yc = lo[1] + (np.arange(T.shape[1]) + 0.5) * dx[1]
    # how far upstream does fuel reach, and how hot does it get there
    has_fuel = np.where(YH2[iu, :].max(axis=1) > 1e-3)[0]
    xfuel = (has_fuel.min() * dx[0]) if has_fuel.size else np.nan
    # subsonic layer at x = 3 cm: fraction of the duct height with M < 1
    i3 = int(3.0 / dx[0])
    sub = (M[i3, :] < 1.0).sum() * dx[1]
    print(f"{label}")
    print(f"   upstream (x<4cm):  max T = {T[iu,:].max():7.1f} K    "
          f"max Y(H2) = {YH2[iu,:].max():.4f}")
    print(f"   fuel (Y_H2>1e-3) reaches upstream to x = {xfuel:.2f} cm")
    print(f"   subsonic layer at x=3cm: {sub*10:.2f} mm of 15 mm "
          f"({100*sub/1.5:.1f} % of the duct)\n")

INK, MUTED = "#1a1a1a", "#6b6b6b"
fig, axes = plt.subplots(2, 2, figsize=(13.2, 5.0), sharex=True, sharey=True,
                         constrained_layout=True)
for col, (label, d, _c) in enumerate(CASES):
    T, YH2, M, lo, dx = D[label]
    xe = lo[0] + np.arange(T.shape[0] + 1) * dx[0]
    ye = lo[1] + np.arange(T.shape[1] + 1) * dx[1]
    mT = axes[0, col].pcolormesh(xe, ye, T.T, cmap="inferno", vmin=300,
                                 vmax=2800, shading="flat", rasterized=True)
    mY = axes[1, col].pcolormesh(xe, ye, YH2.T, cmap="Blues", vmin=0,
                                 vmax=0.06, shading="flat", rasterized=True)
    axes[0, col].set_title(f"{label}   —   t = {tc*1e3:.3f} ms", loc="left",
                           fontsize=10, color=INK, pad=4)
    for r in (0, 1):
        ax = axes[r, col]
        ax.set_aspect("equal")
        ax.axvline(XJET, color="w", lw=0.8, ls=":")
        ax.tick_params(colors=MUTED, labelsize=8, length=3)
        for s in ax.spines.values():
            s.set_color("#cccccc")
    axes[1, col].set_xlabel("x  [cm]", color=INK, fontsize=9)
axes[0, 0].set_ylabel("y  [cm]", color=INK, fontsize=9)
axes[1, 0].set_ylabel("y  [cm]", color=INK, fontsize=9)

cb1 = fig.colorbar(mT, ax=axes[0, :], location="right", shrink=0.9, pad=0.01)
cb1.set_label("T  [K]", color=INK, fontsize=9)
cb2 = fig.colorbar(mY, ax=axes[1, :], location="right", shrink=0.9, pad=0.01)
cb2.set_label(r"$Y_{H_2}$", color=INK, fontsize=9)
for cb in (cb1, cb2):
    cb.ax.tick_params(colors=MUTED, labelsize=8)
    cb.outline.set_edgecolor("#cccccc")

fig.suptitle("Upstream ignition: temperature (top) and H$_2$ (bottom), "
             "mid-plane z = 0.25 cm   —   dotted line: injector",
             fontsize=11.5, color=INK, x=0.01, ha="left")
out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "upstream_flame.png")
fig.savefig(out, dpi=170, facecolor="white")
print("wrote", out)
