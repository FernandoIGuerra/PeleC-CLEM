"""Causality: in the CLEM run, does the duct go subsonic BEFORE fuel and heat
appear upstream of the injector, or after? If choking comes first, the SGS
blockage caused the upstream ignition; if ignition comes first, the heat
release choked the duct (the paper's actual physics)."""
import os
import glob
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_inlet_pressure import read_plotfile_component

CDIR = "/home/pc08/Fernando/CLEMNew/PeleC-CLEM/Exec/RegTests/CrossJet"
RDIR = "/home/pc08/PeleC/PeleC/Exec/RegTests/CrossJet"
XPROBE, XUP = 3.0, 4.0   # probe station, and "upstream" cutoff


def series(pattern, label):
    rows = []
    for d in sorted(glob.glob(pattern)):
        try:
            M, t, lo, dx = read_plotfile_component(d, "MachNumber")
            T, _, _, _ = read_plotfile_component(d, "Temp")
            Y, _, _, _ = read_plotfile_component(d, "Y(H2)")
        except Exception as e:
            print("  skip", os.path.basename(d), e)
            continue
        zc = lo[2] + (np.arange(M.shape[2]) + 0.5) * dx[2]
        k = np.argsort(np.abs(zc - 0.25))[:2]
        i3, iu = int(XPROBE / dx[0]), slice(0, int(XUP / dx[0]))
        m = M[i3, :, k].mean(axis=0)
        rows.append((
            t,
            100.0 * (m < 1.0).sum() / m.size,          # % of duct subsonic
            T[iu, :, k].mean(axis=0).max(),            # max upstream T
            Y[iu, :, k].mean(axis=0).max(),            # max upstream Y_H2
        ))
    rows.sort()
    a = np.array(rows)
    print(f"{label}: {len(a)} plotfiles, t = {a[0,0]:.2e} .. {a[-1,0]:.2e}")
    return a


clem = series(f"{CDIR}/plt_React*", "CLEM")
ref = series(f"{RDIR}/plt_React*", "reference")

INK, MUTED = "#1a1a1a", "#6b6b6b"
CC, CR = "#b2182b", "#2166ac"
fig, ax = plt.subplots(3, 1, figsize=(8.6, 7.0), sharex=True,
                       constrained_layout=True)
for a, lab, c in ((clem, "CLEM  (SGS on)", CC), (ref, "no-model reference", CR)):
    ax[0].plot(a[:, 0] * 1e3, a[:, 1], color=c, lw=2.0, label=lab)
    ax[1].plot(a[:, 0] * 1e3, a[:, 2], color=c, lw=2.0, label=lab)
    ax[2].plot(a[:, 0] * 1e3, a[:, 3], color=c, lw=2.0, label=lab)

# when does the CLEM duct first go >50 % subsonic, and when does fuel arrive?
i_choke = np.argmax(clem[:, 1] > 50.0)
i_fuel = np.argmax(clem[:, 3] > 1e-3)
i_burn = np.argmax(clem[:, 2] > 1300.0)
tch, tf, tb = clem[i_choke, 0] * 1e3, clem[i_fuel, 0] * 1e3, clem[i_burn, 0] * 1e3
print(f"\nCLEM: duct >50% subsonic at t = {tch:.3f} ms")
print(f"CLEM: fuel upstream (Y_H2>1e-3) at t = {tf:.3f} ms")
print(f"CLEM: upstream T > 1300 K at t = {tb:.3f} ms")
for a_ in ax:
    a_.axvline(tch, color=MUTED, lw=0.9, ls="--")
    a_.grid(color="#eeeeee", lw=0.7)
    a_.tick_params(colors=MUTED, labelsize=8, length=3)
    for s in a_.spines.values():
        s.set_color("#cccccc")
ax[0].annotate(f"duct >50 % subsonic\nt = {tch:.3f} ms", (tch, 55),
               xytext=(tch + 0.04, 62), fontsize=8, color=INK,
               arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.8))

ax[0].set_ylabel(f"% of duct height\nsubsonic at x = {XPROBE} cm",
                 color=INK, fontsize=9)
ax[1].set_ylabel("max T upstream\nof injector  [K]", color=INK, fontsize=9)
ax[2].set_ylabel("max $Y_{H_2}$ upstream\nof injector", color=INK, fontsize=9)
ax[2].set_xlabel("t  [ms]", color=INK, fontsize=9)
ax[0].legend(fontsize=8, frameon=False, loc="center right")
fig.suptitle("Upstream of the injector (x < 4 cm): choking precedes fuel, "
             "fuel precedes burning", fontsize=11.5, color=INK, x=0.01,
             ha="left")
out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "choking_history.png")
fig.savefig(out, dpi=170, facecolor="white")
print("wrote", out)
