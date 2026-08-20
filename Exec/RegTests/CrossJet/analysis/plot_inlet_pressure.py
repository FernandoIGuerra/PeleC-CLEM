"""Compare the inlet pressure field of the CLEM run against the no-model
reference, on the spanwise mid-plane z = 0.25 cm, x in [0, 5.5] cm.

Minimal single-level AMReX plotfile reader (max_level = 0), then a two-panel
diverging map of (p - p_inlet)/p_inlet.
"""
import os
import re
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BOX_RE = re.compile(r"\(\((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\)")


def read_plotfile_component(pltdir, varname):
    """Return (data[nx,ny,nz], time, prob_lo, dx) for one component."""
    hdr = open(os.path.join(pltdir, "Header")).read().split("\n")
    nvar = int(hdr[1])
    names = [hdr[2 + i].strip() for i in range(nvar)]
    comp = names.index(varname)
    p = 2 + nvar
    dim = int(hdr[p])
    time = float(hdr[p + 1])
    prob_lo = np.array([float(v) for v in hdr[p + 3].split()])
    dom = BOX_RE.search(hdr[p + 6])
    lo = np.array([int(dom.group(i)) for i in (1, 2, 3)])
    hi = np.array([int(dom.group(i)) for i in (4, 5, 6)])
    dx = np.array([float(v) for v in hdr[p + 8].split()])
    shape = hi - lo + 1

    ch = open(os.path.join(pltdir, "Level_0", "Cell_H")).read().split("\n")
    ncomp = int(ch[2])
    nboxes = int(ch[4].lstrip("(").split()[0])
    boxes = []
    for line in ch[5 : 5 + nboxes]:
        m = BOX_RE.search(line)
        boxes.append(
            (
                np.array([int(m.group(i)) for i in (1, 2, 3)]),
                np.array([int(m.group(i)) for i in (4, 5, 6)]),
            )
        )
    fabs = [l.split()[1:] for l in ch if l.startswith("FabOnDisk:")]
    assert len(fabs) == nboxes, (len(fabs), nboxes)

    out = np.full(shape, np.nan)
    for (blo, bhi), (fname, offset) in zip(boxes, fabs):
        bshape = bhi - blo + 1
        npts = int(np.prod(bshape))
        with open(os.path.join(pltdir, "Level_0", fname), "rb") as f:
            f.seek(int(offset))
            f.readline()  # ASCII "FAB ((8, ...))((lo) (hi) (0,0,0)) ncomp"
            f.seek(comp * npts * 8, os.SEEK_CUR)
            buf = np.frombuffer(f.read(npts * 8), dtype="<f8")
        # AMReX writes each component in Fortran order: i fastest, then j, k
        out[
            blo[0] : bhi[0] + 1, blo[1] : bhi[1] + 1, blo[2] : bhi[2] + 1
        ] = buf.reshape(bshape, order="F")
    assert not np.isnan(out).any(), "plotfile did not cover the domain"
    return out, time, prob_lo, dx


def midplane(pltdir, z_cut=0.25, x_max=5.5):
    p, t, lo, dx = read_plotfile_component(pltdir, "pressure")
    nz = p.shape[2]
    zc = lo[2] + (np.arange(nz) + 0.5) * dx[2]
    # z = 0.25 falls on a cell face; average the two straddling cell centres
    k = np.argsort(np.abs(zc - z_cut))[:2]
    sl = p[:, :, k].mean(axis=2)
    nx = int(round((x_max - lo[0]) / dx[0]))
    xe = lo[0] + np.arange(nx + 1) * dx[0]
    ye = lo[1] + np.arange(sl.shape[1] + 1) * dx[1]
    return sl[:nx, :], xe, ye, t


CLEM = "/home/pc08/Fernando/CLEMNew/PeleC-CLEM/Exec/RegTests/CrossJet/plt_React15631"
REF = "/home/pc08/PeleC/PeleC/Exec/RegTests/CrossJet/plt_React15200"
P0 = 1.01325e6  # inlet pressure, dyn/cm^2

panels = []
for label, d in (("CLEM  (pelec.do_les = 1, Cs = 0.16)", CLEM),
                 ("no-model reference  (pelec.do_les = 0)", REF)):
    sl, xe, ye, t = midplane(d)
    panels.append((label, 100.0 * (sl - P0) / P0, xe, ye, t))
    print(f"{label}: t = {t:.6e} s   "
          f"(p-p0)/p0 range = [{100*(sl.min()-P0)/P0:+.2f}, "
          f"{100*(sl.max()-P0)/P0:+.2f}] %   "
          f"rms = {100*np.std(sl)/P0:.3f} %")

# Symmetric diverging limits, robust to the jet's strong bow shock
lim = np.percentile(np.abs(panels[0][1]), 99.0)
print(f"\ncolour limits: +/- {lim:.2f} %")

INK, MUTED = "#1a1a1a", "#6b6b6b"
fig, axes = plt.subplots(
    2, 1, figsize=(11.5, 5.4), sharex=True, sharey=True,
    constrained_layout=True,
)
for ax, (label, f, xe, ye, t) in zip(axes, panels):
    m = ax.pcolormesh(xe, ye, f.T, cmap="coolwarm", vmin=-lim, vmax=lim,
                      shading="flat", rasterized=True)
    ax.set_aspect("equal")   # so wave angles are geometrically true
    ax.set_ylabel("y  [cm]", color=INK, fontsize=9)
    ax.set_title(f"{label}      t = {t*1e3:.4f} ms",
                 loc="left", fontsize=10, color=INK, pad=4)
    ax.axvline(4.2, color=MUTED, lw=0.8, ls=":")
    for s in ax.spines.values():
        s.set_color("#cccccc")
    ax.tick_params(colors=MUTED, labelsize=8, length=3)

# Mach-angle guide, mu = asin(1/1.5) = 41.8 deg, drawn from the inlet corners
mu = np.arcsin(1.0 / 1.5)
dxdy = 1.0 / np.tan(mu)
for y0, sgn in ((0.0, +1), (1.5, -1)):
    axes[0].plot([0.0, dxdy * 1.5], [y0, y0 + sgn * 1.5],
                 color="#222222", lw=1.0, ls="--", alpha=0.55)
axes[0].text(1.75, 1.30, r"Mach angle $\mu=41.8^\circ$", fontsize=8,
             color="#222222", alpha=0.8)
axes[1].set_xlabel("x  [cm]      (dotted line: H$_2$ injector at x = 4.2 cm)",
                   color=INK, fontsize=9)

cb = fig.colorbar(m, ax=axes, location="right", shrink=0.85, pad=0.015)
cb.set_label(r"$(p-p_\infty)/p_\infty$   [%]", color=INK, fontsize=9)
cb.ax.tick_params(colors=MUTED, labelsize=8)
cb.outline.set_edgecolor("#cccccc")

fig.suptitle(
    "Pressure on the spanwise mid-plane z = 0.25 cm — inlet region",
    fontsize=11.5, color=INK, x=0.01, ha="left",
)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "inlet_pressure_clem_vs_reference.png")
fig.savefig(out, dpi=170, facecolor="white")
print("wrote", out)
