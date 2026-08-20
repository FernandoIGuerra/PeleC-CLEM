"""Minimal single-level AMReX plotfile reader (max_level = 0)."""
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



def tof(d):
    import os
    f = open(os.path.join(d, "Header")).read().split("\n")
    return float(f[2 + int(f[1]) + 1])


def midplane(pltdir, var, z_cut=0.25, x_max=None):
    """Component on the z = z_cut plane, averaged over the two straddling cells."""
    import numpy as np
    a, t, lo, dx = read_plotfile_component(pltdir, var)
    zc = lo[2] + (np.arange(a.shape[2]) + 0.5) * dx[2]
    k = np.argsort(np.abs(zc - z_cut))[:2]
    nx = a.shape[0] if x_max is None else int(round((x_max - lo[0]) / dx[0]))
    return a[:nx, :, k].mean(axis=2), t, lo, dx
