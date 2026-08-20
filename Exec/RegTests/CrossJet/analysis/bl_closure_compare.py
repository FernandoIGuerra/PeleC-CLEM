"""Non-reacting boundary-layer comparison across SGS closures, all CLEM/Maxwell.

Every case is read at the SAME physical time from data that already exists.
Panels:
  (a) incoming BL profile at x = 3.5 cm, with the first cell centre marked
  (b) wall-cell Mach along x   -> where a subsonic channel exists
  (c) wall-cell u/uinf, zoom on the injector -> separation bubble
  (d) eddy viscosity nu_t/nu at x = 3.5 cm, recomputed offline from the
      resolved velocity field with each model's own formula
"""
import os
import glob
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from amrread import read_plotfile_component as rd, tof

C = "/home/pc08/Fernando/CLEMNew/PeleC-CLEM/Exec/RegTests/CrossJet/runs"
R = "/home/pc08/PeleC/PeleC/Exec/RegTests/CrossJet"

# Okabe-Ito (published CVD-safe set); the reference is a neutral baseline.
CASES = [
    ("reference: no model (pristine PeleC)", R + "/plt_React*",                      "#000000", "smag"),
    ("(i) no closure, do_les=0",             C + "/01-nonreacting-CLEM/plt[0-9]*",   "#0072B2", "smag"),
    ("(iii) WALE",                           C + "/04-reacting-CLEM-WALE-nowallmodel/plt_WALE*", "#009E73", "wale"),
    ("(iv) WALE + wall model",               C + "/03-reacting-CLEM-WALE-wallmodel/plt_WM*",     "#CC79A7", "wale"),
    ("(ii) Smagorinsky, Cs=0.16",            C + "/02-reacting-CLEM-smagorinsky/plt_React*",     "#D55E00", "smag"),
]
TT, UINF, XSTA = 1.8e-4, 9.30e4, 3.5
CS, CW = 0.16, 0.325
INK, MUTED = "#1a1a1a", "#6b6b6b"


def near(pat, t):
    c = [p for p in glob.glob(pat) if ".old." not in p]
    return min(c, key=lambda p: abs(tof(p) - t))


def grads(u, v, w, dx):
    """Cell-centred central-difference velocity gradient tensor g[a][b]=du_a/dx_b."""
    f = [u, v, w]
    return [[np.gradient(f[a], dx[b], axis=b, edge_order=2) for b in range(3)]
            for a in range(3)]


def nut_of(model, u, v, w, dx):
    g = grads(u, v, w, dx)
    delta = float(np.prod(dx)) ** (1.0 / 3.0)
    S = [[0.5 * (g[a][b] + g[b][a]) for b in range(3)] for a in range(3)]
    SS = sum(S[a][b] ** 2 for a in range(3) for b in range(3))
    if model == "smag":
        return (CS * delta) ** 2 * np.sqrt(2.0 * SS)
    # WALE: traceless symmetric part of g^2
    g2 = [[sum(g[a][c] * g[c][b] for c in range(3)) for b in range(3)]
          for a in range(3)]
    tr = sum(g2[a][a] for a in range(3)) / 3.0
    Sd = [[0.5 * (g2[a][b] + g2[b][a]) - (tr if a == b else 0.0)
           for b in range(3)] for a in range(3)]
    SdSd = sum(Sd[a][b] ** 2 for a in range(3) for b in range(3))
    den = SS ** 2.5 + SdSd ** 1.25
    return np.where(den > 0.0, (CW * delta) ** 2 * SdSd ** 1.5 / np.maximum(den, 1e-300), 0.0)


def bl(u, y):
    i = int(np.argmax(u >= 0.90 * UINF))
    d90 = np.interp(0.90 * UINF, u[: i + 1], y[: i + 1]) if i > 0 else np.nan
    m = y <= max(4.0 * d90, 3.0 * (y[1] - y[0]))
    f = np.clip(u[m] / UINF, 0.0, None)
    ds = np.trapezoid(1.0 - f, y[m])
    th = np.trapezoid(f * (1.0 - f), y[m])
    return d90, ds, th, (ds / th if th > 0 else np.nan)


fig, ax = plt.subplots(2, 2, figsize=(12.6, 8.0), constrained_layout=True)
(a, b), (c, d) = ax
rows = []

for lab, pat, col, model in CASES:
    p = near(pat, TT)
    u, t, lo, dx = rd(p, "x_velocity")
    v, *_ = rd(p, "y_velocity")
    w, *_ = rd(p, "z_velocity")
    M, *_ = rd(p, "MachNumber")
    T, *_ = rd(p, "Temp")
    rho, *_ = rd(p, "density")
    try:
        yh2o = rd(p, "Y(H2O)")[0].max()
    except ValueError:
        yh2o = 0.0

    xc = lo[0] + (np.arange(u.shape[0]) + 0.5) * dx[0]
    yc = lo[1] + (np.arange(u.shape[1]) + 0.5) * dx[1]
    ix = int(XSTA / dx[0])

    up = u[ix].mean(axis=1)                       # spanwise-averaged profile
    d90, ds, th, H = bl(up, yc)
    Mw, uw = M[:, 0, :].mean(axis=1), u[:, 0, :].mean(axis=1)

    # Sutherland-free nu: use PeleC's own transport is not available offline;
    # nu from mu(T) of air via Sutherland is accurate to a few % for 77% N2.
    mu = 1.716e-4 * (T / 273.15) ** 1.5 * (273.15 + 110.4) / (T + 110.4)
    nu = mu / rho
    nut = nut_of(model, u, v, w, dx)
    ratio = (nut / nu)[ix].mean(axis=1)

    # (i) and (iii) land on top of each other - draw (i) thick underneath and
    # (iii) dashed on top so the coincidence is visible rather than hidden.
    ls, lw = ("--", 2.0) if col == "#000000" else ("-", 2.0)
    if lab.startswith("(i) "):
        lw = 5.0
    if lab.startswith("(iii)"):
        ls, lw = (0, (5, 3)), 1.8
    mk = dict(marker="o", ms=3.2, mfc="white", mew=1.1) if lab.startswith("(i) ") else {}
    a.plot(up / UINF, yc * 10, color=col, lw=lw, ls=ls, label=lab, **mk)
    b.plot(xc, Mw, color=col, lw=lw, ls=ls)
    c.plot(xc, uw / UINF, color=col, lw=lw, ls=ls)
    d.semilogx(np.maximum(ratio, 1e-4), yc * 10, color=col, lw=lw, ls=ls, **mk)

    nsub = int((Mw < 1.0).sum())
    rows.append((lab, t * 1e3, T.max(), yh2o, d90 * 10, ds * 10, H,
                 Mw[ix], 100.0 * nsub / len(xc), ratio[0]))

    if lab.startswith("reference"):
        # Is this layer laminar at all? Freestream state, well above the wall.
        nu_inf = float(nu[ix, -6:, :].mean())
        Rex = UINF * XSTA / nu_inf
        print(f"freestream at x={XSTA} cm: nu = {nu_inf:.3f} cm^2/s, "
              f"Re_x = {Rex:.2e}  ->  {'LAMINAR' if Rex < 5e5 else 'transitional/turbulent'}"
              f"\nBlasius delta99 = {5.0*XSTA/np.sqrt(Rex)*10:.3f} mm  vs  "
              f"cell dy = {dx[1]*10:.4f} mm"
              f"\nwall-model argument R = y1*u1/nu = "
              f"{0.5*dx[1]*u[ix,0,:].mean()/nu_inf:.0f}  (Spalding reads this as a log layer)\n")

y1 = 0.5 * dx[1] * 10
for axx in (a, d):
    axx.axhline(y1, color=MUTED, lw=1.0, ls=":")
    axx.text(0.98, y1, "  first cell centre", transform=axx.get_yaxis_transform(),
             ha="right", va="bottom", fontsize=8, color=MUTED)
    axx.set_ylim(0, 2.2)
    axx.set_ylabel("y [mm]", color=INK, fontsize=9)

a.set_xlabel(r"$u/u_\infty$", color=INK, fontsize=9)
a.set_title(f"(a) incoming boundary layer, x = {XSTA} cm  "
            r"(markers = cell centres)", loc="left", fontsize=10.5, color=INK)
a.legend(fontsize=8.5, frameon=False, loc="center left",
         bbox_to_anchor=(0.02, 0.62))
d.set_xlabel(r"$\nu_t/\nu$   (recomputed offline)", color=INK, fontsize=9)
d.set_title(r"(d) eddy viscosity each closure produces, x = 3.5 cm",
            loc="left", fontsize=10.5, color=INK)
d.axvline(1.0, color=MUTED, lw=1.0, ls="--")
d.text(0.03, 0.97, "for (i) and the reference this $\\nu_t$ acts on SCALARS only\n"
       "(CLEM stirring) — neither adds an SGS momentum stress",
       transform=d.transAxes, va="top", ha="left", fontsize=8, color=MUTED)

b.axhline(1.0, color=MUTED, lw=1.0, ls="--")
b.axvline(4.2, color=MUTED, lw=0.9, ls=":")
b.set_xlabel("x [cm]", color=INK, fontsize=9)
b.set_ylabel("Mach, first cell", color=INK, fontsize=9)
b.set_title("(b) wall-cell Mach — below 1 is a flashback path", loc="left",
            fontsize=10.5, color=INK)

c.axhline(0.0, color=MUTED, lw=1.0, ls="--")
c.axvline(4.2, color=MUTED, lw=0.9, ls=":")
c.set_xlim(2.5, 5.5)
c.set_xlabel("x [cm]   (dotted: injector)", color=INK, fontsize=9)
c.set_ylabel(r"$u/u_\infty$, first cell", color=INK, fontsize=9)
c.set_title("(c) separation upstream of the jet (u < 0)", loc="left",
            fontsize=10.5, color=INK)

for axx in (a, b, c, d):
    axx.grid(color="#ededed", lw=0.7)
    axx.tick_params(colors=MUTED, labelsize=8, length=3)
    for s in axx.spines.values():
        s.set_color("#cccccc")

fig.suptitle("CLEM (Maxwell) non-reacting: what each SGS closure does to the wall, "
             f"t $\\approx$ {TT*1e3:.2f} ms", fontsize=12, color=INK, x=0.005, ha="left")
out = "/home/pc08/Fernando/CLEMNew/PeleC-CLEM/Exec/RegTests/CrossJet/analysis/bl_closure_compare.png"
fig.savefig(out, dpi=160, facecolor="white")
print("wrote", out, "\n")

hdr = (f"{'case':<38}{'t[ms]':>7}{'maxT':>7}{'YH2O':>9}{'d90[mm]':>9}"
       f"{'d*[mm]':>8}{'H':>6}{'M_w':>7}{'M<1 %x':>8}{'nut/nu':>9}")
print(hdr); print("-" * len(hdr))
for r in rows:
    print(f"{r[0]:<38}{r[1]:7.3f}{r[2]:7.0f}{r[3]:9.1e}{r[4]:9.3f}"
          f"{r[5]:8.3f}{r[6]:6.2f}{r[7]:7.3f}{r[8]:8.0f}{r[9]:9.2f}")
