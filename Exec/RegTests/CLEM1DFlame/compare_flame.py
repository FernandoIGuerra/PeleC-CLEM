#!/usr/bin/env python3
"""
Compare the CLEM 1D flame against the reference flame solution.

Overlays the profiles (temperature, velocity, H2 and H2O mass fraction) of a
PeleC/CLEM plotfile line-out against the PREMIX/Cantera reference
(LiDryer_H2_p1_phi0_4000tu0300.dat), after aligning the two flame fronts
(location of max |dT/dx|).  Also prints the key integral metrics:

    * laminar flame speed  S_L     = mean(rho*u) / rho_u
    * thermal thickness    delta_th = (T_b - T_u) / max|dT/dx|
    * burning rate         = rho_u * S_L * (Y_H2,u - Y_H2,b)

Run inside the env that has yt (e.g. `conda activate fpce`):

    python compare_flame.py <plotfile> [reference.dat] [-o out.png]

Example:
    python compare_flame.py CurrentResults/pltFlame_05000
"""
import argparse
import numpy as np

# molar masses [g/mol] for the LiDryer species (used for mole->mass conversion)
MW = {
    "H2": 2.01588, "O2": 31.9988, "H2O": 18.01528, "H": 1.00794,
    "O": 15.9994, "OH": 17.00734, "HO2": 33.00674, "H2O2": 34.01468,
    "N2": 28.0134,
}


# --------------------------------------------------------------------------
# Reference flame (Tecplot POINT format; species columns are MOLE fractions)
# --------------------------------------------------------------------------
def read_reference(path):
    cols, data = None, []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("VARIABLES"):
                cols = [c.strip() for c in s.split("=", 1)[1].split('"')
                        if c.strip() and c.strip() != ","]
                continue
            if s.startswith("ZONE") or not s:
                continue
            try:
                data.append([float(v) for v in s.split()])
            except ValueError:
                continue
    arr = np.array(data)
    d = {name: arr[:, i] for i, name in enumerate(cols)}

    # mean molar mass at every point -> mole to mass fractions
    Mmix = np.zeros_like(d["X"])
    for sp in MW:
        if sp in d:
            Mmix += d[sp] * MW[sp]
    out = {
        "x": d["X"], "T": d["temp"], "u": d["u"], "rho": d["rho"],
        "Y_H2": d["H2"] * MW["H2"] / Mmix,
        "Y_H2O": d["H2O"] * MW["H2O"] / Mmix,
    }
    return out


# --------------------------------------------------------------------------
# CLEM plotfile line-out
# --------------------------------------------------------------------------
def read_plotfile(path):
    import yt
    yt.set_log_level(50)
    ds = yt.load(path)
    c = ds.domain_center
    ray = ds.ortho_ray(0, (c[1], c[2]))
    order = np.argsort(np.array(ray["boxlib", "x"]))

    def fld(*names):
        for n in names:
            try:
                return np.array(ray["boxlib", n][order])
            except Exception:
                continue
        raise KeyError(f"none of {names} present in plotfile")

    return {
        "x": fld("x"),
        "T": fld("Temp"),
        "u": fld("x_velocity"),
        "p":   fld("pressure"),
        "rho": fld("density"),
        "Y_H2": fld("Y(H2)", "y_H2"),
        "Y_H2O": fld("Y(H2O)", "y_H2O"),
        "time": float(ds.current_time),
    }


# --------------------------------------------------------------------------
# metrics
# --------------------------------------------------------------------------
def front_x(p):
    return p["x"][np.argmax(np.abs(np.gradient(p["T"], p["x"])))]


def metrics(p):
    rho_u = p["rho"][0]
    massflux = p["rho"] * p["u"]
    S_L = np.mean(massflux) / rho_u
    dth = (p["T"].max() - p["T"].min()) / np.max(np.abs(np.gradient(p["T"], p["x"])))
    burn = rho_u * S_L * (p["Y_H2"][0] - p["Y_H2"][-1])
    return {
        "S_L": S_L, "delta_th": dth, "burn_rate": burn,
        "massflux_rel_std": np.std(massflux) / np.mean(massflux),
    }


def main():
    ap = argparse.ArgumentParser(description="Compare CLEM 1D flame vs reference.")
    ap.add_argument("plotfile")
    ap.add_argument("reference", nargs="?",
                    default="LiDryer_H2_p1_phi0_4000tu0300.dat")
    ap.add_argument("-o", "--out", default="compare_flame.png")
    args = ap.parse_args()

    sim = read_plotfile(args.plotfile)
    ref = read_reference(args.reference)

    # align reference front onto the simulation front
    shift = front_x(sim) - front_x(ref)
    ref_x = ref["x"] + shift

    ms, mr = metrics(sim), metrics(ref)

    # ---- console report ----
    print(f"\nPlotfile : {args.plotfile}   (t = {sim['time']:.4e} s)")
    print(f"Reference: {args.reference}\n")
    flag = "STEADY" if ms["massflux_rel_std"] < 0.02 else \
        "NOT steady (rho*u varies) - run longer / retune u_in"
    print(f"  mass-flux rel. std = {ms['massflux_rel_std']:.2e}  ({flag})\n")
    print(f"  {'quantity':<26}{'CLEM':>13}{'reference':>13}{'err %':>10}")
    print("  " + "-" * 62)
    for k, lab, unit in (("S_L", "flame speed S_L", "cm/s"),
                         ("delta_th", "thermal thickness", "cm"),
                         ("burn_rate", "burning rate", "g/cm^2/s")):
        e = 100.0 * (ms[k] - mr[k]) / mr[k]
        print(f"  {lab+' ['+unit+']':<26}{ms[k]:>13.5g}{mr[k]:>13.5g}{e:>10.2f}")
    print()
    #print(sim["T"])
    # ---- overlay plots ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    print(f"  rho_u={sim['rho'][0]:.4e}  rho_b={sim['rho'][-1]:.4e}  (ref rho_b≈2.3e-4)")
    print(f"  p:    min={sim['p'].min():.4e}  max={sim['p'].max():.4e}   (1 atm = 1.013e6)")
    print(f"  rho*u: inlet={sim['rho'][0]*sim['u'][0]:.4e}  outlet={sim['rho'][-1]*sim['u'][-1]:.4e}")
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    panels = [
        (ax[0, 0], "T", "Temperature [K]"),
        (ax[0, 1], "u", "x-velocity [cm/s]"),
        (ax[1, 0], "Y_H2", "Y(H2)"),
        (ax[1, 1], "Y_H2O", "Y(H2O)"),
    ]
    for a, key, ylab in panels:
        a.plot(ref_x, ref[key], "k-", lw=2, label="reference")
        a.plot(sim["x"], sim[key], "r--o", ms=3, lw=1.2, label="CLEM")
        a.set_xlabel("x [cm]")
        a.set_ylabel(ylab)
        a.grid(alpha=0.3)
        a.legend()

    # zoom x-window around the flame
    xf = front_x(sim)
    for a, *_ in panels:
        a.set_xlim(xf - 0.15, xf + 0.15)

    fig.suptitle(
        f"CLEM vs reference  |  S_L {ms['S_L']:.2f} vs {mr['S_L']:.2f} cm/s   "
        f"delta_th {ms['delta_th']:.4f} vs {mr['delta_th']:.4f} cm",
        fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(args.out, dpi=130)
    print(f"  figure written to {args.out}\n")


if __name__ == "__main__":
    main()
