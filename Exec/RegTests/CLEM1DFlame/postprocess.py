#!/usr/bin/env python3
"""
Accuracy diagnostics for the CLEM 1D premixed flame.

Extracts a 1D line-out from a PeleC plotfile and reports the three target
quantities, then compares them against the reference flame solution
(LiDryer_H2_p1_phi0_4000tu0300.dat):

  * Displacement / laminar flame speed  S_L  [cm/s]
        measured as the mass flux through the (stationary) flame divided by the
        unburned density:  S_L = mean(rho*u) / rho_u.  In a steady 1D flame the
        mass flux is constant, so this is robust and frame-independent.
  * Thermal flame thickness  delta_th  [cm]
        delta_th = (T_b - T_u) / max|dT/dx|
  * Burning / production rate per unit area  [g/cm^2/s]
        m_dot * (Y_H2,u - Y_H2,b) = rho_u * S_L * Delta Y_H2
        (= integral of the H2 consumption rate across the flame)

Usage:
    python postprocess.py <plotfile> [reference.dat]

Requires: numpy, yt
"""
import sys
import numpy as np

MW = {"H2": 2.01588, "O2": 31.9988, "H2O": 18.01528, "N2": 28.0134}


def read_reference(path):
    """Parse the Tecplot POINT reference flame. Species columns are mole frac."""
    cols = None
    data = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("VARIABLES"):
                cols = [c.strip().strip('"') for c in s.split("=", 1)[1].split('"') if c.strip() and c.strip() != ","]
                continue
            if s.startswith("ZONE") or not s:
                continue
            try:
                data.append([float(v) for v in s.split()])
            except ValueError:
                continue
    arr = np.array(data)
    d = {name: arr[:, i] for i, name in enumerate(cols)}

    x, T, u, rho = d["X"], d["temp"], d["u"], d["rho"]
    # mole -> mass fraction for H2
    Mmix = d["H2"] * MW["H2"] + d["O2"] * MW["O2"] + d["H2O"] * MW["H2O"] + d["N2"] * MW["N2"]
    YH2 = d["H2"] * MW["H2"] / Mmix

    SL = u[0]                                  # cold-boundary velocity = S_L
    dth = (T.max() - T.min()) / np.max(np.abs(np.gradient(T, x)))
    burn = rho[0] * SL * (YH2[0] - YH2[-1])
    return {"S_L": SL, "delta_th": dth, "burn_rate": burn,
            "T_u": T.min(), "T_b": T.max(), "rho_u": rho[0]}


def analyze_plotfile(path):
    import yt
    yt.set_log_level(50)
    ds = yt.load(path)

    # 1D ray along x at the domain mid-line
    c = ds.domain_center
    ray = ds.ortho_ray(0, (c[1], c[2]))
    order = np.argsort(ray["boxlib", "x"])
    x = np.array(ray["boxlib", "x"][order])
    T = np.array(ray["boxlib", "Temp"][order])
    rho = np.array(ray["boxlib", "density"][order])
    u = np.array(ray["boxlib", "x_velocity"][order])
    try:
        YH2 = np.array(ray["boxlib", "Y(H2)"][order])
    except Exception:
        YH2 = np.array(ray["boxlib", "y_H2"][order])

    rho_u = rho[0]
    massflux = rho * u
    SL = np.mean(massflux) / rho_u
    flux_var = np.std(massflux) / np.mean(massflux)  # should be << 1 when steady

    dth = (T.max() - T.min()) / np.max(np.abs(np.gradient(T, x)))
    burn = rho_u * SL * (YH2[0] - YH2[-1])

    # flame front = location of max |dT/dx|
    x_front = x[np.argmax(np.abs(np.gradient(T, x)))]

    return {"S_L": SL, "delta_th": dth, "burn_rate": burn,
            "x_front": x_front, "massflux_rel_std": flux_var,
            "T_u": T.min(), "T_b": T.max(), "rho_u": rho_u}


def err(a, b):
    return 100.0 * (a - b) / b if b else float("nan")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    plt = sys.argv[1]
    ref_path = sys.argv[2] if len(sys.argv) > 2 else "LiDryer_H2_p1_phi0_4000tu0300.dat"

    sim = analyze_plotfile(plt)
    ref = read_reference(ref_path)

    print(f"\nPlotfile : {plt}")
    print(f"Reference: {ref_path}\n")
    print(f"  flame front x         = {sim['x_front']:.4f} cm")
    print(f"  mass-flux rel. std    = {sim['massflux_rel_std']:.2e}  "
          f"({'STEADY' if sim['massflux_rel_std'] < 0.02 else 'NOT steady - keep running / retune u_in'})")
    print(f"  T_u / T_b             = {sim['T_u']:.1f} / {sim['T_b']:.1f} K"
          f"   (ref {ref['T_u']:.1f} / {ref['T_b']:.1f})\n")

    print(f"  {'quantity':<24}{'CLEM':>14}{'reference':>14}{'error %':>12}")
    print("  " + "-" * 62)
    for key, label, unit in (
        ("S_L", "flame speed S_L", "cm/s"),
        ("delta_th", "thermal thickness", "cm"),
        ("burn_rate", "burning rate", "g/cm^2/s"),
    ):
        print(f"  {label+' ['+unit+']':<24}{sim[key]:>14.5g}{ref[key]:>14.5g}{err(sim[key], ref[key]):>12.2f}")
    print()


if __name__ == "__main__":
    main()
