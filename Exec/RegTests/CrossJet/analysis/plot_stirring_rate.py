"""Triplet-map (stirring) activity per timestep, WALE+wall-model vs Smagorinsky.

Reads the CSVs written by stirring_from_log.py. Both runs are identical in the
CLEM block (n_lem = 12, do_stir = 1, Sc_t = 1) and differ only in the LES
closure, so every difference below is the closure acting through D_T = nu_t.

    /home/pc08/miniconda3/envs/felics/bin/python plot_stirring_rate.py
"""
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

A = "/home/pc08/Fernando/CLEMNew/PeleC-CLEM/Exec/RegTests/CrossJet/analysis"
CASES = [("CLEM + Smagorinsky", A + "/stir_smag.csv", "#b2182b"),
         ("CLEM + WALE + wall model", A + "/stir_walewm.csv", "#1b7837")]
INK, MUTED = "#1a1a1a", "#6b6b6b"
WIN = 201                      # rolling-mean window, in steps
MS = 1e3                       # s -> ms on the time axis


def roll(y, n=WIN):
    """Centred rolling mean with the ends blanked. mode='same' zero-pads, which
    drags the first and last n/2 points toward zero and draws a cliff at the
    end of each run that looks like the stirring shutting off."""
    if len(y) < n:
        return y
    out = np.convolve(y, np.ones(n) / n, mode="same")
    out[:n // 2] = np.nan
    out[-(n // 2):] = np.nan
    return out


def main():
    fig, ax = plt.subplots(3, 1, figsize=(9.2, 9.0), sharex=True,
                           constrained_layout=True)
    print(f"{'case':26s} {'steps':>7s} {'maps/step':>10s} {'maps/s':>11s} "
          f"{'%lines/step':>12s} {'rej%':>7s} {'dt [ns]':>9s}")
    for lab, path, c in CASES:
        d = np.genfromtxt(path, delimiter=",", names=True)
        t = d["time"] * MS
        ax[0].plot(t, d["maps"], color=c, lw=0.4, alpha=0.18)
        ax[0].plot(t, roll(d["maps"]), color=c, lw=1.7, label=lab)
        ax[1].plot(t, roll(d["maps_per_sec"]), color=c, lw=1.7, label=lab)
        ax[2].plot(t, roll(d["rejected_pct"]), color=c, lw=1.7, label=lab)
        print(f"{lab:26s} {len(t):7d} {d['maps'].mean():10.2f} "
              f"{d['maps_per_sec'].mean():11.3e} "
              f"{100 * d['maps'].mean() / d['lines'][0]:12.4f} "
              f"{d['rejected_pct'].mean():7.2f} {1e9 * d['dt'].mean():9.3f}")

    ax[0].set_ylabel("triplet maps applied per step")
    ax[0].set_yscale("symlog", linthresh=1.0)
    ax[0].set_title("Stirring activity per timestep: the LES closure sets "
                    "$D_T=\\nu_t$, and $\\nu_t$ sets the eddy event rate",
                    fontsize=10.5, color=INK)
    ax[1].set_ylabel("triplet maps per second\nof physical time")
    ax[1].set_yscale("log")
    ax[1].set_title("dt-normalised: removes the two runs' different step sizes",
                    fontsize=9.5, color=MUTED)
    ax[2].set_ylabel("eddies rejected as\ntoo small [%]")
    ax[2].set_ylim(98.5, 100.05)
    ax[2].set_title("Both runs are starved by $n_{lem}=12$: >99% of eddies "
                    "drawn from $f(\\ell)$ cannot be mapped",
                    fontsize=9.5, color=MUTED)
    ax[2].set_xlabel("time [ms]")
    for a in ax:
        a.grid(alpha=0.25, lw=0.6)
        a.legend(frameon=False, fontsize=9)
        a.tick_params(colors=INK)
    for ext in ("png", "pdf"):
        fig.savefig(f"{A}/stirring_per_step.{ext}", dpi=180)
        print(f"wrote {A}/stirring_per_step.{ext}")


if __name__ == "__main__":
    main()
