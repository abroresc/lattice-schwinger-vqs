"""
Replication: does this code compute what Martinez et al. (2016) computed?

The Hamiltonian is checked term by term elsewhere. Here the check is on the
physics: at their parameters (N = 4, J = w, m = 0.5 w, from the bare vacuum),
the code must reproduce the three theory curves of their Fig. 2(a), read off
the published figure by pixel analysis (schwinger/replication.py says how).

  A. Their pulse order against ours: the same nu at every step.
  B. Their three curves against this code, point by point, with residuals.
  C. What the match pins down: their step length, and their error model.

Run:  python scripts/11_martinez_replication.py
"""
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from schwinger import replication as R
from schwinger.plotstyle import MARKERS, MUTED, PALETTE, annotate, save, use_report_style

import matplotlib.pyplot as plt

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"


def main():
    ks = np.arange(1, 5)
    times = ks * R.STEP
    exact = R.exact_density(times)
    trot = R.trotter_density()
    model = R.error_model_density()
    ideal = np.array(R.FIG2A_IDEAL)

    print("=" * 70)
    print("A. GATE ORDER: their pulse sequence against this code's first-order step")
    print("=" * 70)
    worst = np.max(np.abs(R.martinez_pulse_order_density() - trot))
    print(f"  largest difference in nu over four steps of pi/8: {worst:.1e}")

    print()
    print("=" * 70)
    print("B. MARTINEZ ET AL. FIG. 2(a) AGAINST THIS CODE (reading uncertainty 0.002)")
    print("=" * 70)
    print("  wt       ideal (exact)    discretization    error model     measured")
    print("           this code        theirs   this     theirs  this    (theirs)")
    for k, t in enumerate(times):
        print(f"  {k + 1}pi/8   {exact[k]:.3f}            {R.FIG2A_TROTTER[k]:.3f}   {trot[k]:.3f}    "
              f"{R.FIG2A_MODEL[k]:.3f}   {model[k]:.3f}   {R.FIG2A_DATA[k]:.3f}")
    r_ideal = np.max(np.abs(R.exact_density(ideal[:, 0]) - ideal[:, 1]))
    r_trot = np.max(np.abs(trot - np.array(R.FIG2A_TROTTER)))
    r_model = np.max(np.abs(model - np.array(R.FIG2A_MODEL)))
    print(f"  largest residual: ideal curve (33 points) {r_ideal:.4f};  discretization {r_trot:.4f};"
          f"  error model {r_model:.4f}")
    fine = np.linspace(0, np.pi / 2, 3001)
    nu_fine = R.exact_density(fine)
    k = int(np.argmax(nu_fine))
    print(f"  exact peak: nu = {nu_fine[k]:.4f} at wt = {fine[k]:.3f}  (their green curve peaks at 0.464)")

    print()
    print("=" * 70)
    print("C. WHAT THE MATCH PINS DOWN")
    print("=" * 70)
    for f in (0.9, 1.0, 1.1):
        r = np.max(np.abs(R.trotter_density(step=f * R.STEP) - np.array(R.FIG2A_TROTTER)))
        print(f"  step {f:.1f} x pi/8: largest residual against their discretization curve {r:.3f}")
    r0 = np.max(np.abs(R.error_model_density(p=0.0) - np.array(R.FIG2A_MODEL)))
    print(f"  error model without their dephasing (p = 0): residual {r0:.3f}; with p = 0.038: {r_model:.4f}")
    gap = np.array(R.FIG2A_DATA) - model
    print("  their measured data minus their own model: "
          + ", ".join(f"{g:+.3f}" for g in gap) + "  (data sit above the model late on)")

    use_report_style()
    fig, ax = plt.subplots(figsize=(6.4, 3.9))
    ax.plot(fine, nu_fine, color=MUTED, linewidth=1.8, label="this code: exact")
    ax.plot(ideal[:, 0], ideal[:, 1], linestyle="none", marker=".", markersize=4, color="black",
            label="Martinez et al.: ideal (read off Fig. 2a)")
    t0 = np.concatenate([[0.0], times])
    ax.plot(t0, np.concatenate([[0.0], trot]), color=PALETTE[1], marker=MARKERS[1],
            label=r"this code: first order, $w\delta t = \pi/8$")
    ax.plot(times, R.FIG2A_TROTTER, linestyle="none", marker="^", markersize=9, markerfacecolor="none",
            markeredgecolor=PALETTE[1], markeredgewidth=1.4, label="Martinez et al.: discretization")
    ax.plot(t0, np.concatenate([[0.0], model]), color=PALETTE[3], marker=MARKERS[0],
            label="this code: + their dephasing, p = 0.038")
    ax.plot(times, R.FIG2A_MODEL, linestyle="none", marker="s", markersize=9, markerfacecolor="none",
            markeredgecolor=PALETTE[3], markeredgewidth=1.4, label="Martinez et al.: error model")
    ax.plot(times, R.FIG2A_DATA, linestyle="none", marker="o", markersize=5, color=PALETTE[2],
            label="Martinez et al.: measured")
    ax.set_xticks([0, np.pi / 8, np.pi / 4, 3 * np.pi / 8, np.pi / 2])
    ax.set_xticklabels(["0", r"$\pi/8$", r"$\pi/4$", r"$3\pi/8$", r"$\pi/2$"])
    ax.set_xlabel(r"time $wt$")
    ax.set_ylabel(r"particle density $\nu$")
    ax.set_ylim(-0.01, 0.62)
    ax.set_title("Martinez et al. (2016) Fig. 2(a), recomputed")
    ax.legend(loc="upper left", fontsize=7, ncol=2)
    annotate(ax, f"Same Hamiltonian, step and error model, computed here: largest residual "
                 f"{max(r_ideal, r_trot, r_model):.3f},\nwithin the 0.002 to which their figure can be read. "
                 "Their measured points are not fitted.")
    out = FIG / "11_martinez_replication.png"
    save(fig, out)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
