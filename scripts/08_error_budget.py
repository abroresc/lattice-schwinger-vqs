"""
The error budget: Trotter, gate, readout and shot noise, separated, against step size, with the optimum step size it implies.

  A. The crossover: Trotter bias against shot noise, no hardware error.
  B. Total root-mean-square error of an M-shot measurement of nu(t_trough)
     against dt, for several two-qubit error rates. The minimum is the optimum
     step size, and it moves to larger dt as the hardware gets worse.
  C. The decomposition at the optimum for three observables: the particle
     density and the flux on the outer and middle links.

All biases are exact (density-matrix simulation); shot noise is the exact
standard error of the measured distribution. No sampling anywhere.

Run:  python scripts/08_error_budget.py
"""
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from schwinger.budget import budget, first_crossing, point_circuit
from schwinger.diagnostics import shot_noise_sigma
from schwinger.noise import gate_noise_model
from schwinger.plotstyle import MARKERS, MUTED, PALETTE, save, use_report_style

import matplotlib.pyplot as plt

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
N, M, W, J = 4, 0.5, 1.0, 1.0
T = 1.76                      # first minimum of nu(t): where the signal lives
SHOTS = 1000
READOUT = 0.01
STEPS = np.array([2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 48, 64])
EPS = (0.0, 3e-4, 1e-3, 3e-3)


def main():
    circuits = {int(s): point_circuit(N, M, W, J, T, int(s)) for s in STEPS}
    dts = T / STEPS

    print("=" * 72)
    print(f"A. CROSSOVER at t = {T}: Trotter bias vs shot noise, no hardware error")
    print("=" * 72)
    # Exact shot noise of the exact distribution, Eq. (sigma), and the crossing
    # interpolated in log dt (dts runs from large to small step).
    tb = np.array([abs(budget(circuits[int(s)], N, None, 0.0, 1000).trotter_bias) for s in STEPS])
    for shots in (256, 1000, 4096):
        sig = shot_noise_sigma(circuits[int(STEPS[0])].exact, N, shots)
        cross = first_crossing(dts, tb - sig)
        print(f"  M = {shots:>5}: exact sigma = {sig:.4f}; Trotter bias falls below it at "
              f"dt = {cross:.3f} (interpolated)")

    print()
    print("=" * 72)
    print(f"B. TOTAL RMS ERROR of nu(t={T}) with {SHOTS} shots, readout q = {READOUT}")
    print("=" * 72)
    rows = {}
    for eps in EPS:
        nm = gate_noise_model(eps) if eps > 0 else None
        bs = [budget(circuits[int(s)], N, nm, READOUT, SHOTS) for s in STEPS]
        rms = np.array([b.rms for b in bs])
        k = int(np.argmin(rms))
        rows[eps] = (bs, rms, k)
        b = bs[k]
        print(f"  eps2 = {eps:<7g}: optimum dt = {dts[k]:.3f} ({circuits[int(STEPS[k])].two_qubit_gates:>4} CX)"
              f"  rms = {rms[k]:.4f}  [Trotter {b.trotter_bias:+.4f}, gate {b.gate_bias:+.4f},"
              f" readout {b.readout_bias:+.4f}, shot {b.sigma:.4f}]")

    print()
    print("=" * 72)
    print("C. DECOMPOSITION AT THE OPTIMUM, eps2 = 1e-3 (shares of the summed squares)")
    print("=" * 72)
    k = rows[1e-3][2]
    pc = circuits[int(STEPS[k])]
    shares = {}
    for obs, label in (("particle_density", "particle density"), ("link1", "outer link field"),
                       ("link2", "middle link field")):
        b = budget(pc, N, gate_noise_model(1e-3), READOUT, SHOTS, observable=obs)
        sq = np.array([b.trotter_bias**2, b.gate_bias**2, b.readout_bias**2, b.sigma**2])
        shares[label] = sq / sq.sum()
        print(f"  {label:<18}: total rms {b.rms:.4f}; shares Trotter/gate/readout/shot = "
              + " / ".join(f"{100 * x:.0f}%" for x in shares[label]))

    use_report_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 3.9), gridspec_kw={"width_ratios": [1.4, 1]})
    fig.subplots_adjust(wspace=0.3)
    for i, eps in enumerate(EPS):
        bs, rms, k = rows[eps]
        ax1.loglog(dts, rms, color=PALETTE[i], marker=MARKERS[i],
                   label=rf"$\epsilon_2$ = {eps:g}" if eps else r"$\epsilon_2$ = 0")
        ax1.plot(dts[k], rms[k], marker="*", markersize=13, color=PALETTE[i], markeredgecolor="k")
    tb = np.array([abs(rows[0.0][0][j].trotter_bias) for j in range(len(STEPS))])
    ax1.loglog(dts, tb, color=MUTED, linewidth=1.0, linestyle=(0, (4, 3)), label="Trotter bias alone")
    ax1.axhline(rows[0.0][0][-1].sigma, color=MUTED, linewidth=1.0, linestyle=(0, (1, 2)))
    ax1.text(0.2, rows[0.0][0][-1].sigma * 0.72, f"shot noise alone, {SHOTS} shots",
             fontsize=7.5, color=MUTED)
    ax1.set_xlabel(r"Trotter step $\delta t$")
    ax1.set_ylabel(r"rms error of $\nu(t_{\rm t})$")
    ax1.set_title("(a) Total error of one measurement; stars mark the optimum")
    ax1.legend(loc="lower right", fontsize=7.5)

    labels = list(shares)
    bottom = np.zeros(len(labels))
    names = ("Trotter", "gate", "readout", "shot noise")
    for j, name in enumerate(names):
        vals = np.array([shares[l][j] for l in labels]) * 100
        ax2.bar(labels, vals, bottom=bottom, color=PALETTE[j], label=name, width=0.55)
        bottom += vals
    ax2.set_ylabel("share of summed squared error (%)")
    ax2.set_title(r"(b) Where the error comes from, $\epsilon_2 = 10^{-3}$")
    ax2.legend(loc="upper center", ncol=4, fontsize=7.5, bbox_to_anchor=(0.5, -0.1))
    ax2.set_ylim(0, 100)
    ax2.grid(axis="x", visible=False)
    out = FIG / "08_error_budget.png"
    save(fig, out)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
