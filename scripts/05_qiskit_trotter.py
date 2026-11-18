"""
M7 / M8 prototype: the model as actual quantum circuits, simulated with Qiskit.

Three things, in the order a referee would want them:

  A. Trotter error against step size, log-log, with the SLOPE fitted. The slope
     is the test. An implementation can be accurate and still scale wrongly,
     which means it is wrong for a reason that will bite at larger t.
  B. Particle density measured from finite shots, with shot-noise error bars,
     against the exact curve. This is the first appearance of statistical error.
  C. The crossover: at what step size does Trotter error drop below shot noise?
     That is headline result number one from the proposal, and this is a first
     pass at it with no gate error yet.

Run:  python scripts/05_qiskit_trotter.py
"""
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from schwinger.budget import first_crossing
from schwinger.circuits import (
    circuit_metrics, particle_density_from_counts, run_with_shots, statevector,
    trotter_circuit,
)
from schwinger.diagnostics import shot_noise_sigma
from schwinger.observables import evolve, particle_density
from schwinger.plotstyle import MARKERS, PALETTE, annotate, save, use_report_style
from schwinger.schwinger import bare_vacuum, hamiltonian

import matplotlib.pyplot as plt

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)

N, M, W, J = 4, 0.5, 1.0, 1.0
T = 3.0
SEED = 20260910


def phase_align(a, b):
    return a * np.exp(-1j * np.angle(np.vdot(b, a)))


def main() -> None:
    use_report_style()
    H = hamiltonian(N, M, W, J)
    v0 = bare_vacuum(N)
    exact_state = evolve(v0, H, np.array([T]))[0]
    exact_density = particle_density(exact_state, N)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14.5, 3.9))
    fig.subplots_adjust(wspace=0.32)

    # ---- A. Trotter error scaling -----------------------------------------
    print("=" * 66)
    print("A. TROTTER ERROR SCALING (the slope is the test)")
    print("=" * 66)
    # Error metric: the largest amplitude difference from exp(-iHt)|vac> after the
    # global phase is removed, at t = 3, for dt = 3/10 down to 3/320.
    steps_list = np.array([10, 20, 40, 80, 160, 320])
    print(f"  metric: max |psi_trotter - psi_exact| after phase alignment, t = {T}, "
          f"dt from {T / steps_list[0]:.4f} to {T / steps_list[-1]:.4f}")
    dts = T / steps_list
    for i, order in enumerate((1, 2)):
        errs = []
        for s in steps_list:
            got = phase_align(statevector(trotter_circuit(N, M, W, J, T, int(s), order)),
                              exact_state)
            errs.append(float(np.max(np.abs(got - exact_state))))
        errs = np.array(errs)
        slope, intercept = np.polyfit(np.log(dts), np.log(errs), 1)
        ax1.loglog(dts, errs, color=PALETTE[i], marker=MARKERS[i],
                   label=f"order {order}: slope {slope:.3f}")
        print(f"  order {order}: fitted slope {slope:.4f}   (predicted {float(order)})   "
              f"error at dt={dts[-1]:.4f} is {errs[-1]:.2e}")
    ax1.set_xlabel("Trotter step  dt")
    ax1.set_ylabel("max state error at t = 3")
    ax1.set_title("Trotter error scales as predicted")
    ax1.legend(loc="lower right")
    annotate(ax1, "Slopes fitted, not assumed. First order\n"
                  "goes as dt, second as dt^2. A wrong\n"
                  "implementation can be accurate at one dt\n"
                  "and still scale wrongly.")

    # ---- B. Measured density with shot noise -------------------------------
    print()
    print("=" * 66)
    print("B. PARTICLE DENSITY FROM FINITE SHOTS")
    print("=" * 66)
    times = np.linspace(0.0, 8.0, 33)
    # The exact line is drawn on its own fine grid: through the 33 points above
    # it showed corners the evolution does not have. The sampled points below
    # still use times[::4], so every printed number and seed is unchanged.
    fine = np.linspace(0.0, 8.0, 801)
    exact_curve = np.array([particle_density(s, N)
                            for s in evolve(v0, H, fine)])
    ax2.plot(fine, exact_curve, color="#6b6b6b", linewidth=1.6, label="exact")

    for i, shots in enumerate((256, 4096)):
        sample_t = times[::4]
        est, err = [], []
        for t in sample_t:
            if t == 0.0:
                est.append(0.0); err.append(0.0); continue
            qc = trotter_circuit(N, M, W, J, float(t), 200, order=2)
            e, s = particle_density_from_counts(run_with_shots(qc, shots, seed=SEED), N)
            est.append(e); err.append(s)
        ax2.errorbar(sample_t, est, yerr=err, color=PALETTE[i], marker=MARKERS[i],
                     linestyle="none", capsize=2.5, elinewidth=1.4,
                     label=f"{shots} shots")
        print(f"  {shots:>5} shots: typical uncertainty {np.mean(err[1:]):.4f}")

    ax2.set_ylim(-0.02, 0.58)
    ax2.set_xlabel("time  (units of 1/w)")
    ax2.set_ylabel("particle density")
    ax2.set_title("Measured on the simulator, with error bars")
    ax2.legend(loc="upper left", ncol=3)
    annotate(ax2, "The particle-number operator is diagonal\n"
                  "in Z, so one measurement basis gives it.\n"
                  "Error bars are the standard error over\n"
                  "shots.")

    # ---- C. The crossover --------------------------------------------------
    print()
    print("=" * 66)
    print("C. CROSSOVER: WHERE DOES TROTTER ERROR FALL BELOW SHOT NOISE?")
    print("=" * 66)
    # The noise floor is the EXACT standard error of Eq. (sigma), computed from the
    # exact outcome distribution at t = 3, not estimated from one sampled run (an
    # earlier version did that, and its single-sample sigma was 7 % low at 256
    # shots). The step grid is dense and the crossing is interpolated in log dt,
    # so the quoted crossover is not just the nearest point of a coarse grid.
    steps_c = np.array(list(range(2, 33)) + [40, 48, 64, 96, 128, 192, 256])
    dts_c = T / steps_c
    systematic, twoq = [], []
    for s in steps_c:
        qc = trotter_circuit(N, M, W, J, T, int(s), order=2)
        systematic.append(abs(particle_density(statevector(qc), N) - exact_density))
        twoq.append(circuit_metrics(qc)["two_qubit_gates"])
    systematic = np.array(systematic)

    ax3.loglog(dts_c, np.maximum(systematic, 1e-12), color=PALETTE[2], marker=MARKERS[2],
               markersize=3, label="Trotter error (systematic)")
    exact_probs = np.abs(exact_state) ** 2
    for shots in (256, 1000, 4096):
        sigma = shot_noise_sigma(exact_probs, N, shots)
        cross = first_crossing(dts_c, systematic - sigma)
        k = int(np.argmax(systematic - sigma <= 0))
        print(f"  {shots:>5} shots: exact sigma = {sigma:.5f};  Trotter bias falls below it at "
              f"dt = {cross:.3f} (interpolated between steps {T / dts_c[k - 1]:.0f} and "
              f"{T / dts_c[k]:.0f})")
        # All three shot counts are drawn: the 1000-shot crossover is the one
        # the thesis quotes in the abstract and conclusion.
        colour = PALETTE[{256: 0, 1000: 3, 4096: 1}[shots]]
        ax3.axhline(sigma, color=colour, linewidth=1.4, linestyle=(0, (4, 3)))
        ax3.text(dts_c.min() * 1.15, sigma * 1.2, f"{shots}-shot noise",
                 fontsize=7.5, color=colour)

    ax3.set_xlim(dts_c.min() * 0.55, dts_c.max() * 2.2)
    marks = [int(np.where(steps_c == s)[0][0]) for s in (2, 16, 256)]
    # Offsets keep each label off the curve and off the shot-noise lines: the
    # 512 CX label goes right of the curve, where the plot is empty.
    for k, off in zip(marks, ((10, 4), (14, -7), (8, -9))):
        ax3.annotate(f"{twoq[k]} CX", (dts_c[k], max(systematic[k], 1e-12)),
                     textcoords="offset points", xytext=off, fontsize=7,
                     color="#6b6b6b")
    ax3.set_xlabel("Trotter step  dt")
    ax3.set_ylabel("error on particle density")
    ax3.set_title("Systematic vs statistical, at t = 3")
    ax3.legend(loc="lower right")
    annotate(ax3, "Left of a dashed line, shot noise dominates\n"
                  "and shrinking dt buys nothing. CX counts are\n"
                  "annotated because gate error grows with them.\n"
                  "Gate error is not included in this panel.")

    out = FIG / "05_qiskit_trotter.png"
    save(fig, out)
    print(f"\nwrote {out}")

    print()
    print("=" * 66)
    print("CIRCUIT COST (what a real device would have to run)")
    print("=" * 66)
    print(f"{'steps':>7} {'dt':>8} {'depth':>8} {'2q gates':>10} {'total':>8}")
    for s in (10, 50, 200):
        mtr = circuit_metrics(trotter_circuit(N, M, W, J, T, s, order=2))
        print(f"{s:>7} {T/s:>8.4f} {mtr['depth']:>8} {mtr['two_qubit_gates']:>10} "
              f"{mtr['total_gates']:>8}")


if __name__ == "__main__":
    main()
