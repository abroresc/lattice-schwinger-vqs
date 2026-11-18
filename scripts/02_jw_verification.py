"""
Jordan-Wigner checked two ways, plus a trap in periodic boundaries.

Panel A: the fermionic algebra itself, checked as matrices.
Panel B: the periodic-boundary spectrum, done naively and done correctly. The
naive version is wrong by exactly |2(h - J)| in the ordered phase and right in
the disordered phase, so the error is easy to miss.

Run:  python scripts/02_jw_verification.py
"""
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from schwinger.ising import (
    tfim_hamiltonian, tfim_periodic_spectrum, tfim_periodic_spectrum_naive,
)
from schwinger.jordan_wigner import (
    anticommutation_residuals, hopping_sign_demo, pauli_exclusion_residual,
)
from schwinger.paulis import eigh_sorted
from schwinger.plotstyle import MARKERS, PALETTE, annotate, save, use_report_style

import matplotlib.pyplot as plt

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)


def main() -> None:
    use_report_style()

    print("=" * 66)
    print("FERMIONIC ALGEBRA, checked as matrices on the full Hilbert space")
    print("=" * 66)
    for N in (2, 3, 4, 5, 6):
        r = anticommutation_residuals(N)
        worst = max(r.values())
        print(f"  N = {N}:  worst anticommutator residual = {worst:.3e}"
              f"   c^dag c^dag = {pauli_exclusion_residual(N):.3e}")

    print()
    print("=" * 66)
    print("HOPPING SIGN, checked numerically")
    print("=" * 66)
    for k, v in hopping_sign_demo(4).items():
        print(f"  {k:<40s} residual {v:.3e}")
    print("  The plain Z string gives a minus sign and the alternating (-Z) string")
    print("  removes it. The spectra are identical, so the sign is a relabelling.")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.6))

    sizes = [2, 3, 4, 5, 6, 7]
    worst = [max(anticommutation_residuals(N).values()) for N in sizes]
    ax1.bar([str(N) for N in sizes], np.maximum(worst, 1e-18), color=PALETTE[0], width=0.55)
    ax1.set_yscale("log")
    ax1.set_ylim(1e-18, 1e-10)
    ax1.axhline(1e-12, color="#9a9a9a", linewidth=1.0, linestyle=(0, (4, 3)))
    ax1.text(-0.4, 1.4e-12, "test tolerance", fontsize=7.5, color="#6b6b6b")
    ax1.set_xlabel("number of sites N")
    ax1.set_ylabel("worst  |{c_i, c_j*} - delta_ij|")
    ax1.set_title("The fermion algebra holds exactly")
    annotate(ax1, "Bars sit at or below 1e-18 (plotted floor), i.e. exactly zero. This checks\n"
                  "the transformation directly rather than inferring it from an energy.")

    ratios = np.linspace(0.05, 2.5, 50)
    err_fixed, err_naive = [], []
    N = 6
    for r in ratios:
        ed, _ = eigh_sorted(tfim_hamiltonian(N, 1.0, r, periodic=True))
        ed = np.sort(ed)
        err_fixed.append(np.max(np.abs(ed - tfim_periodic_spectrum(N, 1.0, r))))
        err_naive.append(np.max(np.abs(ed - tfim_periodic_spectrum_naive(N, 1.0, r))))
    print()
    print(f"  periodic TFIM, N = {N}, full spectrum: parity-sector formula residual "
          f"{max(err_fixed):.1e}; naive formula {max(err_naive):.2f}")

    ax2.plot(ratios, np.maximum(err_naive, 1e-17), color=PALETTE[1], marker=MARKERS[1],
             markevery=6, label="naive: sqrt at every k")
    ax2.plot(ratios, np.maximum(err_fixed, 1e-17), color=PALETTE[2], marker=MARKERS[2],
             markevery=6, label="correct: signed at k = 0, pi")
    ax2.plot(ratios[ratios < 1.0], np.abs(2 * (ratios[ratios < 1.0] - 1.0)),
             color="#9a9a9a", linewidth=1.2, linestyle=(0, (4, 3)), label="| 2(h - J) |")
    ax2.set_yscale("log")
    ax2.set_xlabel("h / J")
    ax2.set_ylabel("max |ED  -  closed form|")
    ax2.set_title(f"The periodic-boundary trap, N = {N}")
    ax2.legend(loc="lower left")
    annotate(ax2, "At k = 0 and k = pi the pairing term vanishes, so those modes need the SIGNED\n"
                  "energy 2(h - J cos k), not the positive square root. The naive version is\n"
                  "right for h > J, which is exactly why the bug survives casual testing.")

    out = FIG / "02_jordan_wigner.png"
    save(fig, out)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
