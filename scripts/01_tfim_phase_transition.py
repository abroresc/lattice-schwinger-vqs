"""
M2: the transverse-field Ising model.

Two things at once: the physics (a crossover near h/J = 1 that sharpens into a
phase transition as N grows) and the verification (exact diagonalisation
against the closed-form free-fermion answer, at machine precision).

Run:  python scripts/01_tfim_phase_transition.py
"""
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from schwinger.ising import tfim_exact_spectrum, tfim_hamiltonian
from schwinger.paulis import eigh_sorted
from schwinger.plotstyle import MARKERS, PALETTE, annotate, save, use_report_style

import matplotlib.pyplot as plt

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)


def main() -> None:
    use_report_style()
    ratios = np.linspace(0.0, 2.5, 51)
    sizes = [2, 4, 6, 8]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.6))

    residuals = {}
    for i, N in enumerate(sizes):
        e_per_site = np.empty(len(ratios))
        resid = np.empty(len(ratios))
        for k, r in enumerate(ratios):
            H = tfim_hamiltonian(N, 1.0, r)
            ed, _ = eigh_sorted(H)
            e_per_site[k] = ed[0] / N
            resid[k] = np.max(np.abs(np.sort(ed) - tfim_exact_spectrum(N, 1.0, r)))
        residuals[N] = resid
        ax1.plot(ratios, e_per_site, color=PALETTE[i], marker=MARKERS[i],
                 markevery=6, label=f"N = {N}")

    ax1.axvline(1.0, color="#9a9a9a", linewidth=1.0, linestyle=(0, (4, 3)))
    ax1.text(1.04, ax1.get_ylim()[1] - 0.06, "h/J = 1", fontsize=7.5, color="#6b6b6b")
    ax1.set_xlabel("h / J")
    ax1.set_ylabel("ground energy per site")
    ax1.set_title("Ground state softens as the field wins")
    ax1.legend(loc="lower left")
    annotate(ax1, "At N = 4 this is a crossover, not a sharp transition. The kink only\n"
                  "becomes a true non-analyticity as N grows.")

    for i, N in enumerate(sizes):
        ax2.semilogy(ratios, np.maximum(residuals[N], 1e-17), color=PALETTE[i],
                     marker=MARKERS[i], markevery=6, label=f"N = {N}")
    ax2.axhline(1e-12, color="#9a9a9a", linewidth=1.0, linestyle=(0, (4, 3)))
    ax2.text(0.05, 1.4e-12, "test tolerance", fontsize=7.5, color="#6b6b6b")
    ax2.set_xlabel("h / J")
    ax2.set_ylabel("max |ED  -  free fermion|")
    ax2.set_title("Jordan-Wigner verified against the closed form")
    ax2.legend(loc="upper left", ncol=2)
    annotate(ax2, "Every one of the 2^N levels, not just the ground state. Residuals sit at\n"
                  "double precision, so the transformation is right, not approximately right.")

    out = FIG / "01_tfim.png"
    save(fig, out)
    print(f"wrote {out}")
    print(f"worst residual across all N and h/J: {max(r.max() for r in residuals.values()):.3e}")


if __name__ == "__main__":
    main()
