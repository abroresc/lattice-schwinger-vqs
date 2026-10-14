"""
M5 / Gate 2: the N = 4 Schwinger Hamiltonian, built, diagonalised, and checked.

Prints the Pauli-string table (the Gate 2 deliverable) and every verification
residual, then draws the spectrum and the chiral condensate against the
fermion mass.

Run:  python scripts/03_schwinger_spectrum.py
"""
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from schwinger.analytics import (
    charge_commutator_residual, free_fermion_spectrum, gauss_law_residual,
    hellmann_feynman_residual, strong_coupling_ground_energy,
)
from schwinger.observables import chiral_condensate
from schwinger.plotstyle import MARKERS, PALETTE, annotate, save, use_report_style
from schwinger.schwinger import (
    charge_sectors, construction_residual, fermion_construction_residual,
    format_terms, ground_state, hamiltonian_terms, spectrum,
)

import matplotlib.pyplot as plt

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)
N = 4


def symbolic_table(n: int) -> str:
    """Print each Pauli string with its coefficient in terms of m, w, J."""
    def coeffs(m, w, J):
        return {k: v.real for k, v in hamiltonian_terms(n, m, w, J)}
    cm, cw, cJ = coeffs(1, 0, 0), coeffs(0, 1, 0), coeffs(0, 0, 1)
    keys = sorted(set(list(cm) + list(cw) + list(cJ)),
                  key=lambda s: (sum(c != "I" for c in s), s))
    lines = [f"{'string':<10} {'coefficient':<26} origin", "-" * 60]
    for k in keys:
        a, b, c = cm.get(k, 0.0), cw.get(k, 0.0), cJ.get(k, 0.0)
        parts, origin = [], []
        for val, sym, name in ((a, "m", "mass"), (b, "w", "hopping"), (c, "J", "electric")):
            if abs(val) > 1e-12:
                parts.append(f"{val:+g} {sym}")
                origin.append(name)
        lines.append(f"{k:<10} {' '.join(parts):<26} {', '.join(origin)}")
    lines.append("-" * 60)
    lines.append(f"{len(keys)} terms on {n} qubits, Hilbert space dimension {2**n}")
    return "\n".join(lines)


def main() -> None:
    use_report_style()

    print("=" * 60)
    print(f"GATE 2 ARTIFACT: N = {N} SCHWINGER HAMILTONIAN")
    print("=" * 60)
    print(symbolic_table(N))
    print()
    print("Numerical values at m = w = J = 1:")
    print(format_terms(hamiltonian_terms(N, 1.0, 1.0, 1.0)))

    print()
    print("=" * 60)
    print("VERIFICATION (residuals, not pass/fail flags)")
    print("=" * 60)
    probes = [(1.0, 1.0, 1.0), (0.3, 1.0, 0.7), (2.0, 0.5, 1.5)]
    print(f"  Pauli expansion vs direct construction : "
          f"{max(construction_residual(N, *p) for p in probes):.3e}")
    print(f"  fermion construction vs spin form      : "
          f"{max(fermion_construction_residual(N, *p) for p in probes):.3e}")
    print(f"  Gauss's law is an identity             : {gauss_law_residual(N):.3e}")
    print(f"  [H, Q] = 0                             : "
          f"{max(charge_commutator_residual(N, *p) for p in probes):.3e}")
    sc = max(abs(spectrum(N, m, 0.0, J)[0] - strong_coupling_ground_energy(N, m))
             for m in (0.5, 1.0, 2.0) for J in (0.0, 1.0, 3.0))
    print(f"  strong coupling E0 = -mN/2 (w = 0)     : {sc:.3e}")
    ff = max(np.max(np.abs(spectrum(N, 0.0, w, 0.0) - free_fermion_spectrum(N, w)))
             for w in (1.0, 0.6))
    print(f"  free fermion limit, all 16 levels      : {ff:.3e}")
    print(f"  Hellmann-Feynman on the condensate     : "
          f"{max(hellmann_feynman_residual(N, m, 1.0, 1.0) for m in (0.5, 1.5)):.3e}")
    print()
    print(f"  charge sectors present: {sorted(set(charge_sectors(N).tolist()))}")

    masses = np.linspace(-2.0, 3.0, 121)
    levels = np.array([spectrum(N, m, 1.0, 1.0) for m in masses])
    cond = np.array([chiral_condensate(ground_state(N, m, 1.0, 1.0)[1], N) for m in masses])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.6))

    for k in range(levels.shape[1]):
        ax1.plot(masses, levels[:, k], color="#c2c2c2", linewidth=1.0)
    ax1.plot(masses, levels[:, 0], color=PALETTE[0], linewidth=2.2, label="ground state")
    ax1.plot(masses, levels[:, 1], color=PALETTE[1], linewidth=1.8, label="first excited")
    ax1.plot(masses, -0.5 * masses * N, color="#6b6b6b", linewidth=1.2,
             linestyle=(0, (4, 3)), label="strong coupling, -mN/2")
    ax1.set_xlabel("fermion mass  m / w")
    ax1.set_ylabel("energy")
    ax1.set_title(f"All {2**N} levels, N = {N}, w = J = 1")
    ax1.legend(loc="upper right")
    annotate(ax1, "Grey lines are the remaining levels. The ground state approaches the exact\n"
                  "strong-coupling line at large m, where the vacuum locks into the Dirac sea.")

    ax2.plot(masses, cond, color=PALETTE[2], marker=MARKERS[2], markevery=10,
             label="from the ground state")
    dm = 1e-4
    deriv = np.array([(ground_state(N, m + dm, 1.0, 1.0)[0]
                       - ground_state(N, m - dm, 1.0, 1.0)[0]) / (2 * dm) / N for m in masses])
    ax2.plot(masses, deriv, color=PALETTE[3], linewidth=1.2, linestyle=(0, (4, 3)),
             label="from  (1/N) dE0/dm")
    ax2.axhline(-0.5, color="#9a9a9a", linewidth=1.0)
    ax2.text(-1.95, -0.475, "bare Dirac sea, -1/2", fontsize=7.5, color="#6b6b6b")
    ax2.set_xlabel("fermion mass  m / w")
    ax2.set_ylabel("chiral condensate")
    ax2.set_title("Order parameter, computed two independent ways")
    ax2.legend(loc="upper right")
    annotate(ax2, "Hellmann-Feynman: the operator expectation and the derivative of the ground\n"
                  "energy must agree. They do, to 3e-11, the size of the finite-difference error.")

    out = FIG / "03_schwinger_spectrum.png"
    save(fig, out)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
