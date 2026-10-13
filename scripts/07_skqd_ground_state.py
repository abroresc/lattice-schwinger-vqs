"""
M6 / M7: the ground state by Sample-based Krylov Quantum Diagonalization.

  A. Validation: full-sector SKQD against exact diagonalisation; the
     variational bound on every possible subspace; the exact finite-shot
     expectation against brute-force Monte Carlo.
  B. Which bitstrings matter: ground-state weight of each neutral bitstring.
  C. The SKQD energy and condensate bias against shots per Krylov step,
     computed EXACTLY (no sampling) for noiseless and noisy circuits, for the
     project's full Trotter step and for Rosanowski et al.'s kinetic step.

Run:  python scripts/07_skqd_ground_state.py
"""
import itertools
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from schwinger.noise import gate_noise_model, noisy_probabilities
from schwinger.observables import chiral_condensate
from schwinger.paulis import basis_label
from schwinger.plotstyle import MARKERS, MUTED, PALETTE, save, use_report_style
from schwinger.schwinger import ground_state, hamiltonian
from schwinger.skqd import (
    expected_skqd_energy, krylov_circuits, run_skqd, sector_indices, subspace_ground_state,
)

import matplotlib.pyplot as plt

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
N, M, W, J = 4, 0.5, 1.0, 1.0
SHOTS = np.array([2, 3, 5, 8, 12, 20, 30, 50, 80, 120, 200, 300, 500])
CONFIGS = [  # label, dt, Krylov steps, evolution, eps2
    ("full step, dt=0.3, K=4", 0.3, 4, "full", 0.0),
    ("full step, dt=0.3, K=4, eps2=3e-3", 0.3, 4, "full", 3e-3),
    ("kinetic step, dt=0.77, K=1", 0.77, 1, "kinetic", 0.0),
]


def main():
    ham = hamiltonian(N, M, W, J)
    e0, gs = ground_state(N, M, W, J)
    basis = sector_indices(N)

    print("=" * 70)
    print("A. SKQD VALIDATION (residuals)")
    print("=" * 70)
    e_full, v_full = subspace_ground_state(basis, ham)
    print(f"  full neutral sector vs exact diagonalisation, energy    : {abs(e_full - e0):.1e}")
    print(f"  full neutral sector vs exact, 1 - |overlap|              : {1 - abs(np.vdot(gs, v_full)):.1e}")
    margins = [subspace_ground_state(s, ham)[0] - e0
               for r in range(1, len(basis) + 1) for s in itertools.combinations(basis, r)]
    print(f"  smallest E(S) - E0 over all 63 subspaces (must be >= 0)  : {min(margins):.1e}")
    dists = [noisy_probabilities(c) for c in krylov_circuits(N, M, W, J, 0.3, 4)]
    ex = expected_skqd_energy(dists, 15, N, M, W, J)
    rng = np.random.default_rng(20261012)
    runs = np.array([run_skqd(dists, 15, rng, N, M, W, J)[-1].energy for _ in range(20000)])
    print(f"  exact expected bias at 15 shots                          : {ex.mean_bias:.4e}")
    print(f"  Monte Carlo (20000 runs)                                 : "
          f"{runs.mean() - e0:.4e} +/- {runs.std() / np.sqrt(runs.size):.1e}")

    print()
    print("=" * 70)
    print("B. GROUND-STATE WEIGHT OF EACH NEUTRAL BITSTRING (m/w = 0.5, J = w)")
    print("=" * 70)
    for b in sorted(basis, key=lambda i: -abs(gs[i]) ** 2):
        print(f"  |{basis_label(b, N)}>   {abs(gs[b]) ** 2:.5f}")

    print()
    print("=" * 70)
    print("C. EXPECTED BIAS AGAINST SHOTS PER KRYLOV STEP (exact, no sampling)")
    print("=" * 70)
    curves = {}
    for label, dt, k, evo, eps in CONFIGS:
        nm = gate_noise_model(eps) if eps > 0 else None
        dd = [noisy_probabilities(c, nm) for c in krylov_circuits(N, M, W, J, dt, k, evolution=evo)]
        bias, pfull, cond = [], [], []
        for s in SHOTS:
            ex = expected_skqd_energy(dd, int(s), N, M, W, J)
            bias.append(ex.mean_bias)
            pfull.append(ex.prob_full_support)
            # expected condensate error, weighting each subspace by its probability
            c = 0.0
            for mask, prob in enumerate(ex.subset_probs):
                if mask and prob > 0:
                    idx = ex.basis[[(mask >> i) & 1 == 1 for i in range(ex.basis.size)]]
                    c += prob * (chiral_condensate(subspace_ground_state(idx, ham)[1], N)
                                 - chiral_condensate(gs, N))
            cond.append(c)
        curves[label] = (np.array(bias), np.array(pfull), np.array(cond))
        ok = np.array(bias) / abs(e0) < 1e-3
        need = SHOTS[np.argmax(ok)] if np.any(ok) else None
        print(f"  {label:<36}: bias at 5/50/500 shots = "
              f"{bias[2]:.2e} / {bias[7]:.2e} / {bias[12]:.2e};  shots for 1e-3 relative: {need}")
        print(f"  {'':<36}  condensate bias at 5/50 shots = {cond[2]:+.2e} / {cond[7]:+.2e}")

    use_report_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.0, 3.7))
    fig.subplots_adjust(wspace=0.3)
    for i, (label, *_rest) in enumerate(CONFIGS):
        bias, pfull, cond = curves[label]
        ax1.loglog(SHOTS, np.maximum(bias, 1e-16), color=PALETTE[i], marker=MARKERS[i], label=label)
        ax2.semilogx(SHOTS, np.maximum(1 - pfull, 1e-16), color=PALETTE[i], marker=MARKERS[i],
                     label=label)
    rng = np.random.default_rng(7)
    mc = [np.mean([run_skqd(dists, int(s), rng, N, M, W, J)[-1].energy for _ in range(4000)]) - e0
          for s in SHOTS[:8]]
    ax1.loglog(SHOTS[:8], np.maximum(mc, 1e-16), color=MUTED, marker="x", linestyle="none",
               label="Monte Carlo check (4000 runs)")
    ax1.axhline(1e-3 * abs(e0), color=MUTED, linewidth=0.8, linestyle=(0, (4, 3)))
    ax1.set_ylim(1e-9, 1)
    ax1.set_xlabel("shots per Krylov step")
    ax1.set_ylabel(r"expected $E_{\rm SKQD} - E_0$")
    ax1.set_title("(a) SKQD energy bias, computed exactly")
    ax1.legend(loc="upper right", fontsize=7.5)
    ax2.set_xlabel("shots per Krylov step")
    ax2.set_ylabel("P(a ground-state bitstring is missed)")
    ax2.set_title("(b) Coverage of the six neutral bitstrings")
    ax2.set_yscale("log")
    ax2.set_ylim(1e-6, 1.2)
    out = FIG / "07_skqd.png"
    save(fig, out)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
