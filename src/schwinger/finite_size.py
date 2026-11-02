"""
The same Hamiltonian on longer chains, N = 4 to 12, by exact diagonalisation.

Why this module exists
----------------------
The thesis runs everything at N = 4. A fair question is whether what it
measures, the pair-production rise and recombination dip of nu(t) and the
ground state, is physics of the lattice Schwinger model or an artefact of four
sites. Classical exact diagonalisation answers it on chains three times longer.
Nothing here is a quantum simulation and nothing changes the N = 4 results;
it only checks that the four-site model behaves like a longer one.

Why the neutral sector
----------------------
H commutes with the total charge Q, and the bare vacuum is neutral, so the
evolution never leaves the Q = 0 block: the bitstrings with exactly N/2
occupied sites, binom(N, N/2) of them (924 at N = 12, against 4096 in all).
Diagonalising that block alone is exact, not an approximation, and it is what
makes N = 12 quick. The trap is to forget that the block must be taken from a
matrix that really has no elements between sectors; `neutral_block` checks it.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np


def neutral_indices(n_sites: int) -> np.ndarray:
    """Basis indices of the Q = 0 sector, read off the diagonal of the charge operator."""
    from .schwinger import total_charge_op

    q = np.real(np.diag(total_charge_op(n_sites)))
    return np.where(np.abs(q) < 1e-9)[0]


@lru_cache(maxsize=8)
def _block(n_sites: int, m: float, w: float, J: float):
    from .schwinger import hamiltonian

    h = hamiltonian(n_sites, m, w, J)
    idx = neutral_indices(n_sites)
    rest = np.setdiff1d(np.arange(h.shape[0]), idx)
    leak = float(np.abs(h[np.ix_(idx, rest)]).max()) if rest.size else 0.0
    if leak > 1e-12:
        raise RuntimeError(f"H couples the neutral sector to others ({leak:.1e}); it should not")
    energies, vecs = np.linalg.eigh(h[np.ix_(idx, idx)])
    return idx, energies, vecs


def neutral_block(n_sites: int, m: float, w: float, J: float):
    """(indices, eigenvalues, eigenvectors) of H restricted to Q = 0."""
    return _block(n_sites, float(m), float(w), float(J))


def _embed(n_sites: int, idx: np.ndarray, v: np.ndarray) -> np.ndarray:
    full = np.zeros(2 ** n_sites, dtype=complex)
    full[idx] = v
    return full


def quench_density(n_sites: int, m: float, w: float, J: float, times) -> np.ndarray:
    """nu(t) after a quench from the bare vacuum, exactly, inside the neutral block."""
    from .observables import particle_number_op
    from .schwinger import bare_vacuum

    idx, energies, vecs = neutral_block(n_sites, m, w, J)
    vac = bare_vacuum(n_sites)[idx].astype(complex)
    coeff = vecs.conj().T @ vac
    counts = np.real(np.diag(particle_number_op(n_sites)))[idx]
    out = []
    for t in np.asarray(times, dtype=float):
        psi = vecs @ (np.exp(-1j * energies * t) * coeff)
        out.append(float(np.abs(psi) ** 2 @ counts) / n_sites)
    return np.array(out)


def ground_state(n_sites: int, m: float, w: float, J: float) -> tuple[float, float]:
    """Ground-state energy and chiral condensate per site, from the neutral block."""
    from .observables import chiral_condensate

    idx, energies, vecs = neutral_block(n_sites, m, w, J)
    return float(energies[0]), chiral_condensate(_embed(n_sites, idx, vecs[:, 0]), n_sites)


def first_peak_and_trough(times, nu) -> tuple[float, float, float, float]:
    """(t_peak, nu_peak, t_trough, nu_trough): the first local maximum and the minimum after it."""
    t, y = np.asarray(times), np.asarray(nu)
    k = next(i for i in range(1, len(y) - 1) if y[i] >= y[i - 1] and y[i] > y[i + 1])
    j = next(i for i in range(k + 1, len(y) - 1) if y[i] <= y[i - 1] and y[i] < y[i + 1])
    return float(t[k]), float(y[k]), float(t[j]), float(y[j])
