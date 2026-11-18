"""
Physical observables for the lattice Schwinger model.

Every observable here is an explicit Hermitian matrix, so its expectation
value is exact for a given state. diagnostics.py computes the same quantities
from measured outcome distributions, with their shot-noise errors.
"""

from __future__ import annotations

import numpy as np

from .paulis import expectation, identity, site_op
from .schwinger import charge_density_op, electric_field_op, stagger, total_charge_op


def particle_number_op(n_sites: int) -> np.ndarray:
    """
    Total number of particles plus antiparticles, counted from the bare vacuum.

    An even physics site holds a particle when occupied; an odd site holds an
    antiparticle when empty (a hole in the Dirac sea). Adding both:

        N_part = (1/2) sum_n [ I + (-1)^n Z_n ]

    which is zero on the bare vacuum by construction, and rises as pairs are
    created. It is the observable in the pair-production figures.
    """
    dim = 2 ** n_sites
    out = np.zeros((dim, dim), dtype=complex)
    for j in range(n_sites):
        out += 0.5 * (identity(n_sites) + stagger(j) * site_op("Z", j, n_sites))
    return out


def particle_density(state: np.ndarray, n_sites: int) -> float:
    """Particle number per site, in [0, 1]."""
    return expectation(state, particle_number_op(n_sites)) / n_sites


def chiral_condensate_op(n_sites: int) -> np.ndarray:
    """
    The lattice chiral condensate  <psi-bar psi>  as an operator.

        psi-bar psi = (1/N) sum_n (-1)^n c_n^dag c_n
                    = (1/2N) sum_n (-1)^n Z_n          (the constant drops for even N)

    It is an order parameter: it tells you which phase the vacuum is in.

    It also gives a free correctness check: the staggered mass term of H
    is exactly  m * N * (psi-bar psi), so by the Hellmann-Feynman theorem

        <psi-bar psi>  =  (1/N) dE_0/dm

    `analytics.hellmann_feynman_residual` checks the operator expectation
    against a numerical derivative of the ground energy. What that proves
    is narrow: the theorem holds for any Hamiltonian H0 + m V, so
    agreement shows that this operator is exactly (1/N) dH/dm, i.e. that the
    observable and the mass term of H match in sign and normalisation. It says
    nothing about the hopping or electric terms; other checks cover those.
    """
    dim = 2 ** n_sites
    out = np.zeros((dim, dim), dtype=complex)
    for j in range(n_sites):
        out += stagger(j) * site_op("Z", j, n_sites)
    return out / (2.0 * n_sites)


def chiral_condensate(state: np.ndarray, n_sites: int) -> float:
    return expectation(state, chiral_condensate_op(n_sites))


def electric_field_profile(state: np.ndarray, n_sites: int) -> np.ndarray:
    """<L_n> on each of the N-1 links. The flux string, if there is one."""
    return np.array([expectation(state, electric_field_op(j, n_sites)) for j in range(n_sites - 1)])


def charge_profile(state: np.ndarray, n_sites: int) -> np.ndarray:
    """<rho_n> on each site."""
    return np.array([expectation(state, charge_density_op(j, n_sites)) for j in range(n_sites)])


def total_charge(state: np.ndarray, n_sites: int) -> float:
    return expectation(state, total_charge_op(n_sites))


def charge_variance(state: np.ndarray, n_sites: int) -> float:
    """
    Var(Q) = <Q^2> - <Q>^2.

    The main noise diagnostic. docs/derivation.md section 8.2 explains why it
    replaces a Gauss's-law check. In this encoding Gauss's law is solved and substituted away, so
    it holds identically on every state in the Hilbert space and there is no
    Gauss's law violation left to measure. What survives as a symmetry
    constraint is total charge: [H, Q] = 0, so physical evolution stays in one
    charge sector. Noise does not respect that. Var(Q) is therefore a direct,
    model-independent measure of how far noise has pushed the state out of the
    physical sector, and it plays the role the proposal assigned to Gauss's law.
    """
    Q = total_charge_op(n_sites)
    q1 = expectation(state, Q)
    q2 = expectation(state, Q @ Q)
    return q2 - q1 ** 2


def evolve(state: np.ndarray, ham: np.ndarray, times: np.ndarray) -> np.ndarray:
    """
    Exact real-time evolution, |psi(t)> = exp(-i H t) |psi(0)>, via eigenbasis.

    Done through the spectral decomposition rather than scipy.expm, so the cost
    is one diagonalisation regardless of how many time points are requested.
    The Trotter circuits are compared against this.
    """
    from .paulis import eigh_sorted

    vals, vecs = eigh_sorted(ham)
    coeffs = vecs.conj().T @ state
    out = np.empty((len(times), len(state)), dtype=complex)
    for i, t in enumerate(times):
        out[i] = vecs @ (np.exp(-1j * vals * t) * coeffs)
    return out
