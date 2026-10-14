"""
Analytic limits. These are the known answers the numerics get checked against.

The design principle of the whole project is that no stage is ever trusted
because it looks right. Each function here computes a quantity a second,
independent way, so a disagreement localises the bug instead of announcing it.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np


def strong_coupling_ground_energy(n_sites: int, m: float) -> float:
    """
    Exact ground energy at w = 0, any J >= 0, m >= 0.

    With no hopping the Hamiltonian is diagonal in the computational basis.
    The electric term is a sum of squares, so it is >= 0 and vanishes only when
    every L_n = 0, which means every charge density is zero. That is the bare
    vacuum, and it also minimises the mass term, since

        (m/2) sum_n (-1)^n Z_n  >=  -m N / 2

    with equality exactly when Z_n = -(-1)^n. One configuration minimises both
    terms at once, so there is no competition and the answer is exact:

        E_0(w=0) = -m N / 2
    """
    return -0.5 * m * n_sites


def free_fermion_single_particle(n_sites: int, w: float) -> np.ndarray:
    """
    Single-particle energies of the pure hopping chain, m = J = 0.

    The hopping term alone is a tight-binding chain on N sites with open ends,
    whose standing-wave modes give

        eps_k = 2 w cos( k pi / (N+1) ),   k = 1 .. N

    The Jordan-Wigner sign discussed in jordan_wigner.py flips w, and the set
    {cos(k pi / (N+1))} is symmetric about zero, so the SET of energies is
    unchanged. That is why the sign never shows up in a spectrum.
    """
    k = np.arange(1, n_sites + 1)
    return 2.0 * w * np.cos(k * np.pi / (n_sites + 1))


def free_fermion_spectrum(n_sites: int, w: float) -> np.ndarray:
    """Full many-body spectrum at m = J = 0: every subset of single-particle modes."""
    eps = free_fermion_single_particle(n_sites, w)
    energies = []
    for r in range(n_sites + 1):
        for subset in combinations(range(n_sites), r):
            energies.append(float(sum(eps[i] for i in subset)))
    return np.sort(np.array(energies))


def hellmann_feynman_residual(n_sites: int, m: float, w: float, J: float, dm: float = 1e-5) -> float:
    """
    Check <psi-bar psi> two ways: as an operator expectation, and as (1/N) dE_0/dm.

    The Hellmann-Feynman theorem says these must agree exactly. They are
    computed by different routes (one from the ground-state vector, one from a
    finite difference of eigenvalues). What agreement proves is narrow: the
    theorem dE0/dm = <dH/dm> holds for ANY H0 + m V, so it confirms that the
    condensate operator is exactly (1/N) dH/dm, i.e. that the observable and
    the mass term agree in sign and normalisation. A wrong hopping or electric
    term would pass it. The constructions and the limits cover those terms.

    Returns the absolute difference. Measured 4e-12 at m = 0.5w and 3e-11 at
    m = 1.5w with a central difference at dm = 1e-5; anything at 1e-3 or larger is a bug, not truncation error.
    """
    from .observables import chiral_condensate
    from .schwinger import ground_state

    _, psi = ground_state(n_sites, m, w, J)
    from_operator = chiral_condensate(psi, n_sites)

    e_plus, _ = ground_state(n_sites, m + dm, w, J)
    e_minus, _ = ground_state(n_sites, m - dm, w, J)
    from_derivative = (e_plus - e_minus) / (2 * dm) / n_sites

    return abs(from_operator - from_derivative)


def gauss_law_residual(n_sites: int) -> float:
    """
    Confirm that Gauss's law is an IDENTITY in this encoding, not a constraint.

    Checks  L_n - L_{n-1} - rho_n = 0  as a matrix, on the whole Hilbert space,
    with L_0 = 0. It should be exactly zero everywhere, because the gauge field
    was eliminated by substituting the solved constraint rather than by
    imposing it.

    Read that consequence carefully. It means there is no Gauss's law violation
    to measure in this encoding, and the diagnostic in the proposal has to be
    total charge conservation instead. See observables.charge_variance.
    """
    from .schwinger import charge_density_op, electric_field_op
    from .paulis import identity

    worst = 0.0
    dim = 2 ** n_sites
    for j in range(n_sites):
        L_n = electric_field_op(j, n_sites)
        L_prev = electric_field_op(j - 1, n_sites) if j > 0 else np.zeros((dim, dim), dtype=complex)
        resid = L_n - L_prev - charge_density_op(j, n_sites)
        worst = max(worst, float(np.max(np.abs(resid))))
    return worst


def charge_commutator_residual(n_sites: int, m: float, w: float, J: float) -> float:
    """[H, Q] must vanish: total charge is conserved, so H blocks by charge sector."""
    from .paulis import commutator
    from .schwinger import hamiltonian, total_charge_op

    return float(np.max(np.abs(commutator(hamiltonian(n_sites, m, w, J), total_charge_op(n_sites)))))


def cp_symmetry_residual(n_sites: int, m: float, w: float, J: float) -> tuple[float, float]:
    """
    Check that CP commutes with H inside the neutral sector but not outside it.

    Returns (residual_in_neutral_sector, residual_on_full_space). The first
    must vanish; the second must NOT, and reporting both is the point. A
    symmetry that holds only on a subspace is easy to state incorrectly, so
    the code states it by measuring both sides.
    """
    from .schwinger import charge_sectors, cp_operator, hamiltonian

    ham = hamiltonian(n_sites, m, w, J)
    cp = cp_operator(n_sites)
    full = float(np.max(np.abs(cp @ ham - ham @ cp)))

    idx = np.where(charge_sectors(n_sites) == 0)[0]
    proj = np.zeros((2 ** n_sites, len(idx)), dtype=complex)
    for a, i in enumerate(idx):
        proj[i, a] = 1.0
    h_s = proj.conj().T @ ham @ proj
    cp_s = proj.conj().T @ cp @ proj
    sector = float(np.max(np.abs(cp_s @ h_s - h_s @ cp_s)))
    return sector, full


def link_symmetry_residual(state: np.ndarray, n_sites: int) -> float:
    """
    On a neutral state that CP leaves unchanged,  <L_n> = <L_{N-n}>.

    Neutrality alone is not enough: |0110> is neutral, not CP-symmetric, and
    has links (0, 0, -1). The states this project evolves are CP-symmetric
    (they start from the bare vacuum, and H commutes with CP at Q = 0).
    Returns the worst violation. Nonzero on such a state means it left the
    neutral sector, lost the symmetry, or the field operators are wrong.
    """
    from .observables import electric_field_profile

    prof = electric_field_profile(state, n_sites)
    return float(np.max(np.abs(prof - prof[::-1])))
