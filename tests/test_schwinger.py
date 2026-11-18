"""
Tests of the N = 4 Schwinger Hamiltonian against every known answer.

The Hamiltonian is built three independent ways and each analytic limit is
checked separately, so a disagreement points at which stage is wrong rather
than announcing that something, somewhere, is.
"""
import numpy as np
import pytest

from schwinger.analytics import (
    charge_commutator_residual, free_fermion_spectrum, gauss_law_residual,
    hellmann_feynman_residual, strong_coupling_ground_energy,
)
from schwinger.observables import (
    charge_profile, charge_variance, chiral_condensate, electric_field_profile,
    evolve, particle_density, particle_number_op, total_charge,
)
from schwinger.paulis import eigh_sorted, is_hermitian
from schwinger.schwinger import (
    bare_vacuum, bare_vacuum_label, charge_sectors, construction_residual,
    electric_field_op, fermion_construction_residual, ground_state,
    hamiltonian, hamiltonian_direct, hamiltonian_from_fermions,
    hamiltonian_terms, spectrum, stagger, total_charge_op,
)

TOL = 1e-10
SIZES = [2, 3, 4, 5, 6]
COUPLINGS = [(1.0, 1.0, 1.0), (0.3, 1.0, 0.7), (2.0, 0.5, 1.5), (0.0, 1.0, 1.0), (1.0, 1.0, 0.0)]


# --------------------------------------------------------------------------
# Construction: three routes to the same operator
# --------------------------------------------------------------------------

@pytest.mark.parametrize("n", SIZES)
@pytest.mark.parametrize("m,w,J", COUPLINGS)
def test_pauli_expansion_matches_direct_construction(n, m, w, J):
    """The hand expansion of sum_n L_n^2 into Pauli strings, checked by matrix."""
    assert construction_residual(n, m, w, J) < TOL


@pytest.mark.parametrize("n", SIZES)
@pytest.mark.parametrize("m,w,J", COUPLINGS)
def test_fermion_construction_matches_spin_construction(n, m, w, J):
    """
    The central test. Building from lattice fermions and routing through
    Jordan-Wigner must give the same operator as the printed spin form.
    This is the step the whole derivation exists to justify.
    """
    assert fermion_construction_residual(n, m, w, J) < TOL


@pytest.mark.parametrize("n", SIZES)
@pytest.mark.parametrize("m,w,J", COUPLINGS)
def test_hermitian_and_real_spectrum(n, m, w, J):
    h = hamiltonian(n, m, w, J)
    assert is_hermitian(h)
    assert np.all(np.isreal(np.linalg.eigvalsh(h)))


def test_n2_hamiltonian_by_hand():
    """
    N = 2, m = w = J = 1, worked out by hand and hardcoded.

    Sites: n = 1 (odd, stagger -1), n = 2 (even, stagger +1). One link, L_1.
      mass     = (1/2)[ -Z_1 + Z_2 ]
      hopping  = (1/2)( X_1X_2 + Y_1Y_2 )
      electric = L_1^2 with L_1 = (Z_1 - 1)/2, so L_1^2 = (Z_1^2 - 2Z_1 + 1)/4
                = (2 - 2 Z_1)/4 = (I - Z_1)/2
    Total Z_1 coefficient: -1/2 - 1/2 = -1. Z_2: +1/2. Identity: +1/2.
    """
    terms = dict((l, c.real) for l, c in hamiltonian_terms(2, 1.0, 1.0, 1.0))
    assert terms["ZI"] == pytest.approx(-1.0)
    assert terms["IZ"] == pytest.approx(0.5)
    assert terms["II"] == pytest.approx(0.5)
    assert terms["XX"] == pytest.approx(0.5)
    assert terms["YY"] == pytest.approx(0.5)
    assert len(terms) == 5


def test_term_list_has_no_duplicates_and_no_zeros():
    terms = hamiltonian_terms(4, 1.0, 1.0, 1.0)
    labels = [l for l, _ in terms]
    assert len(labels) == len(set(labels))
    assert all(abs(c) > 1e-14 for _, c in terms)


def test_long_range_zz_coefficients_match_derived_formula():
    """
    Expanding sum_n L_n^2 gives Z_m Z_m' coefficient (J/2)(N - m') for m < m',
    depending only on the right site: sites m and m' appear together in the
    N - m' link fields L_n with n >= m'.
    """
    N, J = 5, 1.0
    terms = dict((l, c.real) for l, c in hamiltonian_terms(N, 0.0, 0.0, J))
    for a in range(1, N + 1):
        for b in range(a + 1, N + 1):
            chars = ["I"] * N
            chars[a - 1] = "Z"
            chars[b - 1] = "Z"
            label = "".join(chars)
            expected = 0.5 * J * (N - b) if b <= N - 1 else 0.0
            assert terms.get(label, 0.0) == pytest.approx(expected)


# --------------------------------------------------------------------------
# Analytic limits
# --------------------------------------------------------------------------

@pytest.mark.parametrize("n", SIZES)
@pytest.mark.parametrize("m", [0.5, 1.0, 2.0])
@pytest.mark.parametrize("J", [0.0, 1.0, 3.0])
def test_strong_coupling_limit(n, m, J):
    """w = 0: H is diagonal, the bare vacuum is the exact ground state, E0 = -mN/2."""
    e0, psi = ground_state(n, m, 0.0, J)
    assert e0 == pytest.approx(strong_coupling_ground_energy(n, m), abs=TOL)
    overlap = abs(np.vdot(bare_vacuum(n), psi))
    assert overlap == pytest.approx(1.0, abs=1e-8)


@pytest.mark.parametrize("n", SIZES)
@pytest.mark.parametrize("w", [1.0, 0.6])
def test_free_fermion_limit(n, w):
    """m = J = 0: every one of the 2^N levels is a subset sum of 2w cos(k pi/(N+1))."""
    assert np.max(np.abs(spectrum(n, 0.0, w, 0.0) - free_fermion_spectrum(n, w))) < TOL


@pytest.mark.parametrize("n", SIZES)
def test_gauss_law_is_an_identity(n):
    """
    L_n - L_{n-1} - rho_n = 0 on the WHOLE Hilbert space, not just physical states.

    Consequence, and it matters for the proposal: there is no Gauss's law
    violation to measure in this encoding. The symmetry diagnostic has to be
    total charge instead. See test_charge_is_conserved.
    """
    assert gauss_law_residual(n) < TOL


@pytest.mark.parametrize("n", SIZES)
@pytest.mark.parametrize("m,w,J", COUPLINGS)
def test_charge_is_conserved(n, m, w, J):
    """[H, Q] = 0. This is the symmetry that noise will break, so it is the diagnostic."""
    assert charge_commutator_residual(n, m, w, J) < TOL


@pytest.mark.parametrize("n", [2, 4, 6])
@pytest.mark.parametrize("m,w,J", [(0.5, 1.0, 1.0), (1.5, 1.0, 0.4)])
def test_hellmann_feynman_on_condensate(n, m, w, J):
    """<psi-bar psi> from the state must equal (1/N) dE0/dm from the spectrum."""
    assert hellmann_feynman_residual(n, m, w, J) < 1e-7


# --------------------------------------------------------------------------
# States and observables
# --------------------------------------------------------------------------

def test_bare_vacuum_is_neutral_and_field_free():
    n = 4
    v = bare_vacuum(n)
    assert bare_vacuum_label(4) == "0101"
    assert np.allclose(charge_profile(v, n), 0.0, atol=TOL)
    assert np.allclose(electric_field_profile(v, n), 0.0, atol=TOL)
    assert particle_density(v, n) == pytest.approx(0.0, abs=TOL)
    assert total_charge(v, n) == pytest.approx(0.0, abs=TOL)
    assert charge_variance(v, n) == pytest.approx(0.0, abs=TOL)


def test_particle_number_counts_pairs_not_occupations():
    """
    Flipping one site of the bare vacuum creates one particle (or antiparticle),
    and moves the total charge by one unit. Flipping a neighbouring pair creates
    two and leaves the charge alone: that is a particle-antiparticle pair.
    """
    n = 4
    Np = particle_number_op(n)
    from schwinger.paulis import basis_index
    one_flip = np.zeros(2 ** n, dtype=complex)
    one_flip[basis_index("1101")] = 1.0
    assert np.vdot(one_flip, Np @ one_flip).real == pytest.approx(1.0)
    assert total_charge(one_flip, n) == pytest.approx(-1.0)

    pair = np.zeros(2 ** n, dtype=complex)
    pair[basis_index("0011")] = 1.0
    assert np.vdot(pair, Np @ pair).real == pytest.approx(2.0)
    assert total_charge(pair, n) == pytest.approx(0.0)


def test_particle_number_is_zero_to_N_and_integer_valued():
    n = 4
    vals = np.linalg.eigvalsh(particle_number_op(n))
    assert np.allclose(vals, np.round(vals), atol=TOL)
    assert vals.min() == pytest.approx(0.0)
    assert vals.max() == pytest.approx(float(n))


def test_charge_sectors_partition_the_basis():
    n = 4
    sectors = charge_sectors(n)
    assert len(sectors) == 2 ** n
    assert sectors.min() <= 0 <= sectors.max()
    # H must be block diagonal: no matrix element between different sectors.
    h = hamiltonian(n, 0.7, 1.0, 1.3)
    for i in range(2 ** n):
        for j in range(2 ** n):
            if sectors[i] != sectors[j]:
                assert abs(h[i, j]) < TOL


def test_ground_state_lies_in_zero_charge_sector():
    n = 4
    _, psi = ground_state(n, 0.7, 1.0, 1.3)
    assert total_charge(psi, n) == pytest.approx(0.0, abs=1e-9)
    assert charge_variance(psi, n) == pytest.approx(0.0, abs=1e-9)


def test_variational_principle_holds_for_bare_vacuum():
    """<bare vacuum|H|bare vacuum> can never be below E0. A free correctness check."""
    n = 4
    for m, w, J in COUPLINGS:
        e0, _ = ground_state(n, m, w, J)
        v = bare_vacuum(n)
        expval = np.vdot(v, hamiltonian(n, m, w, J) @ v).real
        assert expval >= e0 - 1e-12


def test_condensate_goes_to_minus_half_at_large_mass():
    """
    As m grows the vacuum locks into the bare Dirac sea, where
    <psi-bar psi> = (1/2N) sum_n (-1)^n Z_n with Z_n = -(-1)^n, giving -1/2.
    """
    n = 4
    _, psi = ground_state(n, 50.0, 1.0, 1.0)
    assert chiral_condensate(psi, n) == pytest.approx(-0.5, abs=1e-3)


# --------------------------------------------------------------------------
# Real-time evolution
# --------------------------------------------------------------------------

def test_evolution_is_unitary_and_conserves_energy_and_charge():
    n = 4
    h = hamiltonian(n, 0.5, 1.0, 1.0)
    times = np.linspace(0.0, 5.0, 21)
    states = evolve(bare_vacuum(n), h, times)
    e0 = np.vdot(states[0], h @ states[0]).real
    for s in states:
        assert np.vdot(s, s).real == pytest.approx(1.0, abs=1e-10)
        assert np.vdot(s, h @ s).real == pytest.approx(e0, abs=1e-9)
        assert total_charge(s, n) == pytest.approx(0.0, abs=1e-9)
        assert charge_variance(s, n) == pytest.approx(0.0, abs=1e-9)


def test_pairs_are_actually_produced():
    """
    The physics claim of the whole project: start in the bare vacuum, let it
    evolve, and particle density becomes nonzero. Nothing is injected; the
    energy was already there because the bare vacuum is not the ground state.
    """
    n = 4
    h = hamiltonian(n, 0.5, 1.0, 1.0)
    times = np.linspace(0.0, 10.0, 201)
    states = evolve(bare_vacuum(n), h, times)
    dens = np.array([particle_density(s, n) for s in states])
    assert dens[0] == pytest.approx(0.0, abs=1e-12)
    assert dens.max() > 0.05


def test_time_evolution_matches_scipy_expm():
    """Independent check of the eigenbasis propagator against a matrix exponential."""
    from scipy.linalg import expm
    n = 4
    h = hamiltonian(n, 0.5, 1.0, 1.0)
    v = bare_vacuum(n)
    for t in (0.3, 1.7, 4.2):
        ours = evolve(v, h, np.array([t]))[0]
        theirs = expm(-1j * h * t) @ v
        assert np.max(np.abs(ours - theirs)) < 1e-11


# --------------------------------------------------------------------------
# CP symmetry: holds in the neutral sector, fails outside it
# --------------------------------------------------------------------------

@pytest.mark.parametrize("n", [4, 6])
@pytest.mark.parametrize("m,w,J", [(0.7, 1.0, 1.3), (1.5, 0.6, 2.0), (0.0, 1.0, 1.0)])
def test_cp_is_a_symmetry_only_in_the_neutral_sector(n, m, w, J):
    from schwinger.analytics import cp_symmetry_residual
    sector, full = cp_symmetry_residual(n, m, w, J)
    assert sector < TOL, "CP must commute with H inside Q = 0"
    assert full > 1.0, "CP must NOT commute with H on the full space; check the claim"


@pytest.mark.parametrize("n", [4, 6])
def test_cp_operator_properties(n):
    from schwinger.schwinger import cp_operator, site_reversal_op
    cp = cp_operator(n)
    assert np.allclose(cp.conj().T @ cp, np.eye(2 ** n))
    assert np.allclose(cp @ cp, np.eye(2 ** n))
    r = site_reversal_op(n)
    assert np.allclose(r @ r, np.eye(2 ** n))


@pytest.mark.parametrize("n", [4, 6])
def test_cp_maps_densities_and_links_as_derived(n):
    """
    The two operator identities behind the symmetry argument, checked as
    matrices: CP rho_n CP = -rho_{N+1-n}, and CP L_n CP = L_{N-n} - Q. The
    second is why CP holds only at Q = 0. Reflection alone gives neither.
    """
    from schwinger.schwinger import charge_density_op, cp_operator, site_reversal_op
    cp, r, q = cp_operator(n), site_reversal_op(n), total_charge_op(n)
    for j in range(n):
        rho, mirror = charge_density_op(j, n), charge_density_op(n - 1 - j, n)
        assert np.max(np.abs(cp @ rho @ cp + mirror)) < TOL
        assert np.max(np.abs(r @ rho @ r - mirror)) > 0.5        # reflection alone fails
    for j in range(n - 1):                                         # link n = j + 1
        mapped = cp @ electric_field_op(j, n) @ cp
        assert np.max(np.abs(mapped - (electric_field_op(n - 2 - j, n) - q))) < TOL


def test_link_equality_needs_cp_symmetry_not_just_neutrality():
    """
    A counterexample to "CP forces <L_n> = <L_{N-n}> on any neutral state". |0110> is neutral but not CP-symmetric (CP
    maps it to |1001>), and its links are (0, 0, -1): the outer links differ.
    """
    from schwinger.paulis import basis_index
    from schwinger.schwinger import cp_operator
    n = 4
    s = np.zeros(16, dtype=complex)
    s[basis_index("0110")] = 1.0
    assert total_charge(s, n) == pytest.approx(0.0, abs=TOL)
    assert np.allclose(electric_field_profile(s, n), [0.0, 0.0, -1.0], atol=TOL)
    assert abs((cp_operator(n) @ s)[basis_index("1001")]) == pytest.approx(1.0, abs=TOL)


def test_outer_links_carry_equal_flux_on_neutral_states():
    """
    On a neutral state that CP leaves unchanged, <L_n> = <L_{N-n}>. The bare
    vacuum is such a state and H commutes with CP in the neutral sector, so
    every state on the trajectory is too, and so is the ground state. At N = 4
    the two outer links must agree at every time. A free correctness check.
    """
    from schwinger.analytics import link_symmetry_residual
    n = 4
    h = hamiltonian(n, 0.5, 1.0, 1.0)
    states = evolve(bare_vacuum(n), h, np.linspace(0.0, 20.0, 81))
    for s in states:
        assert link_symmetry_residual(s, n) < 1e-12
    _, psi = ground_state(n, 0.7, 1.0, 1.3)
    assert link_symmetry_residual(psi, n) < 1e-12
