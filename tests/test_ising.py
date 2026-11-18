"""
TFIM tests: Jordan-Wigner verified against a closed-form answer.

This is the warm-up whose entire purpose is that an analytic result exists.
"""
import numpy as np
import pytest

from schwinger.ising import (
    tfim_AB, tfim_exact_spectrum, tfim_hamiltonian, tfim_ground_energy_scan,
    tfim_periodic_spectrum, tfim_periodic_spectrum_naive, tfim_terms,
    verify_tfim_jordan_wigner,
)
from schwinger.paulis import eigh_sorted, from_terms, is_hermitian

PARAMS = [(1.0, 0.3), (1.0, 0.5), (1.0, 1.0), (1.0, 2.0), (0.8, 1.4), (0.0, 1.0), (1.0, 0.0)]


@pytest.mark.parametrize("n", [2, 3, 4, 5])
@pytest.mark.parametrize("J,h", PARAMS)
def test_term_list_matches_matrix(n, J, h):
    assert np.allclose(from_terms(tfim_terms(n, J, h), n), tfim_hamiltonian(n, J, h))


@pytest.mark.parametrize("n", [2, 3, 4, 5])
@pytest.mark.parametrize("J,h", PARAMS)
def test_hermitian(n, J, h):
    assert is_hermitian(tfim_hamiltonian(n, J, h))


@pytest.mark.parametrize("n", [2, 3, 4, 5, 6])
@pytest.mark.parametrize("J,h", PARAMS)
def test_jordan_wigner_rebuilds_spin_hamiltonian(n, J, h):
    """
    The (A, B) matrices derived on paper must rebuild the spin Hamiltonian
    entry for entry. A single wrong sign fails here, localised, immediately.
    """
    assert verify_tfim_jordan_wigner(n, J, h) < 1e-12


@pytest.mark.parametrize("n", [2, 3, 4, 5, 6])
@pytest.mark.parametrize("J,h", PARAMS)
def test_free_fermion_spectrum_matches_exact_diagonalisation(n, J, h):
    """
    The main Jordan-Wigner test. The full 2^N many-body spectrum, rebuilt from N
    single-particle Bogoliubov energies, against brute-force diagonalisation.
    """
    ed, _ = eigh_sorted(tfim_hamiltonian(n, J, h))
    assert np.max(np.abs(np.sort(ed) - tfim_exact_spectrum(n, J, h))) < 1e-10


@pytest.mark.parametrize("n", [4, 6])
@pytest.mark.parametrize("J,h", PARAMS)
def test_periodic_closed_form_with_parity_sectors(n, J, h):
    """Ring geometry: correct treatment of the wrapped string and unpaired modes."""
    ed, _ = eigh_sorted(tfim_hamiltonian(n, J, h, periodic=True))
    assert np.max(np.abs(np.sort(ed) - tfim_periodic_spectrum(n, J, h))) < 1e-10


def test_naive_periodic_form_fails_in_ordered_phase_by_known_amount():
    """
    The trap, pinned as a regression test. Using the positive square root for
    the unpaired k = 0 mode is wrong by exactly |2(h - J)| when h < J, and
    right when h > J. Documented as a test so nobody re-introduces it.
    """
    n, J = 4, 1.0
    for h in (0.3, 0.5):
        ed, _ = eigh_sorted(tfim_hamiltonian(n, J, h, periodic=True))
        err = np.max(np.abs(np.sort(ed) - tfim_periodic_spectrum_naive(n, J, h)))
        assert abs(err - abs(2 * (h - J))) < 1e-9
    for h in (1.5, 2.0):
        ed, _ = eigh_sorted(tfim_hamiltonian(n, J, h, periodic=True))
        assert np.max(np.abs(np.sort(ed) - tfim_periodic_spectrum_naive(n, J, h))) < 1e-10


def test_AB_shapes_and_symmetry():
    A, B, const = tfim_AB(5, 1.0, 0.7)
    assert np.allclose(A, A.T)
    assert np.allclose(B, -B.T)
    assert const == pytest.approx(0.7 * 5)


def test_ground_energy_decreases_with_field():
    """Adding a transverse field can only lower the ground energy (variational)."""
    ratios = np.linspace(0.0, 3.0, 13)
    e = tfim_ground_energy_scan(4, 1.0, ratios)
    assert np.all(np.diff(e) < 1e-12)
