"""
Tests for Jordan-Wigner. This is the crux step, so it is tested hardest.

The algebra tests are the ones that matter: they check the transformation
directly rather than inferring correctness from a downstream energy.
"""
import numpy as np
import pytest

from schwinger.jordan_wigner import (
    anticommutation_residuals, bogoliubov_spectrum, hopping_sign_demo,
    jw_annihilation, jw_creation, jw_number, many_body_spectrum_from_modes,
    pauli_exclusion_residual, quadratic_from_AB,
)
from schwinger.paulis import identity, site_op

TOL = 1e-12


@pytest.mark.parametrize("n", [2, 3, 4, 5])
@pytest.mark.parametrize("alternating", [False, True])
def test_fermionic_anticommutation(n, alternating):
    """{c_i, c_j^dag} = delta_ij, {c_i, c_j} = 0. The definition of a fermion."""
    for name, resid in anticommutation_residuals(n, alternating).items():
        assert resid < TOL, f"{name} = {resid:.3e}"


@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_pauli_exclusion_is_algebraic(n):
    """c^dag c^dag = 0 exactly. Exclusion as an identity, not a rule imposed by hand."""
    assert pauli_exclusion_residual(n) < TOL


@pytest.mark.parametrize("n", [2, 4])
def test_number_operator_is_local(n):
    """n_j = c_j^dag c_j = (I + Z_j)/2: the Z string cancels against itself."""
    for j in range(n):
        built = jw_creation(j, n) @ jw_annihilation(j, n)
        assert np.allclose(built, jw_number(j, n))
        assert np.allclose(built, 0.5 * (identity(n) + site_op("Z", j, n)))


def test_number_operator_eigenvalues_are_zero_or_one():
    n = 4
    for j in range(n):
        vals = np.linalg.eigvalsh(jw_number(j, n))
        assert np.allclose(np.sort(np.unique(np.round(vals, 12))), [0.0, 1.0])


def test_hopping_sign_convention():
    """
    The minus sign is real, removable, and does not change the spectrum.
    All three claims checked numerically so nobody has to argue about them.
    """
    d = hopping_sign_demo(4)
    assert d["plain string gives minus sign"] < TOL
    assert d["alternating string gives plus sign"] < TOL
    assert d["spectra agree"] < 1e-10


def test_string_is_actually_nonlocal():
    """c_3 must involve sites 1 and 2. If it does not, the string is missing."""
    n = 4
    c3 = jw_annihilation(2, n)
    local_guess = np.kron(np.kron(identity(2), np.array([[0, 0], [1, 0]], dtype=complex)), identity(1))
    assert not np.allclose(c3, local_guess)
    assert np.max(np.abs(c3)) > 0


def test_quadratic_from_AB_validates_inputs():
    n = 3
    A = np.eye(n)
    B = np.zeros((n, n))
    quadratic_from_AB(A, B, n)
    with pytest.raises(ValueError):
        quadratic_from_AB(np.array([[0.0, 1.0], [2.0, 0.0]]), np.zeros((2, 2)), 2)
    with pytest.raises(ValueError):
        quadratic_from_AB(np.eye(2), np.array([[0.0, 1.0], [1.0, 0.0]]), 2)


def test_bogoliubov_reproduces_free_fermion_chain():
    """
    Sanity case with no pairing: B = 0, so the quasiparticle energies must be
    |eigenvalues of A|, and the many-body spectrum must be all subset sums.
    """
    n = 4
    A = np.zeros((n, n))
    for i in range(n - 1):
        A[i, i + 1] = A[i + 1, i] = 1.0
    lam, offset = bogoliubov_spectrum(A, np.zeros((n, n)))
    assert np.allclose(np.sort(lam), np.sort(np.abs(np.linalg.eigvalsh(A))))
    built = quadratic_from_AB(A, np.zeros((n, n)), n)
    assert np.allclose(
        np.sort(np.linalg.eigvalsh(built)),
        many_body_spectrum_from_modes(lam, offset),
        atol=1e-10,
    )
