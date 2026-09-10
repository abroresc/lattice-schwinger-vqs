"""
Tests for the Pauli-string machinery.

Every test compares against something known independently of this code:
a hand-written matrix, an algebraic identity, or Qiskit.
"""
import numpy as np
import pytest

from schwinger.paulis import (
    I2, X, Y, Z, SP, SM,
    anticommutator, basis_index, basis_label, commutator, eigh_sorted,
    from_terms, identity, is_hermitian, pauli_string, site_op,
    site_to_qiskit_qubit, qiskit_bitstring_to_sites, two_site_op,
)


def test_single_site_matrices_by_hand():
    assert np.allclose(X @ X, I2)
    assert np.allclose(Y @ Y, I2)
    assert np.allclose(Z @ Z, I2)
    assert np.allclose(X @ Y, 1j * Z)
    assert np.allclose(Y @ Z, 1j * X)
    assert np.allclose(Z @ X, 1j * Y)
    assert np.allclose(SP + SM, X)
    assert np.allclose(SP - SM, 1j * Y)
    assert np.allclose(SP @ SM, 0.5 * (I2 + Z))


def test_pauli_string_two_site_by_hand():
    # Z on site 1 of two: diag(+1, +1, -1, -1) since site 1 is the leading factor
    assert np.allclose(pauli_string("ZI"), np.diag([1, 1, -1, -1]))
    # Z on site 2 of two: diag(+1, -1, +1, -1)
    assert np.allclose(pauli_string("IZ"), np.diag([1, -1, 1, -1]))
    assert np.allclose(pauli_string("XI"), np.kron(X, I2))


def test_ordering_is_site_one_leftmost():
    """The convention that everything else depends on. If this flips, so does the physics."""
    N = 4
    for j in range(N):
        chars = ["I"] * N
        chars[j] = "Z"
        assert np.allclose(site_op("Z", j, N), pauli_string("".join(chars)))
    # site 0 is the most significant factor: its Z eigenvalue splits the basis in half
    d = np.diag(site_op("Z", 0, N)).real
    assert np.all(d[: 2 ** (N - 1)] == 1.0)
    assert np.all(d[2 ** (N - 1):] == -1.0)


def test_string_validation():
    with pytest.raises(ValueError):
        pauli_string("ZQ")
    with pytest.raises(ValueError):
        pauli_string("ZI", 3)
    with pytest.raises(ValueError):
        site_op("Z", 4, 4)
    with pytest.raises(ValueError):
        two_site_op("Z", 1, "Z", 1, 4)


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_all_pauli_strings_hermitian_and_unitary(n):
    import itertools
    for combo in itertools.product("IXYZ", repeat=n):
        mat = pauli_string("".join(combo))
        assert is_hermitian(mat)
        assert np.allclose(mat @ mat, identity(n))


def test_pauli_algebra_on_chain():
    N = 4
    # Same-site X and Z anticommute; different-site operators commute.
    assert np.allclose(anticommutator(site_op("X", 1, N), site_op("Z", 1, N)), 0)
    assert np.allclose(commutator(site_op("X", 1, N), site_op("Z", 2, N)), 0)


def test_basis_label_roundtrip():
    for n in (2, 4, 6):
        for i in range(2 ** n):
            assert basis_index(basis_label(i, n)) == i
    assert basis_label(5, 4) == "0101"


def test_from_terms_matches_manual_sum():
    N = 3
    terms = [("ZII", 0.5), ("IXI", -1.25), ("IIY", 2.0), ("ZZI", 0.75)]
    manual = sum(c * pauli_string(l) for l, c in terms)
    assert np.allclose(from_terms(terms, N), manual)


def test_eigh_sorted_rejects_non_hermitian():
    bad = np.array([[0, 1], [0, 0]], dtype=complex)
    with pytest.raises(ValueError):
        eigh_sorted(bad)


def test_qiskit_labels_coincide_and_qubit_numbering_reverses():
    """
    The convention trap, pinned down. Labels agree; qubit numbers reverse.
    Skipped if Qiskit is not installed, since Gate 2 does not require it.
    """
    qi = pytest.importorskip("qiskit.quantum_info")
    N = 4
    for label in ("IIIZ", "ZIII", "IXXI", "ZIZI", "IIYY", "XYZI"):
        assert np.allclose(qi.SparsePauliOp(label).to_matrix(), pauli_string(label))
    for j in range(N):
        q = site_to_qiskit_qubit(j, N)
        theirs = qi.SparsePauliOp.from_sparse_list([("Z", [q], 1.0)], num_qubits=N).to_matrix()
        assert np.allclose(site_op("Z", j, N), theirs)
    assert site_to_qiskit_qubit(0, 4) == 3
    assert qiskit_bitstring_to_sites("0011") == "1100"
