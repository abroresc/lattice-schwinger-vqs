"""
Pauli-string machinery: build any N-qubit operator as an explicit matrix.

Everything downstream (the Ising warm-up, Jordan-Wigner, the Schwinger
Hamiltonian) is a sum of Pauli strings built with `pauli_string`.

Conventions (see CONVENTIONS.md; change both together or neither):

  - Sites are labelled 1..N in the physics text, 0..N-1 in code.
  - A string "XZIY" means X on site 1, Z on site 2, I on site 3, Y on site 4.
  - The matrix is built as  op(site 1) kron op(site 2) kron ... kron op(site N).
    So site 1 is the most significant factor, as in Martinez et al. (2016).
  - Basis index b of a computational basis state is
        b = sum_j s_j * 2^(N-1-j),  s_j in {0,1}, j = 0..N-1
    i.e. the bit string read left to right is site 1, site 2, ..., site N.
  - |0> is the sigma^z = +1 eigenstate ("spin up"). Under the Jordan-Wigner
    convention in jordan_wigner.py this means |0> is an occupied site.
  - Qiskit's Pauli labels coincide with ours (checked numerically, see
    `site_to_qiskit_qubit`), but its qubit numbering is reversed: our site j
    is Qiskit qubit N-1-j. Convert every circuit index with
    `site_to_qiskit_qubit`; never pass a site index to a Qiskit gate.
"""

from __future__ import annotations

from functools import reduce
from typing import Iterable, Sequence

import numpy as np

# Single-site matrices in the {|0>, |1>} basis.
I2 = np.array([[1, 0], [0, 1]], dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

# sigma^+ = |0><1| raises sigma^z, which creates an occupation (convention above)
SP = np.array([[0, 1], [0, 0]], dtype=complex)
# sigma^- = |1><0|
SM = np.array([[0, 0], [1, 0]], dtype=complex)

_SINGLE = {"I": I2, "X": X, "Y": Y, "Z": Z, "+": SP, "-": SM}


def kron_all(mats: Iterable[np.ndarray]) -> np.ndarray:
    """Kronecker product of a sequence, left to right."""
    return reduce(np.kron, mats)


def pauli_string(ops: str, n_sites: int | None = None) -> np.ndarray:
    """
    Build the full 2^N x 2^N matrix for a Pauli string.

    Parameters
    ----------
    ops
        One character per site, from {I, X, Y, Z, +, -}. Character 0 is site 1.
        '+' and '-' are the ladder operators, included because the Schwinger
        hopping term is naturally written with them.
    n_sites
        If given, `ops` must have exactly this length. Used as a guard so a
        typo in a string shows up as an exception rather than a wrong-size
        matrix that silently broadcasts.

    Examples
    --------
    >>> np.allclose(pauli_string("ZI"), np.diag([1, 1, -1, -1]))
    True
    """
    if n_sites is not None and len(ops) != n_sites:
        raise ValueError(f"string {ops!r} has length {len(ops)}, expected {n_sites}")
    bad = set(ops) - set(_SINGLE)
    if bad:
        raise ValueError(f"unknown operator character(s) {sorted(bad)} in {ops!r}")
    return kron_all(_SINGLE[c] for c in ops)


def site_op(op: str, site: int, n_sites: int) -> np.ndarray:
    """
    Single-site operator `op` acting on `site` (0-indexed), identity elsewhere.

    site_op("Z", 0, 4) is Z on site 1 of a 4-site chain.
    """
    if not 0 <= site < n_sites:
        raise ValueError(f"site {site} out of range for {n_sites} sites")
    chars = ["I"] * n_sites
    chars[site] = op
    return pauli_string("".join(chars))


def two_site_op(op_a: str, site_a: int, op_b: str, site_b: int, n_sites: int) -> np.ndarray:
    """Product of two single-site operators on distinct sites."""
    if site_a == site_b:
        raise ValueError("sites must be distinct; multiply single-site ops directly instead")
    for s in (site_a, site_b):
        if not 0 <= s < n_sites:
            raise ValueError(f"site {s} out of range for {n_sites} sites")
    chars = ["I"] * n_sites
    chars[site_a] = op_a
    chars[site_b] = op_b
    return pauli_string("".join(chars))


def from_terms(terms: Sequence[tuple[str, complex]], n_sites: int | None = None) -> np.ndarray:
    """
    Build a Hamiltonian from a list of (pauli_string, coefficient) pairs.

    Every Hamiltonian in the project is stored this way: the list can be
    printed, checked by hand, and passed to Qiskit unchanged.
    """
    if not terms:
        raise ValueError("empty term list")
    n = n_sites if n_sites is not None else len(terms[0][0])
    dim = 2 ** n
    out = np.zeros((dim, dim), dtype=complex)
    for label, coeff in terms:
        out += coeff * pauli_string(label, n)
    return out


def identity(n_sites: int) -> np.ndarray:
    return np.eye(2 ** n_sites, dtype=complex)


def basis_label(index: int, n_sites: int) -> str:
    """Bit string for a basis index, site 1 leftmost. Inverse of `basis_index`."""
    return format(index, f"0{n_sites}b")


def basis_index(label: str) -> int:
    """Basis index for a bit string, site 1 leftmost."""
    return int(label, 2)


def site_to_qiskit_qubit(site: int, n_sites: int) -> int:
    """
    Translate our site index to a Qiskit qubit number.

    The two libraries agree on one thing and differ on another, and the
    difference was found by a numerical check, not from the documentation.

    What agrees: the Pauli label string. Qiskit's Pauli and
    SparsePauliOp read a label left to right as most significant tensor factor
    first, exactly as `pauli_string` does. So

        SparsePauliOp("ZIZI").to_matrix()  ==  pauli_string("ZIZI")

    entry for entry, with no reversal. Reversing the string is wrong: it
    produces a Hamiltonian reflected through the middle of the chain, which
    for the Schwinger model swaps particles and antiparticles.

    What differs: the qubit number. Qiskit calls the rightmost character
    qubit 0, so our site j (0-indexed from the left) is Qiskit qubit
    n_sites - 1 - j. That matters the moment you build a circuit, because
    qc.cx(0, 1) acts on the right end of our chain, not the left. It also
    matters when reading measurement bitstrings back, which Qiskit prints in
    the reverse of our site order.

    Use this function for every circuit index. Never pass a site index to a
    Qiskit gate directly.
    """
    if not 0 <= site < n_sites:
        raise ValueError(f"site {site} out of range for {n_sites} sites")
    return n_sites - 1 - site


def qiskit_bitstring_to_sites(bits: str) -> str:
    """
    Convert a Qiskit measurement bitstring to our site ordering.

    Qiskit prints qubit 0 last, which, combined with the numbering above,
    means its bitstring runs from our site N-1 down to site 0. Reverse it.
    """
    return bits[::-1]


def is_hermitian(mat: np.ndarray, tol: float = 1e-12) -> bool:
    return bool(np.max(np.abs(mat - mat.conj().T)) < tol)


def commutator(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return a @ b - b @ a


def anticommutator(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return a @ b + b @ a


def expectation(state: np.ndarray, op: np.ndarray) -> float:
    """<psi|op|psi> for a normalised state, returned real (asserts small imaginary part)."""
    val = state.conj() @ (op @ state)
    if abs(val.imag) > 1e-10:
        raise ValueError(f"expectation value has imaginary part {val.imag:.3e}; operator not Hermitian?")
    return float(val.real)


def eigh_sorted(mat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Hermitian eigendecomposition with ascending eigenvalues. Returns (evals, evecs)."""
    if not is_hermitian(mat):
        raise ValueError("matrix is not Hermitian to 1e-12")
    vals, vecs = np.linalg.eigh(mat)
    order = np.argsort(vals)
    return vals[order], vecs[:, order]
