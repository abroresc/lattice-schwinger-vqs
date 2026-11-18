"""
Ground-state preparation by sample-based Krylov quantum diagonalization (SKQD).

What SKQD does
--------------
Following Yu et al. (arXiv:2501.09702) and its Schwinger-model application by
Rosanowski et al. (arXiv:2510.26951, Section III A):

  1. Prepare a reference state |psi_0> on the quantum computer.
  2. Evolve it: |psi_k> = U(dt)^k |psi_0>, for k = 1 .. K.
  3. Measure each |psi_k> in the computational basis, many shots each, and keep
     the set of distinct bitstrings seen. Amplitudes are thrown away.
  4. Discard every bitstring with nonzero total charge (classical postselection
     onto the physical sector; it costs nothing because Q is diagonal).
  5. On a classical computer, write H as a matrix in the span of the kept
     bitstrings and diagonalise it. The lowest eigenvalue is the estimate.

Why the answer can be trusted, and where it cannot
--------------------------------------------------
Step 5 is an exact diagonalisation inside a subspace, so the estimate obeys the
variational bound E_SKQD >= E_0 for any set of bitstrings, however noisy the
device was. Noise can only waste shots or add bitstrings; it cannot push the
energy below the truth. That is the noise resilience Rosanowski et al. describe.

The weakness is coverage. A bitstring that is never sampled cannot be in the
subspace, so the estimate is biased upward whenever an important bitstring is
missed (Miura, arXiv:2603.13536). At N = 4 and m/w = 0.5 the ground state puts
71 % of its weight on the bare vacuum |0101> and only 6.5e-4 on |1100>, the
configuration with the longest flux string. Confinement suppresses exactly the
bitstrings SKQD is most likely to miss.

The exact expectation
---------------------
Because the neutral sector at N = 4 has only six bitstrings, the probability of
every possible sampled set can be computed exactly from the outcome
distributions, with no sampling at all. For a set T of sector bitstrings,

    P(sampled set is inside T) = prod_k ( 1 - sum_{b in sector, b not in T} p_k(b) )^M

since every one of the M shots at every step must avoid the bitstrings outside
T (shots that land on charged bitstrings are discarded and so are harmless).
Inclusion-exclusion (a Moebius inversion over subsets) turns these into
P(sampled set = S), and the expected SKQD energy is sum_S P(S) E_0(S). The
Monte Carlo sampler must agree with this number; tests/test_skqd.py checks it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from qiskit import QuantumCircuit

from .circuits import _append_hopping, bare_vacuum_circuit, trotter_step
from .diagnostics import z_basis_table
from .paulis import basis_index, eigh_sorted
from .schwinger import bare_vacuum_label, hamiltonian


# ----------------------------------------------------------------------------
# The classical half: diagonalise in a subspace of bitstrings
# ----------------------------------------------------------------------------


def sector_indices(n_sites: int, charge: int = 0) -> np.ndarray:
    """Basis indices of every bitstring with total charge `charge`."""
    return np.where(np.isclose(z_basis_table(n_sites).charge, charge))[0]


def subspace_ground_state(indices, ham: np.ndarray) -> tuple[float, np.ndarray]:
    """
    Lowest eigenpair of H restricted to span{|b> : b in indices}.

    Returns the energy and the eigenvector embedded back in the full 2^N space
    (zeros outside the subspace), so every observable function in the package
    can be evaluated on it directly.
    """
    idx = np.asarray(sorted(set(int(i) for i in indices)))
    if idx.size == 0:
        raise ValueError("empty subspace")
    vals, vecs = eigh_sorted(ham[np.ix_(idx, idx)])
    full = np.zeros(ham.shape[0], dtype=complex)
    full[idx] = vecs[:, 0]
    return float(vals[0]), full


# ----------------------------------------------------------------------------
# The quantum half: Krylov circuits
# ----------------------------------------------------------------------------


def krylov_circuits(n_sites: int, m: float, w: float, J: float, dt: float, n_krylov: int,
                    evolution: str = "full", order: int = 2) -> list[QuantumCircuit]:
    """
    Circuits preparing U(dt)^k |bare vacuum> for k = 1 .. n_krylov.

    evolution = "full"    : U(dt) is one Trotter step of the full Hamiltonian,
                            the same circuits.trotter_step the dynamics use. This
                            is the reuse that decided SKQD over VQE.
    evolution = "kinetic" : U(dt) contains only the hopping terms, as in
                            Rosanowski et al. The diagonal terms only change
                            phases, which sampling cannot see at the first step,
                            so the kinetic circuit is shallower and samples a
                            similar set of bitstrings.

    The reference state is the bare vacuum, the ground state of the mass and
    electric terms, which is also Rosanowski et al.'s "ground state of the mass
    Hamiltonian" initialisation (their Section V).
    """
    if evolution == "full":
        step = trotter_step(n_sites, m, w, J, dt, order)
    elif evolution == "kinetic":
        step = QuantumCircuit(n_sites, name="kinetic")
        _append_hopping(step, n_sites, w, dt, 0)
        _append_hopping(step, n_sites, w, dt, 1)
    else:
        raise ValueError("evolution must be 'full' or 'kinetic'")
    out = []
    qc = bare_vacuum_circuit(n_sites)
    for _ in range(n_krylov):
        qc = qc.compose(step)
        out.append(qc.copy())
    return out


# ----------------------------------------------------------------------------
# Running SKQD on sampled data
# ----------------------------------------------------------------------------


@dataclass
class SKQDStep:
    """The SKQD estimate after accumulating Krylov steps 1..k."""

    k: int
    energy: float
    dimension: int
    indices: tuple[int, ...]
    vector: np.ndarray


def run_skqd(distributions: list[np.ndarray], shots: int, rng: np.random.Generator,
             n_sites: int, m: float, w: float, J: float, include_reference: bool = True,
             charge: int = 0) -> list[SKQDStep]:
    """
    Sample `shots` bitstrings from each Krylov distribution, postselect onto the
    charge sector, accumulate the set, and diagonalise after every step.

    `distributions[k-1]` is the exact outcome distribution of the k-th Krylov
    circuit (noise-free or noisy, from noise.probabilities). Sampling from it
    is statistically identical to running that circuit `shots` times.

    include_reference adds the bare-vacuum bitstring to the set from the start.
    It is known classically, costs no shots, and guarantees a nonempty subspace.
    """
    ham = hamiltonian(n_sites, m, w, J)
    in_sector = np.isclose(z_basis_table(n_sites).charge, charge)
    found: set[int] = set()
    if include_reference:
        found.add(basis_index(bare_vacuum_label(n_sites)))
    out = []
    for k, p in enumerate(distributions, start=1):
        clean = np.clip(np.asarray(p, dtype=float), 0, None)
        counts = rng.multinomial(shots, clean / clean.sum())
        found.update(int(b) for b in np.nonzero(counts)[0] if in_sector[b])
        if not found:
            continue
        e, vec = subspace_ground_state(found, ham)
        out.append(SKQDStep(k, e, len(found), tuple(sorted(found)), vec))
    return out


# ----------------------------------------------------------------------------
# The exact finite-shot expectation
# ----------------------------------------------------------------------------


@dataclass
class SKQDExpectation:
    mean_energy: float
    std_energy: float
    exact_energy: float
    prob_full_support: float       # P(every bitstring of the ground state was sampled)
    subset_probs: np.ndarray       # P(sampled set = S), indexed by bitmask over `basis`
    subset_energies: np.ndarray    # E_0 restricted to S
    basis: np.ndarray              # sector basis indices; bit i of a mask is basis[i]

    @property
    def mean_bias(self) -> float:
        return self.mean_energy - self.exact_energy


def expected_skqd_energy(distributions: list[np.ndarray], shots: int, n_sites: int,
                         m: float, w: float, J: float, include_reference: bool = True,
                         charge: int = 0, support_tol: float = 1e-12) -> SKQDExpectation:
    """
    Exact mean and spread of the SKQD energy over all possible M-shot runs.

    Enumerates every subset of the sector basis, so it is limited to small
    sectors (2^6 = 64 subsets at N = 4). See the module docstring for the
    formula. With include_reference, sets lacking the reference bitstring have
    probability zero, which the inversion handles automatically.
    """
    ham = hamiltonian(n_sites, m, w, J)
    basis = sector_indices(n_sites, charge)
    nb = basis.size
    if nb > 14:
        raise ValueError("sector too large for exact subset enumeration")
    masks = np.arange(2**nb)
    ref = basis_index(bare_vacuum_label(n_sites))
    ref_bit = int(np.where(basis == ref)[0][0]) if include_reference else None

    # P(sampled set is inside T) for every T
    sector_p = np.array([[p[b] / np.sum(p) for b in basis] for p in distributions])   # (K, nb)
    member = ((masks[:, None] >> np.arange(nb)[None, :]) & 1).astype(bool)          # (2^nb, nb)
    outside = sector_p @ (~member).T                                                # (K, 2^nb)
    inside_T = np.prod(np.clip(1.0 - outside, 0.0, 1.0) ** shots, axis=0)
    if ref_bit is not None:
        inside_T = np.where(member[:, ref_bit], inside_T, 0.0)

    # Moebius inversion over the subset lattice: P(set = S)
    g = inside_T.copy()
    for i in range(nb):
        has = (masks >> i) & 1 == 1
        g[has] -= g[masks[has] ^ (1 << i)]

    energies = np.full(2**nb, np.nan)
    for mask in masks[1:]:
        energies[mask], _ = subspace_ground_state(basis[member[mask]], ham)

    e0, gs = subspace_ground_state(basis, ham)
    support = np.abs(gs[basis]) ** 2 > support_tol
    full_support = np.all(member[:, support], axis=1)

    valid = masks > 0
    weight = g[valid].sum()                   # < 1 only if an empty set is possible
    mean = float(np.sum(g[valid] * energies[valid]) / weight)
    var = float(np.sum(g[valid] * (energies[valid] - mean) ** 2) / weight)
    return SKQDExpectation(mean, np.sqrt(max(var, 0.0)), e0, float(g[full_support].sum()),
                           g, energies, basis)
