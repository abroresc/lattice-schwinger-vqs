"""
Tests for SKQD ground-state preparation.

Known answers:
  - with every sector bitstring, SKQD is exact diagonalisation of the sector;
  - the variational bound E(S) >= E_0 holds for every one of the 63 subsets;
  - Cauchy interlacing: adding a bitstring never raises the energy;
  - the single-bitstring subspace {bare vacuum} gives <vac|H|vac> = -m N / 2,
    the bare-vacuum energy listed in scripts/04 (independent of w and J);
  - the exact finite-shot expectation agrees with brute-force Monte Carlo;
  - the subset probabilities form a probability distribution.
"""
import itertools

import numpy as np
import pytest

pytest.importorskip("qiskit_aer")

from schwinger.noise import gate_noise_model, noisy_probabilities
from schwinger.observables import chiral_condensate
from schwinger.paulis import basis_index
from schwinger.schwinger import bare_vacuum_label, ground_state, hamiltonian
from schwinger.skqd import (
    expected_skqd_energy, krylov_circuits, run_skqd, sector_indices, subspace_ground_state,
)

N, M, W, J = 4, 0.5, 1.0, 1.0


@pytest.fixture(scope="module")
def krylov_dists():
    return [noisy_probabilities(c) for c in krylov_circuits(N, M, W, J, 0.3, 4)]


@pytest.mark.parametrize("m", [0.0, 0.5, 2.0])
def test_full_sector_reproduces_exact_ground_state(m):
    ham = hamiltonian(N, m, W, J)
    e_s, v_s = subspace_ground_state(sector_indices(N), ham)
    e0, v0 = ground_state(N, m, W, J)
    assert e_s == pytest.approx(e0, abs=1e-12)
    assert abs(abs(np.vdot(v0, v_s)) - 1.0) < 1e-10
    assert chiral_condensate(v_s, N) == pytest.approx(chiral_condensate(v0, N), abs=1e-10)


def test_variational_bound_and_interlacing_on_every_subset():
    ham = hamiltonian(N, M, W, J)
    basis = list(sector_indices(N))
    e0, _ = ground_state(N, M, W, J)
    energies = {}
    for r in range(1, len(basis) + 1):
        for subset in itertools.combinations(basis, r):
            e, _ = subspace_ground_state(subset, ham)
            assert e >= e0 - 1e-12
            energies[frozenset(subset)] = e
    for s, e in energies.items():
        for b in basis:
            if b not in s:
                assert energies[s | {b}] <= e + 1e-12


def test_bare_vacuum_alone_gives_bare_vacuum_energy():
    for m, w, j in [(0.5, 1.0, 1.0), (2.0, 0.3, 1.7)]:
        e, _ = subspace_ground_state([basis_index(bare_vacuum_label(N))], hamiltonian(N, m, w, j))
        assert e == pytest.approx(-m * N / 2, abs=1e-12)


def test_subset_probabilities_are_a_distribution(krylov_dists):
    ex = expected_skqd_energy(krylov_dists, 20, N, M, W, J)
    assert ex.subset_probs.min() > -1e-12
    assert ex.subset_probs.sum() == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize("shots", [5, 15, 40])
def test_exact_expectation_matches_monte_carlo(krylov_dists, shots):
    """
    The inclusion-exclusion formula against 3000 simulated SKQD runs: the mean
    energy, and the frequency of every sampled subspace that is common enough
    to test (expected count >= 10), each within four binomial standard errors.

    The spread of the energy is deliberately not compared. It is dominated by
    subspaces with probability ~1e-4 (both long-string bitstrings missed), so
    its sample estimate swings by a factor of four between 3000-run batches.
    A test of it would fail on a correct formula.
    """
    ex = expected_skqd_energy(krylov_dists, shots, N, M, W, J)
    basis = [int(b) for b in ex.basis]
    rng = np.random.default_rng(2026)
    runs = 3000
    energies = np.empty(runs)
    freq = np.zeros(ex.subset_probs.size)
    for i in range(runs):
        last = run_skqd(krylov_dists, shots, rng, N, M, W, J)[-1]
        energies[i] = last.energy
        freq[sum(1 << basis.index(b) for b in last.indices)] += 1
    sem = energies.std() / np.sqrt(runs) + 1e-15
    assert abs(energies.mean() - ex.mean_energy) < 4 * sem + 1e-12
    for mask in np.where(ex.subset_probs * runs >= 10)[0]:
        p = ex.subset_probs[mask]
        assert abs(freq[mask] / runs - p) < 4 * np.sqrt(p * (1 - p) / runs)


def test_many_shots_converge_to_exact_energy(krylov_dists):
    ex = expected_skqd_energy(krylov_dists, 2000, N, M, W, J)
    assert ex.mean_bias == pytest.approx(0.0, abs=1e-12)
    assert ex.prob_full_support == pytest.approx(1.0, abs=1e-12)


def test_noise_cannot_push_energy_below_truth():
    """Noisy samples only change which bitstrings appear; the bound still holds."""
    dists = [noisy_probabilities(c, gate_noise_model(5e-3))
             for c in krylov_circuits(N, M, W, J, 0.3, 3)]
    e0, _ = ground_state(N, M, W, J)
    rng = np.random.default_rng(3)
    for _ in range(50):
        for step in run_skqd(dists, 20, rng, N, M, W, J):
            assert step.energy >= e0 - 1e-12
