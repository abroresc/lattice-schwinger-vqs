"""
Tests for exact diagonalisation on longer chains, inside the neutral sector.

Every check compares the block calculation with an answer it does not share
code with: the full-space evolution, a closed-form limit, or a count.
"""
from math import comb

import numpy as np
import pytest

from schwinger.finite_size import (
    first_peak_and_trough, ground_state, neutral_block, neutral_indices, quench_density,
)
from schwinger.observables import evolve, particle_density
from schwinger.schwinger import bare_vacuum, hamiltonian


@pytest.mark.parametrize("n", [2, 4, 6, 8])
def test_neutral_sector_has_binomial_size(n):
    """Q = 0 means exactly N/2 occupied sites: binom(N, N/2) bitstrings."""
    assert neutral_indices(n).size == comb(n, n // 2)


@pytest.mark.parametrize("n", [4, 6])
def test_block_evolution_matches_full_space_evolution(n):
    """The neutral block alone must give the same nu(t) as the full 2^N space."""
    t = np.linspace(0, 3, 13)
    full = [particle_density(s, n) for s in evolve(bare_vacuum(n), hamiltonian(n, 0.5, 1, 1), t)]
    assert np.allclose(quench_density(n, 0.5, 1.0, 1.0, t), full, atol=1e-12)


def test_block_ground_state_is_the_global_ground_state():
    """The ground state lies in the neutral sector, so the block minimum is E0."""
    e_full = np.linalg.eigvalsh(hamiltonian(6, 0.5, 1.0, 1.0))[0]
    assert ground_state(6, 0.5, 1.0, 1.0)[0] == pytest.approx(e_full, abs=1e-10)


def test_strong_coupling_limit_on_a_long_chain():
    """At w = 0 the bare vacuum is exact: E0 = -mN/2, condensate -1/2."""
    e0, cond = ground_state(10, 0.7, 0.0, 1.0)
    assert e0 == pytest.approx(-0.7 * 10 / 2, abs=1e-10)
    assert cond == pytest.approx(-0.5, abs=1e-10)


def test_free_fermion_limit_on_a_long_chain():
    """At m = J = 0 the ground state fills the N/2 lowest levels 2w cos(k pi/(N+1))."""
    n, w = 8, 1.0
    levels = np.sort(2 * w * np.cos(np.arange(1, n + 1) * np.pi / (n + 1)))
    assert ground_state(n, 0.0, w, 0.0)[0] == pytest.approx(levels[: n // 2].sum(), abs=1e-10)


def test_quench_starts_at_zero_and_stays_physical():
    nu = quench_density(8, 0.5, 1.0, 1.0, np.linspace(0, 4, 41))
    assert nu[0] == pytest.approx(0.0, abs=1e-14)
    assert np.all(nu >= -1e-12) and np.all(nu <= 1 + 1e-12)


def test_peak_and_trough_finder_on_a_known_curve():
    """On -cos(t), the first maximum is at pi and the minimum after it at 2 pi."""
    t = np.linspace(0.1, 7, 6901)
    tp, yp, tt, yt = first_peak_and_trough(t, -np.cos(t))
    assert tp == pytest.approx(np.pi, abs=2e-3) and yp == pytest.approx(1.0, abs=1e-6)
    assert tt == pytest.approx(2 * np.pi, abs=2e-3) and yt == pytest.approx(-1.0, abs=1e-6)


def test_block_is_cached_not_rebuilt():
    assert neutral_block(6, 0.5, 1.0, 1.0) is neutral_block(6, 0.5, 1.0, 1.0)
