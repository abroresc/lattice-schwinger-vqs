"""
Tests for the Z-basis observables and symmetry diagnostics.

Each test compares against something known independently: a value computed
through the operator route in observables.py (a separate construction), a
closed form for the maximally mixed state, or a distribution built to satisfy a
model exactly. None of them checks the code against itself.
"""
from math import comb

import numpy as np
import pytest

from schwinger.diagnostics import (
    fidelity_from_charge_variance, fidelity_from_neutral_fraction, global_depolarize,
    mixed_state_values, neutral_sector_size, observables_from_counts,
    observables_from_distribution, shot_noise_sigma, undo_global_depolarizing, z_basis_table,
)
from schwinger.observables import (
    charge_variance, electric_field_profile, evolve, particle_density, total_charge,
)
from schwinger.schwinger import bare_vacuum, hamiltonian

TOL = 1e-12


def random_state(n, seed):
    rng = np.random.default_rng(seed)
    v = rng.normal(size=2**n) + 1j * rng.normal(size=2**n)
    return v / np.linalg.norm(v)


@pytest.mark.parametrize("n", [2, 3, 4, 5, 6])
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_diagonal_observables_match_operator_route(n, seed):
    """
    The bitstring table and the operator expectations in observables.py are two
    unrelated constructions of the same numbers. On a random state they agree.
    """
    psi = random_state(n, seed)
    obs = observables_from_distribution(np.abs(psi) ** 2, n)
    assert obs["particle_density"] == pytest.approx(particle_density(psi, n), abs=TOL)
    assert obs["charge_mean"] == pytest.approx(total_charge(psi, n), abs=TOL)
    assert obs["charge_variance"] == pytest.approx(charge_variance(psi, n), abs=TOL)
    assert np.allclose(obs["links"], electric_field_profile(psi, n), atol=TOL)


@pytest.mark.parametrize("n", [2, 4, 6])
def test_bare_vacuum_is_empty_neutral_and_field_free(n):
    obs = observables_from_distribution(np.abs(bare_vacuum(n)) ** 2, n)
    assert obs["particle_density"] == pytest.approx(0.0, abs=TOL)
    assert obs["charge_variance"] == pytest.approx(0.0, abs=TOL)
    assert obs["neutral_fraction"] == pytest.approx(1.0, abs=TOL)
    assert np.allclose(obs["links"], 0.0, atol=TOL)


@pytest.mark.parametrize("n", [2, 4, 6, 8])
def test_maximally_mixed_values_match_closed_form(n):
    """Var(Q) = N/4, P(Q=0) = C(N,N/2)/2^N, density 1/2, for even N."""
    vals = mixed_state_values(n)
    assert vals["charge_variance"] == pytest.approx(n / 4, abs=TOL)
    assert vals["neutral_fraction"] == pytest.approx(comb(n, n // 2) / 2**n, abs=TOL)
    assert vals["particle_density"] == pytest.approx(0.5, abs=TOL)
    assert neutral_sector_size(n) == comb(n, n // 2)


def test_neutral_sector_at_four_sites_is_martinez_six_states():
    """Martinez et al. list exactly six physical states at N = 4."""
    t = z_basis_table(4)
    assert int(np.sum(np.isclose(t.charge, 0.0))) == 6


@pytest.mark.parametrize("fid", [1.0, 0.9, 0.5, 0.123, 0.0])
def test_fidelity_estimators_exact_for_global_depolarizing(fid):
    """
    Build a distribution that IS a global depolarizing mixture of a neutral
    state. Both estimators must return the mixing weight exactly. This is the
    model the diagnostic-based threshold prediction assumes, so it is tested
    where it is exact before it is used where it is not.
    """
    n = 4
    psi = evolve(bare_vacuum(n), hamiltonian(n, 0.5, 1.0, 1.0), np.array([1.3]))[0]
    noisy = global_depolarize(np.abs(psi) ** 2, fid)
    obs = observables_from_distribution(noisy, n)
    assert fidelity_from_charge_variance(obs["charge_variance"], n) == pytest.approx(fid, abs=TOL)
    assert fidelity_from_neutral_fraction(obs["neutral_fraction"], n) == pytest.approx(fid, abs=TOL)


def test_global_depolarizing_round_trip():
    n = 4
    p = np.abs(random_state(n, 7)) ** 2
    assert np.allclose(undo_global_depolarizing(global_depolarize(p, 0.37), 0.37), p, atol=TOL)


def test_cp_asymmetry_vanishes_on_exact_neutral_evolution():
    """<L_1> = <L_3> on the CP-symmetric trajectory from the bare vacuum; the table must see that."""
    n = 4
    states = evolve(bare_vacuum(n), hamiltonian(n, 0.5, 1.0, 1.0), np.linspace(0, 6, 25))
    worst = max(abs(observables_from_distribution(np.abs(s) ** 2, n)["cp_asymmetry"])
                for s in states)
    assert worst < 1e-12


def test_postselected_density_ignores_charged_shots():
    n = 4
    t = z_basis_table(n)
    p = np.zeros(16)
    neutral = np.where(np.isclose(t.charge, 0))[0]
    charged = np.where(~np.isclose(t.charge, 0))[0]
    p[neutral[0]] = 0.5
    p[charged[0]] = 0.5
    obs = observables_from_distribution(p, n)
    assert obs["postselected_density"] == pytest.approx(t.particle_density[neutral[0]], abs=TOL)
    assert obs["neutral_fraction"] == pytest.approx(0.5, abs=TOL)


@pytest.mark.parametrize("shots", [256, 2048])
def test_standard_errors_match_empirical_scatter(shots):
    """
    The quoted standard errors must describe the actual run-to-run scatter.
    Draw many independent M-shot data sets from one exact distribution and
    compare the spread of the estimates with the formulas. Tolerance 15 %,
    which is about four standard errors of a standard deviation from 400 runs.
    """
    n, runs = 4, 400
    psi = evolve(bare_vacuum(n), hamiltonian(n, 0.5, 1.0, 1.0), np.array([0.9]))[0]
    p = global_depolarize(np.abs(psi) ** 2, 0.6)        # some charge violation to measure
    rng = np.random.default_rng(11)
    dens, varq, quoted_d, quoted_v = [], [], [], []
    for _ in range(runs):
        obs = observables_from_counts(rng.multinomial(shots, p).astype(float), n)
        dens.append(obs["particle_density"]); quoted_d.append(obs["particle_density_err"])
        varq.append(obs["charge_variance"]); quoted_v.append(obs["charge_variance_err"])
    assert np.std(dens) == pytest.approx(np.mean(quoted_d), rel=0.15)
    assert np.std(varq) == pytest.approx(np.mean(quoted_v), rel=0.15)
    assert np.mean(quoted_d) == pytest.approx(shot_noise_sigma(p, n, shots), rel=0.05)
