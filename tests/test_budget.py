"""
Tests for the error-budget decomposition and the threshold machinery.

Known answers:
  - with no noise, the gate and readout biases are exactly zero;
  - readout bias on the particle density has the closed form -2q (nu - 1/2);
  - the Trotter bias of the second-order formula falls with slope 2 in dt;
  - the mean squared error of simulated M-shot measurements equals
    total_bias^2 + sigma^2, which is what "the budget adds up" means;
  - on data that obey the global-depolarizing model exactly, the predicted dip
    from one noise level reproduces the true dip at every other level;
  - first_crossing recovers a known root.
"""
import numpy as np
import pytest

pytest.importorskip("qiskit_aer")

from schwinger.budget import (
    budget, dip, exact_distribution, first_crossing, point_circuit, predicted_dip,
)
from schwinger.diagnostics import global_depolarize, observables_from_counts
from schwinger.noise import apply_readout_error, gate_noise_model, probabilities, sample_counts

N, M, W, J = 4, 0.5, 1.0, 1.0


@pytest.fixture(scope="module")
def trough():
    return point_circuit(N, M, W, J, 1.76, 18)


def test_noise_free_budget_has_only_trotter_error(trough):
    b = budget(trough, N, None, 0.0, 1000)
    assert b.gate_bias == pytest.approx(0.0, abs=1e-12)
    assert b.readout_bias == pytest.approx(0.0, abs=1e-12)
    assert b.total_bias == pytest.approx(b.trotter_bias, abs=1e-12)


@pytest.mark.parametrize("q", [0.01, 0.03])
def test_readout_bias_closed_form(trough, q):
    b = budget(trough, N, gate_noise_model(1e-3), q, 1000)
    assert b.readout_bias == pytest.approx(-2 * q * (b.gate - 0.5), abs=1e-12)


def test_biases_telescope_to_total(trough):
    b = budget(trough, N, gate_noise_model(2e-3), 0.02, 1000)
    assert b.trotter_bias + b.gate_bias + b.readout_bias == pytest.approx(b.total_bias, abs=1e-14)


def test_trotter_bias_scales_as_dt_squared():
    """
    The bias on the density, not the state error, so this is a separate check
    from the one in test_circuits.py. Small dt only, to stay in the asymptotic
    regime where the leading power dominates.
    """
    t = 1.76
    steps = np.array([20, 40, 80, 160])
    biases = [abs(budget(point_circuit(N, M, W, J, t, int(s)), N, None, 0.0, 1000).trotter_bias)
              for s in steps]
    slope = np.polyfit(np.log(t / steps), np.log(biases), 1)[0]
    assert slope == pytest.approx(2.0, abs=0.05)


def test_mean_squared_error_of_real_runs_matches_budget(trough):
    """
    Simulate 600 independent 500-shot experiments at a fixed noise level and
    check that their mean squared deviation from the exact answer equals
    total_bias^2 + sigma^2. This is the statement that the budget is complete:
    nothing is left over once the four pieces are accounted for.
    """
    shots, runs, q = 500, 600, 0.02
    nm = gate_noise_model(1e-3)
    b = budget(trough, N, nm, q, shots)
    measured = apply_readout_error(probabilities(trough.prepared, nm), q)
    rng = np.random.default_rng(99)
    est = np.array([observables_from_counts(sample_counts(measured, shots, rng), N)["particle_density"]
                    for _ in range(runs)])
    mse = np.mean((est - b.exact) ** 2)
    assert mse == pytest.approx(b.total_bias**2 + b.sigma**2, rel=0.12)


def test_prediction_exact_when_global_model_holds():
    """
    Synthetic data obeying p(eps) = F(eps) p_ideal + (1 - F)/16 with
    F = exp(-kappa eps). Predicting from one eps0 must reproduce D everywhere.
    """
    p_peak = exact_distribution(N, M, W, J, 0.89)
    p_trough = exact_distribution(N, M, W, J, 1.76)
    k_peak, k_trough, eps0 = 300.0, 600.0, 1e-3
    grid = np.logspace(-4, -2, 9)
    pr = predicted_dip(global_depolarize(p_peak, np.exp(-k_peak * eps0)),
                       global_depolarize(p_trough, np.exp(-k_trough * eps0)), eps0, grid, N, 1000)
    for e, d_pred in zip(grid, pr["D"]):
        d_true, _ = dip(global_depolarize(p_peak, np.exp(-k_peak * e)),
                        global_depolarize(p_trough, np.exp(-k_trough * e)), N, 1000)
        assert d_pred == pytest.approx(d_true, abs=1e-12)


def test_first_crossing_finds_known_root():
    x = np.logspace(-4, -1, 31)
    assert first_crossing(x, np.log(3.3e-3 / x)) == pytest.approx(3.3e-3, rel=1e-12)
    assert first_crossing(x, np.ones_like(x)) == float("inf")
