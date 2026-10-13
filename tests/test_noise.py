"""
Tests for the gate and readout error layer (M11).

Known answers used here, none of them produced by the code under test:
  - noise-free density-matrix probabilities equal |statevector|^2;
  - the NumPy step-channel route and the Aer route agree when both are noise free;
  - Qiskit's own average_gate_fidelity of the channel equals 1 - eps;
  - a symmetric readout flip multiplies every <Z_n> by (1 - 2q), in closed form;
  - a fully depolarizing CX makes the register maximally mixed;
  - a global channel applied per step gives F = (1 - p)^steps exactly;
  - logical-frame dephasing commutes with Q, so Var(Q) = 0 exactly;
  - the Trotterised circuit preserves CP, so the asymmetry is zero without noise.
"""
import numpy as np
import pytest

pytest.importorskip("qiskit_aer")

from schwinger.circuits import statevector, trotter_circuit
from schwinger.diagnostics import (
    fidelity_from_charge_variance, mixed_state_values, observables_from_distribution,
)
from schwinger.noise import (
    apply_readout_error, depolarizing_lambda, gate_noise_model, noisy_probabilities,
    sample_counts, step_channel_probabilities,
)

N, M, W, J = 4, 0.5, 1.0, 1.0
TOL = 1e-10


@pytest.mark.parametrize("time,steps", [(0.9, 9), (1.76, 18)])
def test_noise_free_probabilities_equal_statevector(time, steps):
    qc = trotter_circuit(N, M, W, J, time, steps, 2)
    assert np.allclose(noisy_probabilities(qc, None), np.abs(statevector(qc)) ** 2, atol=TOL)


def test_numpy_step_route_matches_aer_route_without_noise():
    """Two unrelated simulators of the same circuit: Aer density matrix, NumPy products."""
    qc = trotter_circuit(N, M, W, J, 1.76, 18, 2)
    a = noisy_probabilities(qc, None)
    b = step_channel_probabilities(N, M, W, J, 1.76, 18, 2, "depolarizing", 0.0)
    assert np.allclose(a, b, atol=TOL)


@pytest.mark.parametrize("n_qubits", [1, 2])
@pytest.mark.parametrize("eps", [1e-4, 3e-3, 2e-2])
def test_infidelity_convention_against_qiskit(n_qubits, eps):
    """eps is the AVERAGE gate infidelity, the number a data sheet quotes."""
    from qiskit.quantum_info import average_gate_fidelity
    from qiskit_aer.noise import depolarizing_error

    channel = depolarizing_error(depolarizing_lambda(eps, n_qubits), n_qubits).to_quantumchannel()
    assert 1.0 - average_gate_fidelity(channel) == pytest.approx(eps, rel=1e-9)


@pytest.mark.parametrize("q", [0.0, 0.01, 0.05, 0.2])
def test_readout_flip_scales_every_z_by_one_minus_2q(q):
    """
    Closed form: a symmetric flip with probability q maps <Z_n> to (1-2q)<Z_n>.
    For the particle density that means nu -> 1/2 + (1-2q)(nu - 1/2).
    """
    qc = trotter_circuit(N, M, W, J, 0.9, 9, 2)
    p = noisy_probabilities(qc, None)
    nu = observables_from_distribution(p, N)["particle_density"]
    nu_q = observables_from_distribution(apply_readout_error(p, q), N)["particle_density"]
    assert nu_q == pytest.approx(0.5 + (1 - 2 * q) * (nu - 0.5), abs=TOL)


def test_readout_error_on_a_single_bitstring_is_a_product_of_bernoullis():
    p = np.zeros(16); p[0] = 1.0                       # |0000>
    out = apply_readout_error(p, 0.1)
    assert out[0] == pytest.approx(0.9**4, abs=TOL)
    assert out[0b1000] == pytest.approx(0.1 * 0.9**3, abs=TOL)
    assert out.sum() == pytest.approx(1.0, abs=TOL)


def test_classical_readout_matches_aer_readout_error_by_sampling():
    """Aer's ReadoutError applied to real shots versus our classical confusion matrix."""
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, ReadoutError
    from schwinger.circuits import to_hardware_basis

    q, shots = 0.08, 200_000
    qc = trotter_circuit(N, M, W, J, 0.9, 9, 2)
    predicted = apply_readout_error(noisy_probabilities(qc, None), q)
    nm = NoiseModel()
    nm.add_all_qubit_readout_error(ReadoutError([[1 - q, q], [q, 1 - q]]))
    meas = to_hardware_basis(qc); meas.measure_all()
    counts = AerSimulator(noise_model=nm).run(meas, shots=shots, seed_simulator=3).result().get_counts()
    observed = np.zeros(16)
    for bits, c in counts.items():
        observed[int(bits, 2)] = c / shots              # Qiskit's order coincides with ours
    sigma = np.sqrt(predicted * (1 - predicted) / shots) + 1e-12
    assert np.max(np.abs(observed - predicted) / sigma) < 5.0


def test_fully_depolarizing_cx_gives_maximally_mixed_register():
    """
    With eps2 at its maximum, every CX outputs I/4 on its pair. The even-bond
    layer touches all four qubits, so afterwards the register is I/16, and
    unitaries cannot undo that.
    """
    qc = trotter_circuit(N, M, W, J, 0.9, 9, 2)
    p = noisy_probabilities(qc, gate_noise_model(eps2=0.75, eps1=0.0))   # lambda = 1
    obs = observables_from_distribution(p, N)
    mixed = mixed_state_values(N)
    assert np.allclose(p, 1 / 16, atol=1e-10)
    assert obs["charge_variance"] == pytest.approx(mixed["charge_variance"], abs=1e-10)
    assert obs["particle_density"] == pytest.approx(0.5, abs=1e-10)


@pytest.mark.parametrize("p", [0.005, 0.02, 0.1])
def test_global_channel_per_step_gives_exact_fidelity(p):
    steps = 18
    probs = step_channel_probabilities(N, M, W, J, 1.76, steps, 2, "global", p)
    var_q = observables_from_distribution(probs, N)["charge_variance"]
    assert fidelity_from_charge_variance(var_q, N) == pytest.approx((1 - p) ** steps, abs=1e-12)


@pytest.mark.parametrize("p", [0.01, 0.05])
def test_logical_dephasing_is_invisible_to_charge_but_not_to_physics(p):
    """
    The blind spot. Z errors on the logical state commute with Q, so charge
    stays exactly conserved and CP stays exact, while the particle density
    moves well away from its noise-free value.
    """
    clean = observables_from_distribution(
        step_channel_probabilities(N, M, W, J, 1.76, 18, 2, "dephasing", 0.0), N)
    noisy = observables_from_distribution(
        step_channel_probabilities(N, M, W, J, 1.76, 18, 2, "dephasing", p), N)
    assert noisy["charge_variance"] == pytest.approx(0.0, abs=1e-12)
    assert noisy["neutral_fraction"] == pytest.approx(1.0, abs=1e-12)
    assert abs(noisy["cp_asymmetry"]) < 1e-12
    assert abs(noisy["particle_density"] - clean["particle_density"]) > 0.02


def test_gate_level_dephasing_is_visible_to_charge():
    """
    The correction to the blind-spot claim. A Z error between the CXs of a
    compiled RXX is an X or Y error in the logical frame, so gate-level
    dephasing DOES break charge conservation.
    """
    qc = trotter_circuit(N, M, W, J, 1.76, 18, 2)
    p = noisy_probabilities(qc, gate_noise_model(eps2=3e-3, kind="dephasing"))
    assert observables_from_distribution(p, N)["charge_variance"] > 0.1


@pytest.mark.parametrize("dt", [0.5, 0.25, 0.1])
def test_trotter_circuit_preserves_cp_without_noise(dt):
    """Any measured CP asymmetry is therefore hardware error, never Trotter error."""
    for t in (0.9, 1.76, 3.0):
        steps = max(1, int(round(t / dt)))
        p = noisy_probabilities(trotter_circuit(N, M, W, J, t, steps, 2), None)
        assert abs(observables_from_distribution(p, N)["cp_asymmetry"]) < 1e-12


def test_gate_noise_breaks_cp_and_charge():
    qc = trotter_circuit(N, M, W, J, 1.76, 18, 2)
    obs = observables_from_distribution(noisy_probabilities(qc, gate_noise_model(1e-3)), N)
    assert obs["charge_variance"] > 0.1
    assert abs(obs["cp_asymmetry"]) > 1e-3


def test_sample_counts_is_seeded_and_conserves_shots():
    p = np.full(16, 1 / 16)
    a = sample_counts(p, 1000, np.random.default_rng(5))
    b = sample_counts(p, 1000, np.random.default_rng(5))
    assert np.array_equal(a, b) and a.sum() == 1000
