"""
Tests for the Qiskit circuit layer.

The circuit layer is never trusted because Qiskit ran without raising. Every
circuit is compared against the exact NumPy result from the rest of the package.
"""
import numpy as np
import pytest

pytest.importorskip("qiskit")

from schwinger.circuits import (
    bare_vacuum_circuit, circuit_metrics, hamiltonian_op,
    particle_density_from_counts, run_with_shots, statevector, trotter_circuit,
    trotter_step,
)
from schwinger.observables import evolve, particle_density
from schwinger.schwinger import bare_vacuum, hamiltonian

TOL = 1e-10


def phase_align(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Remove the global phase difference; it is unobservable."""
    return a * np.exp(-1j * np.angle(np.vdot(b, a)))


@pytest.mark.parametrize("n", [2, 4, 6])
def test_basis_orderings_coincide(n):
    """
    Qiskit indexes with qubit 0 least significant; our site j is qubit N-1-j,
    so site N-1 is least significant for us too. The orderings coincide and no
    permutation is needed. Checked, not assumed.
    """
    from qiskit import QuantumCircuit
    from schwinger.paulis import basis_index, basis_label, site_to_qiskit_qubit
    for i in range(2 ** n):
        label = basis_label(i, n)
        qc = QuantumCircuit(n)
        for j, bit in enumerate(label):
            if bit == "1":
                qc.x(site_to_qiskit_qubit(j, n))
        sv = statevector(qc)
        assert abs(sv[basis_index(label)]) == pytest.approx(1.0)


@pytest.mark.parametrize("n", [2, 4, 6])
def test_bare_vacuum_circuit_matches_numpy(n):
    assert np.allclose(statevector(bare_vacuum_circuit(n)), bare_vacuum(n))


@pytest.mark.parametrize("n", [2, 3, 4, 5])
@pytest.mark.parametrize("m,w,J", [(1.0, 1.0, 1.0), (0.3, 1.0, 0.7), (2.0, 0.5, 1.5)])
def test_sparse_pauli_op_matches_numpy(n, m, w, J):
    assert np.allclose(hamiltonian_op(n, m, w, J).to_matrix(), hamiltonian(n, m, w, J))


@pytest.mark.parametrize("order", [1, 2])
@pytest.mark.parametrize("m,w,J", [(0.5, 1.0, 1.0), (1.0, 1.0, 2.0)])
def test_trotter_converges_to_exact(order, m, w, J):
    """Many small steps must reproduce exp(-iHt) to high accuracy."""
    n, t = 4, 2.0
    exact = evolve(bare_vacuum(n), hamiltonian(n, m, w, J), np.array([t]))[0]
    got = phase_align(statevector(trotter_circuit(n, m, w, J, t, 4000, order)), exact)
    assert np.max(np.abs(got - exact)) < (2e-3 if order == 1 else 1e-5)


@pytest.mark.parametrize("order,expected_slope", [(1, 1.0), (2, 2.0)])
def test_trotter_error_scaling_exponent(order, expected_slope):
    """
    The main circuit test. Fit the slope of log(error) against log(1/steps) and check it
    matches the predicted exponent. The slope is the test, not the size of the
    error: a wrong implementation can be accurate and still scale wrongly.
    """
    n, m, w, J, t = 4, 0.5, 1.0, 1.0, 2.0
    exact = evolve(bare_vacuum(n), hamiltonian(n, m, w, J), np.array([t]))[0]
    steps = np.array([20, 40, 80, 160, 320])
    errs = []
    for s in steps:
        got = phase_align(statevector(trotter_circuit(n, m, w, J, t, int(s), order)), exact)
        errs.append(np.max(np.abs(got - exact)))
    slope = np.polyfit(np.log(1.0 / steps), np.log(np.array(errs)), 1)[0]
    assert slope == pytest.approx(expected_slope, abs=0.05), f"measured slope {slope:.3f}"


def _step_operators(n, m, w, J, dt):
    """exp(-i dt X) for the diagonal, outer-bond and middle-bond parts, from the Pauli terms."""
    from scipy.linalg import expm

    from schwinger.paulis import pauli_string
    from schwinger.schwinger import hamiltonian_terms

    terms = hamiltonian_terms(n, m, w, J)
    part = lambda keep: sum(c * pauli_string(s) for s, c in terms if keep(s))
    diag = part(lambda s: set(s) <= {"I", "Z"})
    outer = part(lambda s: s in ("XXII", "YYII", "IIXX", "IIYY"))
    middle = part(lambda s: s in ("IXXI", "IYYI"))
    return [lambda f, A=A: expm(-1j * A * dt * f) for A in (diag, outer, middle)]


def test_first_order_step_applies_the_diagonal_terms_first():
    """
    The circuit applies the diagonal layer first, so as an operator it is
    e^{-i H_o dt} e^{-i H_e dt} e^{-i H_Z dt}. The reverse product is a
    different operator, and an easy one to write down by mistake.
    """
    from qiskit.quantum_info import Operator
    n, m, w, J, dt = 4, 0.5, 1.0, 1.0, 0.3
    Z, E, O = _step_operators(n, m, w, J, dt)
    got = Operator(trotter_step(n, m, w, J, dt, order=1)).data
    assert np.max(np.abs(phase_align(got, O(1) @ E(1) @ Z(1)) - O(1) @ E(1) @ Z(1))) < TOL
    assert np.max(np.abs(phase_align(got, Z(1) @ E(1) @ O(1)) - Z(1) @ E(1) @ O(1))) > 0.1


def test_second_order_step_is_the_symmetric_product():
    """Eq. (strang): half diagonal, half outer, full middle, half outer, half diagonal."""
    from qiskit.quantum_info import Operator
    n, m, w, J, dt = 4, 0.5, 1.0, 1.0, 0.3
    Z, E, O = _step_operators(n, m, w, J, dt)
    want = Z(.5) @ E(.5) @ O(1) @ E(.5) @ Z(.5)
    got = Operator(trotter_step(n, m, w, J, dt, order=2)).data
    assert np.max(np.abs(phase_align(got, want) - want)) < TOL


def test_single_trotter_step_is_unitary():
    n = 4
    from qiskit.quantum_info import Operator
    op = Operator(trotter_step(n, 0.5, 1.0, 1.0, 0.1, 2)).data
    assert np.allclose(op.conj().T @ op, np.eye(2 ** n), atol=1e-10)


def test_trotter_step_rejects_bad_order():
    with pytest.raises(ValueError):
        trotter_step(4, 0.5, 1.0, 1.0, 0.1, order=3)
    with pytest.raises(ValueError):
        trotter_circuit(4, 0.5, 1.0, 1.0, 1.0, steps=0)


def test_counts_estimator_is_exact_on_a_known_bitstring():
    """
    Feed known bitstrings and check the particle number by hand.
    Bare vacuum 0101 -> 0 particles. 0011 -> one pair, 2 particles.
    """
    n = 4
    d, err = particle_density_from_counts({"0101": 1000}, n)
    assert d == pytest.approx(0.0) and err == pytest.approx(0.0)
    d, err = particle_density_from_counts({"0011": 1000}, n)
    assert d == pytest.approx(2.0 / n) and err == pytest.approx(0.0)
    d, _ = particle_density_from_counts({"0101": 500, "0011": 500}, n)
    assert d == pytest.approx(1.0 / n)


def test_counts_estimator_validates_length():
    with pytest.raises(ValueError):
        particle_density_from_counts({"010": 10}, 4)
    with pytest.raises(ValueError):
        particle_density_from_counts({}, 4)


def test_shot_sampling_agrees_with_exact_within_error_bars():
    """
    Sample the simulator and check the estimate sits within a few sigma of the
    exact value. This is the first appearance of shot noise, which is the
    statistical half of the error budget.
    """
    n, m, w, J, t = 4, 0.5, 1.0, 1.0, 3.0
    exact_state = evolve(bare_vacuum(n), hamiltonian(n, m, w, J), np.array([t]))[0]
    exact = particle_density(exact_state, n)
    qc = trotter_circuit(n, m, w, J, t, 400, order=2)
    counts = run_with_shots(qc, shots=20000, seed=1234)
    est, err = particle_density_from_counts(counts, n)
    assert err > 0
    assert abs(est - exact) < 5 * err, f"est {est:.4f} +- {err:.4f} vs exact {exact:.4f}"


def test_shot_noise_falls_as_one_over_sqrt_shots():
    n, m, w, J, t = 4, 0.5, 1.0, 1.0, 3.0
    qc = trotter_circuit(n, m, w, J, t, 200, order=2)
    _, e1 = particle_density_from_counts(run_with_shots(qc, 2000, seed=7), n)
    _, e2 = particle_density_from_counts(run_with_shots(qc, 32000, seed=7), n)
    assert e1 / e2 == pytest.approx(4.0, rel=0.25)


def test_circuit_metrics_grow_with_step_count():
    """
    Fewer Trotter error means more steps means deeper circuits means more gate
    error. Quantifying that tension is one of the two headline thesis results,
    so the metric it depends on is tested.
    """
    n = 4
    small = circuit_metrics(trotter_circuit(n, 0.5, 1.0, 1.0, 2.0, 10, 2))
    large = circuit_metrics(trotter_circuit(n, 0.5, 1.0, 1.0, 2.0, 100, 2))
    assert large["two_qubit_gates"] > small["two_qubit_gates"]
    assert large["depth"] > small["depth"]
    assert small["two_qubit_gates"] > 0
