"""
Gate and readout error (milestone M11), kept separate from shot noise.

THE DESIGN DECISION THAT MAKES THE ERROR BUDGET CLEAN
------------------------------------------------------
A noisy shot-sampled simulation mixes two errors in every number it returns:
the BIAS that gate noise puts into the state, and the random SCATTER from a
finite number of shots. Averaging many runs to separate them is slow and never
exact. So this module never samples to get a bias. It runs the circuit on Aer's
density-matrix simulator, which evolves the full 16 x 16 density matrix through
every noisy gate and returns the exact outcome probabilities. The gate-error
bias of any diagonal observable is then an exact number, and shot noise is
added afterwards, either analytically (diagnostics.shot_noise_sigma) or by
sampling from the exact distribution (sample_counts) when a realistic data set
is wanted. The two error sources never contaminate each other.

At N = 4 the density matrix has 256 entries, so this costs milliseconds.

THE ERROR MODELS
----------------
Errors attach to the physical gates of circuits.HARDWARE_BASIS:

- depolarizing on CX with average gate infidelity eps2 (the "two-qubit gate
  error" a hardware data sheet quotes, from randomized benchmarking);
- depolarizing on SX and X with infidelity eps1 (RZ is a software frame change
  and carries no error);
- optionally, pure dephasing on CX instead of depolarizing;
- symmetric readout error q, applied classically to the outcome probabilities.

A second, independent route (step_channel_probabilities) applies noise to the
LOGICAL state once per Trotter step, using plain NumPy density matrices. It
reproduces the error model Martinez et al. fitted to their data (dephasing per
qubit per step) and it is the control for the charge diagnostic: logical-frame
dephasing commutes with Q, so Var(Q) stays exactly zero. Gate-level dephasing
does NOT, because a Z error between the CXs of a compiled RXX is an X or Y
error in the logical frame. See the note in diagnostics.py.

Conversion between the depolarizing parameter lambda used by Aer,
    E(rho) = (1 - lambda) rho + lambda I/d,
and the average gate infidelity eps a data sheet reports:
    eps = lambda (d - 1) / d,   so lambda = eps d / (d - 1)
(d = 2 for one qubit, d = 4 for two). Getting this wrong would mis-state every
threshold by a factor of 4/3; tests/test_noise.py checks it against Qiskit's own
average_gate_fidelity rather than trusting the algebra.
"""

from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit

from .circuits import HARDWARE_BASIS, to_hardware_basis


# ----------------------------------------------------------------------------
# Noise models
# ----------------------------------------------------------------------------


def depolarizing_lambda(infidelity: float, n_qubits: int) -> float:
    """Aer depolarizing parameter for a given average gate infidelity."""
    d = 2 ** n_qubits
    lam = infidelity * d / (d - 1)
    if not 0.0 <= lam <= d * d / (d * d - 1):
        raise ValueError(f"infidelity {infidelity} out of range for {n_qubits} qubit(s)")
    return lam


def gate_noise_model(eps2: float, eps1: float | None = None, kind: str = "depolarizing"):
    """
    Aer NoiseModel with errors on CX (eps2) and on SX, X (eps1).

    eps1 defaults to eps2 / 10, roughly the ratio on current superconducting
    devices. kind = "dephasing" replaces the CX error by independent phase
    flips with probability eps2 on each of the two qubits, the control channel
    that charge conservation cannot see. Returns None when there is no error,
    so callers can pass the result straight to the simulator.
    """
    from qiskit_aer.noise import NoiseModel, depolarizing_error, pauli_error

    if eps1 is None:
        eps1 = eps2 / 10.0
    if eps2 == 0 and eps1 == 0:
        return None
    nm = NoiseModel(basis_gates=list(HARDWARE_BASIS))
    if kind == "depolarizing":
        if eps2 > 0:
            nm.add_all_qubit_quantum_error(depolarizing_error(depolarizing_lambda(eps2, 2), 2), ["cx"])
        if eps1 > 0:
            nm.add_all_qubit_quantum_error(depolarizing_error(depolarizing_lambda(eps1, 1), 1), ["sx", "x"])
    elif kind == "dephasing":
        if eps2 > 0:
            flip = pauli_error([("Z", eps2), ("I", 1.0 - eps2)])
            nm.add_all_qubit_quantum_error(flip.tensor(flip), ["cx"])
        if eps1 > 0:
            nm.add_all_qubit_quantum_error(pauli_error([("Z", eps1), ("I", 1.0 - eps1)]), ["sx", "x"])
    else:
        raise ValueError("kind must be 'depolarizing' or 'dephasing'")
    return nm


# ----------------------------------------------------------------------------
# Exact outcome probabilities
# ----------------------------------------------------------------------------


def prepare_for_probabilities(qc: QuantumCircuit) -> QuantumCircuit:
    """
    Compile to the hardware basis and append an instruction that saves the exact
    Z-basis probabilities. Do this once per circuit; then reuse it across a
    sweep of noise strengths with `probabilities`, since compilation is the slow
    part and does not depend on the noise.
    """
    from qiskit_aer.library import SaveProbabilities

    hw = to_hardware_basis(qc)
    hw.append(SaveProbabilities(hw.num_qubits), hw.qubits)
    return hw


def probabilities(prepared: QuantumCircuit, noise_model=None) -> np.ndarray:
    """
    Exact outcome probabilities of a prepared circuit under a noise model.

    Ordering: Aer indexes with qubit 0 as the least significant bit. Our site j
    is qubit N-1-j, so site N-1 is least significant for us as well and the
    index b is the same basis state in both conventions. This is the same
    coincidence circuits.statevector relies on, and it is checked numerically
    in tests/test_noise.py (noise-free probabilities equal |statevector|^2).
    """
    from qiskit_aer import AerSimulator

    sim = AerSimulator(method="density_matrix", noise_model=noise_model)
    data = sim.run(prepared, shots=1).result().data()
    p = np.asarray(data["probabilities"], dtype=float)
    return p / p.sum()


def noisy_probabilities(qc: QuantumCircuit, noise_model=None) -> np.ndarray:
    """Convenience: prepare and run in one call."""
    return probabilities(prepare_for_probabilities(qc), noise_model)


# ----------------------------------------------------------------------------
# Logical-frame noise, once per Trotter step, in plain NumPy
# ----------------------------------------------------------------------------


def _pauli_on(op: str, site: int, n_sites: int) -> np.ndarray:
    from .paulis import site_op
    return site_op(op, site, n_sites)


def apply_step_channel(rho: np.ndarray, channel: str, p: float, n_sites: int) -> np.ndarray:
    """
    One application of a logical-frame error channel to a density matrix.

    "dephasing"    each qubit: rho -> (1-p) rho + p Z rho Z
                   (Martinez et al.'s phenomenological model, p per qubit per step)
    "depolarizing" each qubit: rho -> (1-p) rho + (p/3)(X rho X + Y rho Y + Z rho Z)
    "global"       whole register: rho -> (1-p) rho + p I/2^N
    """
    if p == 0:
        return rho
    dim = 2 ** n_sites
    if channel == "global":
        return (1 - p) * rho + p * np.trace(rho) * np.eye(dim) / dim
    for j in range(n_sites):
        if channel == "dephasing":
            z = _pauli_on("Z", j, n_sites)
            rho = (1 - p) * rho + p * z @ rho @ z
        elif channel == "depolarizing":
            acc = (1 - p) * rho
            for op in ("X", "Y", "Z"):
                s = _pauli_on(op, j, n_sites)
                acc = acc + (p / 3) * s @ rho @ s
            rho = acc
        else:
            raise ValueError(f"unknown channel {channel!r}")
    return rho


def step_channel_probabilities(n_sites: int, m: float, w: float, J: float, time: float,
                               steps: int, order: int, channel: str, p: float) -> np.ndarray:
    """
    Z-basis probabilities after `steps` Trotter steps from the bare vacuum, with
    `channel` applied to the logical state after every step.

    Built from the step unitary of circuits.trotter_step and NumPy matrix
    products only, with no Aer involved, so it is also an independent check on
    the Aer pipeline: at p = 0 the two must agree, and the tests confirm they do.
    The step unitary comes back from Qiskit's Operator in Qiskit's ordering,
    which coincides with ours for the reason given in `probabilities`.
    """
    from qiskit.quantum_info import Operator

    from .circuits import trotter_step
    from .schwinger import bare_vacuum

    u = Operator(trotter_step(n_sites, m, w, J, time / steps, order)).data
    psi = bare_vacuum(n_sites)
    rho = np.outer(psi, psi.conj())
    for _ in range(steps):
        rho = apply_step_channel(u @ rho @ u.conj().T, channel, p, n_sites)
    probs = np.real(np.diag(rho)).copy()
    return probs / probs.sum()


# ----------------------------------------------------------------------------
# Readout error, applied classically
# ----------------------------------------------------------------------------


def apply_readout_error(probs: np.ndarray, p01: float, p10: float | None = None) -> np.ndarray:
    """
    Push an exact distribution through independent per-qubit readout errors.

    p01 = P(read 1 | true 0), p10 = P(read 0 | true 1); symmetric if p10 is
    omitted. Each qubit's 2 x 2 confusion matrix acts on its own axis of the
    distribution reshaped to (2, 2, ..., 2). This is done classically because
    readout error acts only on the classical outcome, after the state is gone;
    doing it here keeps the readout bias exactly separable from the gate bias.
    Checked against Aer's own ReadoutError by sampling in the tests.
    """
    if p10 is None:
        p10 = p01
    p = np.asarray(probs, dtype=float)
    n = int(round(np.log2(p.size)))
    conf = np.array([[1.0 - p01, p10],
                     [p01, 1.0 - p10]])          # conf[measured, true]
    tensor = p.reshape((2,) * n)
    for axis in range(n):
        tensor = np.moveaxis(np.tensordot(conf, tensor, axes=([1], [axis])), 0, axis)
    return tensor.reshape(-1)


# ----------------------------------------------------------------------------
# Shots
# ----------------------------------------------------------------------------


def sample_counts(probs: np.ndarray, shots: int, rng: np.random.Generator) -> np.ndarray:
    """
    Draw `shots` outcomes from an exact distribution: a multinomial sample.

    Sampling from the exact distribution is statistically identical to running
    the noisy circuit `shots` times (each shot of a circuit is an independent
    draw from the same distribution), and it lets the shot noise be switched on
    and off without rerunning the circuit. Seeded for reproducibility.
    """
    p = np.clip(np.asarray(probs, dtype=float), 0.0, None)
    return rng.multinomial(shots, p / p.sum()).astype(float)
