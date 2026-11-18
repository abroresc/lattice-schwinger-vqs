"""
Qiskit layer: build the Schwinger Hamiltonian as circuits and simulate them.

This is the bridge from exact linear algebra (everything else in this package)
to something a quantum computer could run: state preparation, Trotter
evolution, compilation and shots.

Every function here is checked against the exact NumPy result in
tests/test_circuits.py. The circuit layer is not trusted just because Qiskit
ran without an error.

The index trap, restated because circuits are where it matters.
Our site j (0-indexed, leftmost in a Pauli string) is Qiskit qubit N-1-j.
Every gate index below goes through `site_to_qiskit_qubit`. Pass a site index
straight to a Qiskit gate and the chain silently reverses, which for this model
swaps particles and antiparticles.
"""

from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp, Statevector

from .paulis import site_to_qiskit_qubit
from .schwinger import bare_vacuum_label, hamiltonian_terms, stagger


def hamiltonian_op(n_sites: int, m: float, w: float, J: float) -> SparsePauliOp:
    """
    The Hamiltonian as a Qiskit SparsePauliOp.

    Label strings pass through unchanged: Qiskit reads a label left to right as
    most significant tensor factor first, exactly as `pauli_string` does. Only
    the qubit NUMBERING is reversed, and that only matters for gates.
    """
    return SparsePauliOp.from_list([(label, coeff) for label, coeff in
                                    hamiltonian_terms(n_sites, m, w, J)])


def bare_vacuum_circuit(n_sites: int) -> QuantumCircuit:
    """
    Prepare the bare vacuum (filled Dirac sea) with one layer of X gates.

    The target is Z_n = -(-1)^n, which is |0101...> in our site order: odd
    physics sites occupied, even sites empty. Preparing it costs N/2 X gates
    and no entanglement, which is exactly why it is the natural starting state
    for a real-time simulation. The true ground state is entangled; skqd.py finds it.
    """
    qc = QuantumCircuit(n_sites, name="bare_vacuum")
    for j, bit in enumerate(bare_vacuum_label(n_sites)):
        if bit == "1":
            qc.x(site_to_qiskit_qubit(j, n_sites))
    return qc


# ---------------------------------------------------------------------------
# Trotterisation
# ---------------------------------------------------------------------------
#
# Split H into three groups. Within each group all terms commute, so each
# group exponentiates exactly and the only error comes from the splitting.
#
#   H_Z     all Z and ZZ terms. Diagonal, so they all commute with each other.
#   H_even  hopping on bonds (0,1), (2,3), ...   disjoint, so they commute.
#   H_odd   hopping on bonds (1,2), (3,4), ...   disjoint, so they commute.
#
# Splitting the hopping into even and odd bonds is the standard move: adjacent
# bonds share a site and do not commute, so they cannot go in one group.
# ---------------------------------------------------------------------------


def _z_coefficients(n_sites: int, m: float, w: float, J: float
                    ) -> tuple[dict[int, float], dict[tuple[int, int], float]]:
    """Pull the single-Z and ZZ coefficients out of the term list."""
    singles: dict[int, float] = {}
    pairs: dict[tuple[int, int], float] = {}
    for label, coeff in hamiltonian_terms(n_sites, m, w, J):
        sites = [i for i, c in enumerate(label) if c != "I"]
        kinds = {label[i] for i in sites}
        if kinds and kinds <= {"Z"}:
            if len(sites) == 1:
                singles[sites[0]] = singles.get(sites[0], 0.0) + coeff.real
            elif len(sites) == 2:
                key = (sites[0], sites[1])
                pairs[key] = pairs.get(key, 0.0) + coeff.real
            else:
                raise ValueError(f"unexpected Z string of weight {len(sites)}: {label}")
    return singles, pairs


def _append_diagonal(qc: QuantumCircuit, n_sites: int, m: float, w: float,
                     J: float, dt: float) -> None:
    """exp(-i dt H_Z), exactly. RZ(theta) = exp(-i theta/2 Z), so theta = 2 c dt."""
    singles, pairs = _z_coefficients(n_sites, m, w, J)
    for site, c in singles.items():
        qc.rz(2.0 * c * dt, site_to_qiskit_qubit(site, n_sites))
    for (a, b), c in pairs.items():
        qc.rzz(2.0 * c * dt, site_to_qiskit_qubit(a, n_sites),
               site_to_qiskit_qubit(b, n_sites))


def _append_hopping(qc: QuantumCircuit, n_sites: int, w: float, dt: float,
                    parity: int) -> None:
    """
    exp(-i dt H_hop) on the even (parity 0) or odd (parity 1) bonds, exactly.

    Each bond term is (w/2)(X_n X_{n+1} + Y_n Y_{n+1}). RXX(theta) =
    exp(-i theta/2 XX), so theta = w dt reproduces exp(-i dt (w/2) XX). XX and
    YY commute on the same bond, so the two gates together are exact for that
    bond, with no Trotter error inside it.
    """
    for j in range(parity, n_sites - 1, 2):
        qa = site_to_qiskit_qubit(j, n_sites)
        qb = site_to_qiskit_qubit(j + 1, n_sites)
        qc.rxx(w * dt, qa, qb)
        qc.ryy(w * dt, qa, qb)


def trotter_step(n_sites: int, m: float, w: float, J: float, dt: float,
                 order: int = 2) -> QuantumCircuit:
    """
    One Trotter step of exp(-i H dt).

    order = 1:  e^{-i H_odd dt} e^{-i H_even dt} e^{-i H_Z dt}
                as an operator: the diagonal layer is applied first in time, so
                it is the rightmost factor (checked in tests/test_circuits.py).
                Error per step O(dt^2), so total error over fixed time O(dt).

    order = 2:  symmetric (Strang) splitting, half steps on the outside.
                error per step O(dt^3), total error O(dt^2).

    The measured slope on log-log axes is the test of whether the
    implementation matches the theory. See scripts/05_qiskit_trotter.py.
    """
    qc = QuantumCircuit(n_sites, name=f"trotter{order}")
    if order == 1:
        _append_diagonal(qc, n_sites, m, w, J, dt)
        _append_hopping(qc, n_sites, w, dt, 0)
        _append_hopping(qc, n_sites, w, dt, 1)
    elif order == 2:
        _append_diagonal(qc, n_sites, m, w, J, 0.5 * dt)
        _append_hopping(qc, n_sites, w, 0.5 * dt, 0)
        _append_hopping(qc, n_sites, w, dt, 1)
        _append_hopping(qc, n_sites, w, 0.5 * dt, 0)
        _append_diagonal(qc, n_sites, m, w, J, 0.5 * dt)
    else:
        raise ValueError("order must be 1 or 2")
    return qc


def trotter_circuit(n_sites: int, m: float, w: float, J: float, time: float,
                    steps: int, order: int = 2,
                    initial: QuantumCircuit | None = None) -> QuantumCircuit:
    """Full evolution circuit: state preparation followed by `steps` Trotter steps."""
    if steps < 1:
        raise ValueError("steps must be at least 1")
    qc = QuantumCircuit(n_sites)
    qc.compose(initial if initial is not None else bare_vacuum_circuit(n_sites),
               inplace=True)
    step = trotter_step(n_sites, m, w, J, time / steps, order)
    for _ in range(steps):
        qc.compose(step, inplace=True)
    return qc


# ---------------------------------------------------------------------------
# Reading results back
# ---------------------------------------------------------------------------


def statevector(qc: QuantumCircuit) -> np.ndarray:
    """
    Circuit output as a plain NumPy vector in this package's basis ordering.

    Qiskit's Statevector indexes with qubit 0 as the least significant bit.
    Since our site j maps to qubit N-1-j, site N-1 is also least significant
    for us, so the two orderings coincide and no permutation is needed. That
    is verified numerically in tests/test_circuits.py rather than assumed.
    """
    return np.asarray(Statevector(qc).data, dtype=complex)


def particle_density_from_counts(counts: dict[str, int], n_sites: int
                                 ) -> tuple[float, float]:
    """
    Estimate particle density and its shot-noise uncertainty from measured counts.

    Why counts are enough: the particle number operator

        N_part = (1/2) sum_n [ I + (-1)^n Z_n ]

    is diagonal in the computational basis, so a single Z-basis measurement
    gives it. No basis rotations, no separate circuits per term. Energies, as VQE measures them,
    need several measurement bases; this observable needs one.

    The uncertainty is the standard error of the mean over shots, computed from
    the sample variance of the per-shot particle number. That is the right
    estimator here: the per-shot value is an integer 0..N, not a Bernoulli variable,
    so a naive sqrt(p(1-p)/n) would be wrong.

    Qiskit bitstrings come back with qubit 0 last, which is our site N-1 last,
    so the string reads left to right as site 0, site 1, ... and needs no
    reversal. Verified in tests.
    """
    total = sum(counts.values())
    if total == 0:
        raise ValueError("no shots")
    values = []
    weights = []
    for bits, n in counts.items():
        clean = bits.replace(" ", "")
        if len(clean) != n_sites:
            raise ValueError(f"bitstring {bits!r} has length {len(clean)}, expected {n_sites}")
        # site j occupied when its bit is '0' (|0> is the Z = +1 state)
        z = np.array([1.0 if c == "0" else -1.0 for c in clean])
        npart = 0.5 * sum(1.0 + stagger(j) * z[j] for j in range(n_sites))
        values.append(npart)
        weights.append(n)
    values = np.array(values)
    weights = np.array(weights, dtype=float)
    mean = float(np.sum(values * weights) / total)
    var = float(np.sum(weights * (values - mean) ** 2) / total)
    return mean / n_sites, np.sqrt(var / total) / n_sites


def run_with_shots(qc: QuantumCircuit, shots: int, seed: int | None = None,
                   noise_model=None) -> dict[str, int]:
    """
    Sample the circuit on the Aer simulator and return measurement counts.

    `noise_model` is None by default; pass one from noise.gate_noise_model to
    sample a noisy circuit.
    """
    from qiskit_aer import AerSimulator

    meas = qc.copy()
    meas.measure_all()
    backend = AerSimulator(noise_model=noise_model)
    from qiskit import transpile
    result = backend.run(transpile(meas, backend), shots=shots, seed_simulator=seed).result()
    return result.get_counts()


# The gate set of an IBM-style superconducting device. RZ is implemented in
# software (a frame change) and is error-free on such hardware; SX and X are
# physical single-qubit pulses; CX is the only entangling gate. Noise models in
# noise.py attach errors to exactly these gates, so the circuit has to be
# compiled into this basis before noise is applied. Otherwise Aer would run
# RXX/RYY/RZZ natively and the noise model would never fire.
HARDWARE_BASIS: tuple[str, ...] = ("rz", "sx", "x", "cx")


def to_hardware_basis(qc: QuantumCircuit, basis: tuple[str, ...] = HARDWARE_BASIS
                      ) -> QuantumCircuit:
    """
    Compile to the hardware basis. One function, so that the circuit whose gates
    are counted (circuit_metrics) is the same circuit whose gates carry noise
    (noise.py). A fixed seed makes the compilation reproducible.
    """
    from qiskit import transpile

    return transpile(qc, basis_gates=list(basis), optimization_level=1, seed_transpiler=0)


def circuit_metrics(qc: QuantumCircuit, basis: tuple[str, ...] = HARDWARE_BASIS
                    ) -> dict[str, int]:
    """
    Depth and two-qubit gate count after transpiling to a hardware-like basis.

    Gate error grows with the two-qubit gate count, so shrinking the Trotter
    step lowers Trotter error and raises hardware error at the same time.
    budget.py measures that trade-off.
    """
    t = to_hardware_basis(qc, basis)
    ops = t.count_ops()
    return {
        "depth": t.depth(),
        "two_qubit_gates": int(ops.get("cx", 0) + ops.get("ecr", 0) + ops.get("cz", 0)),
        "total_gates": int(sum(ops.values())),
    }
