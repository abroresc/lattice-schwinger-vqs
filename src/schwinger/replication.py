"""
Replication of Martinez et al., Nature 534, 516 (2016), Fig. 2(a).

Why this module exists
----------------------
The tests elsewhere prove that the Hamiltonian is the one Martinez et al. print,
term by term. That is not the same as showing that this code computes the same
physics they did. Here their published curves are the known answer: at their
parameters (N = 4, J = w, m = 0.5 w, start in the bare vacuum) the code must
reproduce, point by point,

  - their "ideal evolution", the exact particle density nu(t);
  - their "discretization errors" curve, the same evolution in first-order
    Trotter steps of their step length;
  - their "experimental error model", that Trotter evolution with their fitted
    dephasing, a phase flip with probability p = 0.038 on every qubit after
    every step.

Where the numbers come from
---------------------------
They were read off Fig. 2(a) of the arXiv v1 copy (arXiv:1605.04570) by pixel
analysis, not by eye: the page rendered at 576 dpi, both axes calibrated on
their tick marks (largest calibration residual 2e-4), each straight segment of a
plotted curve fitted away from the markers and evaluated at the step times, and
the centroids of the small circles taken for the measured data. The two
segments that meet at a plotted point give the same value to 5e-4, so the
reading uncertainty is about 0.002 in nu, set by the line width.

Martinez et al. never state their step length. Their five plotted times are
equally spaced from wt = 0 to pi/2, so it is pi/8; the agreement below is the
confirmation.

The trap: gate order inside a step
----------------------------------
Their pulse sequence (Extended Data Table I) applies the hopping on bonds
(1,2) and (3,4), then on (2,3), then the diagonal terms. The first-order step
here applies the diagonal terms first, then the same hopping order. From the
bare vacuum, which the diagonal terms only multiply by a phase, the two orders
give the same state up to a final diagonal factor, and so the same nu at every
step. `martinez_pulse_order_density` builds their order independently, with
scipy's expm, so the tests can check that equivalence rather than assert it.
"""
from __future__ import annotations

import numpy as np

# Their parameters (Fig. 2 caption: J = w, m = 0.5 w), with w = 1.
N_SITES, M, W, J = 4, 0.5, 1.0, 1.0
STEP = np.pi / 8          # w T, inferred from the plotted times (see above)
P_DEPHASING = 0.038       # Methods, "Experimental error model", nonzero J
READ_UNCERTAINTY = 0.002  # in nu, from the pixel analysis

# nu after k = 1, 2, 3, 4 steps (wt = pi/8, pi/4, 3pi/8, pi/2), read off Fig. 2(a)
FIG2A_TROTTER = (0.199, 0.471, 0.424, 0.225)   # red triangles
FIG2A_MODEL = (0.198, 0.448, 0.433, 0.312)     # purple squares
FIG2A_DATA = (0.184, 0.439, 0.447, 0.337)      # measured, postselected (circles)

# Their ideal curve (green), sampled where no marker covers it: (wt, nu)
FIG2A_IDEAL = (
    (0.0291, 0.0013), (0.0733, 0.0080), (0.1176, 0.0204), (0.1618, 0.0379),
    (0.2060, 0.0606), (0.2503, 0.0877), (0.2945, 0.1179), (0.3387, 0.1511),
    (0.3830, 0.1876), (0.4714, 0.2590), (0.5157, 0.2945), (0.5599, 0.3282),
    (0.6041, 0.3601), (0.6484, 0.3883), (0.6926, 0.4127), (0.7590, 0.4406),
    (0.8253, 0.4578), (0.8696, 0.4636), (0.9138, 0.4619), (0.9580, 0.4570),
    (1.0023, 0.4478), (1.0686, 0.4237), (1.1129, 0.4041), (1.1571, 0.3814),
    (1.2013, 0.3568), (1.2456, 0.3307), (1.2898, 0.3044), (1.3340, 0.2779),
    (1.3783, 0.2518), (1.4225, 0.2275), (1.4667, 0.2053), (1.5110, 0.1857),
    (1.5552, 0.1685),
)


def _density_of(probs: np.ndarray) -> float:
    """nu from Z-basis probabilities: the particle-number diagonal, per site."""
    from .observables import particle_number_op
    return float(probs @ np.real(np.diag(particle_number_op(N_SITES)))) / N_SITES


def exact_density(times) -> np.ndarray:
    """Exact nu(t) from the bare vacuum, by exact diagonalisation."""
    from .observables import evolve, particle_density
    from .schwinger import bare_vacuum, hamiltonian

    states = evolve(bare_vacuum(N_SITES), hamiltonian(N_SITES, M, W, J), np.asarray(times))
    return np.array([particle_density(s, N_SITES) for s in states])


def trotter_density(steps: int = 4, step: float = STEP) -> np.ndarray:
    """nu after 1..steps first-order Trotter steps, using the circuit's own step."""
    from qiskit.quantum_info import Operator

    from .circuits import trotter_step
    from .observables import particle_density
    from .schwinger import bare_vacuum

    u = Operator(trotter_step(N_SITES, M, W, J, step, order=1)).data
    psi = bare_vacuum(N_SITES).astype(complex)
    out = []
    for _ in range(steps):
        psi = u @ psi
        out.append(particle_density(psi, N_SITES))
    return np.array(out)


def martinez_pulse_order_density(steps: int = 4, step: float = STEP) -> np.ndarray:
    """
    nu after 1..steps of their pulse order: bonds (1,2) and (3,4), then (2,3),
    then all diagonal terms. Built with scipy's expm from the Pauli terms, so it
    shares nothing with the circuit except the Hamiltonian's term list.
    """
    from scipy.linalg import expm

    from .observables import particle_density
    from .paulis import pauli_string
    from .schwinger import bare_vacuum, hamiltonian_terms

    terms = hamiltonian_terms(N_SITES, M, W, J)
    diag = sum(c * pauli_string(s) for s, c in terms if set(s) <= {"I", "Z"})
    outer = sum(c * pauli_string(s) for s, c in terms if s in ("XXII", "YYII", "IIXX", "IIYY"))
    middle = sum(c * pauli_string(s) for s, c in terms if s in ("IXXI", "IYYI"))
    # time order outer, middle, diagonal: as an operator, outer acts first (rightmost)
    u = expm(-1j * diag * step) @ expm(-1j * middle * step) @ expm(-1j * outer * step)
    psi = bare_vacuum(N_SITES).astype(complex)
    out = []
    for _ in range(steps):
        psi = u @ psi
        out.append(particle_density(psi, N_SITES))
    return np.array(out)


def error_model_density(steps: int = 4, step: float = STEP, p: float = P_DEPHASING) -> np.ndarray:
    """nu after 1..steps Trotter steps with their dephasing (p per qubit per step)."""
    from .noise import step_channel_probabilities

    return np.array([
        _density_of(step_channel_probabilities(N_SITES, M, W, J, k * step, k, 1, "dephasing", p))
        for k in range(1, steps + 1)
    ])
