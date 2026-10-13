"""
The error budget (M10 to M12): every source separated, on every observable.

THE DECOMPOSITION
-----------------
For a diagonal observable O at time t, four outcome distributions are computed
EXACTLY (no sampling), each differing from the previous one by one source:

    exact     |<b| e^{-iHt} |vac>|^2           exact diagonalisation
    trotter   |<b| U_trotter |vac>|^2           noiseless circuit
    gate      diag( rho_noisy )                 same circuit, noisy gates
    measured  readout confusion applied to gate

The biases then telescope, so they add up exactly with no cross term:

    O_measured - O_exact = (O_trotter - O_exact)      Trotter bias
                         + (O_gate - O_trotter)       gate-error bias
                         + (O_measured - O_gate)      readout bias

and the shot noise is the standard error an M-shot run would have, computed
from the measured distribution:  sigma = sqrt( Var_per_shot / M ).
The root-mean-square error of a single M-shot measurement is then
sqrt( total_bias^2 + sigma^2 ), the quantity that has an optimum in dt.

Note what the ordering means. The gate bias is defined on the Trotterised
circuit, so it is "the extra error from noise, given this circuit". That is the
operationally relevant split: an experimentalist chooses dt first and then
lives with the hardware.

THE FAILURE THRESHOLD
---------------------
Pair production from the bare vacuum shows a rise and then a fall of the
particle density (Martinez et al.: "an initial phase of rapid pair creation is
followed by a reduction ... due to recombination"). Incoherent noise can only
drive the density monotonically toward the maximally mixed value 1/2, so the
FALL, the recombination dip

    D = nu(t_peak) - nu(t_trough),

is a signature noise cannot fake. The signal is lost when D is no longer
resolvable from shot noise, D <= z sigma_D with sigma_D = sqrt(s_peak^2 +
s_trough^2), z = 2. The failure threshold eps* is the two-qubit gate error at
which that first happens.

PREDICTING THE THRESHOLD FROM THE DIAGNOSTIC
--------------------------------------------
Given measured distributions at one noise level eps0 (all an experiment has),
read the fidelity at each time from charge conservation, F = 1 - Var(Q)/(N/4),
assume the global-depolarizing model and F(eps) = exp(-kappa eps), and invert:

    p_ideal = (p_measured - (1-F)/2^N) / F,       kappa = -ln F / eps0.

Re-noising p_ideal at any other eps predicts D and sigma_D there, and hence a
predicted threshold. Nothing in this uses the exact answer. Comparing it with
the threshold found by sweeping eps directly is the test of whether the
symmetry violation predicts where the physics is lost.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .circuits import circuit_metrics, trotter_circuit
from .diagnostics import (
    fidelity_from_charge_variance, global_depolarize, observables_from_distribution,
    shot_noise_sigma, undo_global_depolarizing,
)
from .noise import apply_readout_error, gate_noise_model, prepare_for_probabilities, probabilities
from .observables import evolve
from .schwinger import bare_vacuum, hamiltonian


# ----------------------------------------------------------------------------
# The four distributions and the decomposition
# ----------------------------------------------------------------------------


def exact_distribution(n_sites: int, m: float, w: float, J: float, time: float) -> np.ndarray:
    psi = evolve(bare_vacuum(n_sites), hamiltonian(n_sites, m, w, J), np.array([time]))[0]
    return np.abs(psi) ** 2


@dataclass
class PointCircuit:
    """A Trotter circuit compiled once, reusable across noise strengths."""

    time: float
    steps: int
    order: int
    prepared: object
    two_qubit_gates: int
    exact: np.ndarray
    trotter: np.ndarray


def point_circuit(n_sites: int, m: float, w: float, J: float, time: float, steps: int,
                  order: int = 2) -> PointCircuit:
    qc = trotter_circuit(n_sites, m, w, J, time, steps, order)
    prep = prepare_for_probabilities(qc)
    return PointCircuit(time, steps, order, prep, circuit_metrics(qc)["two_qubit_gates"],
                        exact_distribution(n_sites, m, w, J, time), probabilities(prep, None))


@dataclass
class Budget:
    exact: float
    trotter: float
    gate: float
    measured: float
    sigma: float

    @property
    def trotter_bias(self) -> float:
        return self.trotter - self.exact

    @property
    def gate_bias(self) -> float:
        return self.gate - self.trotter

    @property
    def readout_bias(self) -> float:
        return self.measured - self.gate

    @property
    def total_bias(self) -> float:
        return self.measured - self.exact

    @property
    def rms(self) -> float:
        return float(np.hypot(self.total_bias, self.sigma))


def _value(probs: np.ndarray, n_sites: int, observable: str) -> float:
    obs = observables_from_distribution(probs, n_sites)
    if observable.startswith("link"):
        return float(obs["links"][int(observable[4:]) - 1])
    return float(obs[observable])


def budget(pc: PointCircuit, n_sites: int, noise_model, readout: float, shots: int,
           observable: str = "particle_density") -> Budget:
    """The four-way decomposition of one diagonal observable at one point."""
    gate = probabilities(pc.prepared, noise_model)
    measured = apply_readout_error(gate, readout) if readout > 0 else gate
    return Budget(_value(pc.exact, n_sites, observable),
                  _value(pc.trotter, n_sites, observable),
                  _value(gate, n_sites, observable),
                  _value(measured, n_sites, observable),
                  shot_noise_sigma(measured, n_sites, shots, observable))


# ----------------------------------------------------------------------------
# The recombination dip and the failure threshold
# ----------------------------------------------------------------------------


def dip(p_peak: np.ndarray, p_trough: np.ndarray, n_sites: int, shots: int) -> tuple[float, float]:
    """D = nu(t_peak) - nu(t_trough) and its shot-noise standard error."""
    d = (observables_from_distribution(p_peak, n_sites)["particle_density"]
         - observables_from_distribution(p_trough, n_sites)["particle_density"])
    s = np.hypot(shot_noise_sigma(p_peak, n_sites, shots), shot_noise_sigma(p_trough, n_sites, shots))
    return float(d), float(s)


def first_crossing(x: np.ndarray, y: np.ndarray) -> float:
    """
    Smallest x at which y changes sign from positive to non-positive, by linear
    interpolation in log x. Returns inf if y never crosses on the grid.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    for i in range(1, len(x)):
        if y[i - 1] > 0 >= y[i]:
            lx0, lx1 = np.log(x[i - 1]), np.log(x[i])
            return float(np.exp(lx0 + (lx1 - lx0) * y[i - 1] / (y[i - 1] - y[i])))
    return float("inf")


def sweep_dip(peak: PointCircuit, trough: PointCircuit, n_sites: int, eps_grid: np.ndarray,
              shots: int, kind: str = "depolarizing") -> dict[str, np.ndarray]:
    """
    Run both circuits across a grid of two-qubit error rates and record the dip,
    its noise floor and the diagnostics at both times.
    """
    rows = {k: [] for k in ("D", "sigma_D", "F_peak", "F_trough", "varQ_peak", "varQ_trough",
                            "cp_peak", "cp_trough")}
    dists = []
    for eps in eps_grid:
        nm = gate_noise_model(eps, kind=kind)
        pp, pt = probabilities(peak.prepared, nm), probabilities(trough.prepared, nm)
        d, s = dip(pp, pt, n_sites, shots)
        op, ot = observables_from_distribution(pp, n_sites), observables_from_distribution(pt, n_sites)
        rows["D"].append(d)
        rows["sigma_D"].append(s)
        rows["varQ_peak"].append(op["charge_variance"])
        rows["varQ_trough"].append(ot["charge_variance"])
        rows["F_peak"].append(fidelity_from_charge_variance(op["charge_variance"], n_sites))
        rows["F_trough"].append(fidelity_from_charge_variance(ot["charge_variance"], n_sites))
        rows["cp_peak"].append(op["cp_asymmetry"])
        rows["cp_trough"].append(ot["cp_asymmetry"])
        dists.append((pp, pt))
    out = {k: np.array(v) for k, v in rows.items()}
    out["eps"] = np.asarray(eps_grid, dtype=float)
    out["distributions"] = dists
    return out


def predicted_dip(p_peak_meas: np.ndarray, p_trough_meas: np.ndarray, eps0: float,
                  eps_grid: np.ndarray, n_sites: int, shots: int) -> dict[str, np.ndarray]:
    """
    Predict D(eps) and sigma_D(eps) everywhere from data taken at ONE noise level
    eps0, using only the charge diagnostic and the global-depolarizing model.
    See the module docstring. Works on exact or shot-sampled distributions.

    The re-noised distributions can carry small negative entries when F was
    estimated from shots (see diagnostics.undo_global_depolarizing). The dip is
    computed from them unclipped, since it is linear in the distribution; only
    the shot-noise estimate, which needs a genuine distribution, uses the
    clipped version.
    """
    preds = []
    for p in (p_peak_meas, p_trough_meas):
        var_q = observables_from_distribution(p, n_sites)["charge_variance"]
        f0 = min(max(fidelity_from_charge_variance(var_q, n_sites), 1e-12), 1.0)
        kappa = -np.log(f0) / eps0
        ideal = undo_global_depolarizing(p, f0)
        preds.append([global_depolarize(ideal, np.exp(-kappa * e)) for e in eps_grid])
    d, s = [], []
    for pp, pt in zip(*preds):
        d.append(observables_from_distribution(pp, n_sites)["particle_density"]
                 - observables_from_distribution(pt, n_sites)["particle_density"])
        s.append(dip(np.clip(pp, 0, None), np.clip(pt, 0, None), n_sites, shots)[1])
    return {"eps": np.asarray(eps_grid, dtype=float), "D": np.array(d), "sigma_D": np.array(s)}
