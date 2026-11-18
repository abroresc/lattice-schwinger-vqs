"""
Observables and symmetry diagnostics read from Z-basis measurements.

Why one measurement basis is enough
-----------------------------------
Every quantity the error budget needs is diagonal in the computational basis:

    particle number   N_part = (1/2) sum_n [ 1 + (-1)^n Z_n ]
    charge density    rho_n  = ( Z_n + (-1)^n ) / 2
    electric field    L_n    = sum_{k<=n} rho_k          (Gauss's law, solved)
    total charge      Q      = sum_n rho_n

So each one is a function of the measured bitstring alone. One circuit, measured
once in Z, gives the physics observable and the diagnostics that say whether to
trust it. Nothing extra has to be run to get the diagnostics, which is the
property that makes them usable on a real device.

Everything here takes a probability vector over the 2^N basis states (exact,
from a density matrix) or a vector of counts (from shots). Index b of either
vector is the basis state whose bit string, read left to right, is site 1 ...
site N, exactly as in paulis.basis_label. Bit '0' means Z = +1, which under
this project's Jordan-Wigner convention is an occupied site.

The diagnostics
---------------
Gauss's law cannot be a diagnostic here: it was solved and substituted, so it
holds on every bitstring (see analytics.gauss_law_residual). Two exact
symmetries remain. Most noise breaks both, but not all noise does: see the
blind spot below.

1. Total charge. [H, Q] = 0 and the evolution starts in Q = 0, so ideally every
   shot has Q = 0. Two equivalent readings of a violation:
     - Var(Q), the spread of the measured charge;
     - 1 - P(Q = 0), the "leakage" Martinez et al. reported as a fidelity proxy.
2. CP. On a neutral state that CP leaves unchanged, <L_n> = <L_{N-n}>; the
   bare vacuum is such a state (see schwinger.cp_operator for why neutrality
   alone is not enough). The noiseless Trotter circuit keeps the state
   CP-symmetric exactly (each Trotter factor maps the neutral sector to itself
   and commutes with CP there), so any measured <L_1> - <L_3> is hardware
   error, not Trotter error.

The global-depolarizing reading
-------------------------------
If noise acted as one global depolarizing channel,

    rho_noisy = F rho_ideal + (1 - F) I / 2^N,

then every diagonal observable would interpolate linearly between its ideal value
and its maximally-mixed value with the same weight F. Since the ideal state has
Var(Q) = 0 and P(Q=0) = 1, F can be read off either diagnostic:

    F_Q   = 1 - Var(Q) / Var_mixed(Q),             Var_mixed(Q) = N / 4
    F_P0  = ( P(Q=0) - f_mixed ) / ( 1 - f_mixed ), f_mixed = C(N, N/2) / 2^N

Real gate noise is local, not global, so these are models, and their accuracy is
what the thesis tests. They are exact for the global channel, which is what the
tests in tests/test_diagnostics.py check.

The blind spot
--------------
Z errors commute with Q, so dephasing that acts on the logical state, between
Trotter steps (the error model Martinez et al. fitted to their data), leaves
Var(Q) exactly zero while the physics decays. The charge diagnostic cannot see
it at all. noise.step_channel_probabilities implements that model and the tests
confirm Var(Q) = 0 to machine precision.

Dephasing on the CX gates is different: it is visible (Var(Q) = 0.88 at 3e-3
per gate). Inside a compiled RXX or RYY the qubits sit in a rotated basis between the CXs, so a Z
error there is an X or Y error in the logical frame, and it does break charge
conservation. Whether a dephasing error is visible depends on where in the
compiled circuit it happens, not only on its type.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb

import numpy as np

from .paulis import basis_label
from .schwinger import stagger


# ----------------------------------------------------------------------------
# Per-bitstring values. Computed once per N and cached.
# ----------------------------------------------------------------------------


@dataclass(frozen=True)
class ZBasisTable:
    """
    The value of each diagonal observable on each computational basis state.

    Arrays are indexed by basis index b (length 2^N). `links` has shape
    (2^N, N-1): column n is L_{n+1}, the field on the link to the right of
    0-indexed site n.
    """

    n_sites: int
    particle_number: np.ndarray
    charge: np.ndarray
    links: np.ndarray

    @property
    def particle_density(self) -> np.ndarray:
        return self.particle_number / self.n_sites


_TABLES: dict[int, ZBasisTable] = {}


def z_basis_table(n_sites: int) -> ZBasisTable:
    """Evaluate N_part, Q and every L_n on every bit string, from the bits alone."""
    if n_sites in _TABLES:
        return _TABLES[n_sites]
    dim = 2 ** n_sites
    npart = np.zeros(dim)
    charge = np.zeros(dim)
    links = np.zeros((dim, max(n_sites - 1, 0)))
    s = np.array([stagger(j) for j in range(n_sites)], dtype=float)
    for b in range(dim):
        z = np.array([1.0 if c == "0" else -1.0 for c in basis_label(b, n_sites)])
        rho = 0.5 * (z + s)
        npart[b] = 0.5 * np.sum(1.0 + s * z)
        charge[b] = np.sum(rho)
        links[b] = np.cumsum(rho)[: n_sites - 1]
    table = ZBasisTable(n_sites, npart, charge, links)
    _TABLES[n_sites] = table
    return table


# ----------------------------------------------------------------------------
# Observables from an exact probability vector
# ----------------------------------------------------------------------------


def _check_distribution(probs: np.ndarray, n_sites: int) -> np.ndarray:
    p = np.asarray(probs, dtype=float)
    if p.shape != (2 ** n_sites,):
        raise ValueError(f"expected a length-{2 ** n_sites} vector, got shape {p.shape}")
    total = p.sum()
    if total <= 0:
        raise ValueError("distribution has no weight")
    return p / total


def observables_from_distribution(probs: np.ndarray, n_sites: int) -> dict[str, float | np.ndarray]:
    """
    Every diagonal observable and diagnostic, from a probability vector.

    Returns a dict with
      particle_density      <N_part>/N, the pair-production observable
      links                 <L_n> for n = 1..N-1
      charge_mean           <Q>
      charge_variance       Var(Q), diagnostic 1
      neutral_fraction      P(Q = 0), the postselection yield
      cp_asymmetry          <L_1> - <L_{N-1}>, diagnostic 2 (zero if CP holds)
      postselected_density  particle density conditioned on Q = 0
    """
    p = _check_distribution(probs, n_sites)
    t = z_basis_table(n_sites)
    q = t.charge
    neutral = np.isclose(q, 0.0)
    p0 = float(p[neutral].sum())
    links = p @ t.links if n_sites > 1 else np.zeros(0)
    out: dict[str, float | np.ndarray] = {
        "particle_density": float(p @ t.particle_density),
        "links": links,
        "charge_mean": float(p @ q),
        "charge_variance": float(p @ q**2 - (p @ q) ** 2),
        "neutral_fraction": p0,
        "cp_asymmetry": float(links[0] - links[-1]) if n_sites > 2 else 0.0,
        "postselected_density": (float(p[neutral] @ t.particle_density[neutral] / p0)
                                 if p0 > 0 else float("nan")),
    }
    return out


# ----------------------------------------------------------------------------
# Observables from finite shots, with their statistical uncertainty
# ----------------------------------------------------------------------------


def observables_from_counts(counts: np.ndarray, n_sites: int) -> dict[str, float | np.ndarray]:
    """
    The same observables estimated from shot counts, each with a standard error.

    `counts[b]` is the number of shots that returned basis state b. Means are
    the sample means; uncertainties are standard errors of those means, from
    the sample variance of the per-shot value (Taylor, ch. 4-5). For Var(Q) the
    standard error of a sample variance is used,

        SE[s^2] ~ sqrt( (mu_4 - s^4) / M ),

    with mu_4 the fourth central moment. This is the large-M formula; at the
    shot counts used here (M >= 256) it is adequate, and the tests compare it
    against the empirical scatter over repeated runs rather than trusting it.
    """
    c = np.asarray(counts, dtype=float)
    if c.shape != (2 ** n_sites,):
        raise ValueError(f"expected a length-{2 ** n_sites} vector, got shape {c.shape}")
    shots = c.sum()
    if shots < 2:
        raise ValueError("need at least two shots")
    out = observables_from_distribution(c, n_sites)
    p = c / shots
    t = z_basis_table(n_sites)

    def sem(values: np.ndarray) -> float:
        mean = p @ values
        var = p @ (values - mean) ** 2 * shots / (shots - 1)   # unbiased sample variance
        return float(np.sqrt(var / shots))

    q = t.charge
    dq = q - p @ q
    mu2 = p @ dq**2
    mu4 = p @ dq**4
    out["shots"] = int(shots)
    out["particle_density_err"] = sem(t.particle_density)
    out["links_err"] = np.array([sem(t.links[:, n]) for n in range(t.links.shape[1])])
    out["charge_variance_err"] = float(np.sqrt(max(mu4 - mu2**2, 0.0) / shots))
    p0 = out["neutral_fraction"]
    out["neutral_fraction_err"] = float(np.sqrt(p0 * (1 - p0) / shots))
    out["cp_asymmetry_err"] = (sem(t.links[:, 0] - t.links[:, -1]) if n_sites > 2 else 0.0)
    return out


def shot_noise_sigma(probs: np.ndarray, n_sites: int, shots: int,
                     observable: str = "particle_density") -> float:
    """
    The standard error an M-shot measurement will have, computed exactly from
    the probability distribution rather than estimated from one sample.

    This is what makes the error budget clean: shot noise is a number computed
    from the distribution, not a fluctuation that has to be averaged away.
    sigma = sqrt( Var_per_shot / M ).
    """
    p = _check_distribution(probs, n_sites)
    t = z_basis_table(n_sites)
    if observable == "particle_density":
        values = t.particle_density
    elif observable == "cp_asymmetry":
        values = t.links[:, 0] - t.links[:, -1]
    elif observable.startswith("link"):
        values = t.links[:, int(observable[4:]) - 1]
    else:
        raise ValueError(f"unknown observable {observable!r}")
    mean = p @ values
    return float(np.sqrt(p @ (values - mean) ** 2 / shots))


# ----------------------------------------------------------------------------
# The global-depolarizing reading of the diagnostics
# ----------------------------------------------------------------------------


def mixed_state_values(n_sites: int) -> dict[str, float]:
    """
    Var(Q), P(Q=0) and particle density of the maximally mixed state, computed
    from the table rather than typed in, so they cannot drift from the
    definitions. For even N these equal N/4, C(N,N/2)/2^N and 1/2.
    """
    uniform = np.full(2 ** n_sites, 1.0 / 2 ** n_sites)
    obs = observables_from_distribution(uniform, n_sites)
    return {"charge_variance": obs["charge_variance"],
            "neutral_fraction": obs["neutral_fraction"],
            "particle_density": obs["particle_density"]}


def fidelity_from_charge_variance(charge_variance: float, n_sites: int) -> float:
    """F_Q = 1 - Var(Q)/Var_mixed(Q). Exact for a global depolarizing channel."""
    return 1.0 - charge_variance / mixed_state_values(n_sites)["charge_variance"]


def fidelity_from_neutral_fraction(neutral_fraction: float, n_sites: int) -> float:
    """F_P0 = (P(Q=0) - f_mixed)/(1 - f_mixed). Exact for a global depolarizing channel."""
    f_mixed = mixed_state_values(n_sites)["neutral_fraction"]
    return (neutral_fraction - f_mixed) / (1.0 - f_mixed)


def neutral_sector_size(n_sites: int) -> int:
    """Number of bitstrings with Q = 0. Six at N = 4: the 'physical subspace'."""
    return comb(n_sites, n_sites // 2) if n_sites % 2 == 0 else 0


def global_depolarize(probs: np.ndarray, fidelity: float) -> np.ndarray:
    """F p + (1 - F)/2^N: what a global depolarizing channel does to a diagonal."""
    p = np.asarray(probs, dtype=float)
    return fidelity * p + (1.0 - fidelity) / p.size


def undo_global_depolarizing(probs: np.ndarray, fidelity: float) -> np.ndarray:
    """
    Invert global_depolarize: p_ideal = (p - (1-F)/2^N) / F.

    This is the mitigation of Vovrosh et al. and Urbanek et al. applied to the
    whole distribution. With an estimated F it can return small negative
    entries; they are kept, not clipped, because clipping would bias every
    observable computed from the result.
    """
    p = np.asarray(probs, dtype=float)
    if fidelity <= 0:
        raise ValueError("fidelity must be positive to invert")
    return (p - (1.0 - fidelity) / p.size) / fidelity
