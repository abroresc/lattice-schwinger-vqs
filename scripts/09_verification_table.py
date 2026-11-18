"""
Every row of the thesis verification table (Table 2), recomputed and printed.

Why this script exists. The tests in tests/ check each of these quantities
against a tolerance. That is right for a test and wrong for a table: a tolerance
says how much disagreement would be accepted, not how much there is. This script
prints the measured disagreement for every row, so each number in the table has
one source.

All rows are at the thesis parameters N = 4, m = 0.5 w, J = w unless the label
says otherwise. The experiment behind each row is the same one the matching
test runs (same times, step counts, seeds), so the table and the tests describe
the same thing.

Run:  python scripts/09_verification_table.py      (about a minute)
"""
import itertools
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from schwinger import analytics
from schwinger.budget import budget, dip, exact_distribution, point_circuit, predicted_dip
from schwinger.circuits import hamiltonian_op, statevector, to_hardware_basis, trotter_circuit
from schwinger.diagnostics import (
    fidelity_from_charge_variance, global_depolarize, mixed_state_values,
    observables_from_counts, observables_from_distribution,
)
from schwinger.ising import tfim_exact_spectrum, tfim_hamiltonian
from schwinger.noise import (
    apply_readout_error, depolarizing_lambda, gate_noise_model, noisy_probabilities,
    probabilities, sample_counts, step_channel_probabilities,
)
from schwinger.observables import evolve
from schwinger.paulis import eigh_sorted, expectation
from schwinger.schwinger import (
    bare_vacuum, charge_density_op, construction_residual, cp_operator, electric_field_op,
    fermion_construction_residual, ground_state, hamiltonian, spectrum, total_charge_op,
)
from schwinger.skqd import (
    expected_skqd_energy, krylov_circuits, run_skqd, sector_indices, subspace_ground_state,
)

N, M, W, J = 4, 0.5, 1.0, 1.0


def row(label: str, value: str) -> None:
    print(f"  {label:<78} {value}")


def phase_align(a, b):
    """Remove the global phase of a relative to b, which no measurement sees."""
    return a * np.exp(-1j * np.angle(np.vdot(b, a)))


def hamiltonian_rows():
    print("HAMILTONIAN")
    row("three constructions: Pauli expansion vs direct matrix",
        f"{construction_residual(N, M, W, J):.1e}")
    row("three constructions: fermionic form through Jordan-Wigner vs spin form",
        f"{fermion_construction_residual(N, M, W, J):.1e}")
    e0_w0, _ = ground_state(N, M, 0.0, J)
    row("strong coupling, w = 0: E0 vs -mN/2",
        f"{abs(e0_w0 - analytics.strong_coupling_ground_energy(N, M)):.1e}")
    row("free fermions, m = J = 0: all 16 levels vs subset sums of 2w cos(k pi/5)",
        f"{np.max(np.abs(spectrum(N, 0.0, W, 0.0) - analytics.free_fermion_spectrum(N, W))):.1e}")
    row("Hellmann-Feynman: condensate from the state vs N^-1 dE0/dm (dm = 1e-5)",
        f"{analytics.hellmann_feynman_residual(N, M, W, J):.1e}")
    row("Gauss's law as an operator identity", f"{analytics.gauss_law_residual(N):.1e}")
    row("[H, Q], largest element", f"{analytics.charge_commutator_residual(N, M, W, J):.1e}")
    inside, full = analytics.cp_symmetry_residual(N, M, W, J)
    row("[CP, H] inside Q = 0 (largest element on the full space, which must not vanish)",
        f"{inside:.1e} ({full:.2g})")
    cp, q = cp_operator(N), total_charge_op(N)
    rho_map = max(np.max(np.abs(cp @ charge_density_op(j, N) @ cp + charge_density_op(N - 1 - j, N)))
                  for j in range(N))
    link_map = max(np.max(np.abs(cp @ electric_field_op(j, N) @ cp
                                 - (electric_field_op(N - 2 - j, N) - q))) for j in range(N - 1))
    row("CP maps as operators: rho_n -> -rho_(N+1-n); L_n -> L_(N-n) - Q",
        f"{rho_map:.1e}; {link_map:.1e}")
    worst = 0.0
    for n in (2, 4, 6, 8):
        for r in np.linspace(0.0, 2.5, 51):
            ed, _ = eigh_sorted(tfim_hamiltonian(n, 1.0, r))
            worst = max(worst, float(np.max(np.abs(np.sort(ed) - tfim_exact_spectrum(n, 1.0, r)))))
    row("transverse-field Ising warm-up: full 2^N spectrum, N = 2, 4, 6, 8, 51 field values",
        f"{worst:.1e}")


def circuit_rows():
    print("CIRCUITS")
    row("Qiskit SparsePauliOp from the Pauli labels vs NumPy matrix",
        f"{np.max(np.abs(hamiltonian_op(N, M, W, J).to_matrix() - hamiltonian(N, M, W, J))):.1e}")
    ham = hamiltonian(N, M, W, J)
    exact2 = evolve(bare_vacuum(N), ham, np.array([2.0]))[0]
    got = phase_align(statevector(trotter_circuit(N, M, W, J, 2.0, 4000, 2)), exact2)
    row("circuit state vs exp(-iHt)|vac>, t = 2, 4000 second-order steps",
        f"{np.max(np.abs(got - exact2)):.1e}")
    # Same fit as scripts/05 section A, so the table and the Results quote one number.
    exact3 = evolve(bare_vacuum(N), ham, np.array([3.0]))[0]
    steps = np.array([10, 20, 40, 80, 160, 320])
    slopes = []
    for order in (1, 2):
        errs = [np.max(np.abs(phase_align(statevector(trotter_circuit(N, M, W, J, 3.0, int(s), order)),
                                          exact3) - exact3)) for s in steps]
        slopes.append(np.polyfit(np.log(3.0 / steps), np.log(errs), 1)[0])
    row("fitted exponent of the state error, first / second order, t = 3, dt 0.3 to 0.0094",
        f"{slopes[0]:.3f} / {slopes[1]:.3f}")
    s2 = np.array([20, 40, 80, 160])
    biases = [abs(budget(point_circuit(N, M, W, J, 1.76, int(s)), N, None, 0.0, 1000).trotter_bias)
              for s in s2]
    row("fitted exponent of the second-order bias on nu, t = 1.76, dt 0.088 to 0.011",
        f"{np.polyfit(np.log(1.76 / s2), np.log(biases), 1)[0]:.3f}")
    worst = 0.0
    for dt in (0.5, 0.25, 0.1):
        for t in (0.89, 1.76, 3.0):
            p = noisy_probabilities(trotter_circuit(N, M, W, J, t, max(1, round(t / dt)), 2), None)
            worst = max(worst, abs(observables_from_distribution(p, N)["cp_asymmetry"]))
    row("CP asymmetry of the noiseless circuit, dt = 0.5, 0.25, 0.1 at t = 0.89, 1.76, 3",
        f"{worst:.1e}")


def noise_rows():
    print("NOISE")
    worst = 0.0
    for t, s in ((0.9, 9), (1.76, 18)):
        qc = trotter_circuit(N, M, W, J, t, s, 2)
        worst = max(worst, np.max(np.abs(noisy_probabilities(qc, None) - np.abs(statevector(qc)) ** 2)))
    row("density-matrix simulation vs statevector, no noise", f"{worst:.1e}")
    qc = trotter_circuit(N, M, W, J, 1.76, 18, 2)
    np_route = step_channel_probabilities(N, M, W, J, 1.76, 18, 2, "depolarizing", 0.0)
    row("NumPy step-matrix route vs Aer route, no noise",
        f"{np.max(np.abs(np_route - noisy_probabilities(qc, None))):.1e}")

    from qiskit.quantum_info import average_gate_fidelity
    from qiskit_aer.noise import depolarizing_error
    rel = max(abs((1 - average_gate_fidelity(depolarizing_error(depolarizing_lambda(e, n), n)
                                             .to_quantumchannel())) - e) / e
              for n in (1, 2) for e in (1e-4, 3e-3, 2e-2))
    row("eps vs 1 - (Qiskit average gate fidelity), relative, 1 and 2 qubits", f"{rel:.1e}")

    p9 = noisy_probabilities(trotter_circuit(N, M, W, J, 0.9, 9, 2), None)
    nu = observables_from_distribution(p9, N)["particle_density"]
    worst = max(abs(observables_from_distribution(apply_readout_error(p9, q), N)["particle_density"]
                    - (0.5 + (1 - 2 * q) * (nu - 0.5))) for q in (0.01, 0.05, 0.2))
    row("readout: nu vs closed form 1/2 + (1-2q)(nu - 1/2), q = 0.01, 0.05, 0.2", f"{worst:.1e}")

    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, ReadoutError
    q, shots = 0.08, 200_000
    predicted = apply_readout_error(p9, q)
    nm = NoiseModel()
    nm.add_all_qubit_readout_error(ReadoutError([[1 - q, q], [q, 1 - q]]))
    meas = to_hardware_basis(trotter_circuit(N, M, W, J, 0.9, 9, 2))
    meas.measure_all()
    counts = AerSimulator(noise_model=nm).run(meas, shots=shots, seed_simulator=3).result().get_counts()
    observed = np.zeros(16)
    for bits, c in counts.items():
        observed[int(bits, 2)] = c / shots
    dev = np.max(np.abs(observed - predicted) / np.sqrt(predicted * (1 - predicted) / shots))
    row("readout: classical confusion matrix vs Aer, q = 0.08, 2e5 shots, worst of 16 bins",
        f"{dev:.1f} sigma")

    worst = 0.0
    for p in (0.005, 0.02, 0.1):
        g = step_channel_probabilities(N, M, W, J, 1.76, 18, 2, "global", p)
        fq = fidelity_from_charge_variance(observables_from_distribution(g, N)["charge_variance"], N)
        worst = max(worst, abs(fq - (1 - p) ** 18))
    row("global channel per step, p = 0.005, 0.02, 0.1, 18 steps: F_Q vs (1-p)^18", f"{worst:.1e}")

    var_q, a_cp, shift = 0.0, 0.0, np.inf
    clean = observables_from_distribution(
        step_channel_probabilities(N, M, W, J, 1.76, 18, 2, "dephasing", 0.0), N)["particle_density"]
    for p in (0.01, 0.05):
        d = observables_from_distribution(
            step_channel_probabilities(N, M, W, J, 1.76, 18, 2, "dephasing", p), N)
        var_q = max(var_q, abs(d["charge_variance"]))
        a_cp = max(a_cp, abs(d["cp_asymmetry"]))
        shift = min(shift, abs(d["particle_density"] - clean))
    row("logical-frame dephasing, p = 0.01, 0.05: Var(Q); A_CP (smallest shift of nu it causes)",
        f"{var_q:.1e}; {a_cp:.1e} ({shift:.3f})")

    full = noisy_probabilities(trotter_circuit(N, M, W, J, 0.9, 9, 2), gate_noise_model(0.75, 0.0))
    row("fully depolarising CX (eps2 = 3/4): max |p_b - 1/16|, |Var(Q) - N/4|",
        f"{np.max(np.abs(full - 1 / 16)):.1e}, "
        f"{abs(observables_from_distribution(full, N)['charge_variance'] - mixed_state_values(N)['charge_variance']):.1e}")


def statistics_and_skqd_rows():
    print("STATISTICS, BUDGET AND GROUND STATE")
    psi = evolve(bare_vacuum(N), hamiltonian(N, M, W, J), np.array([0.9]))[0]
    p = global_depolarize(np.abs(psi) ** 2, 0.6)
    ratios = []
    for shots in (256, 2048):
        rng = np.random.default_rng(11)
        dens, varq, qd, qv = [], [], [], []
        for _ in range(400):
            o = observables_from_counts(rng.multinomial(shots, p).astype(float), N)
            dens.append(o["particle_density"]); qd.append(o["particle_density_err"])
            varq.append(o["charge_variance"]); qv.append(o["charge_variance_err"])
        ratios += [np.std(dens) / np.mean(qd), np.std(varq) / np.mean(qv)]
    row("scatter of 400 runs / quoted standard error, nu and Var(Q), 256 and 2048 shots",
        f"{min(ratios):.2f} to {max(ratios):.2f}")

    tr = point_circuit(N, M, W, J, 1.76, 18)
    nm, q, shots = gate_noise_model(1e-3), 0.02, 500
    b = budget(tr, N, nm, q, shots)
    measured = apply_readout_error(probabilities(tr.prepared, nm), q)
    rng = np.random.default_rng(99)
    est = np.array([observables_from_counts(sample_counts(measured, shots, rng), N)["particle_density"]
                    for _ in range(600)])
    row("budget closes: mean squared error of 600 runs / (b^2 + sigma^2), 500 shots",
        f"{np.mean((est - b.exact) ** 2) / (b.total_bias ** 2 + b.sigma ** 2):.3f}")

    worst = 0.0
    p_peak, p_trough = exact_distribution(N, M, W, J, 0.89), exact_distribution(N, M, W, J, 1.76)
    k_peak, k_trough, eps0 = 300.0, 600.0, 1e-3
    grid = np.logspace(-4, -2, 9)
    pr = predicted_dip(global_depolarize(p_peak, np.exp(-k_peak * eps0)),
                       global_depolarize(p_trough, np.exp(-k_trough * eps0)), eps0, grid, N, 1000)
    for e, d_pred in zip(grid, pr["D"]):
        d_true, _ = dip(global_depolarize(p_peak, np.exp(-k_peak * e)),
                        global_depolarize(p_trough, np.exp(-k_trough * e)), N, 1000)
        worst = max(worst, abs(d_pred - d_true))
    row("threshold prediction on data that obey the global model: max |D_pred - D_true|",
        f"{worst:.1e}")

    ham = hamiltonian(N, M, W, J)
    e0, _ = ground_state(N, M, W, J)
    e_full, _ = subspace_ground_state(sector_indices(N), ham)
    row("SKQD on the whole neutral sector vs exact ground energy", f"{abs(e_full - e0):.1e}")
    # Proper subspaces only (r < 6). The full sector is the row above; including
    # it here would just report its rounding error as the "smallest margin".
    margins = [subspace_ground_state(s, ham)[0] - e0
               for r in range(1, 6) for s in itertools.combinations(sector_indices(N), r)]
    row("smallest E(S) - E0 over the 62 proper neutral subspaces (must be positive)",
        f"{min(margins):.1e}")

    dists = [noisy_probabilities(c) for c in krylov_circuits(N, M, W, J, 0.3, 4)]
    ex = expected_skqd_energy(dists, 15, N, M, W, J)
    rng = np.random.default_rng(20261012)
    runs = np.array([run_skqd(dists, 15, rng, N, M, W, J)[-1].energy for _ in range(20000)])
    row("SKQD mean bias, 15 shots per step: exact formula vs 20000 simulated runs",
        f"{ex.mean_bias:.4e} vs {runs.mean() - e0:.4e} +/- {runs.std() / np.sqrt(runs.size):.1e}")


def dip_rows():
    """
    The threshold uses the recombination dip D as its signal on the grounds
    that noise erodes D and never creates it. That is an observation about the
    channels modelled, not a theorem, so it is checked here on every one of
    them: gate depolarizing and gate dephasing (Aer), and the two
    logical-frame step channels the thesis defines (per-qubit dephasing and the
    global depolarizing channel), at the three step sizes of the threshold.
    """
    print("THRESHOLD SIGNAL")
    from schwinger.budget import sweep_dip
    nu = lambda p: observables_from_distribution(p, N)["particle_density"]
    eps = np.logspace(-5, -1.5, 36)
    ps = np.logspace(-4, -0.5, 30)
    worst, lowest = -np.inf, np.inf
    for dt in (0.1, 0.2, 0.3):
        sp, st = max(1, round(0.89 / dt)), max(1, round(1.76 / dt))
        pk, tr = point_circuit(N, M, W, J, 0.89, sp), point_circuit(N, M, W, J, 1.76, st)
        d0 = nu(pk.trotter) - nu(tr.trotter)
        for kind in ("depolarizing", "dephasing"):
            d = sweep_dip(pk, tr, N, eps, 1000, kind=kind)["D"]
            worst, lowest = max(worst, d.max() - d0), min(lowest, d.min())
        for ch in ("dephasing", "global"):
            d = np.array([nu(step_channel_probabilities(N, M, W, J, 0.89, sp, 2, ch, p))
                          - nu(step_channel_probabilities(N, M, W, J, 1.76, st, 2, ch, p)) for p in ps])
            worst, lowest = max(worst, d.max() - d0), min(lowest, d.min())
    row("largest D - D(noiseless), 4 channels, dt = 0.1, 0.2, 0.3 (must be < 0); lowest D",
        f"{worst:.1e}; {lowest:+.2f}")


def energy_rows():
    """
    Where the pair energy comes from. On the bare vacuum the mass term is at its
    minimum (-mN/2) and the electric term is zero, its minimum, so making pairs
    raises both. The total is conserved, so the hopping term must pay for both.
    H is linear in (m, w, J), so the three pieces are H at one coupling each.
    """
    print("ENERGY AFTER THE QUENCH FROM THE BARE VACUUM")
    h_mass, h_hop, h_el = hamiltonian(N, M, 0, 0), hamiltonian(N, 0, W, 0), hamiltonian(N, 0, 0, J)
    total = hamiltonian(N, M, W, J)
    row("|H_mass + H_hop + H_el - H|", f"{np.max(np.abs(h_mass + h_hop + h_el - total)):.1e}")
    for t, s in zip((0.0, 0.89, 1.76), evolve(bare_vacuum(N), total, np.array([0.0, 0.89, 1.76]))):
        row(f"t = {t:<4}: <H_mass>, <H_hop>, <H_el>, <H>",
            ", ".join(f"{expectation(s, op):+.3f}" for op in (h_mass, h_hop, h_el, total)))
    times = np.linspace(0.0, 12.0, 481)
    drift = 0.0
    for m in (0.0, 0.5, 1.0, 2.0):
        H = hamiltonian(N, m, W, J)
        e = [expectation(s, H) for s in evolve(bare_vacuum(N), H, times)]
        drift = max(drift, float(np.ptp(e)))
    row("largest drift of <H> over t in [0, 12], m/w = 0, 0.5, 1, 2", f"{drift:.1e}")


def main():
    hamiltonian_rows()
    circuit_rows()
    noise_rows()
    statistics_and_skqd_rows()
    dip_rows()
    energy_rows()


if __name__ == "__main__":
    main()
