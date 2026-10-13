"""
M11 / M12: gate error, the failure threshold, and whether the symmetry
diagnostic predicts it.

  A. Validation residuals for the noise layer, printed as numbers.
  B. The recombination dip D = nu(t_peak) - nu(t_trough) against the two-qubit
     gate error, for three Trotter step sizes, and the failure threshold eps*
     where D falls to twice its shot-noise uncertainty.
  C. The prediction of eps* from charge-conservation data taken at ONE lower
     noise level, under the global-depolarizing model.
  D. The control: noise the charge diagnostic cannot see.

Every number here comes from exact outcome distributions (density-matrix
simulation); shot noise enters only through the analytic standard error. So
nothing in this script is a random fluctuation, and rerunning it gives the
same output to the last digit.

Run:  python scripts/06_noise_threshold.py
"""
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from schwinger.budget import first_crossing, point_circuit, predicted_dip, sweep_dip
from schwinger.circuits import statevector, trotter_circuit
from schwinger.diagnostics import (
    fidelity_from_charge_variance, observables_from_distribution,
)
from schwinger.noise import (
    gate_noise_model, noisy_probabilities, prepare_for_probabilities, probabilities,
    step_channel_probabilities,
)
from schwinger.observables import evolve, particle_density
from schwinger.plotstyle import MARKERS, MUTED, PALETTE, save, use_report_style
from schwinger.schwinger import bare_vacuum, hamiltonian

import matplotlib.pyplot as plt

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
N, M, W, J = 4, 0.5, 1.0, 1.0
T_PEAK, T_TROUGH = 0.89, 1.76          # first maximum and minimum of the exact nu(t)
SHOTS = (1000, 10000)
DTS = (0.1, 0.2, 0.3)
EPS = np.logspace(-5, -1.5, 36)


def steps_for(t, dt):
    return max(1, int(round(t / dt)))


def part_a():
    print("=" * 70)
    print("A. NOISE-LAYER VALIDATION (residuals, not pass/fail flags)")
    print("=" * 70)
    qc = trotter_circuit(N, M, W, J, T_TROUGH, 18, 2)
    sv = np.abs(statevector(qc)) ** 2
    dm = noisy_probabilities(qc, None)
    np_route = step_channel_probabilities(N, M, W, J, T_TROUGH, 18, 2, "global", 0.0)
    print(f"  density matrix vs statevector, no noise      : {np.max(np.abs(dm - sv)):.1e}")
    print(f"  NumPy step route vs Aer route, no noise       : {np.max(np.abs(np_route - dm)):.1e}")
    g = step_channel_probabilities(N, M, W, J, T_TROUGH, 18, 2, "global", 0.02)
    fq = fidelity_from_charge_variance(observables_from_distribution(g, N)["charge_variance"], N)
    print(f"  global channel: F_Q - (1-p)^18                : {abs(fq - 0.98**18):.1e}")
    d = observables_from_distribution(
        step_channel_probabilities(N, M, W, J, T_TROUGH, 18, 2, "dephasing", 0.02), N)
    print(f"  logical dephasing: Var(Q)                     : {d['charge_variance']:.1e}")
    print(f"  logical dephasing: CP asymmetry               : {abs(d['cp_asymmetry']):.1e}")
    worst = 0.0
    for dt in (0.5, 0.25, 0.1):
        for t in (T_PEAK, T_TROUGH, 3.0):
            p = noisy_probabilities(trotter_circuit(N, M, W, J, t, steps_for(t, dt), 2), None)
            worst = max(worst, abs(observables_from_distribution(p, N)["cp_asymmetry"]))
    print(f"  noiseless Trotter circuit: worst CP asymmetry : {worst:.1e}")
    full = observables_from_distribution(noisy_probabilities(qc, gate_noise_model(0.75, 0.0)), N)
    print(f"  fully depolarizing CX: Var(Q) - N/4           : {abs(full['charge_variance'] - N / 4):.1e}")


def part_bcd():
    print()
    print("=" * 70)
    print("B. FAILURE THRESHOLD   C. PREDICTION FROM THE CHARGE DIAGNOSTIC")
    print("=" * 70)
    results = {}
    for dt in DTS:
        pk = point_circuit(N, M, W, J, T_PEAK, steps_for(T_PEAK, dt))
        tr = point_circuit(N, M, W, J, T_TROUGH, steps_for(T_TROUGH, dt))
        for shots in SHOTS:
            sw = sweep_dip(pk, tr, N, EPS, shots)
            star = first_crossing(EPS, sw["D"] - 2 * sw["sigma_D"])
            star_inf = first_crossing(EPS, sw["D"])
            line = (f"  dt={dt:.1f}  CX(peak,trough)=({pk.two_qubit_gates:>3},{tr.two_qubit_gates:>3})  "
                    f"M={shots:>5}:  D(0)={sw['D'][0]:.4f}  eps*={star:.2e}  "
                    f"eps*xCX={star * tr.two_qubit_gates:.2f}  eps*(M=inf)={star_inf:.2e}")
            print(line)
            ratios = []
            for frac in (0.1, 0.3):
                i = int(np.argmin(np.abs(np.log(EPS / (frac * star)))))
                pp, pt = sw["distributions"][i]
                pr = predicted_dip(pp, pt, EPS[i], EPS, N, shots)
                ratios.append((EPS[i] / star,
                               first_crossing(EPS, pr["D"] - 2 * pr["sigma_D"]) / star))
            print("        predicted eps*/true eps* from data at eps0 = "
                  + ", ".join(f"{a:.2f} eps*: {b:.3f}" for a, b in ratios))
            results[(dt, shots)] = (pk, tr, sw, star)
    return results


def part_d_control():
    print()
    print("=" * 70)
    print("D. THE BLIND SPOT: logical-frame dephasing, once per Trotter step")
    print("=" * 70)
    ps = np.logspace(-4, -0.5, 30)
    dt = 0.1
    rows = []
    for p in ps:
        a = observables_from_distribution(step_channel_probabilities(
            N, M, W, J, T_PEAK, steps_for(T_PEAK, dt), 2, "dephasing", p), N)
        b = observables_from_distribution(step_channel_probabilities(
            N, M, W, J, T_TROUGH, steps_for(T_TROUGH, dt), 2, "dephasing", p), N)
        rows.append((a["particle_density"] - b["particle_density"], b["charge_variance"]))
    rows = np.array(rows)
    zero = first_crossing(ps, rows[:, 0])
    print(f"  dip vanishes at p = {zero:.2e} per qubit per step; largest Var(Q) seen = {rows[:, 1].max():.1e}")
    print("  The charge diagnostic reads exactly zero while the signal disappears.")
    return ps, rows


def figure(results, control):
    use_report_style()
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14.5, 3.9))
    fig.subplots_adjust(wspace=0.3)

    # (a) trajectories at dt = 0.1
    times = np.linspace(0.0, 4.0, 41)
    exact = [particle_density(s, N) for s in evolve(bare_vacuum(N), hamiltonian(N, M, W, J), times)]
    ax1.plot(times, exact, color="#6b6b6b", linewidth=1.6, label="exact")
    preps = [prepare_for_probabilities(trotter_circuit(N, M, W, J, t, steps_for(t, 0.1), 2))
             if t > 0 else None for t in times]
    for i, eps in enumerate((1e-3, 3e-3, 1e-2)):
        nm = gate_noise_model(eps)
        nu = [0.0 if p is None else observables_from_distribution(probabilities(p, nm), N)["particle_density"]
              for p in preps]
        ax1.plot(times, nu, color=PALETTE[i], marker=MARKERS[i], markevery=4,
                 label=rf"$\epsilon_2$ = {eps:g}")
    for t in (T_PEAK, T_TROUGH):
        ax1.axvline(t, color=MUTED, linewidth=0.8, linestyle=(0, (3, 3)))
    ax1.axhline(0.5, color=MUTED, linewidth=0.8)
    ax1.set_xlabel("time  (units of 1/w)")
    ax1.set_ylabel(r"particle density $\nu$")
    ax1.set_title(r"(a) Pair production under gate error, $\delta t$ = 0.1")
    ax1.legend(loc="lower right")

    # (b) dip and noise floor against eps
    for i, dt in enumerate(DTS):
        pk, tr, sw, star = results[(dt, 1000)]
        ax2.semilogx(EPS, sw["D"], color=PALETTE[i], marker=MARKERS[i], markevery=3,
                     label=rf"$\delta t$ = {dt}  ({tr.two_qubit_gates} CX)")
        ax2.semilogx(EPS, 2 * sw["sigma_D"], color=PALETTE[i], linewidth=1.0,
                     linestyle=(0, (4, 3)))
        ax2.axvline(star, color=PALETTE[i], linewidth=0.8, alpha=0.6)
        j = int(np.argmin(np.abs(np.log(EPS / (0.1 * star)))))
        pp, pt = sw["distributions"][j]
        pr = predicted_dip(pp, pt, EPS[j], EPS, N, 1000)
        ax2.semilogx(EPS, pr["D"], color=PALETTE[i], linewidth=1.0, linestyle=(0, (1, 1.5)))
    ax2.set_xlabel(r"two-qubit gate error $\epsilon_2$")
    ax2.set_ylabel(r"recombination dip $D$")
    ax2.set_title("(b) Signal against noise floor, 1000 shots")
    ax2.legend(loc="upper right", fontsize=7.5)
    ax2.set_ylim(-0.05, 0.4)

    # (c) diagnostic against damage
    pk, tr, sw, star = results[(0.1, 1000)]
    atten = sw["D"] / sw["D"][0]
    ax3.plot(sw["F_trough"], atten, color=PALETTE[0], marker=MARKERS[0], linestyle="none",
             label="depolarizing CX")
    deph = sweep_dip(pk, tr, N, EPS, 1000, kind="dephasing")
    ax3.plot(deph["F_trough"], deph["D"] / deph["D"][0], color=PALETTE[1], marker=MARKERS[1],
             linestyle="none", label="dephasing on CX")
    ps, rows = control
    fq = [fidelity_from_charge_variance(v, N) for v in rows[:, 1]]
    ax3.plot(fq, rows[:, 0] / sw["D"][0], color=PALETTE[2], marker=MARKERS[2], linestyle="none",
             label="dephasing per step (logical)")
    ax3.plot([0, 1], [0, 1], color=MUTED, linewidth=1.0)
    ax3.set_xlabel(r"$F_Q = 1 - \mathrm{Var}(Q)/(N/4)$ at $t_{\rm trough}$")
    ax3.set_ylabel(r"signal kept, $D/D_0$")
    ax3.set_title("(c) Does the charge diagnostic see the damage?")
    ax3.legend(loc="upper left", fontsize=7.5)
    ax3.set_xlim(-0.05, 1.05)
    ax3.set_ylim(-0.15, 1.1)

    out = FIG / "06_noise_threshold.png"
    save(fig, out)
    print(f"\nwrote {out}")


def main():
    part_a()
    results = part_bcd()
    control = part_d_control()
    figure(results, control)


if __name__ == "__main__":
    main()
