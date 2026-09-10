"""
Past Gate 2, and cheap: real-time pair production out of the vacuum.

This is the headline physics of the whole project, and at N = 4 it costs one
diagonalisation. It is exact, so it is also the ground truth that Trotterised
evolution gets compared against in M8.

WHERE THE ENERGY COMES FROM, since this is the question you will be asked.
Nothing is injected and no external agent pulls charges apart. The simulation
starts in the BARE vacuum, which is an eigenstate of the mass and electric
terms but not of the full Hamiltonian once w is nonzero. Total energy is
conserved exactly (checked below). What changes is how that energy is
distributed between electric field and rest mass. Pairs appear because the
initial state was never the true ground state.

Run:  python scripts/04_pair_production_preview.py
"""
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from schwinger.observables import (
    charge_variance, electric_field_profile, evolve, particle_density, total_charge,
)
from schwinger.paulis import expectation
from schwinger.plotstyle import MARKERS, PALETTE, annotate, save, use_report_style
from schwinger.schwinger import bare_vacuum, ground_state, hamiltonian

import matplotlib.pyplot as plt

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)
N = 4


def main() -> None:
    use_report_style()
    times = np.linspace(0.0, 12.0, 481)
    masses = [0.0, 0.5, 1.0, 2.0]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.6))

    print("=" * 62)
    print("CONSERVATION CHECKS during evolution (these must not drift)")
    print("=" * 62)

    for i, m in enumerate(masses):
        H = hamiltonian(N, m, 1.0, 1.0)
        v = bare_vacuum(N)
        states = evolve(v, H, times)
        dens = np.array([particle_density(s, N) for s in states])
        ax1.plot(times, dens, color=PALETTE[i], marker=MARKERS[i], markevery=40,
                 label=f"m/w = {m}")

        e = np.array([expectation(s, H) for s in states])
        q = np.array([total_charge(s, N) for s in states])
        vq = np.array([charge_variance(s, N) for s in states])
        nrm = np.array([abs(np.vdot(s, s)) for s in states])
        print(f"  m/w = {m}:  energy drift {np.ptp(e):.2e}   |Q| max {np.max(np.abs(q)):.2e}"
              f"   Var(Q) max {np.max(np.abs(vq)):.2e}   norm drift {np.ptp(nrm):.2e}")

    ax1.set_ylim(-0.02, 0.80)
    ax1.set_xlabel("time  (units of 1/w)")
    ax1.set_ylabel("particle density")
    ax1.set_title("Pairs created out of the bare vacuum, N = 4")
    ax1.legend(loc="upper right", ncol=2)
    annotate(ax1, "Starts at exactly zero by construction. Heavier fermions cost more to make,\n"
                  "so the density oscillates with smaller amplitude. Energy is conserved\n"
                  "throughout: it moves between electric field and rest mass, nothing enters.")

    m = 0.5
    H = hamiltonian(N, m, 1.0, 1.0)
    states = evolve(bare_vacuum(N), H, times)
    profile = np.array([electric_field_profile(s, N) for s in states])
    ax2.plot(times, profile[:, 1], color=PALETTE[1], marker=MARKERS[1], markevery=48,
             label="link 2 (middle)")
    ax2.plot(times, profile[:, 0], color=PALETTE[0], linewidth=3.2, alpha=0.5,
             label="link 1 (outer)")
    ax2.plot(times, profile[:, 2], color=PALETTE[2], linewidth=1.4,
             linestyle=(0, (5, 3)), label="link 3 (outer)")
    ax2.axhline(0.0, color="#9a9a9a", linewidth=1.0)
    ax2.set_xlabel("time  (units of 1/w)")
    ax2.set_ylabel("electric field  <L_n>")
    ax2.set_title(f"Flux on each link, m/w = {m}")
    ax2.legend(loc="upper right")
    worst = float(np.max(np.abs(profile[:, 0] - profile[:, 2])))
    print(f"\n  CP check: max |<L_1> - <L_3>| over the evolution = {worst:.3e}")
    annotate(ax2, "The two outer links lie exactly on top of each other. That is not a\n"
                  "coincidence: CP symmetry forces <L_n> = <L_(N-n)> on any neutral state, so\n"
                  "it is a free correctness check. The field is not independent either, since\n"
                  "Gauss's law fixes it from the charges to its left.")

    out = FIG / "04_pair_production.png"
    save(fig, out)
    print(f"\nwrote {out}")

    print()
    print("=" * 62)
    print("WHY PAIRS APPEAR AT ALL: bare vacuum is not the ground state")
    print("=" * 62)
    print(f"{'m/w':>5} {'E(bare vacuum)':>16} {'E0 (true)':>12} {'gap':>10} {'peak density':>14}")
    for m in masses:
        H = hamiltonian(N, m, 1.0, 1.0)
        v = bare_vacuum(N)
        e_bare = expectation(v, H)
        e0, _ = ground_state(N, m, 1.0, 1.0)
        peak = max(particle_density(s, N) for s in evolve(v, H, times))
        print(f"{m:>5} {e_bare:>16.6f} {e0:>12.6f} {e_bare - e0:>10.6f} {peak:>14.6f}")
    print()
    print("  The gap column is the whole story. A larger mismatch between the state")
    print("  you prepare and the state the Hamiltonian wants gives more pair production.")


if __name__ == "__main__":
    main()
