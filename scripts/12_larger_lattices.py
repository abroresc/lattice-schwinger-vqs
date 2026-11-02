"""
Is the four-site model the physics of the lattice Schwinger model, or an
artefact of four sites? The same Hamiltonian on chains of N = 4 to 12, solved
exactly (classically) in the neutral sector, at the thesis parameters
(J = w, m = 0.5 w, open boundaries, the same lattice spacing ag = 1).

  A. The quench from the bare vacuum: the first peak of nu(t), the
     recombination minimum after it, and the dip between them, which is the
     signal the failure threshold is built on.
  B. The ground state: energy and chiral condensate per site.

This is a cross-check of the model, not a larger quantum simulation: every
quantum-circuit result in the thesis stays at N = 4. It says nothing about the
continuum limit either, since the lattice spacing is fixed.

Run:  python scripts/12_larger_lattices.py      (about a minute; N = 12 dominates)
"""
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from schwinger.finite_size import first_peak_and_trough, ground_state, neutral_block, quench_density
from schwinger.plotstyle import MARKERS, PALETTE, annotate, save, use_report_style

import matplotlib.pyplot as plt

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
M, W, J = 0.5, 1.0, 1.0
SIZES = (4, 6, 8, 10, 12)


def main():
    t = np.linspace(0, 4, 2001)
    curves, rows = {}, []
    for n in SIZES:
        nu = quench_density(n, M, W, J, t)
        tp, yp, tt, yt = first_peak_and_trough(t, nu)
        e0, cond = ground_state(n, M, W, J)
        curves[n] = nu
        rows.append((n, neutral_block(n, M, W, J)[0].size, tp, yp, tt, yt, yp - yt, e0 / n, cond))

    print("=" * 74)
    print("A. QUENCH FROM THE BARE VACUUM   B. GROUND STATE   (J = w, m = 0.5 w)")
    print("=" * 74)
    print("   N  neutral dim   t_peak  nu_peak   t_trough  nu_trough   dip D    E0/N     condensate")
    for n, dim, tp, yp, tt, yt, d, e, c in rows:
        print(f"  {n:2d}  {dim:8d}     {tp:.3f}   {yp:.4f}    {tt:.3f}     {yt:.4f}    {d:.4f}  "
              f"{e:+.4f}    {c:+.4f}")
    ref, four = rows[-1], rows[0]
    rel = lambda a, b: 100 * (a - b) / abs(b)
    print(f"  N = 4 against N = 12: peak time {rel(four[2], ref[2]):+.1f} %, trough time "
          f"{rel(four[4], ref[4]):+.1f} %, dip {rel(four[6], ref[6]):+.1f} %, "
          f"condensate {rel(four[8], ref[8]):+.1f} %")
    steps = [abs(rows[i + 1][6] - rows[i][6]) for i in range(len(rows) - 1)]
    print("  change in the dip from each N to the next: " + ", ".join(f"{s:.4f}" for s in steps))
    # Beyond the first oscillation the four-site chain departs from the longer
    # ones (the pairs reach the ends of so short a chain): say so, with numbers.
    second = []
    for n in SIZES:
        late = t >= rows[SIZES.index(n)][4]
        tp2, yp2, _, _ = first_peak_and_trough(t[late], curves[n][late])
        second.append((n, tp2, yp2))
    print("  second peak, after the trough (not used in the thesis): "
          + ", ".join(f"N = {n}: nu = {y:.3f} at t = {tp:.2f}" for n, tp, y in second))

    use_report_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.0, 3.6), gridspec_kw={"width_ratios": [1.6, 1]})
    fig.subplots_adjust(wspace=0.32)
    for i, n in enumerate(SIZES):
        ax1.plot(t, curves[n], color=PALETTE[i], linewidth=2.2 if n == 4 else 1.4,
                 label=f"N = {n}" + ("  (the thesis)" if n == 4 else ""))
    ax1.axvline(0.89, color="#9a9a9a", linewidth=0.8, linestyle=(0, (4, 3)))
    ax1.axvline(1.76, color="#9a9a9a", linewidth=0.8, linestyle=(0, (4, 3)))
    ax1.set_xlabel("time  (units of 1/w)")
    ax1.set_ylabel(r"particle density $\nu$")
    ax1.set_title("(a) Pair production on longer chains, solved exactly")
    ax1.set_ylim(-0.01, 0.66)
    ax1.legend(loc="upper right", fontsize=7.5, ncol=2)
    ns = np.array([r[0] for r in rows])
    ax2.plot(ns, [r[6] for r in rows], color=PALETTE[0], marker=MARKERS[0], label="recombination dip $D$")
    ax2.plot(ns, [-r[8] for r in rows], color=PALETTE[1], marker=MARKERS[1],
             label=r"$-\langle\bar\psi\psi\rangle$ in the ground state")
    ax2.set_xticks(ns)
    ax2.set_xlabel("chain length N")
    ax2.set_ylim(0.0, 0.45)
    ax2.set_title("(b) Signal and ground state against N")
    ax2.legend(loc="lower right", fontsize=7.5)
    annotate(ax1, "Dashed lines: the N = 4 peak and trough times used in the thesis. The first\n"
                  "rise and recombination dip are there at every length and shift by under 10 %.\n"
                  "Later, the four-site chain departs: its second peak is a quarter lower.")
    out = FIG / "12_larger_lattices.png"
    save(fig, out)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
