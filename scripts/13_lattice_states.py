"""
What each four-qubit bit string means, as a figure (thesis Figure 1).

It shows the six
arrangements with zero total charge (the only ones the dynamics ever visits
from the bare vacuum), and for each one:

  * the bit string the quantum computer stores (|0> = site occupied),
  * the charge on every site: + for a particle, - for an antiparticle,
  * the electric field on every link, from the solved Gauss's law,
  * the energy of the arrangement with hopping switched off.

Nothing is typed in by hand. Charges come from charge_profile, fields from
electric_field_profile, energies from the diagonal of the same Hamiltonian
matrix every other script uses. The script asserts the three facts the figure
is meant to teach, so if a convention ever changes the figure fails loudly
instead of drawing something wrong:

  1. the bare vacuum has no charge anywhere and no field anywhere,
  2. every neutral arrangement has Gauss's law L_n - L_{n-1} = rho_n exactly,
  3. a pair costs 2m for the two particles plus J for every link its field
     crosses (that is confinement: the stretched pair |1100> costs 3J of field).

Hopping, the one term not drawn, is the term that turns one row into another.
Its expectation value is zero in every row (each row is a single bit string),
which is why the energies here are just the diagonal of H.

Output: figures/13_lattice_states.png and .pdf (the thesis includes the PDF).
Run:    .venv/bin/python scripts/13_lattice_states.py
"""
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from schwinger.observables import charge_profile, electric_field_profile
from schwinger.plotstyle import INK, MUTED, PALETTE, save, use_report_style
from schwinger.schwinger import bare_vacuum_label, charge_sectors, hamiltonian

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
N, M, W, J = 4, 0.5, 1.0, 1.0

# Row order is the order a reader needs: nothing, then the three ways one
# neighbouring pair can sit, then the stretched pair and two pairs.
ORDER = ["0101", "1001", "0011", "0110", "1100", "1010"]
NAMES = {
    "0101": "empty space (bare vacuum)",
    "1001": "one pair, left",
    "0011": "one pair, middle",
    "0110": "one pair, right",
    "1100": "one pair, stretched",
    "1010": "two pairs",
}
PLUS, MINUS = PALETTE[0], PALETTE[1]     # particle blue, antiparticle orange
FIELD = "#b9b9b9"


def basis_state(bits: str) -> np.ndarray:
    psi = np.zeros(2 ** N, dtype=complex)
    psi[int(bits, 2)] = 1.0              # site 1 is the most significant bit
    return psi


def main() -> None:
    use_report_style()
    H = hamiltonian(N, M, W, J)

    neutral = [format(i, f"0{N}b") for i, q in enumerate(charge_sectors(N)) if q == 0]
    assert sorted(neutral) == sorted(ORDER), neutral
    assert ORDER[0] == bare_vacuum_label(N)

    rows = []
    print(f"{'state':8s} {'charges':18s} {'link fields':16s} {'energy':>7s}  "
          f"{'2m x pairs + J x links':>22s}")
    for bits in ORDER:
        psi = basis_state(bits)
        rho = np.real(charge_profile(psi, N))
        L = np.real(electric_field_profile(psi, N))
        E = float(np.real(H[int(bits, 2), int(bits, 2)]))
        # Gauss's law on this state, link by link (L_0 = 0 at the open end).
        gauss = np.abs(np.diff(np.concatenate([[0.0], L, [L[-1] + rho[-1]]])) - rho).max()
        n_charges = int(np.abs(rho).sum())
        predicted = -M * N / 2 + M * n_charges + J * float((L ** 2).sum())
        rows.append((bits, rho, L, E))
        print(f"|{bits}>  {np.array2string(rho, precision=0, sign='+'):18s} "
              f"{np.array2string(L, precision=0, sign='+'):16s} {E:+7.3f}  "
              f"{predicted:+22.3f}   Gauss residual {gauss:.0e}")
        assert gauss < 1e-12
        assert abs(E - predicted) < 1e-12
    assert np.abs(rows[0][1]).max() < 1e-12 and np.abs(rows[0][2]).max() < 1e-12

    # ------------------------------------------------------------------ draw
    xs = np.arange(N) * 1.25 + 4.0          # site positions
    r = 0.30
    fig, ax = plt.subplots(figsize=(6.4, 3.7))
    ax.set_xlim(-0.1, 11.2)
    ax.set_ylim(-len(rows) - 0.05, 1.25)
    ax.set_aspect("equal")
    ax.axis("off")

    for j, x in enumerate(xs):
        parity = "odd" if j % 2 == 0 else "even"
        ax.text(x, 0.95, f"site {j + 1}", ha="center", va="bottom", fontsize=8, color=INK)
        ax.text(x, 0.72, parity, ha="center", va="bottom", fontsize=7, color=MUTED)
    ax.text(0.0, 0.72, "qubits", ha="left", va="bottom", fontsize=7, color=MUTED)
    ax.text(9.35, 0.95, "energy,", ha="left", va="bottom", fontsize=8, color=INK)
    ax.text(9.35, 0.72, "hopping off", ha="left", va="bottom", fontsize=7, color=MUTED)

    for i, (bits, rho, L, E) in enumerate(rows):
        y = -i
        ax.text(0.0, y + 0.08, f"|{bits}⟩", ha="left", va="center", fontsize=9,
                color=INK, family="DejaVu Sans Mono")
        ax.text(0.0, y - 0.30, NAMES[bits], ha="left", va="center", fontsize=6.5, color=MUTED)
        for n in range(N - 1):                      # links first, under the sites
            if abs(L[n]) > 0.5:
                x0, x1 = xs[n] + r, xs[n + 1] - r
                ax.add_patch(FancyBboxPatch((x0, y - 0.07), x1 - x0, 0.14,
                                            boxstyle="round,pad=0,rounding_size=0.05",
                                            fc=FIELD, ec="none", zorder=1))
                ax.text((x0 + x1) / 2, y + 0.16, f"{L[n]:+.0f}".replace("-", "\u2212"), ha="center", va="bottom",
                        fontsize=6.5, color=MUTED)
            else:
                ax.plot([xs[n] + r, xs[n + 1] - r], [y, y], color="#e3e3e3", lw=0.8, zorder=1)
        for n, x in enumerate(xs):
            q = round(rho[n])
            if q == 0:
                ax.add_patch(Circle((x, y), r, fc="white", ec="#9a9a9a", lw=1.0, zorder=2))
            else:
                ax.add_patch(Circle((x, y), r, fc=PLUS if q > 0 else MINUS, ec="white",
                                    lw=2.0, zorder=3))
                ax.text(x, y - 0.01, "+" if q > 0 else "−", ha="center", va="center",
                        fontsize=11, color="white", fontweight="bold", zorder=4)
        ax.text(10.75, y, f"{E:+.1f}".replace("-", "\u2212"), ha="right", va="center", fontsize=9, color=INK)

    # Key, under the rows, in words rather than a colour legend.
    yk = -len(rows) + 0.30
    ax.add_patch(Circle((4.0, yk), 0.16, fc=PLUS, ec="white", lw=1.5))
    ax.text(4.25, yk, "particle", va="center", fontsize=7, color=INK)
    ax.add_patch(Circle((5.55, yk), 0.16, fc=MINUS, ec="white", lw=1.5))
    ax.text(5.8, yk, "antiparticle", va="center", fontsize=7, color=INK)
    ax.add_patch(FancyBboxPatch((7.4, yk - 0.06), 0.5, 0.12,
                                boxstyle="round,pad=0,rounding_size=0.04", fc=FIELD, ec="none"))
    ax.text(8.05, yk, "electric field", va="center", fontsize=7, color=INK)

    save(fig, FIG / "13_lattice_states.png")
    print(f"wrote {FIG / '13_lattice_states.png'} (and .pdf)")
    print("energies in units of w at m = 0.5 w, J = w; every assert passed")


if __name__ == "__main__":
    main()
