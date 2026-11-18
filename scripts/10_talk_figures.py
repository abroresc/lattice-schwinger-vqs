"""
Figures for the capstone talk. One message per figure.

The thesis figures (scripts 04 to 08) are built to be read on paper: three
panels, a legend, a caption. A slide gets a few seconds from the back of a
room, so each figure here keeps one idea from one of those panels, labels the
lines directly instead of using a legend, and drops everything else.

Nothing new is computed. Every curve comes from the same tested functions the
thesis scripts call, at the same parameters, and the script prints each number
the talk quotes next to the value in the thesis text. If the talk and the
thesis ever drift apart, it shows up here as a residual, not in front of an
audience.

Output: figures/talk/*.png, saved at twice the size they appear on a
1920 x 1080 slide, with transparent backgrounds so they sit on the slide colour.
Figure sizes below are in slide pixels; fonts are set so that no label is
smaller than 26 px on the slide.

Run:  .venv/bin/python scripts/10_talk_figures.py
"""
import pathlib
import sys

import numpy as np
from scipy.interpolate import CubicSpline

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from schwinger.budget import (
    budget, exact_distribution, first_crossing, point_circuit, sweep_dip,
)
from schwinger.circuits import circuit_metrics, statevector, trotter_circuit
from schwinger.diagnostics import (
    observables_from_distribution, shot_noise_sigma, z_basis_table,
)
from schwinger.noise import (
    gate_noise_model, prepare_for_probabilities, probabilities, sample_counts,
    step_channel_probabilities,
)
from schwinger.finite_size import first_peak_and_trough, quench_density
from schwinger.observables import evolve, particle_density
from schwinger import replication as R
from schwinger.schwinger import bare_vacuum, hamiltonian

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = pathlib.Path(__file__).resolve().parents[1] / "figures" / "talk"
N, M, W, J = 4, 0.5, 1.0, 1.0
T_PEAK, T_TROUGH = 0.89, 1.76          # first maximum and minimum of the exact nu(t)
SHOTS = 1000

# The deck's palette. Kept here rather than in plotstyle.py because it belongs
# to the slides, not the thesis; the thesis figures keep their own palette.
# The deck has white slides, a dark title slide, Avenir Next, and a vivid blue
# and orange: blue means the exact answer, orange means error or cost.
INK = "#111827"
SOFT = "#4B5563"
BLUE = "#2563EB"
ORANGE = "#EA580C"
GREY = "#9CA3AF"
FAINT = "#EEF2F7"
RAMP = ("#FB923C", "#EA580C", "#9A3412")   # light to dark orange: more gate error
FONT = "Avenir Next"
FAILED = "#9B1C1C"                          # a failed gate: dark red, distinct from every curve


def px(n: float) -> float:
    """Font size in points that shows as n pixels on the slide (100 slide px per inch)."""
    return n * 0.72


def style() -> None:
    plt.rcParams.update({
        "font.family": FONT,
        "font.size": px(28),
        "axes.labelsize": px(30),
        "axes.labelcolor": SOFT,
        "axes.edgecolor": SOFT,
        "axes.linewidth": 1.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "xtick.labelsize": px(26),
        "ytick.labelsize": px(26),
        "xtick.color": SOFT,
        "ytick.color": SOFT,
        "xtick.major.width": 1.6,
        "ytick.major.width": 1.6,
        "xtick.major.size": 7,
        "ytick.major.size": 7,
        "text.color": INK,
        "lines.linewidth": 4.0,
        "lines.markersize": 11,
        "savefig.transparent": True,
    })


def slide_fig(w_px: int, h_px: int):
    fig, ax = plt.subplots(figsize=(w_px / 100, h_px / 100))
    return fig, ax


def save(fig, name: str) -> None:
    # A fixed canvas, not bbox="tight": the slide reserves an exact box for
    # each figure, so the saved image must have exactly that aspect ratio.
    fig.savefig(OUT / f"{name}.png", dpi=200)
    plt.close(fig)
    print(f"  wrote figures/talk/{name}.png")


def steps_for(t: float, dt: float) -> int:
    return max(1, int(round(t / dt)))


def check(label: str, got: float, thesis: float) -> None:
    """Print a quoted number beside the thesis value. A residual, not a flag."""
    print(f"  {label:<52} {got:>9.4g}   thesis {thesis:<8g} diff {got - thesis:+.1e}")


# ----------------------------------------------------------------------------
# 1. The signal: pairs appear, then recombine
# ----------------------------------------------------------------------------


def fig_pairs():
    times = np.linspace(0.0, 4.0, 801)
    nu = np.array([particle_density(s, N)
                   for s in evolve(bare_vacuum(N), hamiltonian(N, M, W, J), times)])
    ip = int(np.argmin(np.abs(times - T_PEAK)))
    it = int(np.argmin(np.abs(times - T_TROUGH)))
    check("dip D = nu(t_peak) - nu(t_trough), exact", nu[ip] - nu[it], 0.33)

    # Every canvas is the size of its box on the slide (100 px an inch), so px()
    # is exact.
    fig, ax = slide_fig(1213, 495)
    fig.subplots_adjust(left=0.1, right=0.83, bottom=0.2, top=0.95)
    ax.plot(times, nu, color=INK)
    # The signal is blue, the deck's colour for the exact answer; orange is kept
    # for errors.
    ax.plot([times[ip], times[it]], [nu[ip], nu[it]], "o", color=BLUE, markersize=16, zorder=5)
    # The arrow stands just past t = 4, so it crosses no part of the curve: after
    # the peak the curve stays below 0.46, and after the trough above 0.20.
    xa = 4.12
    for y in (nu[ip], nu[it]):
        ax.plot([times[ip] if y == nu[ip] else times[it], xa], [y, y], color=BLUE,
                linewidth=1.6, linestyle=(0, (2, 3)), clip_on=False)
    ax.annotate("", xy=(xa, nu[ip]), xytext=(xa, nu[it]), annotation_clip=False,
                arrowprops=dict(arrowstyle="<->", color=BLUE, lw=3.5, shrinkA=0, shrinkB=0))
    ax.text(xa + 0.07, (nu[ip] + nu[it]) / 2, "the dip:\nour signal", color=BLUE, fontsize=px(30),
            va="center", clip_on=False)
    ax.text(times[ip] - 0.12, nu[ip] + 0.035, "pairs appear", ha="center", fontsize=px(30))
    ax.text(times[it] + 0.0, nu[it] - 0.075, "pairs recombine", ha="center", fontsize=px(30))
    ax.set_xlim(0, 4)
    ax.set_ylim(-0.02, 0.56)
    ax.set_yticks([0, 0.25, 0.5])
    ax.set_xlabel("time")
    ax.set_ylabel("particles per site")
    save(fig, "pairs")


def fig_revival():
    """Backup: the same exact curve run ten times longer, to t = 40.

    Why it exists: on the main signal slide the swings shrink between t = 0 and 4, and the
    natural question is "if pairs recombine, why does it damp?". It does not.
    The vacuum overlaps five energy eigenstates, so nu(t) is a constant plus ten
    cosines at their energy differences. They start in phase (nu = 0 exactly)
    and drift apart, which looks like damping over a short window. Evolution
    is unitary, so nothing is lost: near t = 25 the curve climbs above the
    first peak. Real, permanent damping appears only with gate noise.
    """
    times = np.linspace(0.0, 40.0, 8001)
    nu = np.array([particle_density(s, N)
                   for s in evolve(bare_vacuum(N), hamiltonian(N, M, W, J), times)])
    i1 = int(np.argmax(np.where(times <= 2.0, nu, -1)))
    i2 = int(np.argmax(np.where(times > 3.0, nu, -1)))
    # The long-time average is the diagonal-ensemble value: the sum over
    # eigenstates of |overlap|^2 times <n|nu|n>, 0.258 from the spectral
    # decomposition. The mean over [0, 40] sits near it, not on it,
    # because 40 is a finite window.
    check("first peak of nu(t), exact", nu[i1], 0.4631)
    check("highest point after t = 3, exact", nu[i2], 0.476)
    check("time of that point", times[i2], 25.2)
    check("mean of nu over 0 <= t <= 40", float(np.mean(nu)), 0.258)

    fig, ax = slide_fig(730, 500)
    fig.subplots_adjust(left=0.14, right=0.85, bottom=0.2, top=0.95)
    # The window the main slide shows, so the audience can find it here.
    ax.axvspan(0, 4, color=FAINT, zorder=0)
    ax.text(0.4, 0.6, "slide 11", ha="left", va="center", color=SOFT, fontsize=px(22))
    # A thinner line than on the main slide: ten times more time in a narrower box.
    ax.plot(times, nu, color=INK, linewidth=2.2, zorder=3)
    ax.axhline(0.258, color=GREY, linewidth=1.6, linestyle=(0, (2, 3)), zorder=2)
    # The label sits in the right margin: the curve crosses the line everywhere
    # inside the axes, so any spot there would cover part of it.
    ax.text(40.8, 0.258, "average\n0.26", ha="left", va="center", color=SOFT, fontsize=px(22),
            clip_on=False)
    ax.plot([times[i1], times[i2]], [nu[i1], nu[i2]], "o", color=BLUE, markersize=13, zorder=5)
    ax.text(times[i1] + 1.2, nu[i1] + 0.01, f"first peak {nu[i1]:.2f}", color=BLUE,
            va="bottom", fontsize=px(24))
    ax.text(times[i2], nu[i2] + 0.03, f"back to {nu[i2]:.2f}", color=BLUE, ha="center",
            va="bottom", fontsize=px(24))
    ax.set_xlim(0, 40)
    ax.set_ylim(-0.02, 0.63)
    ax.set_xticks([0, 10, 20, 30, 40])
    ax.set_yticks([0, 0.25, 0.5])
    ax.set_xlabel("time")
    ax.set_ylabel("particles per site")
    save(fig, "revival")


# ----------------------------------------------------------------------------
# 2. Error one: chopping time into steps
# ----------------------------------------------------------------------------


def fig_steps():
    times = np.linspace(0.0, 4.0, 401)
    exact = np.array([particle_density(s, N)
                      for s in evolve(bare_vacuum(N), hamiltonian(N, M, W, J), times)])
    fig, ax = slide_fig(1213, 475)
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.2, top=0.95)
    ax.plot(times, exact, color=INK, linewidth=3.0)
    # Small steps as dots only, so the black exact curve shows between them;
    # the label points at the black line itself.
    ax.annotate("exact", xy=(3.3, float(np.interp(3.3, times, exact)) - 0.01), xytext=(3.45, 0.06),
                ha="center", va="center", fontsize=px(30), color=INK,
                arrowprops=dict(arrowstyle="-", color=INK, lw=1.6, shrinkA=10, shrinkB=0))
    for dt, color, label, xy, line in ((0.8, ORANGE, "big steps", (2.45, 0.53), 2.6),
                                       (0.2, BLUE, "small steps", (1.95, 0.04), 0)):
        ks = np.arange(0, int(round(4.0 / dt)) + 1)
        nu = [0.0] + [particle_density(statevector(trotter_circuit(N, M, W, J, k * dt, int(k), 2)), N)
                      for k in ks[1:]]
        ax.plot(ks * dt, nu, color=color, marker="o", linewidth=line, markersize=13, clip_on=False, zorder=4)
        ax.text(*xy, label, color=color, fontsize=px(32))
        k_tr = int(np.argmin(np.abs(ks * dt - T_TROUGH)))
        i_ex = int(np.argmin(np.abs(times - ks[k_tr] * dt)))
        print(f"  dt = {dt}: |Trotter - exact| at t = {ks[k_tr] * dt:.1f} is {abs(nu[k_tr] - exact[i_ex]):.3f}")
    ax.set_xlim(0, 4)
    ax.set_ylim(-0.02, 0.62)
    ax.set_yticks([0, 0.25, 0.5])
    ax.set_xlabel("time")
    ax.set_ylabel("particles per site")
    save(fig, "steps")


# ----------------------------------------------------------------------------
# 3. Error two: finite shots
# ----------------------------------------------------------------------------


def fig_shots():
    p = exact_distribution(N, M, W, J, T_TROUGH)
    truth = observables_from_distribution(p, N)["particle_density"]
    rng = np.random.default_rng(20261025)
    fig, ax = slide_fig(740, 495)
    fig.subplots_adjust(left=0.17, right=0.98, bottom=0.13, top=0.92)
    ax.axhline(truth, color=INK, linewidth=2.4, linestyle=(0, (5, 4)))
    ax.text(2.98, truth + 0.004, "true value", fontsize=px(28), color=INK, va="bottom", ha="right")
    for i, shots in enumerate((250, 1000, 4000)):
        est = np.array([observables_from_distribution(sample_counts(p, shots, rng) / shots, N)
                        ["particle_density"] for _ in range(60)])
        x = i + rng.uniform(-0.22, 0.22, size=est.size)
        ax.plot(x, est, "o", color=ORANGE, alpha=0.7, markersize=12, markeredgewidth=0)   # an error: orange
        sigma = shot_noise_sigma(p, N, shots)
        # Four decimals, so the halving shows on screen (0.0146, 0.0073, 0.0036).
        ax.text(i, truth + 0.075, f"± {sigma:.4f}", ha="center", fontsize=px(28), color=ORANGE)
        print(f"  {shots:>5} shots: exact sigma {sigma:.4f}, scatter of 60 runs {est.std(ddof=1):.4f}")
    ax.set_xticks([0, 1, 2], ["250 shots", "1000 shots", "4000 shots"])
    ax.tick_params(axis="x", length=0, labelsize=px(30), labelcolor=INK)
    ax.set_xlim(-0.5, 3.0)
    ax.set_ylim(truth - 0.07, truth + 0.10)
    ax.set_yticks([0.1, 0.15, 0.2])
    ax.set_ylabel("measured particles per site")
    ax.spines["bottom"].set_visible(False)
    save(fig, "shots")


# ----------------------------------------------------------------------------
# 4. Error three: every two-qubit gate can fail
# ----------------------------------------------------------------------------


def fig_gate_grid():
    """
    The trough circuit at dt = 0.1: 18 Trotter steps of 32 CX gates each, one
    row per step, drawn at 0.1 % per gate so it explains the blue curve beside it
    on the Error 3 slide. A two-qubit depolarizing channel with average
    infidelity r applies a non-identity Pauli with probability 1.25 r (thesis
    Discussion), so a run expects 1.25 x 0.001 x 576 = 0.72 failed gates and
    about half the runs come through clean (the trough circuit's fidelity is
    0.54, and the dip keeps 51 % of its size). One square is highlighted, at a
    seeded random place: one run that was hit.
    """
    pc = point_circuit(N, M, W, J, T_TROUGH, steps_for(T_TROUGH, 0.1))
    rows, cols = steps_for(T_TROUGH, 0.1), 32
    check("CX gates in the trough circuit, dt = 0.1", pc.two_qubit_gates, 576)
    assert rows * cols == pc.two_qubit_gates
    expected = 1.25 * 0.001 * rows * cols
    print(f"  P(no failed gate in a run) = exp(-{expected:.2f}) = {np.exp(-expected):.2f}")
    fails = np.zeros(rows * cols, dtype=bool)
    fails[np.random.default_rng(11).choice(rows * cols, int(round(expected)), replace=False)] = True
    print(f"  highlighted failures: {fails.sum()} (expected at 0.1 %: {expected:.2f})")
    fig, ax = slide_fig(1380, 780)
    fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
    for k in range(rows * cols):
        r, c = divmod(k, cols)
        ax.add_patch(plt.Rectangle((c, rows - 1 - r), 0.8, 0.8,
                                   color=FAILED if fails[k] else FAINT, linewidth=0))
    # A ring around each failed gate: one small square in 576 is easy to miss.
    for k in np.flatnonzero(fails):
        r, c = divmod(int(k), cols)
        ax.add_patch(plt.Circle((c + 0.4, rows - 1 - r + 0.4), 1.6, fill=False, color=FAILED, linewidth=4,
                                clip_on=False))
    ax.set_xlim(-0.1, cols)
    ax.set_ylim(-0.1, rows)
    ax.set_aspect("equal")
    ax.axis("off")
    save(fig, "gate_grid")


# ----------------------------------------------------------------------------
# 5. Result: gate error erases the dip
# ----------------------------------------------------------------------------


def fig_noise():
    # The exact curve on a fine grid like the other plots; the noisy curves at
    # every time step of 0.05 (81 circuits), which is smooth at slide size.
    fine = np.linspace(0.0, 4.0, 401)
    times = np.linspace(0.0, 4.0, 41)    # k steps of exactly 0.1: the gate count grows by 32 a point
    exact = [particle_density(s, N) for s in evolve(bare_vacuum(N), hamiltonian(N, M, W, J), fine)]
    preps = [prepare_for_probabilities(trotter_circuit(N, M, W, J, t, steps_for(t, 0.1), 2))
             if t > 0 else None for t in times]
    # Narrow, to sit beside the gate grid on the Error 3 slide, with the labels
    # at the line ends, where there is room.
    fig, ax = slide_fig(840, 500)
    fig.subplots_adjust(left=0.13, right=0.73, bottom=0.17, top=0.95)
    ax.axhline(0.5, color=GREY, linewidth=2.0, linestyle=(0, (5, 4)))
    ax.text(0.05, 0.525, "fully scrambled", color=SOFT, fontsize=px(26), va="bottom")
    ax.plot(fine, exact, color=INK)
    ax.set_xticks(np.arange(0, 4.01, 0.5))
    ax.text(4.1, exact[-1], "no gate error", color=INK, fontsize=px(28), va="center", clip_on=False)
    # Both noisy curves in the orange ramp, darker = more gate error, as on
    # the budget plot; blue stays for the exact answer.
    for eps, color, label, dy in ((1e-3, RAMP[0], "0.1 % per gate", -0.015),
                                  (1e-2, RAMP[2], "1 % per gate", 0.04)):
        nm = gate_noise_model(eps)
        nu = [0.0 if p is None else observables_from_distribution(probabilities(p, nm), N)["particle_density"]
              for p in preps]
        # A cubic through the computed points: between them the step count is
        # not defined, and straight segments drew kinks the data do not have.
        ax.plot(fine, CubicSpline(times, nu)(fine), color=color)
        ax.text(4.1, nu[-1] + dy, label, color=color, fontsize=px(28), va="center", clip_on=False)
    ax.set_xlim(0, 4)
    ax.set_ylim(-0.02, 0.6)
    ax.set_yticks([0, 0.25, 0.5])
    ax.set_xlabel("time")
    ax.set_ylabel("particles per site")
    save(fig, "noise")


# ----------------------------------------------------------------------------
# 6. Result: the crossover, Trotter error against shot noise
# ----------------------------------------------------------------------------


def fig_crossover():
    # Measured at the dip's trough (time 1.76), where every other talk result
    # lives; the abstract's time-3 crossover (step 0.182, about 17 steps) is in the notes.
    T = T_TROUGH
    exact_state = evolve(bare_vacuum(N), hamiltonian(N, M, W, J), np.array([T]))[0]
    exact_nu = particle_density(exact_state, N)
    steps = np.array(list(range(2, 33)) + [40, 48, 64, 96, 128, 192, 256])
    bias, gates = [], []
    for s in steps:
        qc = trotter_circuit(N, M, W, J, T, int(s), order=2)
        bias.append(abs(particle_density(statevector(qc), N) - exact_nu))
        gates.append(circuit_metrics(qc)["two_qubit_gates"])
    bias, gates = np.array(bias), np.array(gates)
    sigma = shot_noise_sigma(np.abs(exact_state) ** 2, N, SHOTS)
    dt_x = first_crossing(T / steps, bias - sigma)
    check("crossover dt at the trough, 1000 shots", dt_x, 0.36)
    s_x = T / dt_x

    fig, ax = slide_fig(1213, 485)
    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.2, top=0.95)
    ax.axvspan(s_x, 160, color=FAINT, alpha=0.55, linewidth=0)
    ax.loglog(steps, np.maximum(bias, 1e-12), color=INK)    # an error, so not blue
    ax.axhline(sigma, color=ORANGE, linewidth=3.2, linestyle=(0, (6, 4)))
    ax.text(22.0, sigma * 1.4, "shot noise, 1,000 shots", color=ORANGE, fontsize=px(28))
    ax.text(22.0, 1.3e-3, "time-step error", color=INK, fontsize=px(28))
    ax.axvline(s_x, color=INK, linewidth=2.0)
    ax.text(s_x * 1.3, 0.3, "extra steps here:\ninvisible gain,\nmore gates",
            fontsize=px(28), color=INK, va="top")
    # A dot and a label at three step counts, so "gates" is never mistaken for
    # steps: 32 two-qubit gates a step. 17 is the step count nearest the crossing.
    for s_mark, off in ((5, (14, -74)),):
        k = int(np.where(steps == s_mark)[0][0])
        ax.plot(steps[k], bias[k], "o", color=INK, markersize=11, zorder=5)
        ax.annotate(f"{s_mark} steps: {gates[k]:,} gates", (steps[k], bias[k]), textcoords="offset points",
                    xytext=off, fontsize=px(24), color=SOFT)
    ax.set_xlim(1.8, 160)
    ax.set_ylim(1e-5, 0.6)
    ax.set_yticks([1e-4, 1e-3, 1e-2, 1e-1], ["0.0001", "0.001", "0.01", "0.1"])
    ax.set_xticks([2, 5, 10, 20, 50, 100], ["2", "5", "10", "20", "50", "100"])
    ax.minorticks_off()
    ax.set_xlabel("number of time steps")
    ax.set_ylabel("error in particles per site")
    save(fig, "crossover")
    return s_x


# ----------------------------------------------------------------------------
# 7. Result: the failure threshold is set by the expected number of gate errors
# ----------------------------------------------------------------------------


def threshold_sweeps():
    eps = np.logspace(-5, -1.5, 36)
    out = {}
    for dt in (0.1, 0.2, 0.3):
        pk = point_circuit(N, M, W, J, T_PEAK, steps_for(T_PEAK, dt))
        tr = point_circuit(N, M, W, J, T_TROUGH, steps_for(T_TROUGH, dt))
        sw = sweep_dip(pk, tr, N, eps, SHOTS)
        star = first_crossing(eps, sw["D"] - 2 * sw["sigma_D"])
        out[dt] = (pk, tr, sw, star)
    for dt, thesis in ((0.1, 1.92), (0.2, 1.98), (0.3, 2.03)):
        check(f"eps* x n_CX, dt = {dt}, 1000 shots", out[dt][3] * out[dt][1].two_qubit_gates, thesis)
    return eps, out


def fig_collapse(eps, sweeps):
    """
    Two panels, same y axis. Left, against the error rate alone: three
    separate curves, each losing the signal at its own rate (0.33, 0.69 and
    1.06 %). Right, against error rate x number of two-qubit gates: one curve,
    crossing at about 2. The audience sees the "before" next to the collapse.
    """
    fig, (a, b) = plt.subplots(1, 2, figsize=(12.13, 4.85), sharey=True)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.2, top=0.9, wspace=0.12)
    colors = (RAMP[2], RAMP[1], RAMP[0])   # darker = more gates, so more damage
    floor = np.mean([sw["sigma_D"] for (_, _, sw, _) in sweeps.values()], axis=0) * 2
    for (dt, (pk, tr, sw, star)), color in zip(sweeps.items(), colors):
        n = tr.two_qubit_gates
        a.semilogx(eps, sw["D"], color=color, linewidth=3.4)
        a.plot(star, 2 * np.interp(star, eps, sw["sigma_D"]), "o", color=color, markersize=11, zorder=5)
        b.semilogx(eps * n, sw["D"], color=color, linewidth=3.4)
    # Labels where each curve has reached half of its starting dip.
    # Labels as a coloured list in the empty lower left; the dots mark where
    # each circuit loses the signal (0.33, 0.69 and 1.06 % per gate).
    for (dt, (pk, tr, sw, star)), color, y in zip(sweeps.items(), colors, (0.19, 0.14, 0.09)):
        a.text(1.2e-4, y, f"{tr.two_qubit_gates} gates", color=color, fontsize=px(26))
    a.semilogx(eps, floor, color=GREY, linewidth=2.4, linestyle=(0, (5, 4)))
    b.semilogx(eps * sweeps[0.2][1].two_qubit_gates, floor, color=GREY, linewidth=2.4,
               linestyle=(0, (5, 4)))
    b.text(0.012, 0.035, "2 × shot noise", color=SOFT, fontsize=px(24))
    b.axvline(2.0, color=INK, linewidth=2.0)
    # "About 2" is the summed CX infidelity, eps x n_CX. For depolarizing noise
    # that is about 2.5 actual Pauli errors (lambda = 4 eps / 3, and 15/16 of
    # lambda is a non-identity Pauli), so the label says the product, not "errors".
    b.text(2.25, 0.29, "≈ 2", fontsize=px(44), color=INK)
    a.set_title("against error rate", fontsize=px(26), color=SOFT, loc="left")
    b.set_title("against error rate × gates", fontsize=px(26), color=SOFT, loc="left")
    a.set_xlim(1e-4, 3e-2)
    a.set_xticks([1e-4, 1e-3, 1e-2], ["0.01 %", "0.1 %", "1 %"])
    b.set_xlim(0.01, 10)
    b.set_xticks([0.01, 0.1, 1, 2, 10], ["0.01", "0.1", "1", "2", "10"])
    a.set_ylim(-0.03, 0.4)
    a.set_yticks([0, 0.1, 0.2, 0.3])
    for ax in (a, b):
        ax.minorticks_off()
    a.set_xlabel("error per two-qubit gate")
    b.set_xlabel("error per gate × two-qubit gates")
    a.set_ylabel("size of the dip")
    save(fig, "collapse")


# ----------------------------------------------------------------------------
# 8. Result: the best step size depends on the hardware
# ----------------------------------------------------------------------------


def fig_budget():
    steps = np.array([2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 48, 64])
    circuits = {int(s): point_circuit(N, M, W, J, T_TROUGH, int(s)) for s in steps}
    fig, ax = slide_fig(930, 470)
    fig.subplots_adjust(left=0.14, right=0.97, bottom=0.2, top=0.95)
    # (eps, colour, end label, label height). Heights are set by hand because
    # the two worst curves end almost together.
    lines = ((0.0, INK, "perfect gates", None), (3e-4, RAMP[0], "0.03 %", None),
             (1e-3, RAMP[1], "0.1 %", 0.27), (3e-3, RAMP[2], "0.3 % per gate", 0.42))
    best = {}
    for eps, color, label, y_label in lines:
        nm = gate_noise_model(eps) if eps > 0 else None
        rms = np.array([budget(circuits[int(s)], N, nm, 0.01, SHOTS).rms for s in steps])
        k = int(np.argmin(rms))
        best[eps] = (steps[k], T_TROUGH / steps[k], circuits[int(steps[k])].two_qubit_gates)
        ax.loglog(steps, rms, color=color, linewidth=3.6)
        if eps > 0:
            ax.plot(steps[k], rms[k], marker="*", markersize=34, color=color,
                    markeredgecolor=INK, markeredgewidth=1.5, zorder=5)
        # Labels in ink at each line's end (the light-orange "0.03 %" was 2.2:1 on
        # white, too faint to read); the line beside it carries the colour.
        ax.text(steps[-1] * 1.08, y_label or rms[-1], label, color=SOFT if eps else INK, fontsize=px(28),
                va="center")
    for eps, thesis in ((3e-4, 0.35), (1e-3, 0.44), (3e-3, 0.59)):
        check(f"optimal dt, eps2 = {eps:g}", best[eps][1], thesis)
        print(f"      = {best[eps][0]} steps, {best[eps][2]} CX")
    ax.set_xlim(1.8, 64 * 3.2)
    ax.set_ylim(5e-3, 0.6)
    ax.set_yticks([0.01, 0.03, 0.1, 0.3], ["0.01", "0.03", "0.1", "0.3"])
    ax.set_xticks([2, 3, 4, 5, 10, 20, 50], ["2", "3", "4", "5", "10", "20", "50"])
    ax.minorticks_off()
    ax.set_xlabel("number of time steps")
    ax.set_ylabel("error in particles per site")
    save(fig, "budget")
    return best


# ----------------------------------------------------------------------------
# 9. The warning light: total charge per shot
# ----------------------------------------------------------------------------


def fig_charge(sweeps):
    pk, tr, sw, star = sweeps[0.1]
    table = z_basis_table(N)
    qs = np.arange(-2, 3)
    clean = tr.trotter
    noisy = probabilities(tr.prepared, gate_noise_model(3e-3))
    fig, axes = plt.subplots(1, 2, figsize=(12.13, 3.6), sharey=True)
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.25, top=0.85, wspace=0.12)
    for ax, p, title, color in ((axes[0], clean, "perfect machine", BLUE),
                                (axes[1], noisy, "noisy: 0.3 % per gate, 576 gates", ORANGE)):
        share = np.array([p[table.charge == q].sum() for q in qs])
        ax.bar(qs, 100 * share, color=color, width=0.7)
        ax.set_title(title, fontsize=px(28), color=color, loc="left")
        ax.set_xticks(qs, ["−2", "−1", "0", "+1", "+2"])
        ax.set_xlabel("total charge read in one shot")
        ax.tick_params(axis="x", labelsize=px(26))
        for q, s in zip(qs, share):
            if s > 0.005:
                ax.text(q, 100 * s + 3, f"{100 * s:.0f} %", ha="center", fontsize=px(24), color=INK)
        print(f"  {title:<32}: P(Q = 0) = {share[2]:.4f}, Var(Q) = "
              f"{observables_from_distribution(p, N)['charge_variance']:.3f}")
    axes[0].set_ylabel("share of shots")
    axes[0].set_ylim(0, 112)
    axes[0].set_yticks([0, 50, 100], ["0", "50 %", "100 %"])
    save(fig, "charge")


# ----------------------------------------------------------------------------
# 10. The blind spot: noise the warning light cannot see
# ----------------------------------------------------------------------------


def fig_blind(eps, sweeps):
    """
    The blind spot in two panels the audience already knows: the dip and the
    share of shots that read the wrong total charge. Left, gate noise that flips qubits: as the dip shrinks,
    more shots read the wrong charge, so the warning light comes on. Right,
    phase noise between steps: the dip shrinks the same way, but every shot
    still reads charge 0, so the light never comes on.

    Data: the 1,000-shot sweep at dt = 0.1 and the "dephasing" step channel.
    """
    pk, tr, sw, star = sweeps[0.1]
    table = z_basis_table(N)
    d0 = sw["D"][0]
    zero = table.charge == 0
    # Left: the dip and the wrong-charge share at the trough, gate noise.
    wrong_gate = np.array([1 - probabilities(tr.prepared, gate_noise_model(e))[zero].sum() for e in eps])
    # Right: the same two numbers for phase kicks between steps.
    ps = np.logspace(-4, -0.5, 30)
    dip_ph, wrong_ph = [], []
    for p in ps:
        a = step_channel_probabilities(N, M, W, J, T_PEAK, steps_for(T_PEAK, 0.1), 2, "dephasing", p)
        b = step_channel_probabilities(N, M, W, J, T_TROUGH, steps_for(T_TROUGH, 0.1), 2, "dephasing", p)
        dip_ph.append(observables_from_distribution(a, N)["particle_density"]
                      - observables_from_distribution(b, N)["particle_density"])
        wrong_ph.append(1 - b[zero].sum())
    dip_ph, wrong_ph = np.array(dip_ph), np.array(wrong_ph)
    # Checks: the talk's noisy charge bars (0.3 % per gate) say 44 % of shots read 0;
    # phase noise never moves a shot off charge 0.
    p3 = probabilities(tr.prepared, gate_noise_model(3e-3))
    print(f"  gate noise 0.3 %: shots reading charge 0 = {p3[zero].sum():.4f} (talk's charge bar: 44 %)")
    print(f"  phase noise: largest share of shots with charge != 0 = {wrong_ph.max():.1e}, "
          f"dip {dip_ph[0] / d0:.3f} -> {dip_ph[-1] / d0:+.3f} of its noiseless size")
    print(f"  gate noise: wrong-charge share {wrong_gate[0]:.1e} -> {wrong_gate[-1]:.3f}, "
          f"dip {sw['D'][0] / d0:.3f} -> {sw['D'][-1] / d0:+.3f}")

    fig, (l, r) = plt.subplots(1, 2, figsize=(8.2, 4.85), sharey=True)
    fig.subplots_adjust(left=0.11, right=0.98, bottom=0.2, top=0.88, wspace=0.1)
    pct = lambda v: 100 * np.clip(v, 0, None)
    for ax, x, dip, wrong, title, tcol, xlabel in (
            (l, eps, sw["D"] / d0, wrong_gate, "noise that flips qubits", BLUE, "error per two-qubit gate"),
            (r, ps, dip_ph / d0, wrong_ph, "phase noise", ORANGE, "phase kicks per qubit per step")):
        ax.semilogx(x, pct(dip), color=INK, linewidth=3.4)
        ax.semilogx(x, pct(wrong), color=ORANGE, linewidth=3.4)
        ax.set_title(title, color=tcol, fontsize=px(28), loc="left", fontweight="bold")
        ax.set_xlabel(xlabel, fontsize=px(24))
        lo, hi = x[0], x[-1]
        ticks = [t for t in (1e-4, 1e-3, 1e-2, 1e-1) if lo <= t <= hi]
        ax.set_xticks(ticks, [f"{100 * t:g} %" for t in ticks])
        ax.tick_params(axis="x", labelsize=px(22))
    l.text(1.3e-5, 88, "dip left", color=INK, fontsize=px(26))
    l.text(1.6e-3, 67, "shots with\nthe wrong\ncharge", color=ORANGE, fontsize=px(26), va="bottom")
    r.text(1.3e-4, 88, "dip left", color=INK, fontsize=px(26))
    r.text(1.2e-4, 6, "every shot\nreads charge 0", color=ORANGE, fontsize=px(26), va="bottom")
    l.set_ylim(-5, 108)
    l.set_yticks([0, 50, 100], ["0", "50 %", "100 %"])
    save(fig, "blind")


# ----------------------------------------------------------------------------
# Backup slide: the Hamiltonian, thesis Eq. (eq:H), typeset for projection
# ----------------------------------------------------------------------------


def fig_hamiltonian():
    """Typed from theory.tex Eq. (eq:H); the slide labels the three terms underneath."""
    fig = plt.figure(figsize=(15.6, 4.0))
    with plt.rc_context({"mathtext.fontset": "cm"}):
        fig.text(0.5, 0.64, r"$H = \frac{m}{2}\sum_{n=1}^{N}(-1)^n Z_n"
                 r" + w\sum_{n=1}^{N-1}(\sigma^{+}_{n}\,\sigma^{-}_{n+1} + \mathrm{h.c.})"
                 r" + J\sum_{n=1}^{N-1} L_n^2$", ha="center", va="center", fontsize=px(64), color=INK)
        fig.text(0.5, 0.15, r"$L_n = \frac{1}{2}\sum_{k \leq n}\left[Z_k + (-1)^k\right]$",
                 ha="center", va="center", fontsize=px(52), color=SOFT)
    save(fig, "hamiltonian")


# ----------------------------------------------------------------------------
# 11. Check: four sites against longer chains (talk version of script 12)
# ----------------------------------------------------------------------------


def fig_chains():
    """
    Script 12's left panel, cut to one message: the first rise and dip are the
    same on every chain length, and only later does four sites fall behind.
    Twelve sites is the longest chain; 6, 8 and 10 are drawn faint because
    they sit between 4 and 12 and would only add a legend. The shaded span is
    the stretch the thesis uses, from time 0 to the first trough.
    """
    times = np.linspace(0.0, 4.0, 801)
    curves = {n: quench_density(n, M, W, J, times) for n in (4, 6, 8, 10, 12)}
    thesis_dip = {4: 0.330, 6: 0.354, 8: 0.358, 10: 0.360, 12: 0.361}
    for n, nu in curves.items():
        _, yp, _, yt = first_peak_and_trough(times, nu)
        check(f"dip on a chain of {n} sites", yp - yt, thesis_dip[n])
    k2 = times > 2.0
    print(f"  second peak: 4 sites {curves[4][k2].max():.3f}, 12 sites {curves[12][k2].max():.3f}")

    fig, ax = slide_fig(860, 500)
    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.17, top=0.95)
    ax.axvspan(0, T_TROUGH, color=FAINT, alpha=0.55, linewidth=0)
    ax.text(T_TROUGH / 2, 0.585, "the part I use", ha="center", fontsize=px(28), color=SOFT)
    for n in (6, 8, 10):
        ax.plot(times, curves[n], color=GREY, linewidth=2.0)
    ax.plot(times, curves[12], color=BLUE, linewidth=3.4)
    ax.plot(times, curves[4], color=INK, linewidth=4.0)
    ax.text(2.62, 0.495, "12 sites", color=BLUE, fontsize=px(30), ha="center")
    # Below every curve, with a leader up to the four-site second peak: the
    # 12-site curve is 0.11 higher there, so the leader cannot be misread.
    i2 = int(np.argmax(np.where(k2, curves[4], -1)))
    ax.annotate("4 sites (mine)", xy=(times[i2], curves[4][i2] - 0.012), xytext=(times[i2], 0.045),
                ha="center", va="center", color=INK, fontsize=px(30),
                arrowprops=dict(arrowstyle="-", color=INK, lw=1.6, shrinkA=14, shrinkB=0))
    ax.set_xlim(0, 4)
    ax.set_ylim(-0.02, 0.63)
    ax.set_yticks([0, 0.25, 0.5])
    ax.set_xlabel("time")
    ax.set_ylabel("particles per site")
    save(fig, "chains")


# ----------------------------------------------------------------------------
# 12. Check: Martinez et al. (2016) Fig. 2(a), recomputed (talk version of 11)
# ----------------------------------------------------------------------------


def fig_replica():
    """
    Script 11's figure with the legend replaced by labels on the lines. Lines
    are this code; open symbols are the theory curves read off their Fig. 2(a);
    grey dots are their measured points, which nothing here is fitted to. The
    x axis counts their time steps of w*dt = pi/8, which is where their points sit.
    """
    ks = np.arange(1, 5)
    trot = R.trotter_density()
    model = R.error_model_density()
    ideal = np.array(R.FIG2A_IDEAL)
    fine = np.linspace(0, np.pi / 2, 601)
    exact = R.exact_density(fine)
    worst = max(np.max(np.abs(R.exact_density(ideal[:, 0]) - ideal[:, 1])),
                np.max(np.abs(trot - np.array(R.FIG2A_TROTTER))),
                np.max(np.abs(model - np.array(R.FIG2A_MODEL))))
    check("largest gap to their three theory curves", worst, 0.001)

    step = lambda t: t / R.STEP
    fig, ax = slide_fig(760, 520)
    fig.subplots_adjust(left=0.13, right=0.70, bottom=0.17, top=0.95)
    ax.plot(step(fine), exact, color=INK, linewidth=3.0)
    ax.plot(step(ideal[:, 0]), ideal[:, 1], "o", color=INK, markersize=4.5)
    ax.plot(np.r_[0, ks], np.r_[0, trot], color=ORANGE, linewidth=3.0)
    ax.plot(ks, R.FIG2A_TROTTER, "^", markersize=17, markerfacecolor="none",
            markeredgecolor=ORANGE, markeredgewidth=2.6, clip_on=False)
    ax.plot(np.r_[0, ks], np.r_[0, model], color=BLUE, linewidth=3.0)
    ax.plot(ks, R.FIG2A_MODEL, "s", markersize=16, markerfacecolor="none",
            markeredgecolor=BLUE, markeredgewidth=2.6, clip_on=False)
    # On top, smaller than the open symbols: where they agree the dot sits inside.
    ax.plot(ks, R.FIG2A_DATA, "o", color=GREY, markersize=10, zorder=5, clip_on=False)
    ax.text(3.95, 0.015, "lines: my code\nmarkers: read off their figure", color=SOFT, fontsize=px(24),
            ha="right", va="bottom", linespacing=1.4)
    for y, label, color in ((exact[-1] - 0.01, "exact", INK),
                            (trot[-1] - 0.005, "in time steps", ORANGE),
                            (model[-1] - 0.025, "+ their noise", BLUE),
                            (R.FIG2A_DATA[-1] + 0.03, "their experiment", SOFT)):
        ax.text(4.3, y, label, color=color, fontsize=px(28), va="center", clip_on=False)
    ax.set_xlim(0, 4.05)
    ax.set_ylim(-0.02, 0.56)
    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_yticks([0, 0.25, 0.5])
    ax.set_xlabel("time (their steps)")
    ax.set_ylabel("particles per site")
    save(fig, "replica")


# ----------------------------------------------------------------------------
# 0. Title slide: the talk's question, drawn with the talk's own numbers
# ----------------------------------------------------------------------------


def fig_title():
    """
    The exact signal (white line) and what a noisy machine measures (dots):
    1,000 shots a point, 0.3 % error per two-qubit gate, steps of 0.1. Early
    on the dots sit on the line; by the second peak they have drifted toward
    0.5, the fully scrambled value. Drawn for the dark title slide, so the
    colours are light. Seeded, so every run draws the same dots.
    """
    fine = np.linspace(0.0, 4.0, 401)
    exact = [particle_density(s, N) for s in evolve(bare_vacuum(N), hamiltonian(N, M, W, J), fine)]
    times = np.round(np.arange(0.1, 4.001, 0.1), 10)
    nm = gate_noise_model(3e-3)
    rng = np.random.default_rng(2026)
    measured = []
    for t in times:
        probs = probabilities(prepare_for_probabilities(trotter_circuit(N, M, W, J, t, steps_for(t, 0.1), 2)), nm)
        counts = sample_counts(probs, SHOTS, rng)
        measured.append(observables_from_distribution(counts / counts.sum(), N)["particle_density"])
    fig, ax = slide_fig(640, 190)
    fig.subplots_adjust(left=0.0, right=1.0, bottom=0.04, top=0.96)
    ax.axhline(0.5, color="#475569", linewidth=1.6, linestyle=(0, (4, 4)))
    ax.plot(fine, exact, color="#FFFFFF", linewidth=3.2)
    ax.plot(times, measured, "o", color="#FB923C", markersize=7, markeredgewidth=0)
    ax.set_xlim(-0.05, 4.05)
    ax.set_ylim(-0.03, 0.62)
    ax.axis("off")
    save(fig, "title_signal")
    print(f"  title: measured at t = 4 is {measured[-1]:.3f} against exact {exact[-1]:.3f}")


FIGURES = {
    "title": fig_title,
    "pairs": fig_pairs, "revival": fig_revival, "steps": fig_steps, "shots": fig_shots, "gate_grid": fig_gate_grid,
    "noise": fig_noise, "crossover": fig_crossover, "budget": fig_budget,
    "hamiltonian": fig_hamiltonian, "chains": fig_chains, "replica": fig_replica,
}


def main(only=()) -> None:
    """All figures, or only the named ones (collapse, charge and blind need the sweeps, so they run only in a full pass)."""
    OUT.mkdir(parents=True, exist_ok=True)
    style()
    print("=" * 78)
    print("TALK FIGURES (each quoted number beside the thesis value)")
    print("=" * 78)
    if only:
        for name in only:
            FIGURES[name]()
        return
    for name in ("pairs", "revival", "steps", "shots", "gate_grid", "noise", "crossover"):
        FIGURES[name]()
    eps, sweeps = threshold_sweeps()
    fig_collapse(eps, sweeps)
    fig_budget()
    fig_charge(sweeps)
    fig_blind(eps, sweeps)
    fig_hamiltonian()
    fig_chains()
    fig_replica()


if __name__ == "__main__":
    main(sys.argv[1:])
