"""
Two short animations for the front page of the repository.

The thesis figures are built to be read slowly on paper. A visitor to the
GitHub page gives the first screen a few seconds, so these two GIFs each show
ONE thing moving, and loop:

  figures/readme/pairs_from_empty_space.gif
      The four sites after the quench from the bare vacuum, computed exactly.
      Each site's disc fills with the probability that it holds a particle
      (even sites, blue) or an antiparticle (odd sites, orange); each link
      darkens with the field energy on it, <L_n^2>. Beside it the particle
      density nu(t) draws itself, and four bars keep the energy books: the
      mass and field terms rise, the hopping term falls by the same amount,
      and the total never moves. That answers "where do the pairs get their
      energy?" without a word of text.

  figures/readme/signal_dies.gif
      The same density measured through the second-order Trotter circuit at
      dt = 0.1, while the two-qubit gate error climbs from zero to 1.5 %. The
      recombination dip D (first peak minus first trough, the signal the thesis
      tracks) shrinks until it falls below twice its 1000-shot uncertainty, at
      the failure threshold of the thesis's threshold table.

Nothing is drawn by hand. Every number on screen comes from the tested
package: exact evolution from observables.evolve, noisy probabilities from
noise.probabilities on the compiled circuits (the same Aer density-matrix
route as scripts/06), the dip and its noise floor from budget.dip. The script
checks the numbers the animations show against the values the thesis quotes
and stops with an AssertionError if they drift.

Images are opaque white so they read on GitHub's light and dark themes alike.

Run:  .venv/bin/python scripts/14_readme_animations.py      (about a minute)
"""
import pathlib
import sys

import numpy as np
from PIL import Image
from scipy.interpolate import CubicSpline

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from schwinger.budget import dip, first_crossing, point_circuit
from schwinger.circuits import trotter_circuit
from schwinger.diagnostics import observables_from_distribution
from schwinger.noise import gate_noise_model, prepare_for_probabilities, probabilities
from schwinger.observables import charge_profile, evolve, particle_density
from schwinger.paulis import expectation
from schwinger.plotstyle import INK, MUTED, PALETTE, use_report_style
from schwinger.schwinger import bare_vacuum, electric_field_op, hamiltonian

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch

OUT = pathlib.Path(__file__).resolve().parents[1] / "figures" / "readme"
N, M, W, J = 4, 0.5, 1.0, 1.0
T_PEAK, T_TROUGH = 0.89, 1.76          # first maximum and minimum of the exact nu(t)
SHOTS = 1000
BLUE, ORANGE = PALETTE[0], PALETTE[1]
LIGHT = "#e6e6e6"
DPI = 100


# ----------------------------------------------------------------------------
# GIF writing: one shared palette for every frame, so colours never flicker
# ----------------------------------------------------------------------------

def grab(fig) -> Image.Image:
    fig.canvas.draw()
    return Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[..., :3].copy())


def write_gif(frames: list[Image.Image], durations: list[int], path: pathlib.Path) -> None:
    # Build the palette from a strip of evenly spaced frames, then map every
    # frame onto it without dithering (dither noise would shimmer between frames).
    sample = frames[:: max(1, len(frames) // 8)]
    w, h = frames[0].size
    strip = Image.new("RGB", (w, h * len(sample)))
    for i, f in enumerate(sample):
        strip.paste(f, (0, i * h))
    pal = strip.quantize(colors=128, method=Image.Quantize.MEDIANCUT)
    q = [f.quantize(palette=pal, dither=Image.Dither.NONE) for f in frames]
    q[0].save(path, save_all=True, append_images=q[1:], duration=durations, loop=0,
              optimize=False, disposal=1)
    print(f"wrote {path.relative_to(OUT.parents[1])}  ({len(frames)} frames, "
          f"{path.stat().st_size / 1e6:.2f} MB)")


# ----------------------------------------------------------------------------
# 1. Pairs out of empty space
# ----------------------------------------------------------------------------

def pairs_from_empty_space() -> None:
    H = hamiltonian(N, M, W, J)
    parts = {"mass": hamiltonian(N, M, 0.0, 0.0),
             "field": hamiltonian(N, 0.0, 0.0, J),
             "hopping": hamiltonian(N, 0.0, W, 0.0)}
    assert np.abs(sum(parts.values()) - H).max() < 1e-12, "the three parts must add up to H"
    L2 = [electric_field_op(j, N) @ electric_field_op(j, N) for j in range(N - 1)]

    times = np.linspace(0.0, 5.5, 111)
    fine = np.linspace(0.0, 5.5, 551)
    states = evolve(bare_vacuum(N), H, times)
    nu_fine = np.array([particle_density(s, N) for s in evolve(bare_vacuum(N), H, fine)])
    prob = np.array([np.abs(np.real(charge_profile(s, N))) for s in states])
    field = np.array([[np.real(expectation(s, op)) for op in L2] for s in states])
    energy = {k: np.array([np.real(expectation(s, op)) for s in states]) for k, op in parts.items()}
    total = sum(energy.values())

    # The numbers the thesis quotes for this quench (Theory, Observables).
    i_peak = int(np.argmax(nu_fine[fine < 1.3]))
    i_trough = int(np.argmin(np.where((fine > 1.3) & (fine < 2.2), nu_fine, 9)))
    print(f"first peak  nu = {nu_fine[i_peak]:.3f} at t = {fine[i_peak]:.2f}   (thesis: 0.463 at 0.89)")
    print(f"first trough nu = {nu_fine[i_trough]:.3f} at t = {fine[i_trough]:.2f}   (thesis: 0.133 at 1.76)")
    assert abs(nu_fine[i_peak] - 0.463) < 1e-3 and abs(fine[i_peak] - 0.89) < 0.011
    assert abs(nu_fine[i_trough] - 0.133) < 1e-3 and abs(fine[i_trough] - 1.76) < 0.011
    s89 = evolve(bare_vacuum(N), H, np.array([T_PEAK]))[0]
    gains = {k: np.real(expectation(s89, op)) - energy[k][0] for k, op in parts.items()}
    print("energy change at t = 0.89: " + ", ".join(f"{k} {v:+.2f}" for k, v in gains.items())
          + "   (thesis: mass +0.93, field +0.97, hopping -1.90)")
    assert abs(gains["mass"] - 0.93) < 0.006 and abs(gains["field"] - 0.97) < 0.006
    assert abs(gains["hopping"] + 1.90) < 0.006
    print(f"total energy drift over the animation: {np.abs(total - total[0]).max():.1e}")
    assert np.abs(total - total[0]).max() < 1e-12

    use_report_style()
    plt.rcParams.update({"font.size": 11, "axes.labelsize": 11})
    fig = plt.figure(figsize=(8.4, 3.9), dpi=DPI, facecolor="white")
    ax_l = fig.add_axes([0.02, 0.36, 0.44, 0.52])           # the lattice
    ax_e = fig.add_axes([0.10, 0.06, 0.34, 0.22])           # energy books
    ax_c = fig.add_axes([0.56, 0.17, 0.42, 0.68])           # nu(t)
    xs = np.arange(N) * 1.0
    frames, durations = [], []

    for k, t in enumerate(times):
        for ax in (ax_l, ax_e, ax_c):
            ax.cla()
        # ---- lattice
        ax_l.set_xlim(-0.6, 3.6)
        ax_l.set_ylim(-0.75, 0.95)
        ax_l.set_aspect("equal")
        ax_l.axis("off")
        r = 0.32
        for n in range(N - 1):
            g = float(np.clip(field[k, n], 0, 1))
            shade = tuple(1 - 0.62 * g for _ in range(3))
            ax_l.add_patch(FancyBboxPatch((xs[n] + r, -0.07), xs[n + 1] - xs[n] - 2 * r, 0.14,
                                          boxstyle="round,pad=0,rounding_size=0.05",
                                          fc=shade, ec="#cfcfcf", lw=0.6))
        for n, x in enumerate(xs):
            color = ORANGE if n % 2 == 0 else BLUE     # odd physics site = code site 0, 2
            ax_l.add_patch(Circle((x, 0), r, fc="white", ec="#9a9a9a", lw=1.2))
            p = float(prob[k, n])
            if p > 1e-3:
                ax_l.add_patch(Circle((x, 0), r * np.sqrt(p), fc=color, ec="none"))
            ax_l.text(x, 0.47, f"site {n + 1}", ha="center", fontsize=9, color=MUTED)
            ax_l.text(x, -0.62, f"{p:.2f}", ha="center", fontsize=9, color=INK)
        ax_l.text(-0.55, 0.80, "probability of a particle (blue) or antiparticle (orange)",
                  fontsize=9, color=INK)
        ax_l.text(1.5, -0.88, "grey links: field energy", ha="center", fontsize=8, color=MUTED)
        # ---- energy books: change since t = 0
        names = ["mass", "field", "hopping", "total"]
        vals = [energy["mass"][k] - energy["mass"][0], energy["field"][k] - energy["field"][0],
                energy["hopping"][k] - energy["hopping"][0], total[k] - total[0]]
        cols = [MUTED, MUTED, MUTED, INK]
        ax_e.barh(range(4), vals, color=cols, height=0.62)
        for j, v in enumerate(vals):
            label = f"{v:+.2f}".replace("\u2212", "-")
            label = "0.00" if abs(v) < 0.005 else label.replace("-", "\u2212")
            ax_e.text(v + (0.08 if v >= 0 else -0.08), j, label, va="center",
                      ha="left" if v >= 0 else "right", fontsize=8, color=INK)
        ax_e.axvline(0, color=INK, lw=0.8)
        ax_e.set_xlim(-2.6, 2.6)
        ax_e.set_ylim(3.6, -0.6)
        ax_e.set_yticks(range(4), names, fontsize=9)
        ax_e.tick_params(axis="x", labelsize=8)
        ax_e.grid(False)
        ax_e.spines[["left", "top", "right"]].set_visible(False)
        ax_e.set_title("energy gained since the start", fontsize=9, color=INK, loc="left", pad=2)
        # ---- nu(t)
        ax_c.plot(fine, nu_fine, color=LIGHT, lw=2)
        upto = fine <= t + 1e-9
        ax_c.plot(fine[upto], nu_fine[upto], color=INK, lw=2.2)
        nu_t = particle_density(states[k], N)
        ax_c.plot([t], [nu_t], "o", color=INK, ms=7, mec="white", mew=1.5)
        ax_c.axhline(0.5, color=MUTED, lw=0.8, ls=(0, (4, 3)))
        ax_c.text(5.45, 0.515, "one pair on average", ha="right", va="bottom", fontsize=8, color=MUTED)
        ax_c.set_xlim(0, 5.5)
        ax_c.set_ylim(0, 0.62)
        ax_c.set_xlabel("time (units of 1/w)")
        ax_c.set_title("particles per site", fontsize=11, color=INK, loc="left")
        ax_c.text(0.15, 0.585, f"t = {t:.2f}", fontsize=10, color=INK, va="top")
        fig.suptitle("Pairs appear out of empty space: four sites, solved exactly",
                     x=0.02, ha="left", y=0.985, fontsize=12, color=INK, fontweight="bold")
        frames.append(grab(fig))
        durations.append(1200 if k == 0 else (2200 if k == len(times) - 1 else 70))
    plt.close(fig)
    write_gif(frames, durations, OUT / "pairs_from_empty_space.gif")


# ----------------------------------------------------------------------------
# 2. The signal dies as the gates get worse
# ----------------------------------------------------------------------------

def signal_dies() -> None:
    dt = 0.1
    ks = np.arange(0, 26)                      # t = 0.0, 0.1, ..., 2.5 in steps of exactly dt
    times = ks * dt
    preps = [None] + [prepare_for_probabilities(trotter_circuit(N, M, W, J, k * dt, int(k), 2))
                      for k in ks[1:]]
    # The two circuits the thesis's threshold uses at dt = 0.1: 9 and 18 steps.
    peak = point_circuit(N, M, W, J, T_PEAK, 9)
    trough = point_circuit(N, M, W, J, T_TROUGH, 18)
    print(f"trough circuit: {trough.two_qubit_gates} CX   (thesis: 576)")
    assert trough.two_qubit_gates == 576

    fine = np.linspace(0, 2.5, 251)
    exact = np.array([particle_density(s, N) for s in evolve(bare_vacuum(N), hamiltonian(N, M, W, J), fine)])
    eps_grid = np.concatenate([[0.0], np.logspace(-4, np.log10(1.5e-2), 34)])

    rows = []
    for eps in eps_grid:
        nm = gate_noise_model(eps)
        nu = [0.0] + [observables_from_distribution(probabilities(p, nm), N)["particle_density"]
                      for p in preps[1:]]
        pp, pt = probabilities(peak.prepared, nm), probabilities(trough.prepared, nm)
        D, sD = dip(pp, pt, N, SHOTS)
        nu_p = observables_from_distribution(pp, N)["particle_density"]
        nu_t = observables_from_distribution(pt, N)["particle_density"]
        rows.append((eps, np.array(nu), D, sD, nu_p, nu_t))

    D0 = rows[0][2]
    e = np.array([r[0] for r in rows[1:]])
    margin = np.array([r[2] - 2 * r[3] for r in rows[1:]])
    eps_star = float(first_crossing(e, margin))      # interpolates in log(eps) itself
    nm3 = gate_noise_model(1e-3)
    kept = dip(probabilities(peak.prepared, nm3), probabilities(trough.prepared, nm3), N, SHOTS)[0] / D0
    print(f"noiseless dip D = {D0:.3f}")
    print(f"threshold from this grid: eps* = {eps_star:.2e}   (thesis: 3.3e-3 at 1000 shots)")
    print(f"dip kept at eps = 1e-3: {kept:.2f}   (thesis: 0.51)")
    assert abs(eps_star / 3.3e-3 - 1) < 0.03
    assert abs(kept - 0.51) < 0.006

    use_report_style()
    plt.rcParams.update({"font.size": 11, "axes.labelsize": 11})
    fig = plt.figure(figsize=(8.4, 4.2), dpi=DPI, facecolor="white")
    ax = fig.add_axes([0.08, 0.34, 0.56, 0.49])
    ax_s = fig.add_axes([0.08, 0.10, 0.56, 0.06])         # the error-rate slider
    frames, durations = [], []
    for i, (eps, nu, D, sD, nu_p, nu_t) in enumerate(rows):
        ax.cla()
        ax_s.cla()
        lost = D <= 2 * sD
        mark = ORANGE if lost else BLUE
        ax.plot(fine, exact, color=LIGHT, lw=2.5)
        ax.plot(fine, CubicSpline(times, nu)(fine), color=INK, lw=2.2)
        ax.axhline(0.5, color=MUTED, lw=0.8, ls=(0, (4, 3)))
        ax.text(0.04, 0.515, "fully scrambled", va="bottom", fontsize=8, color=MUTED)
        # the dip bracket
        xb = 1.98
        ax.plot([T_PEAK, xb], [nu_p, nu_p], color=MUTED, lw=0.7, ls=":")
        ax.plot([T_TROUGH, xb], [nu_t, nu_t], color=MUTED, lw=0.7, ls=":")
        if nu_p - nu_t > 0.03:
            ax.annotate("", xy=(xb, nu_t), xytext=(xb, nu_p),
                        arrowprops=dict(arrowstyle="<->", color=mark, lw=1.8))
        else:
            ax.plot([xb, xb], [nu_t, nu_p], color=mark, lw=2.5)
        ax.text(xb + 0.05, (nu_p + nu_t) / 2, "dip", va="center", fontsize=10, color=mark)
        for tt, vv in ((T_PEAK, nu_p), (T_TROUGH, nu_t)):
            ax.plot([tt], [vv], "o", color=mark, ms=7, mec="white", mew=1.5)
        ax.set_xlim(0, 2.5)
        ax.set_ylim(0, 0.62)
        ax.set_xlabel("time (units of 1/w)")
        ax.set_title("particles per site, through the circuit (dt = 0.1)", fontsize=11,
                     color=INK, loc="left")
        # slider
        lo, hi = 1e-4, 1.5e-2
        ax_s.set_xscale("log")
        ax_s.set_xlim(lo, hi)
        ax_s.set_ylim(0, 1)
        ax_s.set_yticks([])
        ax_s.grid(False)
        ax_s.spines[["left", "top", "right"]].set_visible(False)
        ax_s.axvspan(eps_star, hi, color="#fbe3d2", lw=0)
        ax_s.axvline(eps_star, color=ORANGE, lw=1.2)
        ax_s.text(eps_star * 1.08, 1.25, "signal lost", color=ORANGE, fontsize=9, va="bottom")
        ax_s.set_xticks([1e-4, 1e-3, 1e-2], ["0.01 %", "0.1 %", "1 %"], fontsize=9)
        ax_s.set_xlabel("error per two-qubit gate", fontsize=10)
        if eps > 0:
            ax_s.plot([eps], [0.5], "o", color=INK, ms=9, mec="white", mew=1.5, clip_on=False)
        # readout panel on the right
        fig.texts.clear()
        fig.text(0.02, 0.965, "The signal dies as the gates get worse", fontsize=12,
                 color=INK, fontweight="bold", va="top")
        x0 = 0.69
        fig.text(x0, 0.84, "error per two-qubit gate", fontsize=9, color=MUTED)
        fig.text(x0, 0.77, f"{100 * eps:.3f} %" if eps < 1e-3 else f"{100 * eps:.2f} %",
                 fontsize=16, color=INK)
        fig.text(x0, 0.66, "error rate \u00d7 CX gates in the\ndeepest circuit (576)", fontsize=9,
                 color=MUTED, va="bottom", linespacing=1.2)
        fig.text(x0, 0.59, f"{eps * 576:.2f}", fontsize=13, color=INK)
        fig.text(x0, 0.48, "dip  vs  noise floor (2σ, 1000 shots)", fontsize=9, color=MUTED)
        fig.text(x0, 0.41, f"{D:.3f}  vs  {2 * sD:.3f}", fontsize=13, color=INK)
        fig.text(x0, 0.28, "dip lost in the noise" if lost else "dip still visible",
                 fontsize=12, color=mark, fontweight="bold")
        fig.text(x0, 0.12, f"threshold: {100 * eps_star:.2f} % per gate,\n"
                 f"error rate × 576 gates ≈ {eps_star * 576:.1f}", fontsize=9, color=MUTED)
        frames.append(grab(fig))
        durations.append(1500 if i == 0 else (2500 if i == len(rows) - 1 else 260))
    plt.close(fig)
    write_gif(frames, durations, OUT / "signal_dies.gif")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pairs_from_empty_space()
    signal_dies()
    print("every assert passed: the animations show the thesis numbers")


if __name__ == "__main__":
    main()
