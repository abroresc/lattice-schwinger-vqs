# schwinger-n4

Quantum simulation of the lattice Schwinger model on four staggered sites:
exact diagonalisation, verified against every analytic limit available.

PHY 420 capstone, Gettysburg College, Fall 2026. Abror Niyazmetov, advised by
Dr. Kurt Andresen.

---

## What this is

Quantum electrodynamics in one space and one time dimension, put on a
four-site lattice, with the gauge field eliminated by solving Gauss's law, and
mapped to four qubits by Jordan-Wigner. The result is a 16 by 16 Hermitian
matrix that confines charges and produces particle-antiparticle pairs out of
the vacuum.

This repository covers milestones M2 to M12 of the semester plan in
simulation: the verified Hamiltonian (Gate 2), Trotterised circuits, SKQD
ground-state preparation, gate and readout error models, the symmetry
diagnostics, and the four-way error budget. The thesis now has a full first
draft, including the discussion against the literature; what remains is
revision.

The design principle throughout: **every stage has a known answer attached**,
so a wrong result is caught in minutes rather than weeks. Nothing is trusted
because it looks right.

## Start here

New to the toolchain? [`SETUP.md`](SETUP.md) is the step-by-step version,
including how to push this to GitHub.

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python -m pytest tests/ -q                    # 496 tests, all against known answers
python scripts/03_schwinger_spectrum.py       # the Gate 2 artifact and its residuals
python scripts/05_qiskit_trotter.py           # circuits, Trotter scaling, shot noise
```

**[`notebooks/schwinger_walkthrough.ipynb`](notebooks/schwinger_walkthrough.ipynb)
is the guided tour**, with all output already committed, so it reads on GitHub
without running anything.

Then read [`docs/derivation.md`](docs/derivation.md), which derives the whole
thing from the continuum Lagrangian for a reader with one quantum mechanics
course. [`CONVENTIONS.md`](CONVENTIONS.md) is the sign and ordering reference;
check it first when a spectrum comes out wrong.

## The Hamiltonian

```
    H  =  (m/2) sum_{n=1..N}   (-1)^n Z_n
       +   w    sum_{n=1..N-1} ( sigma^+_n sigma^-_{n+1} + h.c. )
       +   J    sum_{n=1..N-1} L_n^2,     L_n = (1/2) sum_{m<=n} ( Z_m + (-1)^m )
```

The spin Hamiltonian of Martinez et al., Nature **534**, 516 (2016), Methods Eq. (3)
in arXiv:1605.04570v1 (Eq. (2) there is Gauss's law), who run exactly
this system at N = 4 on a trapped-ion device. Fourteen Pauli strings on four
qubits; `scripts/03_schwinger_spectrum.py` prints the table.

## Verification

Built three independent ways and checked against six analytic results.
Residuals, not pass/fail flags, and measured rather than the tests'
tolerances. `scripts/09_verification_table.py` prints every row of the
thesis verification table at the thesis parameters (m = 0.5w, J = w); the
first block below is `scripts/03`, which takes the worst case over three
parameter sets.

| Check | Residual |
|---|---|
| Pauli expansion vs direct matrix construction | 9e-16 |
| Fermion construction through Jordan-Wigner vs spin form | 9e-16 |
| Gauss's law is an identity on the full Hilbert space | 0 |
| `[H, Q] = 0` | 0 |
| Strong coupling `E0 = -mN/2` at `w = 0` | 0 |
| Free fermion limit, all 16 levels | 9e-16 |
| Hellmann-Feynman on the chiral condensate (checks the observable against the mass term only) | 3e-11 |
| CP symmetry inside the neutral sector (`scripts/09`) | 0 |
| Open TFIM: full 2^N spectrum vs closed-form free fermions, N = 2, 4, 6, 8 | 3e-14 |
| Periodic TFIM with parity sectors, N = 6 (`scripts/02`) | 2e-14 |
| Qiskit `SparsePauliOp` independent rebuild | 0 |
| Circuit statevector vs exact evolution, t = 2, 4000 second-order steps | 7e-8 |
| Trotter exponent, first order (predicted 1) | measured 1.036 |
| Trotter exponent, second order (predicted 2) | measured 2.003 |

## Layout

```
src/schwinger/
  paulis.py          Pauli-string machinery, conventions, Qiskit adapter
  jordan_wigner.py   JW operators, algebra checks, Bogoliubov solver
  ising.py           TFIM and its closed-form spectrum (the warm-up)
  schwinger.py       the Schwinger Hamiltonian, Gauss's law, CP, states
  observables.py     particle density, condensate, field profile, evolution
  analytics.py       the analytic limits everything is checked against
  circuits.py        Qiskit layer: state prep, Trotter steps, shots, gate counts
  noise.py           gate and readout error, exact noisy probabilities (M11)
  diagnostics.py     observables and symmetry diagnostics from Z-basis data
  skqd.py            SKQD ground state with exact finite-shot expectation (M6/M7)
  budget.py          four-way error budget, failure threshold, prediction (M10/M12)
  plotstyle.py       shared figure style
scripts/             runnable, each regenerates one figure and prints residuals
notebooks/           the guided walkthrough, output committed
tests/               496 tests, every one against an independently known answer
docs/derivation.md   the full derivation, from the continuum Lagrangian
docs/error-budget.md noise, diagnostics, error budget and SKQD, explained
thesis/              LaTeX skeleton for the final write-up, compiles today
figures/             regenerated by the scripts, PNG and PDF, no manual steps
```

Figures are written as both PNG and PDF. The PDFs are vector, so they go
straight into the LaTeX thesis without regenerating anything in December.

## Figures

| File | What it shows |
|---|---|
| `01_tfim.png` | TFIM crossover, and JW verified against the closed form |
| `02_jordan_wigner.png` | The fermion algebra, and the periodic-boundary trap |
| `03_schwinger_spectrum.png` | All 16 levels vs mass; condensate two ways |
| `04_pair_production.png` | Pairs created out of the vacuum; flux on each link |
| `05_qiskit_trotter.png` | Trotter scaling with fitted slopes; measured density with shot noise; the systematic-vs-statistical crossover |
| `06_noise_threshold.png` | Pair production under gate error; the failure threshold; whether Var(Q) sees the damage |
| `07_skqd.png` | SKQD energy bias against shots, computed exactly |
| `08_error_budget.png` | Total error against step size with the optimum; the four-way decomposition |

## Things found by checking rather than reading

Recorded because they were not in the reading list and the next person will
hit them too.

- **Qiskit endianness.** Pauli label strings coincide with ours; the qubit
  *numbering* is reversed. Reversing the label, which is the natural guess, is
  wrong. Found by a numerical cross-check. See `CONVENTIONS.md` section 5.
- **A conditioning bug in the standard Bogoliubov recipe.** Forming
  `(A-B)(A+B)` squares the condition number. The quasiparticle energies are
  the singular values of `A+B`, so take an SVD and never form the product.
  Worst residual at N = 8, `h/J = 0.05`: 7e-10 improved to 1e-14.
- **The periodic-boundary Jordan-Wigner trap.** At `k = 0` and `k = pi` the
  pairing term vanishes and the energy is the signed `2(h - J cos k)`, not the
  positive square root. Wrong by exactly `|2(h-J)|`, and only when `h < J`.
- **There is no Gauss's law violation to measure in this encoding.** The
  constraint is substituted, not imposed, so it holds identically on every
  state. The diagnostic has to be total charge variance instead. This changes
  a piece of the proposal; see `docs/derivation.md` section 8.2.
- **A CP symmetry that holds only in the neutral sector**. On the
  CP-symmetric states this project evolves (everything that starts from the
  bare vacuum) it forces `<L_n> = <L_{N-n}>`, which gives a second free
  diagnostic. Neutrality alone is not enough: `|0110>` is a counterexample.
  See section 7.

## First numbers toward the thesis

At N = 4, `m/w = 0.5`, `w = J = 1`, evolved to `t = 3` with second-order Trotter:

| shots | exact shot-noise sigma on particle density | Trotter bias falls below it at |
|---|---|---|
| 256 | 0.0190 | `dt` = 0.254 (between 11 and 12 steps) |
| 1000 | 0.0096 | `dt` = 0.182 (between 16 and 17 steps) |
| 4096 | 0.0047 | `dt` = 0.128 (between 23 and 24 steps) |

Sigma is computed exactly from the outcome distribution and the crossing is
interpolated in log `dt` (`scripts/05`). An earlier version of this table took
sigma from one sampled run and read the crossing off a factor-of-two grid,
which gave 0.19 and 0.094.

Below those step sizes, shrinking `dt` buys nothing the statistics can see, and
costs circuit depth: 320 two-qubit gates at `dt = 0.3` against 6400 at
`dt = 0.015`. Adding gate error turns that two-way trade into the three-way one
the thesis is actually about.

## The write-up

`thesis/` holds the thesis in LaTeX. The Introduction, Theory and Materials
and Methods are drafted (`thesis/intro_methods.pdf` builds them alone), and the
Results, Discussion, Conclusion and Abstract have a first draft
(`thesis/main.pdf` is the whole thesis). Figures are included
from `../figures/*.pdf`, so regenerating a figure updates the thesis on the next
compile with no manual step.

```bash
cd thesis && latexmk -pdf main.tex
```

It also builds on Overleaf: upload `thesis/` and `figures/` and press Recompile.

## Statement on AI assistance

AI tooling was used for implementation scaffolding, debugging (particularly
sign conventions in Jordan-Wigner strings), literature navigation and drafting.
The physics derivations, experimental design and analysis are mine, and every
claim an AI tool made about the literature was verified against the source
before it entered this repository.

Two errors in AI-drafted material were caught by the verification design and
are documented in place rather than quietly corrected: the Qiskit endianness
adapter and the Bogoliubov conditioning issue, both listed above. That is the
verification structure working as intended, and it is the reason the structure
exists.
