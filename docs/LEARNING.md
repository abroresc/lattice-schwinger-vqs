# Learning guide: topic by topic, with references back into this project

Why this file exists: you can get any of these topics explained well by
Gemini or another general tool, but a general explanation doesn't know which
section of `derivation.md` it feeds into, which paper in `papers/` backs it
with a citable result, or which test in `tests/` checks it against a known
answer. This file is the connective tissue — a short "why this matters here"
per topic, plus exact pointers. Read the explanation anywhere you like; come
back here to see where it plugs in.

Not a substitute for `derivation.md` — that document *is* the derivation, in
the order the thesis needs it, written so you can defend it on a whiteboard.
This file is the index that sits one level above it: what to learn before
each section, and what in the literature backs it up.

## 1. Second-quantized fermion operators

**Why it's here:** the whole project starts by representing "a particle is
present at this site" as an operator algebra, not a wavefunction. Everything
downstream — staggering, Jordan-Wigner, the Hamiltonian itself — is built on
creation/annihilation operators and their anticommutation relations.

- `derivation.md` section 2 (the fermion-doubling problem) assumes this.
- Proposal section 6, item 1: the two-site Hamiltonian by hand is the
  concrete check for this topic specifically.

## 2. Lattice discretization and fermion doubling

**Why it's here:** naively putting a relativistic fermion on a grid conjures
extra particle species that aren't physical. Staggering (next topic) is the
fix, but you have to see the problem first to believe the fix is necessary.

- `derivation.md` section 2.1.
- Proposal section 4.2.

## 3. Staggered fermions (Kogut-Susskind)

**Why it's here:** this is the specific trick that fixes doubling for this
project — spreading the components of one fermion across even/odd sites. It
is also where the "vacuum = odd sites filled" picture comes from, which is
the reference state every later circuit starts from.

- `derivation.md` sections 2.2-2.3.
- `papers/hamiltonian-formulations-lqcd.pdf` for general background.
- Primary citation: Kogut and Susskind, PRD 11, 395 (1975) — not on arXiv,
  get it through the library (see `papers/README.md`).

## 4. Gauss's law and eliminating the gauge field

**Why it's here:** this is the single structural decision that makes N=4 cost
4 qubits instead of 4 qubits plus 3 link registers. Understanding *why* the
field can be solved away in 1+1D (not just that it can) is load-bearing —
it's also why there's no "Gauss's law violation" to measure in this encoding
(the constraint holds identically), which is a correction CLAUDE.md flags
explicitly.

- `derivation.md` section 3.
- `papers/martinez-2016-realtime-dynamics.pdf` — the primary encoding
  reference, Methods Eq. (3) in arXiv:1605.04570v1.

## 5. Jordan-Wigner transformation

**Why it's here:** the known failure mode for this entire class of project
(per the proposal's own risk section) is a sign error here. The project
deliberately learns this on the transverse-field Ising model first, where a
closed-form answer exists, before touching the unsolved Schwinger case.

- `derivation.md` section 4 and 4.1 (the Ising warm-up, which is not
  optional).
- CLAUDE.md's "Jordan-Wigner sign" convention note — the minus sign in
  `c_j^dag c_{j+1} = -sigma^+_j sigma^-_{j+1}` is real, not a bug. Run
  `hopping_sign_demo()` before touching any coefficient to make a sign work
  out.

## 6. Building the N=4 Hamiltonian as a Pauli string

**Why it's here:** this is where sections 1-5 combine into the actual object
the rest of the project operates on — a 16x16 matrix, sum of Pauli strings.

- `derivation.md` section 5 (5.1 expanding the electric term, 5.2 the
  result).
- `derivation.md` section 6 for what each check on this Hamiltonian actually
  proves (strong-coupling limit, mass limits).

## 7. Trotterization / product formulas

**Why it's here:** this is how `e^{-iHt}` gets turned into a quantum circuit.
The scaling exponent (slope 1 for first-order, slope 2 for second-order on a
log-log error plot) is a falsifiable prediction, not a sanity check — a wrong
slope means the circuit is built wrong even if the numbers look plausible.

- Proposal section 5.2 and 7.1.
- `papers/childs-2021-theory-of-trotter-error.pdf` — the theory behind the
  bound. Read before Shaw et al.
- `papers/shaw-2020-quantum-algorithms-lattice-schwinger.pdf` — Trotter and
  shot-noise bounds on the same system, with published numbers to compare
  against.

## 8. The variational principle and VQE

**Why it's here:** VQE was the originally proposed method for ground-state
prep (proposal section 5.1, schedule weeks of Oct 13/20). As of 2026-10-01
this project uses SKQD instead (see topic 9) — but the variational principle
itself is still worth understanding, partly because its free correctness
check (`<H> >= E0` always, so underestimating the ground energy is a bug, not
a discovery) is a good model for what a *cheap* check looks like, and partly
because the Methods section should explain why VQE was considered and set
aside.

- Proposal section 5.1.
- `papers/kokail-2019-self-verifying-variational.pdf` — variational
  simulation on trapped ions; closest published relative of the Gauss's-law
  diagnostic idea.

## 9. Krylov subspace methods and SKQD

**Why it's here:** this is the method actually adopted for ground-state prep,
decided 2026-10-01 after reading Rosanowski et al. The core idea — sample
bitstrings from a time-evolved reference state, build a subspace from them,
diagonalize the Hamiltonian classically inside that subspace — reuses the
same Trotter circuit infrastructure from topic 7, rather than requiring a
separate classical optimization loop like VQE does.

- `papers/rosanowski-2025-sample-based-krylov.pdf` — read in full
  2026-10-01. Section III A is the SKQD protocol itself; Section IV A is an
  N=4 demonstration on a trapped-ion processor that is close to a direct
  template for this project (same zero-charge-sector postselection, same
  kinetic-only Trotter circuit, validated against exact diagonalization).
- `papers/miura-2026-active-sampling-sample-based.pdf` — shot count as a
  *structural bias* source in sample-based diagonalization, distinct from
  ordinary shot noise. Relevant once SKQD's own shot count becomes part of
  the error budget (M10).

## 10. The error budget: Trotter error, shot noise, gate/readout error

**Why it's here:** this is the actual thesis, not a supporting topic. The
three sources scale differently (Trotter error is systematic and
step-size-dependent; shot noise falls as 1/sqrt(shots); gate/readout error is
a hardware noise-model parameter) and have to be separated, not lumped,
because the headline results — the crossover point and the failure threshold
— are statements about where one source overtakes another.

- Proposal section 7, all subsections. This is the most important single
  section of the proposal to re-read before M10-M12.
- `papers/shaw-2020-quantum-algorithms-lattice-schwinger.pdf` and
  `papers/childs-2021-theory-of-trotter-error.pdf` for Trotter error.
  `papers/miura-2026-active-sampling-sample-based.pdf` for shot-count bias.
- Not yet identified: standard references on quantum noise channels for M11
  (per `papers/README.md`) — pull these from citations inside the papers
  above once M11 starts.

## 10b. Where this is implemented

- `docs/error-budget.md` is the explanation of everything in this topic, with a
  twelve-question self-test, written 2026-10-13. Read it before the
  Results/Discussion draft.
- Code: `noise.py`, `diagnostics.py`, `budget.py`, `skqd.py`; figures 06-08.

## 11. Symmetries as diagnostics (charge conservation, CP)

**Why it's here:** in this encoding there is no Gauss's-law violation to
measure (the constraint is substituted, not imposed — it holds identically).
The actual noise diagnostics are total charge variance (`[H,Q] = 0` exactly
in the clean theory) and a CP symmetry that only commutes with H inside the
neutral sector. Both are exact on paper and both break under real noise,
which is what makes them useful probes.

- `derivation.md` section 7 ("a symmetry that was not in the plan") and
  section 8.2 (why the Gauss's-law diagnostic doesn't exist in this
  encoding — a correction to the proposal's own text).
- `papers/nguyen-2022-symmetry-protection-trapped-ions.pdf` — the model for
  this diagnostic.

## How to use this with an external tool (Gemini, etc.)

Fine to get the conceptual explanation anywhere. Before trusting it into the
thesis, check it against the specific section and paper listed above — if a
general explanation of, say, Jordan-Wigner doesn't mention the string-of-Z
nonlocality or match the sign convention in CLAUDE.md, it's explaining the
textbook version, not this project's version. The self-test in
`derivation.md` section 9 is the twelve-question check for whether you can
actually trace the whole chain, which is the standard that matters by
December.
