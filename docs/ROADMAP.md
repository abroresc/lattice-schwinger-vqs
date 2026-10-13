# Roadmap: dates, deliverables, status

One file, merging two calendars that otherwise live in separate PDFs: the
PHY 420 course deadlines (`course/ProjectDetails2026.pdf`), the week-by-week
research schedule from the proposal (`course/Schwinger_Capstone_Proposal.pdf`,
Section 8). Regenerate this by hand
when a date changes; it is a snapshot, not a live document.

Last updated: 2026-10-13.

## Status right now

- **Gate 1** (Jordan-Wigner verified against the Ising model) and **Gate 2**
  (N=4 Schwinger Hamiltonian built, diagonalized, limits reproduced): done,
  three weeks ahead of the Sep 29 target.
- **M2/M3/M5** (Hamiltonian, verification, Qiskit circuit layer) done. 394
  tests passing against exact diagonalization.
- **SQD vs VQE decision: closed 2026-10-01.** Rosanowski et al.
  (arXiv:2510.26951) read. Decision: use SKQD (Sample-Based Krylov Quantum
  Diagonalization) for ground-state preparation, not VQE. Reasoning: SKQD
  reuses the Trotter circuit already in `circuits.py`, needs no classical
  optimizer loop, is noise-resilient by construction, and the paper's own N=4
  demonstration is close to a direct template for this project's
  ground-state step.
- **2026-10-12/13: M6/M7, M10, M11 and a first pass at M12 done in simulation.**
  `noise.py` (depolarizing/dephasing gate error on the compiled circuit, exact
  density-matrix probabilities, classical readout error), `diagnostics.py`
  (Var(Q), leakage, CP asymmetry, global-depolarizing fidelity), `skqd.py`
  (SKQD with an exact finite-shot expectation), `budget.py` (four-way
  decomposition, failure threshold, diagnostic prediction). 473 tests.
  Scripts 06-08 and figures 06-08. Headline numbers: failure threshold
  eps* x n_CX ~ 2 at 1000 shots; charge diagnostic predicts eps* to within
  10-17 % (always optimistic) but is blind to logical-frame dephasing; optimal
  dt grows from 0.35 to 0.59 as eps_2 goes from 3e-4 to 3e-3. See
  `docs/error-budget.md`.
- **Writing:** Introduction (2nd draft), Theory, and Materials and Methods
  drafted in `thesis/sections/`; `thesis/intro_methods.pdf` is the partial
  submission. Results drafted from scripts 04-08; Discussion, Conclusion and
  Abstract still TODO.
- **Not done:** Discussion against the literature; hardware run (optional);
  coherent/correlated noise models; amplitude damping. Decide with Andresen
  whether any of these is in scope.

## This week (Oct 1 - Oct 6)

- Introduction, 1st draft, due **Oct 6**. The skeleton already exists in
  `thesis/main.tex` (Introduction section with a "what this thesis
  contributes" stub) — this is filling it in, not starting from zero. Needs
  roughly half of the eventual 30+ references.

## Merged calendar

| Date | Course deliverable | Research milestone | Notes |
|---|---|---|---|
| Sep 8 | Annotated bibliography 1 | Pauli string machinery, 4-site Ising diagonalized | done |
| Sep 15 | Annotated bibliography 2 | Gate 1: Jordan-Wigner vs analytic spectrum | done |
| Sep 22 | Annotated bibliography 3 | Lattice theory reading, no code (GRE week) | done |
| Sep 29 | Career exercises | Gate 2: Schwinger Hamiltonian built & verified | done 3 weeks early |
| **Oct 6** | **Introduction, 1st draft** | Buffer; SQD/VQE decision | decision closed today |
| Oct 13 | Introduction, 2nd draft | SKQD on a known system (replaces "VQE on known system") | |
| Oct 20 | Research check-in, intro talk | SKQD on the N=4 Schwinger Hamiltonian (replaces "VQE on Schwinger model") | this is M6/M7 |
| Oct 27 | Theory and Methods draft | Trotter evolution, scaling exponent matches prediction | M8 |
| Nov 3 | Practice talks | Physics observables measured, figures produced | |
| Nov 10 | Results and Discussion draft | Shot noise propagated onto every observable | |
| Nov 17 | Capstone talks | Hardware noise; M11 gate/readout error models, three sources separated | |
| Nov 21 | **Full draft** | Gauss's law / symmetry diagnostic tested against the failure threshold; M10-M12 synthesis | |
| Dec 1 | Colloquium talk | Analysis closed, repository cleaned | |
| Dec 8 | **Final draft** | — | |
| ~Dec 11 | Colloquium | — | per CLAUDE.md status |

## If a gate is missed

The proposal has a named fallback (Section 9): pivot to testing whether
published quantum hardware calibration data predicts measured circuit
performance, using the same statistical machinery already built for the error
budget. Not needed yet — flagged here so it isn't forgotten under deadline
pressure.
