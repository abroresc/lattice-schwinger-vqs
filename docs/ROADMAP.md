# Timeline: how the project ran

A record of what was done when, and of the decisions that changed the plan on
the way. Dates come from the commit history and the weekly meeting notes. The
semester plan was organised around *gates*: checkpoints with a known answer
attached, where the project would stop and fix things rather than build on
something unverified.

## What happened when

| Date (2026) | Milestone |
|---|---|
| Aug 25 to Sep 22 | Reading: three annotated bibliographies, five papers each. |
| Sep 10 | *Gate 1 and Gate 2 met*, three weeks before the Sep 29 target. Jordan-Wigner verified against the exactly solvable Ising chain; the four-site Schwinger Hamiltonian built three independent ways and checked against every analytic limit available; the Qiskit circuit layer and the first crossover numbers. |
| Oct 1 | *SKQD chosen over VQE* for the ground state, after reading Rosanowski et al. (arXiv:2510.26951). The proposal had named VQE. |
| Oct 9 | Introduction, Theory and Methods submitted as a partial draft. |
| Oct 13 | *The second half of the research plan done in simulation:* gate and readout error, the symmetry diagnostics, SKQD and the four-way error budget. First headline numbers: the threshold rule, the 10 to 17 % prediction, the optimal step. |
| Oct 14 | Thesis physics corrected and every quoted number recomputed from the scripts. |
| Oct 16 | Sections 1 to 3 rewritten from 7,900 to about 5,500 words, to fit the whole thesis in the 8,000-word limit and to be readable by a physics undergraduate. |
| Oct 20 | Introductory talk and research check-in. |
| Nov 2 | Discussion, Conclusion and Abstract drafted; every statement about a paper checked against the paper, which corrected four claims. From the advisor's feedback: Fig. 2(a) of Martinez et al. reproduced, and the model solved exactly on up to 12 sites. |
| Nov 17 | Final capstone talk. Thesis brought to full-draft form; this repository's tour, animations and state diagram added; a fresh run of the whole pipeline reproduced the Nov 2 log line for line. |
| Late Nov, early Dec | Full draft, then the final draft. |

## Decisions that changed the plan

*SKQD instead of VQE (Oct 1).* SKQD reuses the Trotter circuits already
built, needs no classical optimiser loop, and only cares which bitstrings are
sampled, not their exact amplitudes, so it tolerates noise by design.
Rosanowski et al.'s own four-qubit demonstration was close to a template for
this project. Reasoning in [`LEARNING.md`](LEARNING.md), topic 9.

*No Gauss's-law diagnostic.* The proposal planned to measure Gauss's-law
violation as the warning light. In this encoding Gauss's law is solved and
substituted, so it holds on every state and cannot be violated. The charge
variance and a CP symmetry replaced it. See [`derivation.md`](derivation.md)
section 8.2 and [tour chapter 6](tour/06-warning-light.md).

*Klco et al. is not the reference for this encoding.* They keep the
electric field as qubits with a cutoff; this project follows Martinez et al.,
who eliminate it. [`derivation.md`](derivation.md) section 8.1.

*Phase noise on the gates turned out to be visible (Oct 13).* The first
version of the noise code assumed dephasing on the compiled gates would be
invisible to the charge check. The simulation showed otherwise; only dephasing
between time steps is invisible. That blind spot became one of the results.

*Scope held fixed.* Four sites, open boundaries. Larger lattices, periodic
boundaries and two space dimensions were named in the proposal as temptations
to resist until the error budget worked. The only extension, the exact
twelve-site check, was added to answer whether four sites are enough.

## Not done

Coherent and correlated noise, amplitude damping, and a run on real hardware.
[Tour chapter 9](tour/09-where-this-stops.md) says what each would change.
