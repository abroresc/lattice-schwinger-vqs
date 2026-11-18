[← Previous](08-how-we-know-it-is-right.md) · [Tour](README.md) · *9 of 9*

# 9. Where this stops

What the results cover, and what they do not.

## What was found

At four sites, $m = 0.5\,w$ and $J = w$, with every stage checked against the
exact answer:

1. *The crossover.* At 1000 shots, the Trotter error on the particle density
   falls below the shot noise at a step of 0.18 (time 3) or 0.36 (the first
   trough). Smaller steps after that add gates and nothing you can measure.
2. *The failure threshold.* The pair-production signal is lost when the
   two-qubit error rate times the number of CX gates in the deepest circuit
   reaches about two (1.92, 1.98 and 2.03 for three step sizes), or about one
   per qubit.
3. *The best step.* With gate error, there is an optimal step size, and it
   grows as the hardware gets worse: 0.35, 0.44 and 0.59 for two-qubit errors
   of 0.03 %, 0.1 % and 0.3 %.
4. *The warning light.* The charge variance predicts the threshold without
   the exact answer, 10 to 17 % too optimistically for depolarising noise and
   27 to 42 % for dephasing on the gates. It is blind to dephasing between
   time steps, the noise model fitted to the 2016 experiment.

And the ground state comes out of SKQD with a relative energy error of
$10^{-3}$ from 12 shots per circuit, 15 with gate error.

## What these numbers do not cover

*The noise models are the main limitation.* Depolarising and dephasing noise
act on each gate independently and forget the past. Real devices also have:

- *coherent errors*, such as a gate that always over-rotates by a little,
  which could shift the dip or even imitate one;
- *correlated errors*, where a fault on one qubit makes faults on its
  neighbours more likely;
- *amplitude damping*, which drives every qubit toward $|0\rangle$. That ends
  in $|0000\rangle$, which has total charge 2 with zero spread, so the charge
  variance reads zero and $F_Q$ would report perfect fidelity on a state with
  no signal left.

Farrell et al. saw systematic errors that grow with the rotation angle of the
gates, which no depolarising model produces and which would favour smaller
steps than the ones found here.

*Other choices.* Every pair of qubits is assumed to interact directly; a
chip wired only to nearest neighbours would need extra gates for the field
term between sites 1 and 3. Readout error is left out of the threshold.
"Lost" means $D \le 2\sigma_D$, a reasonable criterion and still a choice.

*Scope.* Everything is at one parameter point on four sites (twelve for the
exact classical check), on a coarse lattice far from the continuum limit, in
simulation with assumed noise. At this size a laptop beats any quantum
computer. Nothing here is a claim about QCD.

*Not yet tested on hardware.* The published pair-production experiments do
not quote their gate errors (Martinez et al., both Farrell et al. papers), so
the threshold rule cannot yet be compared with real data.

## The next question

Whether these numbers survive coherent and correlated errors, and real
hardware. Much of the machinery is ready for it: the error budget in
[`budget.py`](../../src/schwinger/budget.py) accepts any noise model Qiskit
Aer can express, and every new result would still have an exact answer to be
checked against.

---

For the full argument with every reference, read
the [thesis](../../thesis/main.pdf). For the derivation of the Hamiltonian
from the continuum, read [`docs/derivation.md`](../derivation.md). To check a
number, run the script named next to it.

[← Previous](08-how-we-know-it-is-right.md) · [Tour](README.md) · [Back to the front page](../../README.md)
