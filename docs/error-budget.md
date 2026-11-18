# The error budget, explained from the ground up

Abror Niyazmetov, PHY 420, Fall 2026. Written 2026-10-13.

Companion to `derivation.md`. That document gets you to the Hamiltonian. This
one covers everything after it: noise models, the symmetry diagnostics, the
error budget, the failure threshold, and SKQD, in the order of thesis
Sections 3.4 to 3.8. Section 10 is a self-test.

Same rule as before: **trust the number, not the prose.** Every claim below is
printed by a script (`scripts/06`, `07`, `08`) or checked by a test.

---

## 1. Why compute probabilities instead of sampling

A real device gives you shots. A shot-sampled simulation of a noisy device mixes
two things in every number: the bias that gate noise puts into the state,
and the random scatter from having finitely many shots. You could average
many runs to separate them, but that is slow and never exact.

So the code never samples to get a bias. Aer's density-matrix simulator carries
the full 16 x 16 density matrix through every noisy gate and hands back the
exact probability of each of the 16 outcomes. From that:

- the bias of any observable is an exact number;
- the shot noise an M-shot run *would* have is an exact number too,
  `sigma = sqrt( sum_b p_b (x_b - xbar)^2 / M )`.

Shots only appear when you want a realistic data set, and then they are drawn
from the exact distribution (multinomial sampling), which is statistically
identical to running the circuit M times.

*The formula to know:* standard error of a mean = sqrt(variance per shot
/ M). The particle number per shot is 0, 1, 2, 3 or 4, so do not use the
Bernoulli `sqrt(p(1-p)/M)`. Checked: real scatter / quoted error = 0.97 to 1.01.

## 2. Everything comes from one measurement basis

Particle number, charge density, link fields and total charge are all built from
`Z_n` only. So each is a function of the bitstring alone, and a single Z-basis
measurement gives the physics and the diagnostics, with no extra circuits.
Most VQE-style observables need several measurement bases; these do not.

Useful identity: `nu = 1/2 + <psibar psi>`. The particle density and the chiral
condensate are the same operator up to a constant.

## 3. The noise models, and the factor of 4/3 that is easy to get wrong

Depolarizing channel on d dimensions: `E(rho) = (1 - lambda) rho + lambda I/d`.

Hardware data sheets quote the *average gate infidelity* epsilon (from
randomized benchmarking), not lambda. The relation is

    epsilon = lambda (d - 1) / d        so for a CX (d = 4): lambda = 4 epsilon / 3

Get this wrong and every threshold in the thesis is off by 4/3. The code checks
it against Qiskit's own `average_gate_fidelity` (agreement 4e-12).

Where the errors go: on CX (epsilon_2), on SX and X (epsilon_1 = epsilon_2/10, an
assumption; changing it to 0 or 0.3 moves the threshold +8 % / -13 %). RZ is a
software frame change on IBM-style hardware, so it gets no error. One
second-order Trotter step compiles to 32 CX, 24 SX and 53 RZ.

Why compile first: if you hand Aer the circuit with RXX/RYY/RZZ gates, it runs
them natively and a noise model attached to "cx" never fires. Compile to
{RZ, SX, X, CX}, then add noise. The same compiled circuit is used for gate
counts.

*Readout error* is a classical bit flip with probability q, applied to the
outcome probabilities. Closed form: `<Z> -> (1 - 2q) <Z>`, so
`nu -> 1/2 + (1 - 2q)(nu - 1/2)`. Code matches to 6e-17.

## 4. The diagnostics, and why a violation can be turned into a number

Ideal evolution conserves Q, so ideally every shot has Q = 0. Three diagnostics:

- `Var(Q)`, the spread of the measured charge;
- leakage `1 - P(Q = 0)` (what Martinez et al. reported as a fidelity proxy);
- CP asymmetry `<L_1> - <L_3>`.

*CP survives Trotterization.* Each Trotter factor (diagonal part, even bonds,
odd bond) maps the neutral sector to itself and commutes with CP there. Even
bonds map to even bonds under reflection at N = 4. So the noiseless circuit has
`<L_1> = <L_3>` exactly at any dt (the argument is general; the check, at
dt = 0.5, 0.25 and 0.1 at three times, gives 5e-16). A measured CP asymmetry is
hardware error, never Trotter error.

*From violation to fidelity.* Assume global depolarizing noise,
`rho -> F rho + (1 - F) I/16`. Every diagonal expectation is then linear in F:
`<x> = F <x>_ideal + (1 - F) <x>_mixed`. For Q:

    <Q> = 0 in both, so Var(Q) = <Q^2> = F * 0 + (1 - F) * (N/4)

because in the maximally mixed state each Z_n is an independent +-1 coin, Q is
half their sum (plus a constant), and the variance of half a sum of 4 coins is
4/4 = 1 = N/4. Solve: F_Q = 1 - Var(Q)/(N/4). Exact for the global channel
(7e-16).

## 5. The blind spot, and the correction found by checking

Z errors commute with Q. So dephasing applied to the logical state, between
Trotter steps, never changes the charge: Var(Q) stays exactly 0 (1e-33) while
the particle density drifts badly. The diagnostic is blind to it. Martinez et
al.'s fitted error model was exactly this kind of dephasing.

Dephasing on the CX gates is different: the simulation gives Var(Q) = 0.88.
Why: a compiled RXX is a basis change,
two CXs with a rotation between them, then the basis change undone. Between the
CXs the qubits are in a rotated basis, and a Z error there is an X or Y error on
the logical state, which does move charge. Whether an error is visible depends
on where in the compiled circuit it happens, not only on its type.

## 6. The error budget

Four exact distributions, each differing from the last by one source:

    exact  ->  Trotter (noiseless circuit)  ->  gate (noisy)  ->  measured (readout)

The biases are differences of consecutive ones, so they telescope and add up
exactly: `total = b_Trot + b_gate + b_read`, no cross terms. Shot noise sigma
comes from the last distribution. One measurement's rms error is
`sqrt(total_bias^2 + sigma^2)`. "Complete" was tested: 600 simulated 500-shot
experiments have mean squared error 0.999 times `b^2 + sigma^2`.

Results at t = 1.76, with 1000 shots and 1 % readout error:

| eps_2 | best dt | CX | rms error |
|---|---|---|---|
| 3e-4 | 0.35 | 160 | 0.032 |
| 1e-3 | 0.44 | 128 | 0.063 |
| 3e-3 | 0.59 | 96 | 0.119 |

Worse hardware means a coarser best step. Refining dt past the optimum makes
the answer worse.

## 7. The failure threshold, and why it is about two errors

The exact density rises to 0.463 at t = 0.89 and falls to 0.133 at t = 1.76 as
pairs recombine. Each channel modelled here, applied to the bare vacuum with no
dynamics, can only push nu up toward 1/2; it cannot make nu fall. Under these
channels noise therefore erodes a fall rather than creating one, so the
recombination dip `D = nu(0.89) - nu(1.76)` is the part of the signal this
noise cannot fake. That is a heuristic, not a theorem: in the circuit the noise
acts between coherent steps, and a coherent error (an over-rotated RXX, say) is
a different unitary with its own oscillations, which could shift or mimic the
dip. Coherent errors are not modelled.

Threshold eps*: the smallest two-qubit error at which `D <= 2 sigma_D`.

| dt | CX (trough) | eps* (1000 shots) | eps* x CX |
|---|---|---|---|
| 0.1 | 576 | 3.3e-3 | 1.92 |
| 0.2 | 288 | 6.9e-3 | 1.98 |
| 0.3 | 192 | 1.06e-2 | 2.03 |

*Why the product is constant.* eps* x n_CX ~ 2 says the signal is lost once
the summed infidelity of the CX gates in the deepest circuit reaches about two.
For depolarizing noise that is about 2.5 Pauli errors on average, since a
two-qubit depolarizing channel makes a non-identity Pauli error with probability
(15/16) lambda = 1.25 eps. At small error D shrinks roughly in proportion to the
fidelity of the trough circuit (at eps2 = 1e-3: D/D0 = 0.51, 0.72, 0.81 against
F_Q = 0.54, 0.74, 0.82 for dt = 0.1, 0.2, 0.3). Near the threshold D falls
faster than F, and F_Q at eps* is still 0.11 to 0.13. D then reaches zero at a
finite eps, eps*(M = inf), only about 1.37 times the 1000-shot threshold. That
zero crossing is why ten times the shots moves eps* only ~20 %. A pure
exponential, D ~ F ~ exp(-c eps n_CX) all the way down, would give +44 % for
ten times the shots and an infinite threshold without shot noise; scripts/06
shows neither.

## 8. Does the diagnostic predict the threshold?

The test uses only data an experiment has: distributions at one noise level
eps_0. Read F at each time from Var(Q). Assume `F(eps) = exp(-kappa eps)`, so
`kappa = -ln F / eps_0`. Undo the global channel to estimate the ideal
distribution, `p_ideal = (p - (1 - F)/16) / F`, then re-noise it at every other
eps to predict D and sigma_D. Where they cross is the predicted threshold.

Result: for depolarizing gate noise, predicted / true = 1.10 to 1.17 at 1000
shots when Var(Q) at eps_0 comes from the full, unsampled distribution. Every
case is an overestimate (optimistic): real local noise erodes the dip slightly
faster than a global channel with the same charge violation. A real experiment
has shot-sampled data at eps_0, so `scripts/06` part E repeats the prediction
on 200 data sets of 1000 shots sampled from the eps_0 distribution: the median
ratio is 1.09 to 1.18, the 16-84 % range has a half-width of 0.11 to 0.16, and
8 to 20 % of runs come out pessimistic. One run is optimistic 80 to 92 % of the
time, but its scatter is as large as the bias, so it cannot measure its own
bias. For dephasing on the CX gates (part G) the prediction is much worse:
predicted / true = 1.27 to 1.42. How well Var(Q) predicts the threshold depends
on the kind of noise. For logical dephasing the prediction is infinite, because
Var(Q) never moves.

## 9. SKQD in one paragraph, and the formula that checks it

Sample bitstrings from `U^k |vacuum>` for k = 1..K, keep the distinct ones with
Q = 0, diagonalize H restricted to their span. Exact diagonalization in a
subspace can never go below E_0 (variational bound, checked on all 63
subsets), so noise can only hurt by wasting shots or missing bitstrings.

At N = 4 the neutral sector has 6 bitstrings, so the bias is computable exactly.
For the sampled set to lie inside T, every shot at every step must avoid the
neutral bitstrings outside T:

    P(sampled set inside T) = prod_k ( 1 - sum_{b neutral, b not in T} p_k(b) )^M

Inclusion-exclusion over the 64 subsets turns "inside T" into "exactly S", and
the expected energy is the probability-weighted average of the subspace
energies. Checked against 20,000 simulated runs (1.020e-3 vs 1.002(40)e-3).

The physics: the ground state puts only 6.5e-4 of its weight on |1100>, the
longest flux string. Confinement suppresses exactly the bitstring SKQD is most
likely to miss. Noise makes it slightly worse (15 shots per step instead of 12
for 1e-3 accuracy), because noisy shots land on charged bitstrings and are
thrown away.

---

## 10. Self-test

1. Why is a density-matrix simulation better than averaging many shot-sampled
   runs for measuring gate-error bias?
2. A data sheet says the CX error is 3e-3. What lambda goes into Aer, and why?
3. Why must the circuit be compiled before the noise model is attached?
4. Derive Var(Q) = N/4 for the maximally mixed state in two lines.
5. Why does the CP asymmetry vanish for the noiseless Trotter circuit at any dt?
6. Give the physical reason logical dephasing is invisible to Var(Q), and the
   reason CX-level dephasing is not.
7. Why are the four biases guaranteed to add up with no cross terms?
8. Why is the recombination dip, not the peak, the signal the modelled noise
   cannot fake? What kind of error could fake it?
9. Explain eps* x n_CX ~ 2 without a formula.
10. For depolarizing noise, with Var(Q) taken from the unsampled distribution,
    the predicted threshold is always too high. What does that say about real
    local noise compared with a global channel, why can a single 1000-shot run
    not measure it, and why is the prediction worse for CX dephasing?
11. Why can SKQD never report an energy below E_0, even on a broken device?
12. Why does gate noise make SKQD slightly worse at N = 4 instead of better?
