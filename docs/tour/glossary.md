[Tour](README.md) · **Glossary**

# Glossary

Each term links to the chapter where it first matters.

**Amplitude.** One complex number in the list that makes up a quantum state.
Its squared magnitude is the probability of that outcome. Four sites have 16.
([1](01-why-simulate-empty-space.md))

**Antiparticle.** The partner of a particle with the opposite charge; here the
positron, charge $-1$, which sits on odd sites. ([2](02-four-sites-four-qubits.md))

**Bare vacuum.** Empty space as the mass and field terms define it: no
particles, no field. Stored as $|0101\rangle$. It is the starting state of
every simulation and is not the true ground state once hopping is on.
([3](03-why-pairs-appear.md))

**Bitstring.** The result of measuring all four qubits once, such as `1001`.
Site 1 is the leftmost digit; `0` means the site is occupied by a fermion.
([2](02-four-sites-four-qubits.md))

**Charge variance, Var(Q).** The spread of the total charge across shots. Zero
without noise, 1 when the state is fully scrambled. The main warning light.
([6](06-warning-light.md))

**Chiral condensate.** The ground-state average, per site, of the operator the
mass multiplies. On this lattice it equals the particle density minus
one half. ([7](07-ground-state-from-bitstrings.md))

**Confinement.** In one dimension the energy between two opposite charges
grows in proportion to their separation, so they can never be pulled apart;
eventually a new pair forms instead. ([1](01-why-simulate-empty-space.md))

**CP.** Swapping every particle with its antiparticle (C) and reflecting the
chain end to end (P). An exact symmetry of this model among states with zero
total charge. ([6](06-warning-light.md))

**Crossover.** The step size below which Trotter error is smaller than the
shot noise, so further refinement cannot be seen. 0.18 at 1000 shots and
time 3. ([4](04-time-steps-and-shots.md))

**CX.** The controlled-NOT, a two-qubit gate that flips one qubit when the
other is $|1\rangle$. The gate whose count sets how fast errors pile up.
32 per time step here. ([4](04-time-steps-and-shots.md))

**Density matrix.** The description of a quantum state that is uncertain in
the classical sense too: a mixture of states with probabilities. Needed once
noise acts. Here 16 by 16. ([5](05-gate-errors.md))

**Dephasing.** Noise that scrambles the relative phase of $|0\rangle$ and
$|1\rangle$ without flipping either. Invisible to the charge check when it
acts between time steps. ([6](06-warning-light.md))

**Depolarising noise.** With some probability, a gate's qubits are replaced by
a completely random state. The main gate-noise model.
([5](05-gate-errors.md))

**Dip.** The drop in particle density from its first peak ($t = 0.89$) to its
first trough ($t = 1.76$), 0.33 in the exact evolution. The signal the
project tries to protect. ([3](03-why-pairs-appear.md))

**Exact diagonalisation.** Finding every energy level and state of the
Hamiltonian by computer, possible because the matrix is only 16 by 16. The
answer key for everything. ([8](08-how-we-know-it-is-right.md))

**Failure threshold, $\epsilon^\ast$.** The two-qubit error rate at which the
dip is no bigger than twice its shot noise. ([5](05-gate-errors.md))

**Gate.** One operation on one or two qubits. A circuit is a sequence of
gates followed by a measurement. ([1](01-why-simulate-empty-space.md))

**Gauss's law.** The rule that ties the electric field to the charges. In one
dimension it fixes the field completely, which is why the field needs no
qubits. ([2](02-four-sites-four-qubits.md))

**Hamiltonian.** The energy operator. Here a 16-by-16 matrix made of three
terms: mass, hopping and field. ([2](02-four-sites-four-qubits.md))

**Hopping.** The term that moves a fermion to a neighbouring site; on this
lattice every hop creates or removes a particle-antiparticle pair. It is the
lattice version of kinetic energy. ([2](02-four-sites-four-qubits.md))

**Infidelity, $\epsilon$.** One minus the average overlap between what a gate
does and what it should do. The error rate hardware papers quote.
([5](05-gate-errors.md))

**Jordan-Wigner transformation.** The recipe that writes fermion operators as
qubit operators, adding strings of $\hat Z$ to get the fermion signs right.
([2](02-four-sites-four-qubits.md))

**Link.** The gap between two neighbouring sites, where the electric field
lives. Four sites have three links. ([2](02-four-sites-four-qubits.md))

**Particle density, $\nu$.** The average number of particles plus
antiparticles per site. $\nu = 0.5$ means one pair on average at four sites.
([3](03-why-pairs-appear.md))

**Pauli string.** A product of one Pauli matrix ($X$, $Y$, $Z$) or the
identity per qubit, such as `ZIZI`. The Hamiltonian is a sum of 14.
([2](02-four-sites-four-qubits.md))

**Quench.** Changing the Hamiltonian suddenly under a state, here by switching
hopping on under the bare vacuum. ([3](03-why-pairs-appear.md))

**Qubit.** A two-level quantum system: one per site here, recording whether
the site is occupied. ([1](01-why-simulate-empty-space.md))

**Readout error.** A measured bit reported wrong. Correctable on average by
inverting a 2-by-2 matrix per qubit. ([5](05-gate-errors.md))

**Residual.** The size of the disagreement between two independent
calculations of the same thing. Reported instead of "pass" or "fail".
([8](08-how-we-know-it-is-right.md))

**Shot.** One run of a circuit ending in one measured bitstring. Averages need
many; their uncertainty falls as one over the square root of the number.
([4](04-time-steps-and-shots.md))

**Sign problem.** Why lattice Monte Carlo fails in real time: the weights of
configurations become complex phases that cancel.
([1](01-why-simulate-empty-space.md))

**Site.** One of the evenly spaced points that replace the continuous line.
An approximation, like a grid in a finite-difference calculation.
([2](02-four-sites-four-qubits.md))

**SKQD.** Sample-based Krylov quantum diagonalisation: run time-step circuits,
keep which bitstrings appear, and diagonalise the Hamiltonian in their span on
a classical computer. ([7](07-ground-state-from-bitstrings.md))

**Staggered fermions.** Putting the particle and antiparticle parts of the
electron on alternating sites, to avoid a fake extra particle the lattice
would otherwise create. ([2](02-four-sites-four-qubits.md))

**Trotter step.** One short slice of time evolution, built from pieces of the
Hamiltonian applied one after another. Its error shrinks as the step shrinks.
([4](04-time-steps-and-shots.md))

**VQE.** The variational quantum eigensolver: tune a circuit's angles with a
classical optimiser to minimise the measured energy. Not used here.
([7](07-ground-state-from-bitstrings.md))
