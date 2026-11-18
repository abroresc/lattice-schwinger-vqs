# Reading guide, topic by topic

Each topic below is something the project uses. For each one: why it is
needed here, where it appears in this repository, and the paper or book to
read for it. The [tour](tour/README.md) tells the story in order;
[`derivation.md`](derivation.md) does the full derivation. This file is the
index between them and the literature.

## 1. Second-quantised fermion operators

The project represents "a particle is present at this site" as an operator
algebra, not a wavefunction. Staggering, Jordan-Wigner and the Hamiltonian are
all built from creation and annihilation operators and their anticommutation
relations.

- In this repository: `derivation.md` section 2; `jordan_wigner.py`.
- Any graduate quantum mechanics text covers the algebra.

## 2. Lattice discretisation and fermion doubling

Putting a relativistic fermion on a grid naively creates extra particle
species that are not physical. Staggering (topic 3) is the fix, and the
problem has to be seen first for the fix to make sense.

- In this repository: `derivation.md` section 2.1; thesis Section 2.2.
- Nielsen and Ninomiya, Nucl. Phys. B 185, 20 (1981),
  [doi:10.1016/0550-3213(81)90361-8](https://doi.org/10.1016/0550-3213(81)90361-8),
  for why doublers cannot be removed for free.

## 3. Staggered fermions

The trick that fixes doubling here: the components of one fermion are spread
across even and odd sites. It is also where the picture of the vacuum as
"odd sites filled" comes from, and that state is where every circuit starts.

- In this repository: `derivation.md` sections 2.2 and 2.3;
  [tour chapter 2](tour/02-four-sites-four-qubits.md).
- Kogut and Susskind, Phys. Rev. D 11, 395 (1975),
  [doi:10.1103/PhysRevD.11.395](https://doi.org/10.1103/PhysRevD.11.395).

## 4. Gauss's law and eliminating the gauge field

In one dimension Gauss's law fixes the electric field completely, so four
sites cost four qubits rather than four qubits plus three field registers.
The same fact means there is no Gauss's-law violation to measure in this
encoding, because the constraint holds on every state.

- In this repository: `derivation.md` section 3; `schwinger.py`;
  `analytics.gauss_law_residual`.
- Martinez et al., Nature 534, 516 (2016),
  [doi:10.1038/nature18318](https://doi.org/10.1038/nature18318), the encoding
  this project follows (Methods Eq. (3) in arXiv:1605.04570v1).
- Hamer et al., Phys. Rev. D 56, 55 (1997),
  [doi:10.1103/PhysRevD.56.55](https://doi.org/10.1103/PhysRevD.56.55), for
  eliminating the field.

## 5. Jordan-Wigner transformation

A sign error here is the classic way projects like this go wrong. The project
learns the transformation on the transverse-field Ising model first, where a
closed-form answer exists, before applying it to the Schwinger model.

- In this repository: `derivation.md` sections 4 and 4.1; `jordan_wigner.py`
  (`hopping_sign_demo` settles the hopping sign numerically); `ising.py`;
  [`CONVENTIONS.md`](../CONVENTIONS.md) section 4.
- Jordan and Wigner, Z. Phys. 47, 631 (1928),
  [doi:10.1007/BF01331938](https://doi.org/10.1007/BF01331938).
- Sachdev, *Quantum Phase Transitions*, 2nd ed. (2011), which treats
  Jordan-Wigner and the Ising chain together.

## 6. The four-site Hamiltonian as Pauli strings

Topics 1 to 5 combine into the object the rest of the project works on: a
16-by-16 matrix written as a sum of 14 Pauli strings.

- In this repository: `derivation.md` sections 5 and 6;
  `scripts/03_schwinger_spectrum.py` prints the table; thesis Table 1.

## 7. Trotter product formulas

This is how `e^{-iHt}` becomes a circuit. The scaling exponent (slope 1 for
first order, 2 for second order on a log-log error plot) is a prediction that
can fail, which makes it a real test of the circuit.

- In this repository: `circuits.py`; `scripts/05_qiskit_trotter.py`;
  [tour chapter 4](tour/04-time-steps-and-shots.md).
- Childs et al., Phys. Rev. X 11, 011020 (2021),
  [doi:10.1103/PhysRevX.11.011020](https://doi.org/10.1103/PhysRevX.11.011020),
  for Trotter error bounds.
- Shaw et al., Quantum 4, 306 (2020),
  [doi:10.22331/q-2020-08-10-306](https://doi.org/10.22331/q-2020-08-10-306),
  for Trotter and shot-noise bounds on this model.

## 8. The variational principle and VQE

The project proposal planned to use VQE for the ground state. SKQD replaced
it (topic 9), but the variational principle still matters: `<H> >= E_0` for
every state, so an energy below the true ground energy is a bug, not a
discovery. That is the kind of cheap check the project relies on.

- Peruzzo et al., Nat. Commun. 5, 4213 (2014),
  [doi:10.1038/ncomms5213](https://doi.org/10.1038/ncomms5213), the original
  VQE.
- Kokail et al., Nature 569, 355 (2019),
  [doi:10.1038/s41586-019-1177-4](https://doi.org/10.1038/s41586-019-1177-4),
  variational ground states of this model on up to 20 trapped-ion qubits.
- Wang et al., Nat. Commun. 12, 6961 (2021),
  [doi:10.1038/s41467-021-27045-6](https://doi.org/10.1038/s41467-021-27045-6),
  on how noise flattens the cost function until its gradients vanish.

## 9. Krylov subspace methods and SKQD

The method used for the ground state, chosen on 1 October 2026. Sample
bitstrings from time-evolved states, build a subspace from them, and
diagonalise the Hamiltonian inside it on a classical computer. It reuses the
Trotter circuits from topic 7 and needs no optimiser loop.

- In this repository: `skqd.py`; `scripts/07_skqd_ground_state.py`;
  [tour chapter 7](tour/07-ground-state-from-bitstrings.md).
- Yu et al., arXiv:2501.09702 (2025), the method.
- Rosanowski et al., arXiv:2510.26951 (2025), the method applied to the
  Schwinger model. Its four-qubit demonstration (Section IV A) was close to a
  template for this project.
- Miura (2026), [doi:10.1109/QCNC69040.2026.00163](https://doi.org/10.1109/QCNC69040.2026.00163),
  on unsampled bitstrings as a source of bias.

## 10. The error budget

The three error sources scale differently. Trotter error is systematic and
depends on the step size; shot noise falls as one over the square root of the
number of shots; gate and readout error are properties of the hardware. The
main results, the crossover and the failure threshold, are statements about
where one source overtakes another, so the sources have to be separated.

- In this repository: `noise.py`, `diagnostics.py`, `budget.py`;
  `scripts/06` and `scripts/08`; [`error-budget.md`](error-budget.md);
  [tour chapter 5](tour/05-gate-errors.md).
- Nielsen and Chuang, *Quantum Computation and Quantum Information* (2010),
  for the depolarising channel and its link to average gate infidelity.
- Magesan et al., Phys. Rev. Lett. 106, 180504 (2011),
  [doi:10.1103/PhysRevLett.106.180504](https://doi.org/10.1103/PhysRevLett.106.180504),
  for infidelity as the randomised-benchmarking number hardware quotes.
- McKay et al., Phys. Rev. A 96, 022330 (2017),
  [doi:10.1103/PhysRevA.96.022330](https://doi.org/10.1103/PhysRevA.96.022330),
  for why RZ gates carry no error on superconducting hardware.
- Bravyi et al., Phys. Rev. A 103, 042605 (2021),
  [doi:10.1103/PhysRevA.103.042605](https://doi.org/10.1103/PhysRevA.103.042605),
  for correcting readout error by inverting the confusion matrix.

## 11. Symmetries as diagnostics

With no Gauss's-law violation to measure, the noise diagnostics are the
variance of the total charge (`[H, Q] = 0` exactly) and a CP symmetry that
commutes with `H` only inside the neutral sector. Most noise breaks both.
Dephasing on the logical state between Trotter steps breaks neither, which is
the blind spot.

- In this repository: `derivation.md` sections 7 and 8.2;
  `diagnostics.py`; [tour chapter 6](tour/06-warning-light.md).
- Nguyen et al., PRX Quantum 3, 020324 (2022),
  [doi:10.1103/PRXQuantum.3.020324](https://doi.org/10.1103/PRXQuantum.3.020324),
  the same encoding on trapped ions, where keeping only zero-charge outcomes
  helped more than their symmetry-protection scheme.
