# Conventions

Every convention this project depends on, in one place, with the date it was
fixed. Half of all "my spectrum is wrong" bugs are a convention that drifted
between two files. When you adopt a new one, write it here first.

Last reviewed: 2026-09-10

---

## 1. Site indexing

- Physics text uses sites `n = 1 .. N`. Code uses `j = 0 .. N-1`, with `n = j+1`.
- The staggering factor `(-1)^n` is computed in exactly one place,
  `schwinger.schwinger.stagger(j)`, which returns `-1` for even `j`
  (physics site 1, 3, 5, odd `n`) and `+1` for odd `j`.
- Getting this backwards swaps particles and antiparticles and moves the
  vacuum. It does not throw an error. It gives a plausible wrong answer.

## 2. Tensor product ordering

- A Pauli string `"XZIY"` means X on site 1, Z on site 2, I on site 3, Y on site 4.
- The matrix is built as `op(site 1) kron op(site 2) kron ... kron op(site N)`,
  so **site 1 is the most significant tensor factor**.
- Basis index: `b = sum_j s_j 2^(N-1-j)`. The bit string read left to right is
  site 1, site 2, ..., site N. `basis_label(5, 4) == "0101"`.

## 3. Spin and occupation

- `|0>` is the `Z = +1` eigenstate ("spin up").
- `sigma^+ = |0><1|`, `sigma^- = |1><0|`, so `sigma^+ sigma^- = (I + Z)/2`.
- Under the Jordan-Wigner convention below this means **`|0>` is OCCUPIED**.
- Charge density: `rho_n = (Z_n + (-1)^n)/2`. On an odd physics site the
  reference occupation is 1 (the Dirac sea is full there); on an even site it
  is 0.
- Bare vacuum: `Z_n = -(-1)^n`, i.e. `|0101...>` for N even. Odd physics sites
  occupied, even sites empty, every `rho_n = 0`, every `L_n = 0`.

## 4. Jordan-Wigner

- Default: `c_j = ( Z_1 Z_2 ... Z_{j-1} ) sigma^-_j`.
- Consequence: `c_j^dag c_{j+1} = - sigma^+_j sigma^-_{j+1}`. **There is a
  minus sign.** It is removable by `c_j -> (-1)^j c_j`, and
  `jw_annihilation(..., alternating_string=True)` implements that version.
- Both conventions give identical spectra. `hopping_sign_demo()` proves all
  three claims numerically. Do not "fix" a coefficient until you have run it.
- `hamiltonian_from_fermions` uses the alternating string so the printed spin
  Hamiltonian carries a plus sign on `w`, matching the published form.

## 5. Qiskit interoperability

Established by numerical check on 2026-09-10, not by reading documentation.

- **Pauli label strings coincide.** `SparsePauliOp("ZIZI").to_matrix()` equals
  `pauli_string("ZIZI")` entry for entry. Do **not** reverse the string.
  Reversing it reflects the Hamiltonian through the middle of the chain, which
  here swaps particles and antiparticles.
- **Qubit numbering is reversed.** Our site `j` is Qiskit qubit `N-1-j`. Use
  `site_to_qiskit_qubit(j, N)` for every circuit index. Never pass a site
  index to a Qiskit gate directly.
- Measurement bitstrings come back in Qiskit's order; convert with
  `qiskit_bitstring_to_sites`.

## 6. Hamiltonian parameters

`H(N, m, w, J)` with

- `m` fermion mass,
- `w` hopping strength (`w = 1/(2a)` in Martinez et al.),
- `J` electric coupling (`J = g^2 a / 2` in Martinez et al.).

Physics depends on ratios, so the tests fix `w = 1` and vary `m/w`, `J/w`.

## 7. Boundary conditions

Open, with zero background electric field (`L_0 = 0`). This is what makes
Gauss's law solvable rather than merely a constraint, and it is what lets N
sites fit in N qubits. Periodic boundaries would reintroduce a global
constraint and are out of scope; the periodic TFIM in `ising.py` exists only
as a Jordan-Wigner exercise.

## 8. Noise conventions (fixed 2026-10-12)

- **Gate error is quoted as average gate infidelity** epsilon, the number from
  randomized benchmarking. Aer's depolarizing parameter is
  `lambda = epsilon d / (d - 1)`, so `lambda = 4 epsilon / 3` for a CX. Always
  go through `noise.depolarizing_lambda`; never pass epsilon to Aer directly.
- **Noise attaches to the compiled circuit.** Compile with
  `circuits.to_hardware_basis` to {RZ, SX, X, CX} first. RZ carries no error.
  The same compiled circuit is used for `circuit_metrics`.
- **Default single-qubit error** is `epsilon_1 = epsilon_2 / 10`. It is an
  assumption; state it whenever a threshold is quoted.
- **Readout error** is applied classically with `noise.apply_readout_error`,
  never inside the density-matrix run, so its bias stays separable.
- **Probabilities from Aer** are indexed like our basis (site N-1 least
  significant), the same coincidence `circuits.statevector` relies on. Checked
  in `tests/test_noise.py`.
