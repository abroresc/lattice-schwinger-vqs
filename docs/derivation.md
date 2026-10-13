# The N = 4 lattice Schwinger model, derived from scratch

Abror Niyazmetov, PHY 420, Fall 2026. Last revised 2026-09-10.

This document does one thing: it gets you from "quantum electrodynamics in one
space dimension" to an explicit list of fourteen Pauli strings on four qubits,
with every step reproducible on a whiteboard and every claim checkable against
a number. It assumes one quantum mechanics course and no field theory.

Read section 0 first even if you skip everything else. It is the reason the
project is safe.

---

## 0. Why this is not as hard as it sounds

The words "quantum electrodynamics" and "lattice gauge theory" suggest a PhD's
worth of prerequisites. They are not required here, because of one fact:

> Once you put a field theory on a finite lattice, it stops being a field
> theory and becomes a finite-dimensional matrix problem.

The four-site Schwinger model is a 16 by 16 Hermitian matrix. You already know
how to find eigenvalues of a Hermitian matrix. Everything below is machinery
for constructing that matrix.

What is actually hard is bookkeeping. Minus signs, operator ordering, index
conventions. That is not a comprehension failure, it is the normal experience,
and it is why every stage below ends with a number that has a known value.

**The rule for this whole document: trust the spectrum, not the prose.** If
something written here disagrees with an eigenvalue, the eigenvalue is right.

---

## 1. The continuum theory, in one page

The Schwinger model is QED with two of the three spatial dimensions deleted.
Its Lagrangian is

```
    L  =  psi-bar ( i gamma^mu D_mu  -  m ) psi   -   (1/4) F_mu-nu F^mu-nu
```

Three pieces, and you need to know what each is for rather than how to compute
with it.

- `psi` is the matter field, a two-component spinor. Two components because in
  1+1 dimensions the Dirac equation needs two, one for the particle and one
  for the antiparticle. There is no spin up and spin down here; the two
  components are a particle-antiparticle pair, not a spin doublet.
- `D_mu = partial_mu + i g A_mu` is the covariant derivative. The gauge field
  `A_mu` is there because the theory is required to be unchanged under a local
  phase rotation of `psi`. That requirement is what forces the photon to
  exist. This is the only thing you need from gauge invariance.
- `F_mu-nu` is the field strength, the photon's kinetic term.

**The simplification that makes everything work.** The magnetic field is
`F_12`, and there is no second spatial direction, so there is no `F_12`. The
electromagnetic field has exactly one component, the scalar electric field `E`.
No polarisation, no propagating photon, one number per point.

**What survives.** In one dimension the electric field between two opposite
charges cannot spread out. Field lines have nowhere to go. So the field
between them is constant, independent of separation, and the potential energy
grows linearly with distance:

```
    V(r)  =  (g^2 / 2) r
```

Pull them apart and the energy keeps climbing. It never becomes cheaper to
separate them further, so isolated charges do not exist. That is confinement,
and it is the property that makes this model the standard toy version of QCD.

Keep pulling and eventually the stored field energy exceeds `2m`, the cost of
making a fresh particle-antiparticle pair out of nothing. At that point the
vacuum makes one, the new pair screens the original charges, and the flux
string snaps. That is string breaking.

**Why anyone simulates this on a quantum computer.** The workhorse of lattice
gauge theory is Monte Carlo sampling in imaginary time. Continue those results
to real time and the sampling weights become complex, so the statistical error
grows exponentially with the time you want to reach. That is the sign problem.
A quantum processor evolves in real time natively and never forms those
weights. This is one of the few places where the argument for quantum
simulation is a physics argument rather than a sales pitch.

---

## 2. Putting fermions on a lattice, and the problem that creates

Discretise space into `N` sites separated by lattice spacing `a`. Replace the
derivative in the Dirac equation by a finite difference. Something goes wrong
immediately.

### 2.1 Fermion doubling

The naive symmetric difference `(psi_{n+1} - psi_{n-1}) / 2a` has momentum
space representation `(i/a) sin(pa)`, which vanishes **twice** in the Brillouin
zone: at `p = 0`, which is the physical particle, and at `p = pi/a`, which is
not. You asked for one fermion species and the lattice gave you two. In `d`
dimensions it gives you `2^d`. These extra species are called doublers and
they are an artifact of the discretisation, not physics.

This is not a small technical annoyance. It is a theorem (Nielsen and
Ninomiya) that you cannot remove doublers while keeping locality, chiral
symmetry and the correct continuum limit all at once. Something has to give.

### 2.2 Kogut and Susskind's answer: stagger the components

Give up on keeping both spinor components at the same point. Instead put one
component on even sites and the other on odd sites. The two-component spinor
becomes a single-component field on a lattice of twice the spacing, and the
doubling is absorbed into the staggering rather than added on top of it.

The consequence you will actually use:

> **Even physics sites hold particles. Odd physics sites hold antiparticles,
> and an antiparticle is a HOLE, not an occupation.**

So the vacuum is not "everything empty". It is the filled Dirac sea: every odd
site occupied, every even site empty. That is why the charge density is

```
    rho_n  =  c_n^dag c_n  -  (1 - (-1)^n) / 2
```

The subtracted term is the reference occupation: 1 on odd sites, 0 on even
sites. Subtracting it is exactly what makes the Dirac sea neutral.

**Attribution, since your bibliography annotation depends on it.** Kogut and
Susskind, Phys. Rev. D 11, 395 (1975) is the correct citation for the
Hamiltonian lattice formulation and the staggered construction you are using.
Susskind, Phys. Rev. D 16, 3031 (1977) is the later extension to higher
dimensions. For a 1+1 dimensional Hamiltonian project, cite the 1975 paper.
Your existing annotation is right; the flag in the corrections file can be
closed.

### 2.3 The Kogut-Susskind Hamiltonian

Fermions live on sites, the gauge field lives on the links between them, as a
phase `U_n = exp(i a g A_n)` on the link from site `n` to site `n+1`, with
conjugate electric field `E_n`. The Hamiltonian is

```
    H  =  w sum_n ( c_n^dag U_n c_{n+1} + h.c. )        hopping
       +  m sum_n (-1)^n c_n^dag c_n                    staggered mass
       +  (g^2 a / 2) sum_n E_n^2                       electric energy
```

The staggered `(-1)^n` on the mass is the price of staggering: the two spinor
components enter the mass term with opposite sign.

**Where pair production comes from, mechanically.** The hopping term moves a
fermion from site `n+1` to site `n`. If `n` is even and `n+1` is odd, that
takes an occupied odd site (part of the sea) and an empty even site, and
produces an occupied even site and an empty odd site. In the language of
particles: it has created a particle on the even site and an antiparticle
(a hole) on the odd site. **The hopping term IS the pair creation operator.**
That is the single most useful thing to notice in this section.

---

## 3. Gauss's law, and why there is no gauge field left

This is the structural heart of the project.

Gauss's law in the continuum is `dE/dx = rho`. On the lattice it becomes a
difference equation:

```
    E_n  -  E_{n-1}  =  rho_n
```

In three dimensions this is a constraint: it restricts which states are
physical, but it does not tell you what `E` is, because the field can circulate.

In one dimension there is nowhere to circulate. With open boundaries and no
background field, `E_0 = 0`, and the equation can be **solved**:

```
    E_n  =  L_n  =  sum_{k <= n} rho_k
```

Every link's field is completely determined by the total charge to its left.
The gauge field is not an independent degree of freedom at all. So instead of
imposing Gauss's law on a bigger Hilbert space, substitute it in and delete
the gauge field entirely.

**What this buys.** Four sites need four qubits. Not four qubits plus a
truncated bosonic register on every link. No cutoff on the field, no
truncation error, no unphysical sector to project out.

**What it costs.** `L_n` is a sum over everything to the left, so `L_n^2` is
long-ranged in the spins. Expanding the square gives `Z_m Z_m'` couplings
between every pair of sites. You have traded locality for qubits.

**And that long-range coupling is the confinement.** Work out the coefficient
(section 5) and you get `(J/2)(N - m')` for the pair `(m, m')`: it depends only
on the right-hand site, and it grows linearly with distance from the right
edge. A pair separated by more links costs more energy, linearly, because the
flux between them is stored on every link in between. The linear potential of
section 1 is sitting in the coefficient table.

**A consequence that catches people out, and it affects your proposal.**
Because Gauss's law is substituted rather than imposed, it holds identically
on **every** state in the Hilbert space, physical or not. There is no Gauss's
law violation left to measure. `analytics.gauss_law_residual` confirms this:
it is exactly zero on the whole space, not just on physical states.

The proposal's diagnostic has to change. What survives as a symmetry is total
charge, `Q = sum_n rho_n`, which commutes with `H` exactly, so physical
evolution never leaves its charge sector. A noisy device does not respect that.
So `Var(Q)` is the model-independent measure of how far noise has pushed the
state out of the physical sector, and it plays the role the proposal assigned
to Gauss's law. There is a second one, described in section 7.

---

## 4. Jordan-Wigner, the crux

Now turn a fermion problem into a qubit problem.

**The problem.** Fermion operators on different sites anticommute. Spin
operators on different sites commute, because they act on different tensor
factors. So "spin up equals occupied" cannot be right on its own: it gets the
statistics wrong.

**The fix.** Attach a string of `Z` operators:

```
    c_j  =  ( Z_1 Z_2 ... Z_{j-1} )  sigma^-_j
```

Each `Z` contributes `-1` for every occupied site to the left of `j`. That is
exactly the minus sign fermionic statistics demands. The price is that the
operator is nonlocal: an operator on site 5 touches sites 1 through 4.

**Check it, do not believe it.** Take `c_1 = sigma^-_1` and
`c_2^dag = Z_1 sigma^+_2`. Using `sigma^- Z = sigma^-` and
`Z sigma^- = -sigma^-`:

```
    c_1 c_2^dag  =  + sigma^-_1 sigma^+_2
    c_2^dag c_1  =  - sigma^-_1 sigma^+_2
```

They sum to zero. The operators anticommute, which is what was wanted.
`jordan_wigner.anticommutation_residuals` does this for every pair on the full
Hilbert space and returns exactly zero.

**Two useful consequences.**

The number operator is local, because the string cancels against itself:

```
    n_j  =  c_j^dag c_j  =  sigma^+_j sigma^-_j  =  ( I + Z_j ) / 2
```

With `sigma^+ = |0><1|` this means `|0>` is the OCCUPIED state. Combined with
the staggered reference occupation from section 2, that gives

```
    rho_n  =  ( Z_n + (-1)^n ) / 2
```

and the bare vacuum is `Z_n = -(-1)^n`, which for N = 4 is the bit string
`0101`.

**The sign that will cost you an afternoon.** For nearest neighbours the
strings almost cancel, leaving one `Z`:

```
    c_j^dag c_{j+1}  =  sigma^+_j Z_j sigma^-_{j+1}  =  - sigma^+_j sigma^-_{j+1}
```

There is a minus sign. Published spin Hamiltonians are usually printed with a
plus. Both are right. The two differ by the relabelling `c_j -> (-1)^j c_j`,
which changes no spectrum. `hopping_sign_demo()` verifies all three claims
numerically: that the plain string gives a minus, that the alternating string
`prod (-Z_i)` gives a plus, and that the spectra are identical.

**Do not adjust a coefficient to make a sign come out. Run that function.**

### 4.1 Why the warm-up on the Ising model is not optional

The transverse-field Ising model, `H = -J sum X_n X_{n+1} - h sum Z_n`, is
exactly solvable. Jordan-Wigner maps it to free fermions, a Bogoliubov
transformation diagonalises those, and the full spectrum comes out in closed
form. So you can learn the technique where an analytic answer exists, verify
that every one of the `2^N` levels matches, and only then point it at the
Schwinger model, where no analytic answer exists.

Two things were found by doing this rather than reading about it.

**The ring trap.** On a periodic chain the string wraps and the boundary bond
picks up a factor `-P`, with `P` the total fermion parity. So the problem
splits into two sectors, antiperiodic momenta for even parity and periodic for
odd. Worse, at `k = 0` and `k = pi` the pairing term vanishes, so those modes
are already diagonal and their energy is the **signed** `2(h - J cos k)`, not
the positive square root. Using the square root there is wrong by exactly
`|2(h - J)|`, and only in the ordered phase `h < J`. The disordered phase looks
perfect, which is what makes the bug shippable. `figures/02_jordan_wigner.png`
shows both versions against exact diagonalisation.

**A conditioning trap that is not about physics at all.** The textbook route
to the Bogoliubov energies is to form `(A - B)(A + B)`, take its eigenvalues
and square-root them. Since `A` is symmetric and `B` antisymmetric,
`(A + B)^T = A - B`, so that product is `M^T M`. Forming it squares the
condition number, and the square root then amplifies the error on the small
eigenvalues, which are exactly the ones that matter in the ordered phase. The
`Lambda_k` are by definition the singular values of `M`, so take them from an
SVD and never form the product. Measured effect at N = 8, `h/J = 0.05`: worst
residual improved from 7e-10 to 1e-14.

Neither of these is in the reading list. Both were found because there was a
known answer to disagree with.

---

## 5. The N = 4 Hamiltonian

Apply Jordan-Wigner to the Kogut-Susskind Hamiltonian with the gauge field
eliminated, and everything collapses to a spin chain:

```
    H  =  (m/2) sum_{n=1..N}   (-1)^n Z_n
       +   w    sum_{n=1..N-1} ( sigma^+_n sigma^-_{n+1} + h.c. )
       +   J    sum_{n=1..N-1} L_n^2

    L_n  =  (1/2) sum_{m=1..n} ( Z_m + (-1)^m )
```

with `w = 1/(2a)` and `J = g^2 a / 2`.

**This is the spin Hamiltonian of Martinez et al., Nature 534, 516 (2016)**
(Methods Eq. (3) in arXiv:1605.04570v1; Eq. (2) there is Gauss's law), who run
exactly this system at N = 4 with four qubits on a trapped-ion device, open
boundaries, zero background field. It is the reference to check your
derivation against, and it is a printed peer-reviewed equation for precisely
the object Gate 2 asks for.

**It is not Klco et al. (2018).** See section 8.

### 5.1 Expanding the electric term

Write `L_n = (S_n + s_n)/2` with `S_n = sum_{m<=n} Z_m` and
`s_n = sum_{m<=n} (-1)^m`, which is 0 for even `n` and -1 for odd `n`. Then

```
    L_n^2  =  (1/4) [ S_n^2  +  2 s_n S_n  +  s_n^2 ]
    S_n^2  =  n I  +  2 sum_{m < m' <= n} Z_m Z_m'
```

Summing over links `n = 1 .. N-1`:

| term | coefficient |
|---|---|
| `Z_m Z_m'`, `m < m'` | `(J/2)(N - m')` for `m' <= N-1`, else 0 |
| `Z_m` | `(J/2) sum_{n=m}^{N-1} s_n` |
| identity | `(J/4) sum_{n=1}^{N-1} ( n + s_n^2 )` |

The hopping term in Pauli form is
`sigma^+_n sigma^-_{n+1} + h.c. = ( X_n X_{n+1} + Y_n Y_{n+1} ) / 2`.

`schwinger.hamiltonian_terms` implements this expansion and
`schwinger.hamiltonian_direct` builds the same operator without it, by
squaring `L_n` as a matrix. They agree to 1e-16, so the algebra above is
verified rather than asserted.

### 5.2 The result

**N = 4, the Gate 2 artifact.**

| Pauli string | coefficient | origin |
|---|---|---|
| `IIII` | `+2 J` | electric |
| `IIIZ` | `+0.5 m` | mass |
| `IIZI` | `-0.5 m - 0.5 J` | mass, electric |
| `IZII` | `+0.5 m - 0.5 J` | mass, electric |
| `ZIII` | `-0.5 m - 1 J` | mass, electric |
| `IIXX` | `+0.5 w` | hopping |
| `IIYY` | `+0.5 w` | hopping |
| `IXXI` | `+0.5 w` | hopping |
| `IYYI` | `+0.5 w` | hopping |
| `XXII` | `+0.5 w` | hopping |
| `YYII` | `+0.5 w` | hopping |
| `IZZI` | `+0.5 J` | electric |
| `ZIZI` | `+0.5 J` | electric |
| `ZZII` | `+1 J` | electric |

Fourteen terms on four qubits. Three things in that table are worth being able
to explain without notes.

1. **Site 4 has no `ZZ` coupling.** With four sites there are three links, and
   `L_3` involves sites 1, 2, 3 only. The rightmost site never enters the
   electric term. That is a direct consequence of open boundaries and the
   left-counting definition of `L_n`.
2. **`ZZII` has twice the coefficient of `ZIZI` and `IZZI`.** From
   `(J/2)(N - m')`: the pair (1,2) has `m' = 2` and gets `J`, the pairs (1,3)
   and (2,3) have `m' = 3` and get `J/2`. Sites closer to the left edge are
   separated by more links, so they cost more. Confinement, in a table.
3. **At `m = J` the `IZII` coefficient vanishes.** An accidental cancellation
   between mass and electric contributions at that one parameter point, which
   is why the printed list at `m = w = J = 1` has thirteen terms rather than
   fourteen. Worth knowing before someone asks why the count changed.

---

## 6. The checks, and what each one actually proves

Every number below is produced by `scripts/03_schwinger_spectrum.py`. None of
them is a pass or fail flag; they are residuals, so you can see how much room
is left.

| Check | Residual | What it proves |
|---|---|---|
| Pauli expansion vs direct matrix construction | 9e-16 | The algebra in 5.1 is right |
| Fermion construction vs spin form | 9e-16 | **Jordan-Wigner was applied correctly. This is Gate 2.** |
| Gauss's law is an identity | 0 | The gauge field really was eliminated, not just constrained |
| `[H, Q] = 0` | 0 | Charge bookkeeping is consistent |
| Strong coupling `E0 = -mN/2` at `w = 0` | 0 | Mass and electric terms |
| Free fermion limit, all 16 levels | 9e-16 | Hopping term and the JW string |
| Hellmann-Feynman on the condensate | 3e-11 | Hamiltonian and observable agree |
| CP symmetry inside `Q = 0` | 4e-16 | Structure of the electric term |
| Qiskit `SparsePauliOp` rebuild | 0 | The implementation, independent of the derivation |

Three of these deserve a sentence.

**Strong coupling, `w = 0`.** With no hopping, `H` is diagonal. The electric
term is a sum of squares so it is at least zero, and it vanishes only when
every `L_n = 0`, which forces the bare vacuum. That same configuration also
minimises the mass term, since `(m/2) sum (-1)^n Z_n >= -mN/2` with equality
exactly when `Z_n = -(-1)^n`. One configuration minimises both terms at once,
so there is no competition and `E0 = -mN/2` exactly, for any `J`.

**Free fermions, `m = J = 0`.** Only the hopping survives, which is a
tight-binding chain on `N` sites with open ends. Standing waves give
`eps_k = 2w cos( k pi / (N+1) )` for `k = 1..N`, and the many-body spectrum is
every subset sum. All 16 levels, not just the ground state.

**Hellmann-Feynman.** The staggered mass term is exactly `m N` times the
chiral condensate operator, so `<psi-bar psi> = (1/N) dE0/dm`. Computing the
condensate from the ground-state vector and from a finite difference of
eigenvalues are completely different calculations. They agree to 3e-11. That
is the check that would catch a wrong Hamiltonian and a matching wrong
observable, which no single check can.

---

## 7. A symmetry that was not in the plan

While checking the electric field profile it turned out that `<L_1>` and
`<L_3>` are equal to machine precision at every time. That is not an accident.

Define `CP = R . X^(x)N`, where `R` reflects `n -> N+1-n` and `X^(x)N` flips
every spin (charge conjugation, since it exchanges occupied and empty).
Neither factor alone commutes with `H`. The product does, but **only inside
the neutral sector**, and the reason is the same asymmetry that made the
electric term long-ranged. `L_n` counts charge to the left only, so under
reflection

```
    L_n  ->  Q - L_{N-n}
```

On states with `Q = 0` that is just `-L_{N-n}`, so `L_n^2 -> L_{N-n}^2` and
`sum_n L_n^2` is invariant. On charged states the cross term survives and the
symmetry is broken. Measured: the commutator is 4e-16 inside `Q = 0` and
about 6.5 on the full space.

Two payoffs.

- **A free correctness check.** On any neutral state, `<L_n> = <L_{N-n}>`. At
  N = 4 the two outer links must always carry equal flux. Any simulation of
  the neutral sector that violates this is wrong, and it costs nothing to
  measure.
- **A second noise diagnostic.** CP joins total charge as a symmetry a noisy
  device will break, which matters for M11 and M12 now that Gauss's law is
  unavailable for that job.

This is a small result, but it is the project's own, it was found by checking
rather than reading.

---

## 8. Corrections to the project documents

### 8.1 Klco et al. (2018) is not the reference for this encoding

The proposal and the study map both treat Klco et al., Phys. Rev. A 98, 032331
(2018) as the primary methodological reference and imply the four-site
Hamiltonian can follow their construction. It cannot, and the reason is more
fundamental than a site count.

Klco et al. do **not** eliminate the gauge field. They keep explicit flux-link
registers with an energy cutoff `|l_n| < 1`, and then reduce the qubit count
with symmetry projections. Their two-spatial-site system is 4 fermion sites
plus 4 flux links, 1296 states in 12 qubits, cut to 2 or 3 qubits by Gauss's
law, momentum and parity projections. That is a different encoding with a
truncation error your project does not have.

The correction file already noted the two-site versus four-site issue and
already identified Martinez et al. as the cleaner reference for gauge
elimination. Both are right, and the second point is the important one. Promote
Martinez et al. (2016) to primary encoding reference. Keep Klco et al. for the
symmetry-projection method and as an early hardware precedent, and rewrite any
sentence that implies you are following their construction.

### 8.2 The Gauss's law diagnostic does not exist in this encoding

Section 3 above. The proposal calls this "the piece of the design I most want
scrutinised", so it is worth catching now rather than in November. The
replacement is total charge variance, plus the CP check from section 7. Both
are exact, both are cheap, and both break under noise in the same way Gauss's
law would have.

This is a strengthening rather than a loss: `Var(Q)` and the CP asymmetry are
two independent diagnostics instead of one.

### 8.3 Staggered fermion attribution

Kogut and Susskind (1975) is correct for what you are using. Susskind (1977)
is the higher-dimensional extension. The flag in the corrections file can be
closed. Section 2.2 has the wording to defend it.

---

## 9. Self-test

If you cannot answer these without notes, the corresponding section did not
land.

1. Why does the Schwinger model confine, in one sentence, without writing a
   formula?
2. What is fermion doubling, and what do staggered fermions trade away to fix
   it?
3. Why can the gauge field be eliminated in 1+1 dimensions but not in 3+1?
4. What exactly is bought and what is paid for by that elimination?
5. Where does the `(-1)^n` in the mass term come from?
6. What goes wrong if you set "spin up equals occupied" with no Z string?
7. Your hopping term has a plus sign; the plain Jordan-Wigner string gives a
   minus. Which is wrong?
8. The bare vacuum has zero particles and zero electric field. Pairs appear
   under time evolution. Where does the energy come from?
9. Why is `ZZII` twice `ZIZI`, and what physical statement is that?
10. Why is there no Gauss's law violation to measure here, and what replaces it?
11. You claim your Hamiltonian is right. Name three independent checks and say
    what each would fail to catch on its own.
12. Why is the transverse-field Ising warm-up worth the time it costs?

Answer to 8, since it is the one most often got wrong and the proposal has it
framed loosely: nothing is injected and no external agent pulls charges apart.
The simulation starts in the bare vacuum, which is not an eigenstate of the
full Hamiltonian once `w` is nonzero. Total energy is conserved exactly. What
changes is how it is distributed between electric field and rest mass. The
size of the effect tracks the gap between the bare vacuum energy and the true
ground state energy, which
`scripts/04_pair_production_preview.py` tabulates.

---

## 10. References

Martinez, E. A. et al., *Real-time dynamics of lattice gauge theories with a
few-qubit quantum computer*, Nature **534**, 516 (2016).
**Primary encoding reference.** N = 4, four qubits, open boundaries, gauge
field eliminated. Their Methods Eq. (3) (arXiv v1 numbering) is the Hamiltonian in section 5.

Kogut, J. and Susskind, L., *Hamiltonian formulation of Wilson's lattice gauge
theories*, Phys. Rev. D **11**, 395 (1975). Staggered fermions and the
Hamiltonian lattice formulation.

Susskind, L., *Lattice fermions*, Phys. Rev. D **16**, 3031 (1977). The
higher-dimensional extension. Cited for completeness, not used here.

Klco, N. et al., *Quantum-classical computation of Schwinger model dynamics
using quantum computers*, Phys. Rev. A **98**, 032331 (2018). Flux-link
encoding with a cutoff plus symmetry projections. A different encoding from
this one. Cite for the projection method and as hardware precedent.

Shaw, A. F. et al., *Quantum algorithms for simulating the lattice Schwinger
model*, Quantum **4**, 306 (2020). Source for the `x = 1/(ag)^2`,
`mu = 2m/(ag^2)` dimensionless conventions and for the Gauss's law constraint
in the form `E_r - E_{r-1} - rho_r`.

Sachdev, S., *Quantum Phase Transitions*, 2nd ed., ch. 4. Jordan-Wigner and
the transverse-field Ising model together, which is the pairing you want.

Bañuls, M. C. et al., *Simulating lattice gauge theories within quantum
technologies*, Eur. Phys. J. D **74**, 165 (2020). Orientation review.

---

## Statement on AI assistance

This document and the accompanying code were drafted with AI assistance. The
derivations were checked against exact diagonalisation at every stage, and
every claim about a paper was verified against the paper before it entered
this document. Two errors in the AI-drafted material were caught by those
checks and are documented in place rather than quietly fixed: the Qiskit
endianness adapter in section 5 of `CONVENTIONS.md`, and the Bogoliubov
conditioning issue in section 4.1. The corrections in section 8 were found the
same way.
