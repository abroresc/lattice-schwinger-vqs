"""
The lattice Schwinger model on N staggered sites, open boundaries, gauge field
eliminated.

The Hamiltonian
---------------
After Kogut-Susskind staggering, solving Gauss's law with open boundary
conditions and zero background field, and applying Jordan-Wigner, the model
becomes a pure spin chain on N qubits with no gauge register at all:

    H = (m/2) sum_{n=1..N}   (-1)^n Z_n                      [staggered mass]
      +  w    sum_{n=1..N-1} ( sig^+_n sig^-_{n+1} + h.c. )  [pair creation / hopping]
      +  J    sum_{n=1..N-1} L_n^2                           [electric field energy]

    L_n = (1/2) sum_{m=1..n} ( Z_m + (-1)^m )

This is the spin Hamiltonian of Martinez et al., Nature 534, 516 (2016)
(Methods Eq. (3) in arXiv:1605.04570v1; Eq. (2) there is Gauss's law), which runs
exactly this system at N = 4 on a trapped-ion device. That paper is the
reference to check the derivation against. Klco et al. (2018) is not: they
keep explicit flux-link registers with an energy cutoff and reduce with
symmetry projections, which is a different encoding.

Why there is no gauge field
---------------------------
Gauss's law in 1+1 dimensions is the difference equation

    E_n - E_{n-1} = rho_n,      rho_n = ( Z_n + (-1)^n ) / 2

With open boundaries and no background field, E_0 = 0, so the equation does
not merely constrain the electric field, it determines it:

    E_n = L_n = sum_{m<=n} rho_m

Every link value is fixed by the charges to its left. So the gauge field can
be substituted away rather than stored. The price is that the electric energy
sum_n L_n^2 becomes long-range in the spins: expanding the square produces
Z_m Z_m' couplings between every pair of sites. The encoding trades locality
for qubits.

Parameters
----------
m : fermion mass
w : hopping strength, w = 1/(2a) in Martinez et al.
J : electric coupling, J = g^2 a / 2 in Martinez et al.

The physics depends on the ratios, so the tests fix w = 1 and vary m/w and J/w.
"""

from __future__ import annotations

import numpy as np

from .paulis import (
    eigh_sorted,
    from_terms,
    identity,
    is_hermitian,
    site_op,
    two_site_op,
)

# ----------------------------------------------------------------------------
# Staggering helper. In code sites are 0-indexed: code site j is physics site
# n = j + 1. So (-1)^n = (-1)^(j+1) = -(-1)^j. Getting this backwards swaps
# particles and antiparticles and moves the vacuum, so it lives in one place.
# ----------------------------------------------------------------------------


def stagger(j: int) -> int:
    """(-1)^n for physics site n = j + 1, given 0-indexed code site j."""
    return -1 if j % 2 == 0 else +1


# ----------------------------------------------------------------------------
# Gauss's law
# ----------------------------------------------------------------------------


def charge_density_op(j: int, n_sites: int) -> np.ndarray:
    """rho_n = ( Z_n + (-1)^n ) / 2 on 0-indexed site j."""
    return 0.5 * (site_op("Z", j, n_sites) + stagger(j) * identity(n_sites))


def electric_field_op(j: int, n_sites: int) -> np.ndarray:
    """
    L_n on the link to the right of 0-indexed site j, i.e. L_n = sum_{m<=n} rho_m.

    This is the solved Gauss's law, not an independent degree of freedom.
    """
    dim = 2 ** n_sites
    out = np.zeros((dim, dim), dtype=complex)
    for m in range(j + 1):
        out += charge_density_op(m, n_sites)
    return out


def total_charge_op(n_sites: int) -> np.ndarray:
    """Q = sum_n rho_n. Commutes with H, so it labels superselection sectors."""
    dim = 2 ** n_sites
    out = np.zeros((dim, dim), dtype=complex)
    for j in range(n_sites):
        out += charge_density_op(j, n_sites)
    return out


# ----------------------------------------------------------------------------
# The Hamiltonian, built two independent ways
# ----------------------------------------------------------------------------


def hamiltonian_direct(n_sites: int, m: float, w: float, J: float) -> np.ndarray:
    """
    Build H by assembling operators directly: mass, hopping, then sum_n L_n^2
    with L_n squared as a matrix.

    This is the transparent construction. It follows the printed equation term
    for term and does no algebra, so it is the one to trust when the two
    constructions disagree.
    """
    dim = 2 ** n_sites
    ham = np.zeros((dim, dim), dtype=complex)

    for j in range(n_sites):
        ham += 0.5 * m * stagger(j) * site_op("Z", j, n_sites)

    for j in range(n_sites - 1):
        hop = two_site_op("+", j, "-", j + 1, n_sites)
        ham += w * (hop + hop.conj().T)

    for j in range(n_sites - 1):
        L = electric_field_op(j, n_sites)
        ham += J * (L @ L)

    return ham


def hamiltonian_terms(n_sites: int, m: float, w: float, J: float) -> list[tuple[str, complex]]:
    """
    Build H as an explicit list of (Pauli string, coefficient).

    This is the form printed in the thesis (Table 1) and the form Qiskit
    accepts.

    The electric term is expanded by hand here. Writing L_n = (S_n + s_n)/2
    with S_n = sum_{m<=n} Z_m and s_n = sum_{m<=n} (-1)^m,

        L_n^2 = (1/4) [ S_n^2 + 2 s_n S_n + s_n^2 ]
        S_n^2 = n I + 2 sum_{m < m' <= n} Z_m Z_m'

    so summing over links n = 1..N-1 gives

        Z_m Z_m' coefficient  =  (J/2) (N - m')        for m < m' <= N-1
        Z_m      coefficient  =  (J/2) sum_{n=m}^{N-1} s_n   plus the mass term
        identity coefficient  =  (J/4) sum_{n=1}^{N-1} ( n + s_n^2 )

    The ZZ coefficient depends only on the right-hand site m', because sites m
    and m' appear together in every L_n with n >= m', and there are N - m' such
    links. That pattern comes from counting the field from the left end. The
    physics does not depend on the choice: a neutral pair costs J for every link
    between its charges, which is the linear potential.

    `hamiltonian_terms` and `hamiltonian_direct` are checked against each other
    in the test suite, so the algebra above is never taken on trust.
    """
    terms: list[tuple[str, complex]] = []

    # s_n = sum_{m=1..n} (-1)^m, indexed by physics link n = 1..N-1
    s = [0] * (n_sites + 1)
    for n in range(1, n_sites + 1):
        s[n] = s[n - 1] + (1 if n % 2 == 0 else -1)

    z_coeff = np.zeros(n_sites)
    const = 0.0

    # staggered mass
    for j in range(n_sites):
        z_coeff[j] += 0.5 * m * stagger(j)

    # electric energy, expanded
    for n in range(1, n_sites):                      # physics link index
        const += 0.25 * J * (n + s[n] ** 2)
        for mm in range(1, n + 1):                   # physics site index
            z_coeff[mm - 1] += 0.5 * J * s[n]
        for a in range(1, n + 1):
            for b in range(a + 1, n + 1):
                chars = ["I"] * n_sites
                chars[a - 1] = "Z"
                chars[b - 1] = "Z"
                terms.append(("".join(chars), 0.5 * J))

    # hopping: w ( sig^+_n sig^-_{n+1} + h.c. ) = (w/2)( X_n X_{n+1} + Y_n Y_{n+1} )
    for j in range(n_sites - 1):
        for p in ("X", "Y"):
            chars = ["I"] * n_sites
            chars[j] = p
            chars[j + 1] = p
            terms.append(("".join(chars), 0.5 * w))

    for j in range(n_sites):
        if abs(z_coeff[j]) > 0:
            chars = ["I"] * n_sites
            chars[j] = "Z"
            terms.append(("".join(chars), z_coeff[j]))

    if abs(const) > 0:
        terms.append(("I" * n_sites, const))

    return collect_terms(terms)


def collect_terms(terms: list[tuple[str, complex]]) -> list[tuple[str, complex]]:
    """Merge duplicate Pauli strings and drop numerically zero coefficients."""
    acc: dict[str, complex] = {}
    for label, coeff in terms:
        acc[label] = acc.get(label, 0.0) + coeff
    out = [(k, v) for k, v in acc.items() if abs(v) > 1e-14]
    # Sort by weight then alphabetically, so printed tables are stable.
    out.sort(key=lambda kv: (sum(c != "I" for c in kv[0]), kv[0]))
    return out


def hamiltonian(n_sites: int, m: float, w: float, J: float) -> np.ndarray:
    """H from the Pauli-string list. The default constructor for the project."""
    return from_terms(hamiltonian_terms(n_sites, m, w, J), n_sites)


def construction_residual(n_sites: int, m: float, w: float, J: float) -> float:
    """Largest entrywise difference between the two independent constructions."""
    return float(np.max(np.abs(hamiltonian(n_sites, m, w, J) - hamiltonian_direct(n_sites, m, w, J))))


# ----------------------------------------------------------------------------
# States
# ----------------------------------------------------------------------------


def bare_vacuum(n_sites: int) -> np.ndarray:
    """
    The strong-coupling (bare) vacuum: the filled Dirac sea.

    Every charge density is zero, so every link field is zero. That requires
    Z_n = -(-1)^n, i.e. odd physics sites occupied and even sites empty. In
    bit-string form with |0> = occupied that is |0101...>.

    This is the state the real-time simulation starts from. It is not the true
    ground state of H at finite w, which is why pairs appear.
    Nothing is injected and the total energy stays at -mN/2. On this state the
    mass term sits at its minimum (-mN/2) and the electric term at its minimum
    (0), so making pairs raises both. The hopping term, zero here, goes negative
    and pays for both (scripts/09 prints the split: at t = 0.89, the first
    maximum of the density, mass -0.07, electric +0.97, hopping -1.90).
    """
    bits = "".join("0" if j % 2 == 0 else "1" for j in range(n_sites))
    state = np.zeros(2 ** n_sites, dtype=complex)
    state[int(bits, 2)] = 1.0
    return state


def bare_vacuum_label(n_sites: int) -> str:
    return "".join("0" if j % 2 == 0 else "1" for j in range(n_sites))


def ground_state(n_sites: int, m: float, w: float, J: float) -> tuple[float, np.ndarray]:
    """Exact ground state by dense diagonalisation; every other method is checked against it."""
    vals, vecs = eigh_sorted(hamiltonian(n_sites, m, w, J))
    return float(vals[0]), vecs[:, 0]


def spectrum(n_sites: int, m: float, w: float, J: float) -> np.ndarray:
    vals, _ = eigh_sorted(hamiltonian(n_sites, m, w, J))
    return vals


def charge_sectors(n_sites: int) -> np.ndarray:
    """Total charge eigenvalue of each computational basis state."""
    return np.round(np.real(np.diag(total_charge_op(n_sites)))).astype(int)


def format_terms(terms: list[tuple[str, complex]], m: str = "m", w: str = "w", J: str = "J") -> str:
    """Render a term list as a printable table."""
    lines = [f"{'Pauli string':<12} {'coefficient':>16}", "-" * 30]
    for label, coeff in terms:
        c = coeff.real if abs(coeff.imag) < 1e-14 else coeff
        lines.append(f"{label:<12} {c:>16.6f}")
    lines.append("-" * 30)
    lines.append(f"{len(terms)} terms")
    return "\n".join(lines)


def hamiltonian_from_fermions(n_sites: int, m: float, w: float, J: float) -> np.ndarray:
    """
    Build the same Hamiltonian from the fermionic Kogut-Susskind form, routed
    through Jordan-Wigner, as an independent construction.

    The lattice fermion Hamiltonian, after Gauss's law has been solved but
    before anything is written in spins, is

        H = w sum_n ( c_n^dag c_{n+1} + h.c. )
          + m sum_n (-1)^n c_n^dag c_n
          + J sum_n L_n^2,     L_n = sum_{k<=n} rho_k,
                               rho_k = c_k^dag c_k - (1 - (-1)^k)/2

    Note the charge definition: on an odd site the reference occupation is 1
    (the Dirac sea is full there), on an even site it is 0. Subtracting it is
    what makes the bare vacuum neutral.

    This routine uses the alternating Jordan-Wigner string, prod_{i<j} (-Z_i),
    so that c_n^dag c_{n+1} maps to +sig^+_n sig^-_{n+1} and the printed spin
    Hamiltonian comes out with a plus sign on w. With the plain Z string the
    hopping picks up a minus, which is a relabelling and not physics; see
    jordan_wigner.hopping_sign_demo.

    Agreement between this and `hamiltonian` means the spin form was derived
    correctly from the fermion form.
    """
    from .jordan_wigner import jw_annihilation, jw_number

    dim = 2 ** n_sites
    ham = np.zeros((dim, dim), dtype=complex)
    c = [jw_annihilation(j, n_sites, alternating_string=True) for j in range(n_sites)]
    cd = [op.conj().T for op in c]
    num = [jw_number(j, n_sites) for j in range(n_sites)]
    ident = identity(n_sites)

    for j in range(n_sites - 1):
        hop = cd[j] @ c[j + 1]
        ham += w * (hop + hop.conj().T)

    for j in range(n_sites):
        ham += m * stagger(j) * num[j]

    rho = [num[j] - 0.5 * (1 - stagger(j)) * ident for j in range(n_sites)]
    running = np.zeros((dim, dim), dtype=complex)
    for j in range(n_sites - 1):
        running = running + rho[j]
        ham += J * (running @ running)

    return ham


def fermion_construction_residual(n_sites: int, m: float, w: float, J: float) -> float:
    """
    Difference between the spin-form Hamiltonian and the fermion-form one.

    For even N these agree exactly. For odd N the mass terms differ by the
    constant (m/2) sum_n (-1)^n = -m/2, because the spin form drops a constant
    that only vanishes when N is even. The constant shifts every level equally
    and changes no physics, so it is subtracted here rather than hidden.
    """
    offset = 0.0
    if n_sites % 2 == 1:
        offset = 0.5 * m * sum(stagger(j) for j in range(n_sites))
    a = hamiltonian(n_sites, m, w, J)
    b = hamiltonian_from_fermions(n_sites, m, w, J) - offset * identity(n_sites)
    return float(np.max(np.abs(a - b)))


def site_reversal_op(n_sites: int) -> np.ndarray:
    """Permutation operator for the spatial reflection n -> N + 1 - n."""
    from .paulis import basis_index, basis_label

    dim = 2 ** n_sites
    out = np.zeros((dim, dim), dtype=complex)
    for i in range(dim):
        out[basis_index(basis_label(i, n_sites)[::-1]), i] = 1.0
    return out


def cp_operator(n_sites: int) -> np.ndarray:
    """
    CP: spatial reflection composed with charge conjugation, CP = R . X^(x)N.

    Neither reflection nor charge conjugation is a symmetry on its own. Their
    product is, but only inside the neutral sector, for the same reason the
    electric term is long-range.

    Under CP the charge density maps as

        rho_n  ->  -rho_{N+1-n}

    Reflection moves it to the mirror site and charge conjugation flips its
    sign. (Reflection alone does not even map rho_n to a density, because at
    even N it swaps the staggering; X^(x)N swaps it back.) L_n counts the
    charge to the left of link n only, so summing the mapped densities gives

        L_n  ->  L_{N-n} - Q

    where Q is the total charge. On states with Q = 0 that is just L_{N-n},
    so the electric energy sum_{n=1}^{N-1} L_n^2 is invariant. On states with
    Q != 0 the cross term -2 Q L_{N-n} survives and the symmetry is broken.
    The mass and hopping terms are invariant in every sector. Both maps are
    checked as operator identities in tests/test_schwinger.py.

    Two consequences:
      - On a neutral state that CP leaves unchanged, <L_n> = <L_{N-n}>. The
        bare vacuum is one, so is every state evolved from it by H (H and CP
        commute in that sector), and so is the ground state. At N = 4 the two
        outer links of these states carry equal flux, visible in
        figures/04_pair_production.png. A neutral state without the symmetry
        need not: |0110> has links (0, 0, -1).
      - CP joins total charge as a symmetry that most hardware noise breaks,
        so <L_1> - <L_3> is a second noise diagnostic. Not all noise
        breaks it: logical-frame dephasing respects both (see noise.py).
    """
    from .paulis import pauli_string

    return site_reversal_op(n_sites) @ pauli_string("X" * n_sites)
