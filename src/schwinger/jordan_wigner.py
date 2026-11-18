"""
The Jordan-Wigner transformation, built explicitly and verified algebraically.

The problem it solves
---------------------
Fermion operators on different sites anticommute: {c_i, c_j^dag} = delta_ij.
Spin operators on different sites commute, because they act on different
tensor factors. So "spin up = occupied" cannot be right on its own; it gets
the statistics wrong.

Jordan-Wigner fixes this by attaching a string of Z operators:

    c_j = ( Z_1 Z_2 ... Z_{j-1} ) sigma^-_j

The Z string contributes exactly one factor of (-1) for every occupied site
to the left of j. That is the minus sign fermionic statistics demands, and it
is why the transformation is nonlocal.

Convention used here
--------------------
    sigma^+ = |0><1|,  sigma^- = |1><0|,  Z|0> = +|0>
    n_j = c_j^dag c_j = sigma^+_j sigma^-_j = (I + Z_j)/2

so |0> (spin up, Z = +1) is the occupied state and |1> is empty. This is the
convention that makes the Schwinger charge come out as rho_n = (Z_n + (-1)^n)/2,
matching Martinez et al. (2016) exactly. If you flip it you flip the sign of
every staggered mass term and the vacuum moves.

The hopping sign
----------------
With the plain Z string above,

    c_j^dag c_{j+1} = - sigma^+_j sigma^-_{j+1}

There is a minus sign. It is a convention artifact, removable by the on-site
phase redefinition c_j -> (-1)^j c_j, and it does not change any spectrum.
`hopping_sign_demo` in this module proves both statements numerically. Do not
"fix" it by editing a coefficient until you have run that function.
"""

from __future__ import annotations

import numpy as np

from .paulis import SM, SP, Z, anticommutator, identity, kron_all, site_op, I2


def jw_annihilation(site: int, n_sites: int, alternating_string: bool = False) -> np.ndarray:
    """
    The Jordan-Wigner annihilation operator c_j for `site` (0-indexed).

    Parameters
    ----------
    alternating_string
        If False (default), use the textbook string  prod_{i<j} Z_i.
        If True, use  prod_{i<j} (-Z_i), which is the same transformation
        composed with c_j -> (-1)^j c_j. The alternating version makes
        c_j^dag c_{j+1} = + sigma^+_j sigma^-_{j+1}, matching the hopping term
        as it is usually printed. Both give identical spectra; see
        `hopping_sign_demo`.
    """
    if not 0 <= site < n_sites:
        raise ValueError(f"site {site} out of range for {n_sites} sites")
    string_factor = -Z if alternating_string else Z
    mats = [string_factor] * site + [SM] + [I2] * (n_sites - site - 1)
    return kron_all(mats)


def jw_creation(site: int, n_sites: int, alternating_string: bool = False) -> np.ndarray:
    return jw_annihilation(site, n_sites, alternating_string).conj().T


def jw_number(site: int, n_sites: int) -> np.ndarray:
    """n_j = c_j^dag c_j = (I + Z_j)/2. The string cancels, so it is local."""
    return 0.5 * (identity(n_sites) + site_op("Z", site, n_sites))


def anticommutation_residuals(n_sites: int, alternating_string: bool = False) -> dict[str, float]:
    """
    Check the full fermionic algebra to machine precision.

    Returns the largest absolute deviation for each of the three relations:
        {c_i, c_j^dag} = delta_ij * I
        {c_i, c_j}     = 0
        {c_i^dag, c_j^dag} = 0

    This is the check that catches a broken Z string directly, rather than
    letting it surface later as a wrong energy that could have a dozen causes.
    """
    c = [jw_annihilation(j, n_sites, alternating_string) for j in range(n_sites)]
    cd = [op.conj().T for op in c]
    ident = identity(n_sites)

    r_mixed = 0.0
    r_cc = 0.0
    r_dd = 0.0
    for i in range(n_sites):
        for j in range(n_sites):
            target = ident if i == j else 0.0 * ident
            r_mixed = max(r_mixed, float(np.max(np.abs(anticommutator(c[i], cd[j]) - target))))
            r_cc = max(r_cc, float(np.max(np.abs(anticommutator(c[i], c[j])))))
            r_dd = max(r_dd, float(np.max(np.abs(anticommutator(cd[i], cd[j])))))
    return {
        "{c_i, c_j^dag} - delta_ij": r_mixed,
        "{c_i, c_j}": r_cc,
        "{c_i^dag, c_j^dag}": r_dd,
    }


def pauli_exclusion_residual(n_sites: int) -> float:
    """c_j^dag c_j^dag must be exactly zero. Pauli exclusion as an algebraic identity."""
    worst = 0.0
    for j in range(n_sites):
        cd = jw_creation(j, n_sites)
        worst = max(worst, float(np.max(np.abs(cd @ cd))))
    return worst


def hopping_sign_demo(n_sites: int = 4) -> dict[str, float]:
    """
    Demonstrate the hopping sign convention numerically instead of arguing about it.

    Returns three residuals:

      "plain string: c_j^dag c_{j+1} + sigma^+_j sigma^-_{j+1}"
          should be 0, i.e. the plain Z string gives a minus sign.

      "alternating string: c_j^dag c_{j+1} - sigma^+_j sigma^-_{j+1}"
          should be 0, i.e. the (-Z) string gives a plus sign.

      "spectra agree"
          largest difference between the full eigenvalue spectra of the
          nearest-neighbour hopping Hamiltonian built both ways. Zero, because
          the sign is a relabelling and not physics.
    """
    from .paulis import eigh_sorted, two_site_op

    plain = 0.0
    alt = 0.0
    for j in range(n_sites - 1):
        target = two_site_op("+", j, "-", j + 1, n_sites)
        cd_c_plain = jw_creation(j, n_sites) @ jw_annihilation(j + 1, n_sites)
        cd_c_alt = jw_creation(j, n_sites, True) @ jw_annihilation(j + 1, n_sites, True)
        plain = max(plain, float(np.max(np.abs(cd_c_plain + target))))
        alt = max(alt, float(np.max(np.abs(cd_c_alt - target))))

    def hop(alternating: bool) -> np.ndarray:
        dim = 2 ** n_sites
        h = np.zeros((dim, dim), dtype=complex)
        for j in range(n_sites - 1):
            t = jw_creation(j, n_sites, alternating) @ jw_annihilation(j + 1, n_sites, alternating)
            h += t + t.conj().T
        return h

    ev_plain, _ = eigh_sorted(hop(False))
    ev_alt, _ = eigh_sorted(hop(True))

    return {
        "plain string gives minus sign": plain,
        "alternating string gives plus sign": alt,
        "spectra agree": float(np.max(np.abs(ev_plain - ev_alt))),
    }


def quadratic_from_AB(
    A: np.ndarray, B: np.ndarray, n_sites: int, alternating_string: bool = False
) -> np.ndarray:
    """
    Build  H = sum_ij A_ij c_i^dag c_j + (1/2) sum_ij ( B_ij c_i^dag c_j^dag + h.c. )
    as an explicit matrix, using the Jordan-Wigner operators.

    A must be Hermitian, B antisymmetric. Used to check that a hand-derived
    (A, B) pair really does reproduce a given spin Hamiltonian, which is the
    step where sign errors normally hide.
    """
    if not np.allclose(A, A.conj().T):
        raise ValueError("A must be Hermitian")
    if not np.allclose(B, -B.T):
        raise ValueError("B must be antisymmetric")
    c = [jw_annihilation(j, n_sites, alternating_string) for j in range(n_sites)]
    cd = [op.conj().T for op in c]
    dim = 2 ** n_sites
    h = np.zeros((dim, dim), dtype=complex)
    for i in range(n_sites):
        for j in range(n_sites):
            if A[i, j] != 0:
                h += A[i, j] * (cd[i] @ c[j])
            if B[i, j] != 0:
                pair = B[i, j] * (cd[i] @ cd[j])
                h += 0.5 * (pair + pair.conj().T)
    return h


def bogoliubov_spectrum(A: np.ndarray, B: np.ndarray) -> tuple[np.ndarray, float]:
    """
    Diagonalise a quadratic fermion Hamiltonian by the Lieb-Schultz-Mattis route.

    For  H = sum A_ij c_i^dag c_j + (1/2) sum ( B_ij c_i^dag c_j^dag + h.c. )
    with A real symmetric and B real antisymmetric, there exist Bogoliubov
    quasiparticles eta_k with

        H = sum_k Lambda_k ( eta_k^dag eta_k - 1/2 ) + (1/2) Tr A

    and the Lambda_k^2 are the eigenvalues of (A - B)(A + B).

    Numerics. The textbook route is to form the product (A - B)(A + B), get its
    eigenvalues, and take square roots. That is a bad idea in floating point.
    Because A is symmetric and B antisymmetric, (A + B)^T = A - B, so the
    product is exactly M^T M with M = A + B. Forming M^T M squares the
    condition number, and the square root then amplifies the relative error on
    the small eigenvalues, which are precisely the ones that matter in the
    ordered phase where modes become near-degenerate.

    Since Lambda_k are by definition the singular values of M, take them
    directly from an SVD and never form the product at all. Measured effect on
    the N = 8 open TFIM at h/J = 0.05: worst spectrum residual improves from
    about 7e-10 to about 1e-14.

    Returns
    -------
    (Lambda, offset)
        Lambda: non-negative single-particle energies, ascending.
        offset: (1/2) Tr A - (1/2) sum_k Lambda_k, the ground energy of the
        quadratic part.
    """
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    if not np.allclose(A, A.T, atol=1e-12):
        raise ValueError("A must be real symmetric")
    if not np.allclose(B, -B.T, atol=1e-12):
        raise ValueError("B must be real antisymmetric")
    lam = np.sort(np.linalg.svd(A + B, compute_uv=False))
    offset = 0.5 * float(np.trace(A)) - 0.5 * float(np.sum(lam))
    return lam, offset


def many_body_spectrum_from_modes(lam: np.ndarray, offset: float) -> np.ndarray:
    """
    Rebuild the full 2^N many-body spectrum from single-particle energies.

    Every subset of occupied quasiparticle modes gives one eigenvalue:
        E = offset + sum_k n_k Lambda_k,  n_k in {0, 1}
    """
    n = len(lam)
    energies = np.empty(2 ** n)
    for occ in range(2 ** n):
        bits = np.array([(occ >> k) & 1 for k in range(n)], dtype=float)
        energies[occ] = offset + float(bits @ lam)
    return np.sort(energies)
