"""
The transverse-field Ising model: the training ground.

    H = -J sum_{n} X_n X_{n+1}  -  h sum_n Z_n

Two reasons this model is here, both about safety rather than content:

 1. It is the simplest nontrivial spin Hamiltonian, so it is where the
    Pauli-string machinery gets debugged.
 2. It is exactly solvable. Jordan-Wigner maps it to free fermions, which
    diagonalise by Bogoliubov transformation. So the transformation can be
    learned somewhere an analytic answer exists, and verified, BEFORE it is
    pointed at the Schwinger model where no analytic answer exists.

The model has a quantum phase transition at h/J = 1 in the thermodynamic
limit. At N = 4 it is a smooth crossover, not a sharp transition, and the
code says so rather than pretending otherwise.

Note on convention: this file uses the XX + Z form, not ZZ + X. The two are
related by a global Hadamard rotation and describe the same physics, but only
the XX + Z form Jordan-Wigner transforms into the standard (A, B) shape
without an extra basis change, so it is the one to derive on paper.
"""

from __future__ import annotations

import numpy as np

from .paulis import identity, site_op, two_site_op


def tfim_terms(n_sites: int, J: float, h: float, periodic: bool = False) -> list[tuple[str, complex]]:
    """The TFIM as a printed list of (Pauli string, coefficient)."""
    terms: list[tuple[str, complex]] = []
    bonds = n_sites if periodic else n_sites - 1
    for n in range(bonds):
        m = (n + 1) % n_sites
        chars = ["I"] * n_sites
        chars[n] = "X"
        chars[m] = "X"
        terms.append(("".join(chars), -J))
    for n in range(n_sites):
        chars = ["I"] * n_sites
        chars[n] = "Z"
        terms.append(("".join(chars), -h))
    return terms


def tfim_hamiltonian(n_sites: int, J: float, h: float, periodic: bool = False) -> np.ndarray:
    """Explicit 2^N x 2^N matrix for the TFIM."""
    dim = 2 ** n_sites
    ham = np.zeros((dim, dim), dtype=complex)
    bonds = n_sites if periodic else n_sites - 1
    for n in range(bonds):
        m = (n + 1) % n_sites
        ham += -J * two_site_op("X", n, "X", m, n_sites)
    for n in range(n_sites):
        ham += -h * site_op("Z", n, n_sites)
    return ham


def tfim_AB(n_sites: int, J: float, h: float) -> tuple[np.ndarray, np.ndarray, float]:
    """
    The (A, B) matrices and constant for the open-boundary TFIM after Jordan-Wigner.

    Derivation, which you should reproduce on paper:

      X_n X_{n+1} = -( c_n^dag - c_n )( c_{n+1}^dag + c_{n+1} )
      Z_n         = 2 c_n^dag c_n - 1

    so

      H = J sum ( c_n^dag c_{n+1} + c_{n+1}^dag c_n )
        + J sum ( c_n^dag c_{n+1}^dag - c_n c_{n+1} )
        - 2h sum c_n^dag c_n
        + h N

    Matching against  H = sum A_ij c_i^dag c_j + (1/2) sum ( B_ij c_i^dag c_j^dag + h.c. ) + const:

      A_nn = -2h,  A_{n,n+1} = A_{n+1,n} = J
      B_{n,n+1} = J,  B_{n+1,n} = -J
      const = h N

    Returns (A, B, const). `verify_tfim_jordan_wigner` checks this numerically,
    so you never have to take the algebra above on trust.
    """
    A = np.zeros((n_sites, n_sites))
    B = np.zeros((n_sites, n_sites))
    for n in range(n_sites):
        A[n, n] = -2.0 * h
    for n in range(n_sites - 1):
        A[n, n + 1] = J
        A[n + 1, n] = J
        B[n, n + 1] = J
        B[n + 1, n] = -J
    return A, B, h * n_sites


def verify_tfim_jordan_wigner(n_sites: int, J: float, h: float) -> float:
    """
    Check that the hand-derived (A, B) really do rebuild the spin Hamiltonian.

    This is the tightest possible test of the derivation: it compares two
    2^N x 2^N matrices entry by entry, so a single wrong sign anywhere shows
    up immediately and localised, instead of as a wrong energy later.

    Returns the largest absolute entry of the difference.
    """
    from .jordan_wigner import quadratic_from_AB

    A, B, const = tfim_AB(n_sites, J, h)
    rebuilt = quadratic_from_AB(A, B, n_sites) + const * identity(n_sites)
    return float(np.max(np.abs(rebuilt - tfim_hamiltonian(n_sites, J, h))))


def tfim_exact_spectrum(n_sites: int, J: float, h: float) -> np.ndarray:
    """
    Full many-body spectrum of the open-boundary TFIM from free fermions.

    Route: (A, B) -> Bogoliubov single-particle energies Lambda_k -> every
    subset of occupied modes. This never touches the 2^N Hilbert space
    except to enumerate occupations, so agreeing with brute-force
    diagonalisation is a genuine check and not a tautology.
    """
    from .jordan_wigner import bogoliubov_spectrum, many_body_spectrum_from_modes

    A, B, const = tfim_AB(n_sites, J, h)
    lam, offset = bogoliubov_spectrum(A, B)
    return many_body_spectrum_from_modes(lam, offset + const)


def tfim_periodic_spectrum(n_sites: int, J: float, h: float) -> np.ndarray:
    """
    Full many-body spectrum of the PERIODIC TFIM in closed form, done properly.

    This is the classic Jordan-Wigner trap, and it is worth doing once because
    getting it right is what separates owning the technique from owning the
    formula.

    Trap 1: the wrapped string. On a ring the JW string closes on itself, and
    the boundary bond X_N X_1 picks up a factor of -P with P = prod_n Z_n the
    total fermion parity. So this is not one fermion problem but two:

        P = +1 (even parity): ANTIperiodic fermions, k = (2j+1) pi / N
        P = -1 (odd parity):  periodic fermions,     k = 2 pi j / N

    and in each sector only the states with matching parity are physical.

    Trap 2, the one that actually bites. The textbook dispersion

        eps(k) = 2 sqrt( J^2 + h^2 - 2 J h cos k )

    comes from a Bogoliubov rotation that mixes the modes k and -k. At k = 0
    and k = pi we have sin k = 0, the pairing term vanishes, and those modes
    are ALREADY diagonal. Their energy is the SIGNED

        eps(k) = 2 ( h - J cos k )

    which is negative at k = 0 whenever h < J. Taking the positive square root
    there silently flips that one mode, shifting the whole ordered-phase
    spectrum by exactly |2(h - J)|. The disordered phase h > J looks fine,
    which is what makes the bug so easy to ship.

    `scripts/02_jw_verification.py` reproduces both the broken and the fixed
    version side by side.
    """
    energies: list[float] = []
    for parity in (+1, -1):
        if parity == +1:
            ks = np.array([(2 * j + 1) * np.pi / n_sites for j in range(n_sites)])
        else:
            ks = np.array([2 * np.pi * j / n_sites for j in range(n_sites)])

        eps = np.empty(n_sites)
        for i, k in enumerate(ks):
            if abs(np.sin(k)) < 1e-12:
                # Unpaired mode: no Bogoliubov rotation, keep the sign.
                eps[i] = 2.0 * (h - J * np.cos(k))
            else:
                eps[i] = 2.0 * np.sqrt(J ** 2 + h ** 2 - 2 * J * h * np.cos(k))

        base = -0.5 * float(np.sum(eps))
        for occ in range(2 ** n_sites):
            bits = np.array([(occ >> m) & 1 for m in range(n_sites)])
            if (-1) ** int(bits.sum()) != parity:
                continue
            energies.append(base + float(bits @ eps))

    out = np.sort(np.array(energies))
    if len(out) != 2 ** n_sites:
        raise ValueError(f"kept {len(out)} states, expected {2 ** n_sites}")
    return out


def tfim_periodic_spectrum_naive(n_sites: int, J: float, h: float) -> np.ndarray:
    """
    The WRONG version, kept deliberately so the failure can be shown, not asserted.

    Identical to `tfim_periodic_spectrum` except that it uses the positive
    square root for every k, including the unpaired k = 0 and k = pi modes.
    Agrees with exact diagonalisation for h > J and is off by exactly
    |2(h - J)| for h < J.
    """
    energies: list[float] = []
    for parity in (+1, -1):
        if parity == +1:
            ks = np.array([(2 * j + 1) * np.pi / n_sites for j in range(n_sites)])
        else:
            ks = np.array([2 * np.pi * j / n_sites for j in range(n_sites)])
        eps = 2.0 * np.sqrt(J ** 2 + h ** 2 - 2 * J * h * np.cos(ks))
        base = -0.5 * float(np.sum(eps))
        for occ in range(2 ** n_sites):
            bits = np.array([(occ >> m) & 1 for m in range(n_sites)])
            if (-1) ** int(bits.sum()) != parity:
                continue
            energies.append(base + float(bits @ eps))
    return np.sort(np.array(energies))


def tfim_ground_energy_scan(n_sites: int, J: float, h_over_J: np.ndarray) -> np.ndarray:
    """Ground-state energy per site as h/J is swept. Used for the crossover figure."""
    from .paulis import eigh_sorted

    out = np.empty(len(h_over_J))
    for i, r in enumerate(h_over_J):
        vals, _ = eigh_sorted(tfim_hamiltonian(n_sites, J, r * J))
        out[i] = vals[0] / n_sites
    return out
