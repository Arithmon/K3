#!/usr/bin/env python3
"""
invariant_quotient_ring.py — base standard C[Z]/⟨Q_0,Q_1,Q_2⟩ Hermitiens
                             invariants Z_2^3.

Design (GPT 5.6 + Brieuc 07-12) : la route quotient devient obligatoire
after the sampler sweep found r_eff = 216 strictly constant against
N ∈ {500,1000,2000} × 3 seeds, tol 1e-10) et le frozen manifest test V_≤3→V_≤4
(nothing frozen; the min-max was expected, but on ill-defined spaces).

This module produces a quotient analysis basis COMPATIBLE with
signature {"ij", "kl", "type"} du moteur spectral (basis_values,
basis_chart_derivs: only the starting list of monomials changes
retient uniquement les holomorphiques STANDARD hors ⟨Z_p²⟩_{p∈PIVOT}).

Base standard :
  monome Z^I holo de degré d avec I = (i_1, ..., i_d) trié
  ↔ multiplicité α_i = # {k : i_k = i} ∈ {0, 1} pour i ∈ PIVOT.
    In other words each pivot appears at most once in I.

Prédiction : dim quotient invariant Hermitien (Z_2^3, bidegré (d,d))
  d=0: 1  d=1: 10  d=2: 58  d=3: 218  d=4: 610  d=5: 1402
Cumul V_≤d^{G-inv,Herm} sans constante :
  V_≤1: 10  V_≤2: 68  V_≤3: 286  V_≤4: 896  V_≤5: 2298

Usage :
    from .invariant_quotient_ring import (
        basis_at_deg_quotient, constant_coordinates_exact,
        inclusion_matrix_exact,
    )
    basis = basis_at_deg_quotient(d)
    inclusion = inclusion_matrix_exact(d)   # shape (dim V_{d+1}, dim V_d)
    c_next = inclusion.matvec_exact(c_d)
    c_constant = constant_coordinates_exact(d)
"""
from __future__ import annotations

import io
import sys
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations_with_replacement
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

from .metric_oracle import COORD_CHARS, MU, enumerate_metric_basis

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

PIVOT = (0, 1, 2)                                # S — figé (choix stable)
NONPIVOT = tuple(i for i in range(6) if i not in PIVOT)


def is_standard_holo(I: Tuple[int, ...], pivot=PIVOT) -> bool:
    """True ssi le monôme Z^I (I trié) a chaque pivot ≤ 1 fois.

    Équivalent : pour tout p ∈ PIVOT, I.count(p) ≤ 1.
    Complète Z^I = Π Z_{I[k]} pour I ∈ combinations_with_replacement."""
    for p in pivot:
        if I.count(p) > 1:
            return False
    return True


def enumerate_standard_holo(d: int, pivot=PIVOT) -> List[Tuple[int, ...]]:
    """List of sorted multi-indices I of length d with alpha_p in {0, 1} for p in pivot.

    Format identique à combinations_with_replacement(range(6), d), filtré."""
    return [I for I in combinations_with_replacement(range(6), d)
            if is_standard_holo(I, pivot)]


def char_of(I: Tuple[int, ...]) -> Tuple[int, int, int]:
    """Caractère Z_2^3 du monôme Z^I (I trié) : Σ COORD_CHARS[i] mod 2."""
    c = np.zeros(3, dtype=int)
    for i in I:
        c += COORD_CHARS[i]
    return tuple((c % 2).tolist())


def enumerate_sector_quotient(d: int, pivot=PIVOT) -> List[dict]:
    """Enumerate the Z_2^3 invariant Hermitian elements of bidegree (d,d) on the
    base STANDARD hors ⟨Z_p²⟩. Format identique à enumerate_sector du moteur.

    Each element is {"type": ..., "ij": I, "kl": K} with type in
    {"self", "real_pair", "imag_pair"} et char(I) == char(K) (invariance).
    Ordering convention I < K to avoid double counting, and a pair
    (real_pair, imag_pair) chacune contribuant à la dim réelle."""
    holo = enumerate_standard_holo(d, pivot)
    out = []
    for I in holo:
        for K in holo:
            if char_of(I) != char_of(K):
                continue
            if I < K:
                out.append({"type": "real_pair", "ij": I, "kl": K})
                out.append({"type": "imag_pair", "ij": I, "kl": K})
            elif I == K:
                out.append({"type": "self", "ij": I, "kl": K})
    return out


def basis_upto_quotient(d_max: int, pivot=PIVOT) -> List[dict]:
    """DEPRECATED — voir basis_at_deg_quotient. Cumul (direct sum) : conserve
    the semantics of basis_upto (raw), but REDUNDANT on the surface when V_e is
    included through multiplication by s^{d-e} in V_d. To be used only for
    comparison with the raw basis."""
    b = []
    for dg in range(1, d_max + 1):
        b += enumerate_sector_quotient(dg, pivot)
    return b


def basis_at_deg_quotient(d: int, pivot=PIVOT) -> List[dict]:
    """Base réelle Hermitienne invariante Z_2^3 au SEUL bidegré (d, d), mod
    ⟨Q_0, Q_1, Q_2⟩. Dim = Σ_χ m²_{d,χ} = 10, 58, 218, 610, 1402 pour d=1..5.

    Holds because V_{<=d} is V_d: functions of bidegree at most d are
    naturally included in bidegree (d, d) through multiplication by s^k
    au numérateur ET dénominateur ; s > 0 partout). Élimine la redondance
    de basis_upto_quotient (68 fn superflues à V_≤3, 286 à V_≤4)."""
    return enumerate_sector_quotient(d, pivot)


# ================================================================
#  Inclusions exactes V_d -> V_{d+1} par multiplication par s
# ================================================================
@dataclass(frozen=True)
class ExactSparseMatrix:
    """Matrice rationnelle creuse, stockee par colonnes.

    ``entries_by_column[j]`` contient les couples ``(i, a_ij)`` tries par
    row, with ``a_ij`` a non-zero :class:`fractions.Fraction`.  This
    representation est volontairement independante de SciPy : elle reste
    exacte lors de la construction et de la serialisation, tout en offrant
    une conversion CSC explicite au consommateur numerique.
    """

    nrows: int
    ncols: int
    entries_by_column: Tuple[Tuple[Tuple[int, Fraction], ...], ...]

    def __post_init__(self):
        if len(self.entries_by_column) != self.ncols:
            raise ValueError("one sparse column is required per matrix column")
        for column in self.entries_by_column:
            rows = [row for row, value in column]
            if rows != sorted(set(rows)):
                raise ValueError("sparse row indices must be sorted and unique")
            if any(row < 0 or row >= self.nrows for row in rows):
                raise ValueError("sparse row index outside matrix shape")
            if any(not isinstance(value, Fraction) or value == 0
                   for _, value in column):
                raise ValueError("sparse entries must be nonzero Fractions")

    @property
    def shape(self) -> Tuple[int, int]:
        return self.nrows, self.ncols

    @property
    def nnz(self) -> int:
        return sum(len(column) for column in self.entries_by_column)

    def column(self, j: int) -> Dict[int, Fraction]:
        """Return column ``j`` as a fresh row-to-rational dictionary."""
        return dict(self.entries_by_column[j])

    def matvec_exact(self, vector: Sequence) -> Tuple[Fraction, ...]:
        """Compute ``self @ vector`` over Q, without float conversion."""
        if len(vector) != self.ncols:
            raise ValueError(f"expected vector of length {self.ncols}")
        out = [Fraction(0) for _ in range(self.nrows)]
        for j, x in enumerate(vector):
            xq = x if isinstance(x, Fraction) else Fraction(x)
            if xq == 0:
                continue
            for i, aij in self.entries_by_column[j]:
                out[i] += aij * xq
        return tuple(out)

    def to_scipy_csc(self, dtype=np.float64):
        """Convert explicitly to a numeric SciPy CSC matrix.

        The exact object is kept as the source of truth; this method is only
        the boundary used by QR/SVD feature-level code.
        """
        from scipy.sparse import csc_matrix

        data, indices, indptr = [], [], [0]
        for column in self.entries_by_column:
            for row, value in column:
                indices.append(row)
                data.append(value)
            indptr.append(len(data))
        return csc_matrix((np.asarray(data, dtype=dtype),
                           np.asarray(indices, dtype=np.int64),
                           np.asarray(indptr, dtype=np.int64)),
                          shape=self.shape)

    def to_jsonable(self) -> dict:
        """Lossless JSON payload using ``[row, numerator, denominator]``."""
        return {
            "shape": [self.nrows, self.ncols],
            "orientation": "target_rows_by_source_columns",
            "entries_by_column": [
                [[row, value.numerator, value.denominator]
                 for row, value in column]
                for column in self.entries_by_column
            ],
        }


def basis_character(element: dict) -> Tuple[int, int, int]:
    """Common holomorphic character of an invariant Hermitian element."""
    char_ij = char_of(element["ij"])
    char_kl = char_of(element["kl"])
    if char_ij != char_kl:
        raise ValueError("basis element is not Z_2^3 invariant")
    return char_ij


def _nonpivot(pivot: Tuple[int, ...]) -> Tuple[int, ...]:
    return tuple(i for i in range(6) if i not in pivot)


def _mu_as_fraction(i: int) -> Fraction:
    """Read the frozen integral MU convention without inheriting float error."""
    value = float(MU[i])
    if not value.is_integer():
        raise ValueError("exact quotient reduction requires integral MU entries")
    return Fraction(int(value))


def _rat_matinv3(matrix: Sequence[Sequence[Fraction]]):
    """Exact inverse of a 3x3 rational matrix via its adjugate."""
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    if det == 0:
        raise ValueError("singular Vandermonde pivot")
    return (
        ((e * i - f * h) / det, -(b * i - c * h) / det,
         (b * f - c * e) / det),
        (-(d * i - f * g) / det, (a * i - c * g) / det,
         -(a * f - c * d) / det),
        ((d * h - e * g) / det, -(a * h - b * g) / det,
         (a * e - b * d) / det),
    )


@lru_cache(maxsize=None)
def vandermonde_reduction_coeffs_exact(
        pivot: Tuple[int, ...] = PIVOT) -> Tuple[Tuple[Fraction, ...], ...]:
    """Return exact ``Z_p^2 = sum_t A[p,t] Z_t^2`` coefficients.

    The rows are ordered by ``pivot`` and the columns by its complement.
    For the frozen pivot ``(0,1,2)`` all entries happen to be integers, but
    the API retains :class:`Fraction` so another valid pivot stays exact.
    """
    pivot = tuple(pivot)
    nonpivot = _nonpivot(pivot)
    if len(pivot) != 3 or len(set(pivot)) != 3:
        raise ValueError("the CI(2,2,2) reduction requires three pivots")
    mu_s = tuple(_mu_as_fraction(i) for i in pivot)
    mu_t = tuple(_mu_as_fraction(i) for i in nonpivot)
    vand_s = ((Fraction(1),) * 3, mu_s, tuple(x * x for x in mu_s))
    vand_t = ((Fraction(1),) * 3, mu_t, tuple(x * x for x in mu_t))
    vand_s_inv = _rat_matinv3(vand_s)
    return tuple(tuple(sum(-vand_s_inv[i][k] * vand_t[k][j]
                           for k in range(3))
                       for j in range(3))
                 for i in range(3))


def reduce_holo_times_coordinate_exact(
        monomial: Tuple[int, ...], coordinate: int,
        pivot: Tuple[int, ...] = PIVOT) -> Dict[Tuple[int, ...], Fraction]:
    """Reduce ``Z_coordinate * Z^monomial`` in the standard quotient basis.

    ``monomial`` must already be standard.  A one-step multiplication can
    create only one forbidden square ``Z_p^2``.  Replacing that square with
    the exact Vandermonde relation therefore completes the reduction without
    a numerical Groebner calculation or a truncation.
    """
    pivot = tuple(pivot)
    monomial = tuple(monomial)
    if coordinate not in range(6):
        raise ValueError("coordinate must lie in range(6)")
    if tuple(sorted(monomial)) != monomial or not is_standard_holo(monomial, pivot):
        raise ValueError("input monomial must be sorted and standard")
    product = tuple(sorted(monomial + (coordinate,)))
    if is_standard_holo(product, pivot):
        return {product: Fraction(1)}

    if coordinate not in pivot or product.count(coordinate) != 2:
        raise AssertionError("one-step standard reduction found an unexpected term")
    remainder = list(product)
    remainder.remove(coordinate)
    remainder.remove(coordinate)
    row = pivot.index(coordinate)
    out = {}
    for target_coordinate, coefficient in zip(
            _nonpivot(pivot), vandermonde_reduction_coeffs_exact(pivot)[row]):
        if coefficient == 0:
            continue
        target = tuple(sorted(tuple(remainder)
                              + (target_coordinate, target_coordinate)))
        if not is_standard_holo(target, pivot):
            raise AssertionError("Vandermonde reduction did not yield a standard monomial")
        out[target] = out.get(target, Fraction(0)) + coefficient
    return {target: value for target, value in out.items() if value}


@lru_cache(maxsize=None)
def reduce_holomorphic_monomial_exact(
        monomial: Tuple[int, ...],
        pivot: Tuple[int, ...] = PIVOT) -> Tuple[Tuple[Tuple[int, ...], Fraction], ...]:
    """Reduce an arbitrary monomial to the standard quotient basis over Q."""
    pivot = tuple(pivot)
    monomial = tuple(sorted(monomial))
    offending = next((coordinate for coordinate in pivot
                      if monomial.count(coordinate) > 1), None)
    if offending is None:
        return ((monomial, Fraction(1)),)
    remainder = list(monomial)
    remainder.remove(offending)
    remainder.remove(offending)
    row = pivot.index(offending)
    accumulated: Dict[Tuple[int, ...], Fraction] = {}
    for target_coordinate, coefficient in zip(
            _nonpivot(pivot), vandermonde_reduction_coeffs_exact(pivot)[row]):
        if coefficient == 0:
            continue
        substituted = tuple(sorted(
            tuple(remainder) + (target_coordinate, target_coordinate)))
        for target, reduced_coefficient in reduce_holomorphic_monomial_exact(
                substituted, pivot):
            accumulated[target] = (accumulated.get(target, Fraction(0))
                                   + coefficient * reduced_coefficient)
    return tuple(sorted((target, value) for target, value in accumulated.items()
                        if value))


@lru_cache(maxsize=None)
def raw_to_quotient_matrix_exact(
        d: int, pivot: Tuple[int, ...] = PIVOT) -> ExactSparseMatrix:
    """Exact coordinate map from the raw invariant degree-d basis to V_d.

    The source ordering is :func:`enumerate_metric_basis`; the target ordering
    is :func:`basis_at_deg_quotient`.  Thus archived raw coefficients can be
    audited without a numerical projection or a Gram cutoff.
    """
    pivot = tuple(pivot)
    source_basis = enumerate_metric_basis(d)
    target_holo = enumerate_standard_holo(d, pivot)
    target_index = {item: index for index, item in enumerate(target_holo)}
    target_basis = basis_at_deg_quotient(d, pivot)
    target_basis_index = {
        (item["type"], item["ij"], item["kl"]): index
        for index, item in enumerate(target_basis)
    }
    columns = []
    for element in source_basis:
        typ = element["type"]
        left = {target_index[item]: coefficient for item, coefficient
                in reduce_holomorphic_monomial_exact(element["ij"], pivot)}
        right = {target_index[item]: coefficient for item, coefficient
                 in reduce_holomorphic_monomial_exact(element["kl"], pivot)}
        support = sorted(set(left) | set(right))
        column: Dict[int, Fraction] = {}
        if typ == "self":
            for position, a in enumerate(support):
                ca = left.get(a, Fraction(0))
                if ca == 0:
                    continue
                for b in support[position:]:
                    cb = left.get(b, Fraction(0))
                    if cb == 0:
                        continue
                    target_type = "self" if a == b else "real_pair"
                    key = (target_type, target_holo[a], target_holo[b])
                    _add_sparse_entry(column, target_basis_index[key], ca * cb)
        elif typ == "real_pair":
            for position, a in enumerate(support):
                ua, va = left.get(a, Fraction(0)), right.get(a, Fraction(0))
                if 2 * ua * va:
                    key = ("self", target_holo[a], target_holo[a])
                    _add_sparse_entry(column, target_basis_index[key], 2 * ua * va)
                for b in support[position + 1:]:
                    value = (ua * right.get(b, Fraction(0))
                             + va * left.get(b, Fraction(0)))
                    if value:
                        key = ("real_pair", target_holo[a], target_holo[b])
                        _add_sparse_entry(column, target_basis_index[key], value)
        elif typ == "imag_pair":
            for position, a in enumerate(support):
                ua, va = left.get(a, Fraction(0)), right.get(a, Fraction(0))
                for b in support[position + 1:]:
                    value = (ua * right.get(b, Fraction(0))
                             - va * left.get(b, Fraction(0)))
                    if value:
                        key = ("imag_pair", target_holo[a], target_holo[b])
                        _add_sparse_entry(column, target_basis_index[key], value)
        else:
            raise ValueError(f"unknown Hermitian basis type: {typ}")
        columns.append(tuple(sorted(column.items())))
    return ExactSparseMatrix(len(target_basis), len(source_basis), tuple(columns))


@lru_cache(maxsize=None)
def holomorphic_multiplication_matrix_exact(
        d: int, coordinate: int,
        pivot: Tuple[int, ...] = PIVOT) -> ExactSparseMatrix:
    """Exact matrix of ``Z_coordinate : R_d -> R_{d+1}`` in standard bases."""
    pivot = tuple(pivot)
    if d < 0:
        raise ValueError("degree must be nonnegative")
    source = enumerate_standard_holo(d, pivot)
    target = enumerate_standard_holo(d + 1, pivot)
    target_index = {monomial: i for i, monomial in enumerate(target)}
    columns = []
    for monomial in source:
        reduced = reduce_holo_times_coordinate_exact(monomial, coordinate, pivot)
        columns.append(tuple(sorted((target_index[item], coefficient)
                                    for item, coefficient in reduced.items()
                                    if coefficient)))
    return ExactSparseMatrix(len(target), len(source), tuple(columns))


def _add_sparse_entry(column: Dict[int, Fraction], row: int, value: Fraction):
    if value:
        column[row] = column.get(row, Fraction(0)) + value
        if column[row] == 0:
            del column[row]


@lru_cache(maxsize=None)
def inclusion_matrix_exact(
        d: int, pivot: Tuple[int, ...] = PIVOT) -> ExactSparseMatrix:
    """Exact inclusion ``V_d -> V_{d+1}``, oriented target x source.

    On Hermitian numerators the map is

    ``H |-> sum_j T_j H T_j^T``,

    where ``T_j`` is exact multiplication by ``Z_j`` in the holomorphic
    quotient.  Hence a source coefficient vector ``c`` is transported as
    ``c_next = inclusion_matrix_exact(d) @ c``.  No Gram matrix, tolerance,
    SVD cutoff or floating-point coefficient enters this construction.
    """
    pivot = tuple(pivot)
    if d < 0:
        raise ValueError("degree must be nonnegative")
    source_holo = enumerate_standard_holo(d, pivot)
    target_holo = enumerate_standard_holo(d + 1, pivot)
    source_holo_index = {item: i for i, item in enumerate(source_holo)}
    source_basis = basis_at_deg_quotient(d, pivot)
    target_basis = basis_at_deg_quotient(d + 1, pivot)
    target_basis_index = {
        (item["type"], item["ij"], item["kl"]): i
        for i, item in enumerate(target_basis)
    }
    multipliers = tuple(holomorphic_multiplication_matrix_exact(d, j, pivot)
                        for j in range(6))

    columns = []
    for element in source_basis:
        typ = element["type"]
        i_source = source_holo_index[element["ij"]]
        k_source = source_holo_index[element["kl"]]
        column = {}
        for multiplier in multipliers:
            u = multiplier.column(i_source)
            v = multiplier.column(k_source)
            support = sorted(set(u) | set(v))
            if typ == "self":
                for pos, a in enumerate(support):
                    ca = u.get(a, Fraction(0))
                    if ca == 0:
                        continue
                    for b in support[pos:]:
                        cb = u.get(b, Fraction(0))
                        if cb == 0:
                            continue
                        if a == b:
                            key = ("self", target_holo[a], target_holo[a])
                        else:
                            key = ("real_pair", target_holo[a], target_holo[b])
                        _add_sparse_entry(column, target_basis_index[key], ca * cb)
            elif typ == "real_pair":
                for pos, a in enumerate(support):
                    ua, va = u.get(a, Fraction(0)), v.get(a, Fraction(0))
                    diagonal = 2 * ua * va
                    if diagonal:
                        key = ("self", target_holo[a], target_holo[a])
                        _add_sparse_entry(column, target_basis_index[key], diagonal)
                    for b in support[pos + 1:]:
                        value = ua * v.get(b, Fraction(0)) + va * u.get(b, Fraction(0))
                        if value:
                            key = ("real_pair", target_holo[a], target_holo[b])
                            _add_sparse_entry(column, target_basis_index[key], value)
            elif typ == "imag_pair":
                for pos, a in enumerate(support):
                    ua, va = u.get(a, Fraction(0)), v.get(a, Fraction(0))
                    for b in support[pos + 1:]:
                        value = ua * v.get(b, Fraction(0)) - va * u.get(b, Fraction(0))
                        if value:
                            key = ("imag_pair", target_holo[a], target_holo[b])
                            _add_sparse_entry(column, target_basis_index[key], value)
            else:
                raise ValueError(f"unknown Hermitian basis type: {typ}")
        columns.append(tuple(sorted(column.items())))
    return ExactSparseMatrix(len(target_basis), len(source_basis), tuple(columns))


@lru_cache(maxsize=None)
def constant_coordinates_exact(
        d: int, pivot: Tuple[int, ...] = PIVOT) -> Tuple[Fraction, ...]:
    """Coordinates of the constant function ``1 = s^d/s^d`` in ``V_d``.

    The degree-zero basis is the singleton ``|1|^2``.  Repeated exact
    inclusions therefore produce the reduced numerator ``s^d`` and provide
    the structural constant direction required by the spectral pipeline.
    """
    pivot = tuple(pivot)
    if d < 0:
        raise ValueError("degree must be nonnegative")
    coordinates = (Fraction(1),)
    for degree in range(d):
        coordinates = inclusion_matrix_exact(degree, pivot).matvec_exact(coordinates)
    return coordinates


def exact_rank_mod_prime(matrix: ExactSparseMatrix,
                         prime: int = 2_147_483_647) -> int:
    """Compute sparse column rank over F_prime.

    Full column rank modulo a prime (not dividing any denominator) proves
    full column rank over Q.  This is used only by the self-check; it never
    defines or truncates the analysis space.
    """
    pivots = {}
    for column in matrix.entries_by_column:
        vector = {}
        for row, value in column:
            denominator = value.denominator % prime
            if denominator == 0:
                raise ValueError("rank prime divides a matrix denominator")
            residue = (value.numerator % prime) * pow(denominator, -1, prime) % prime
            if residue:
                vector[row] = residue
        while vector:
            pivot_row = min(vector)
            pivot_column = pivots.get(pivot_row)
            if pivot_column is None:
                scale = pow(vector[pivot_row], -1, prime)
                vector = {row: value * scale % prime
                          for row, value in vector.items() if value % prime}
                pivots[pivot_row] = vector
                break
            factor = vector[pivot_row]
            for row, value in pivot_column.items():
                reduced = (vector.get(row, 0) - factor * value) % prime
                if reduced:
                    vector[row] = reduced
                else:
                    vector.pop(row, None)
    return len(pivots)


def full_column_rank_proof_prime(
        matrix: ExactSparseMatrix,
        primes: Sequence[int] = (2_147_483_647, 2_147_483_629,
                                 2_147_483_587)) -> Tuple[int, int] | None:
    """Return a prime witnessing full column rank over Q, if one is found.

    Full rank after reduction modulo one admissible prime proves full rank over
    Q.  Rank loss modulo a prime proves nothing over Q, so this helper tries a
    short deterministic list and returns ``None`` instead of making a false
    non-injectivity claim.
    """
    for prime in primes:
        try:
            rank = exact_rank_mod_prime(matrix, prime)
        except ValueError:
            continue
        if rank == matrix.ncols:
            return prime, rank
    return None


def validate_vandermonde_reduction_exact(
        pivot: Tuple[int, ...] = PIVOT) -> None:
    """Check independently that ``V_S A + V_T`` vanishes over Q."""
    pivot = tuple(pivot)
    nonpivot = _nonpivot(pivot)
    mu_s = tuple(_mu_as_fraction(i) for i in pivot)
    mu_t = tuple(_mu_as_fraction(i) for i in nonpivot)
    vand_s = ((Fraction(1),) * 3, mu_s, tuple(x * x for x in mu_s))
    vand_t = ((Fraction(1),) * 3, mu_t, tuple(x * x for x in mu_t))
    reduction = vandermonde_reduction_coeffs_exact(pivot)
    for equation in range(3):
        for target in range(3):
            residual = (sum(vand_s[equation][source]
                            * reduction[source][target]
                            for source in range(3))
                        + vand_t[equation][target])
            if residual != 0:
                raise AssertionError("V_S A + V_T is not exactly zero")


def _evaluate_basis_on_points(basis: Sequence[dict], points: np.ndarray) -> np.ndarray:
    """Small independent evaluator used by the non-persistent self-check."""
    out = np.empty((len(points), len(basis)), dtype=np.longdouble)
    for column, element in enumerate(basis):
        m_i = np.ones(len(points), dtype=np.clongdouble)
        m_k = np.ones(len(points), dtype=np.clongdouble)
        for coordinate in element["ij"]:
            m_i *= points[:, coordinate]
        for coordinate in element["kl"]:
            m_k *= points[:, coordinate]
        value = m_i * np.conj(m_k)
        if element["type"] == "self":
            numerator = value.real
        elif element["type"] == "real_pair":
            numerator = 2 * value.real
        else:
            numerator = -2 * value.imag
        degree = len(element["ij"])
        scale = np.sum(np.abs(points) ** 2, axis=1) ** degree
        out[:, column] = numerator / scale
    return out


def _self_check_points(seed: int = 20260712, count: int = 5) -> np.ndarray:
    """Generate deterministic floating evaluations on the exact radical chart."""
    rng = np.random.default_rng(seed)
    nonpivot = _nonpivot(PIVOT)
    reduction = np.asarray(vandermonde_reduction_coeffs_exact(PIVOT),
                           dtype=np.clongdouble)
    points = np.zeros((count, 6), dtype=np.clongdouble)
    free = ((rng.uniform(-0.45, 0.45, (count, 3))
             + 1j * rng.uniform(-0.45, 0.45, (count, 3)))
            .astype(np.clongdouble))
    free[:, 0] += 1
    points[:, list(nonpivot)] = free
    pivot_squares = (reduction @ (free * free).T).T
    points[:, list(PIVOT)] = np.sqrt(pivot_squares)
    return points


def validate_raw_reductions_exact(
        degrees: Iterable[int] = (2, 3),
        evaluation_tolerance: float = 2e-10) -> List[dict]:
    """Validate the exact raw-to-quotient maps without a numerical fit."""
    points = _self_check_points(seed=20260713, count=12)
    reports = []
    for degree in degrees:
        raw_basis = enumerate_metric_basis(degree)
        quotient_basis = basis_at_deg_quotient(degree)
        reduction = raw_to_quotient_matrix_exact(degree)
        if reduction.shape != (len(quotient_basis), len(raw_basis)):
            raise AssertionError("raw reduction shape does not match its bases")
        rank = exact_rank_mod_prime(reduction)
        if rank != len(quotient_basis):
            raise AssertionError("raw reduction is not onto over Q")
        raw_values = _evaluate_basis_on_points(raw_basis, points)
        quotient_values = _evaluate_basis_on_points(quotient_basis, points)
        transported = quotient_values @ reduction.to_scipy_csc(dtype=np.longdouble)
        function_error = float(np.max(np.abs(raw_values - transported)))
        if function_error > evaluation_tolerance:
            raise AssertionError(
                f"raw degree-{degree} reduction error {function_error:.3e}")
        reports.append({
            "degree": degree,
            "shape": reduction.shape,
            "nnz": reduction.nnz,
            "rank_mod_2147483647": rank,
            "functional_max_abs_error": function_error,
        })
    return reports


def validate_exact_inclusions(degrees: Iterable[int] = (3, 4),
                              evaluation_tolerance: float = 2e-10) -> List[dict]:
    """Validate dimensions, characters, injectivity and function identities."""
    validate_vandermonde_reduction_exact()
    points = _self_check_points()
    reports = []
    for d in degrees:
        source_basis = basis_at_deg_quotient(d)
        target_basis = basis_at_deg_quotient(d + 1)
        inclusion = inclusion_matrix_exact(d)
        if inclusion.shape != (len(target_basis), len(source_basis)):
            raise AssertionError("inclusion shape does not match quotient bases")

        # Each T_j shifts the holomorphic character by char(Z_j); therefore
        # every resulting Hermitian entry still has equal left/right chars.
        source_holo = enumerate_standard_holo(d)
        target_holo = enumerate_standard_holo(d + 1)
        for coordinate in range(6):
            multiplier = holomorphic_multiplication_matrix_exact(d, coordinate)
            shift = tuple(int(x) for x in COORD_CHARS[coordinate])
            for source_index, column in enumerate(multiplier.entries_by_column):
                expected = tuple((a + b) % 2
                                 for a, b in zip(char_of(source_holo[source_index]), shift))
                if any(char_of(target_holo[row]) != expected for row, _ in column):
                    raise AssertionError("holomorphic character shift is inconsistent")
        if any(basis_character(element) != char_of(element["ij"])
               for element in target_basis):
            raise AssertionError("target basis character labelling failed")

        rank_proof = full_column_rank_proof_prime(inclusion)
        if rank_proof is None:
            raise AssertionError(
                f"inclusion {d}->{d + 1} rank over Q is inconclusive: "
                "no tested prime witnessed full column rank")
        rank_prime, rank = rank_proof

        source_values = _evaluate_basis_on_points(source_basis, points)
        target_values = _evaluate_basis_on_points(target_basis, points)
        transported = np.empty_like(source_values)
        for column, entries in enumerate(inclusion.entries_by_column):
            transported[:, column] = sum(
                (np.longdouble(value.numerator) / np.longdouble(value.denominator))
                * target_values[:, row]
                for row, value in entries
            )
        function_error = float(np.max(np.abs(transported - source_values)))
        if function_error > evaluation_tolerance:
            raise AssertionError(f"functional inclusion error {function_error:.3e}")

        constant_source = constant_coordinates_exact(d)
        constant_target = constant_coordinates_exact(d + 1)
        if inclusion.matvec_exact(constant_source) != constant_target:
            raise AssertionError("constant coordinates do not commute with inclusion")
        constant_values = target_values @ np.asarray(
            [np.longdouble(value.numerator) / np.longdouble(value.denominator)
             for value in constant_target], dtype=np.longdouble)
        constant_error = float(np.max(np.abs(constant_values - 1)))
        if constant_error > evaluation_tolerance:
            raise AssertionError(f"constant evaluation error {constant_error:.3e}")

        reports.append({
            "degree": d,
            "shape": inclusion.shape,
            "nnz": inclusion.nnz,
            "rank_proof_prime": rank_prime,
            "rank_mod_proof_prime": rank,
            "functional_max_abs_error": function_error,
            "constant_max_abs_error": constant_error,
            "constant_nonzeros_source": sum(value != 0 for value in constant_source),
            "constant_nonzeros_target": sum(value != 0 for value in constant_target),
        })
    return reports


# ================================================================
#  Vérifs auto-consistantes (invocable en script)
# ================================================================
if __name__ == "__main__":
    # counts prédits par quotient_basis (Σ m² avec m = # std monoms
    # of char χ à degré d, PIVOT = (0,1,2)) :
    expected_H = {1: 10, 2: 58, 3: 218, 4: 610, 5: 1402}
    expected_cumul = {1: 10, 2: 68, 3: 286, 4: 896, 5: 2298}
    print("=" * 68)
    print(f"K3 CAP quotient engine — self-check  |  PIVOT = {PIVOT}")
    print("=" * 68)
    print(f"{'d':>2} {'std_holo':>10} {'H_d^{G,Herm}':>16} {'expected':>10} "
          f"{'V_leq d':>10} {'expected':>10}")
    cumul = 0
    for d in range(1, 6):
        holo = enumerate_standard_holo(d)
        sec = enumerate_sector_quotient(d)
        cumul += len(sec)
        ok_H = "✓" if len(sec) == expected_H[d] else "✗"
        ok_c = "✓" if cumul == expected_cumul[d] else "✗"
        print(f"{d:>2} {len(holo):>10} {len(sec):>16} {expected_H[d]:>10} "
              f"{cumul:>10} {expected_cumul[d]:>10}   {ok_H}{ok_c}")
    # Total à V_≤3 : 286  (sans constante)
    base = basis_upto_quotient(3)
    print(f"\nbasis_upto_quotient(3)  →  {len(base)} fn  "
          f"({'✓' if len(base) == 286 else '✗ attendu 286'})")
    base4 = basis_upto_quotient(4)
    print(f"basis_upto_quotient(4)  →  {len(base4)} fn  "
          f"({'✓' if len(base4) == 896 else '✗ attendu 896'})")

    print("\nExact raw invariant reductions:")
    for report in validate_raw_reductions_exact():
        print(
            f"  raw d={report['degree']} -> quotient: shape={report['shape']}, "
            f"nnz={report['nnz']}, rank_Q={report['rank_mod_2147483647']}, "
            f"function_err={report['functional_max_abs_error']:.3e}  ✓"
        )

    print("\nExact inclusions q -> s*q/s:")
    for report in validate_exact_inclusions():
        d = report["degree"]
        print(
            f"  V_{d} -> V_{d + 1}: shape={report['shape']}, "
            f"nnz={report['nnz']}, rank_Q={report['rank_mod_proof_prime']} "
            f"(mod {report['rank_proof_prime']}), "
            f"function_err={report['functional_max_abs_error']:.3e}, "
            f"constant_err={report['constant_max_abs_error']:.3e}  ✓"
        )
