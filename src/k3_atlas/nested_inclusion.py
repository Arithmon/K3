#!/usr/bin/env python3
"""
k3_cap_nested_inclusion.py — inclusion exacte E_{d→d+1} entre bases
quotientées (module importable).

Extraction de k3_cap_ritz_diag_nested_34.py (qui s'exécute à l'import)
pour réutilisation : batterie gate du moteur Kähler cohérent, refit,
futur emboîtement V₄→V₅. Prescription commune GPT §4 / Aristotle Q2 :
E est construite par multiplication du numérateur par s + réduction
quotient exacte (Fractions), puis vérifiée pointwise Q_d = Q_{d+1}·E.

Convention colonnes : E (1+nb_hi, 1+nb_lo), colonne 0 = constante ;
Q_d(x) = [1, q̃_1..q̃_nb] au format des assembleurs (constante en tête).
"""
from __future__ import annotations

import io
import sys
from fractions import Fraction

import numpy as np

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

from .invariant_quotient_ring import (
    basis_at_deg_quotient, char_of, PIVOT, NONPIVOT,
)
from .spectral_basis import minor_inv_times_T_float


def exact_pivot_relations():
    """A[p][t] Fraction : Z_p² = Σ_{t∈NONPIVOT} A[p][t] Z_t² mod ⟨Q_0,Q_1,Q_2⟩.
    Système : Σ_i μ_i^a Z_i² = 0, a = 0,1,2 ⟹ VS·x_S = −VT·x_T."""
    mu = [Fraction(1), Fraction(2), Fraction(3),
          Fraction(5), Fraction(7), Fraction(11)]
    VS = [[mu[p] ** a for p in PIVOT] for a in range(3)]
    VT = [[mu[t] ** a for t in NONPIVOT] for a in range(3)]
    (a, b, c), (d, e, f), (g, h, i) = VS
    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    inv = [[(e * i - f * h) / det, (c * h - b * i) / det, (b * f - c * e) / det],
           [(f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det],
           [(d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det]]
    A = {}
    for r, p in enumerate(PIVOT):
        A[p] = {}
        for cix, t in enumerate(NONPIVOT):
            A[p][t] = -sum(inv[r][k] * VT[k][cix] for k in range(3))
    return A


A_EXACT = exact_pivot_relations()
# contrôle vs float engine (au chargement — bon marché et bloquant)
_Afl = minor_inv_times_T_float(PIVOT, NONPIVOT)
for _r, _p in enumerate(PIVOT):
    for _c, _t in enumerate(NONPIVOT):
        assert abs(float(A_EXACT[_p][_t]) - _Afl[_r, _c]) < 1e-12


def raise_holo(I, p):
    """z^I · z_p en monômes STANDARD de degré |I|+1 : [(Fraction, J), ...]."""
    if p not in PIVOT or I.count(p) == 0:
        return [(Fraction(1), tuple(sorted(I + (p,))))]
    J0 = list(I)
    J0.remove(p)
    return [(A_EXACT[p][t], tuple(sorted(J0 + [t, t]))) for t in NONPIVOT]


def herm_pairs(be):
    """Décomposition C = Σ coef·E_IK de l'élément réel (coef ∈ {1, ±i})."""
    I, K, typ = be["ij"], be["kl"], be["type"]
    if typ == "self":
        return [((Fraction(1), Fraction(0)), I, K)]
    if typ == "real_pair":
        return [((Fraction(1), Fraction(0)), I, K),
                ((Fraction(1), Fraction(0)), K, I)]
    # imag_pair : i·E_IK − i·E_KI
    return [((Fraction(0), Fraction(1)), I, K),
            ((Fraction(0), Fraction(-1)), K, I)]


def build_E(d_lo=3):
    """E (1+nb_hi, 1+nb_lo) : colonnes = éléments V_lo exprimés dans V_hi.
    Exact (Fractions), converti float en sortie. Retourne (E, b_lo, b_hi)."""
    b_lo = basis_at_deg_quotient(d_lo)
    b_hi = basis_at_deg_quotient(d_lo + 1)
    idx_hi = {}
    for j, be in enumerate(b_hi):
        idx_hi[(be["ij"], be["kl"], be["type"])] = j
    nb_lo, nb_hi = len(b_lo), len(b_hi)
    E = np.zeros((nb_hi + 1, nb_lo + 1))
    E[0, 0] = 1.0                                     # constante → constante
    for e, be in enumerate(b_lo):
        C4 = {}                                       # (J,L) → [Fr_re, Fr_im]
        for (cre, cim), I, K in herm_pairs(be):
            for p in range(6):
                for aI, J in raise_holo(I, p):
                    for aK, L in raise_holo(K, p):
                        w = aI * aK                   # réel exact
                        cur = C4.setdefault((J, L), [Fraction(0), Fraction(0)])
                        cur[0] += cre * w
                        cur[1] += cim * w
        col = np.zeros(nb_hi)
        seen = set()
        for (J, L), (vre, vim) in C4.items():
            if (J, L) in seen:
                continue
            assert char_of(J) == char_of(L), "caractère non invariant"
            if J == L:
                assert vim == 0, "coeff diagonal non réel"
                col[idx_hi[(J, L, "self")]] += float(vre)
                seen.add((J, L))
            else:
                Jl, Ll = (J, L) if J < L else (L, J)
                v = C4.get((Jl, Ll), [Fraction(0), Fraction(0)])
                col[idx_hi[(Jl, Ll, "real_pair")]] += float(v[0])
                col[idx_hi[(Jl, Ll, "imag_pair")]] += float(v[1])
                seen.add((Jl, Ll))
                seen.add((Ll, Jl))
        E[1:, 1 + e] = col
    return E, b_lo, b_hi
