#!/usr/bin/env python3
"""
THE TAYLOR MODEL, in the proper sense, on the field
Q_gamma: a polynomial of degree at most N in the box symbols in [-1,1]^4
PLUS a remainder bounded by the truncation TAIL, not by a derivative
evaluated in interval arithmetic on the box.

Why this layer exists (measured, not assumed):

  * the earlier layers are "exact jets at the centre plus a remainder
    D^{n}(BOÎTE) enclos en intervalle ». Quatre ablations (C48, T4b,
    showed that the central jet always passes and that 100 percent of the
    wall is the remainder: the interval evaluation of the derivative on the
    box adds the moduli of the 218 elements where the true derivative
    cancels them (a factor 1e3 to 6e9 depending on the rank).
  * the AFFINE spike measured the other
    end: the part LINEAR in delta is transported exactly (determinant:
    7.5e-5 at 4e-3, 3.2e-4 at 1.7e-2, the order of the true float span, about 9 percent
    of the margin to target) but its RADIUS explodes (3.4 to 1.4e3): at
    first order the giant intermediates cancel in the value,
    never in the radius.

  So the synthesis is needed: keep the polynomial EXACTLY up to
  degree N (the cancellations survive to that rank) and bound the
  remainder by the tail of the truncated product, whose scale is that
  of the local COEFFICIENTS (tame), not that of an enclosed derivative.

Algèbre (standard TM, toutes les bornes extérieures) :
  x = P(ε) + I,  P = Σ_{|α| ≤ N} p_α ε^α,  |I| ≤ rem
  · produit : P·Q tronqué au degré N ; queue = Σ_{i+j>N} A_i B_j
    (A_i/B_j = normes par degré) ; + ‖P‖·rem_Q + ‖Q‖·rem_P + rem·rem
  · inverse : z = (x − p₀)/p₀ (constante nulle) ⟹ 1/x =
    (1/p₀)(1 − z + z² − z³) + queue, |queue| ≤ |1/p₀|·q⁴/(1−q)
  · racine : √x = √p₀·(1 + z/2 − z²/8 + z³/16) + queue,
    |queue| ≤ |√p₀|·q⁴/(8(1−q))   (|C(½,k)| ≤ ⅛ pour k ≥ 2)
  · enclosure : p₀ + Σ_{|α|≥1} p_α·rng(α) ± rem, avec rng(α) = [0,1]
    if every exponent is EVEN (a free lesson here and
    liante partout), [−1,1] sinon.

Garde de branche : identique t2/t3/t4 (`civ_sqrt_principal` sur p₀)
PLUS the range guard (the whole range must avoid the cut).

Ordre : `K3_TM_ORDER` (défaut 3). N=2 → 15 monômes, N=3 → 35, N=4 → 70.

Self-test (négatifs inclus) :
  M1 polynomial exactness: a polynomial of degree at most N is enclosed
     EXACTEMENT (rem = 0) ; NÉGATIF : degré N+1 ⟹ rem > 0
  M2 identités : f·inv(f) ∋ 1, (√f)² ∋ f (bornes)
  M3 soundness réelle : 400 points float ⊆ enclosure (4 composantes
     ET det) sur B@{8e-3, 1.7e-2}
  M4 dégénéré h=0 ≡ moteur float (rel < 5e-12)
  M5 NÉGATIFS : w = −1 échoue ; t_bad det.hi < 0
  M6 garde de branche : adverse C63 refusé, boîte latérale acceptée
  M7 order monotonicity: the enclosure at N=3 is tighter than at N=2
     on the hard cell (the mechanism does what it announces)

Usage : k3_cap_tm_kernel.py --selftest
"""
from __future__ import annotations

import io
import itertools
import os
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

_argv = sys.argv
sys.argv = [sys.argv[0]]
from mpmath import iv, mp                                          # noqa: E402
from .witness_parametrisation import B3, B3_MULTIS, B3_IDX               # noqa: E402
from .interval_arithmetic import (                            # noqa: E402
    CIV, CZERO, CONE, MU_INT, NG, BranchCutError, civ_sqrt_principal,
    iv_bounds, minor_inv_times_T_exact, riv)
sys.argv = _argv

IV0 = iv.mpf(0)
IV1 = iv.mpf(1)
IVPM = iv.mpf([-1, 1])
IV01 = iv.mpf([0, 1])

TM_ORDER = int(os.environ.get("K3_TM_ORDER", "3"))
# C74 / a reviewer §4 (les DEUX reviews, indépendamment) : `TM_ORDER` ne fixait
# than the polynomial BASE: inverse and root were truncated at z^3 IN HARD CODE, hence
# N-independent, and their tail q^4/(1-q) CAPPED the remainder (about 2.6e-3 at
# h = 1.7e-2: a naive N=6 would have gained nothing). Both parameters are
# désormais séparés et sérialisables :
#   POLY_DEG = degré de la base monomiale (convolution tronquée)
#   UNARY_SERIES_DEG = profondeur des séries analytiques inv/√
UNARY_SERIES_DEG = int(os.environ.get("K3_TM_SERIES", str(TM_ORDER)))
if UNARY_SERIES_DEG < 1:
    raise ValueError("K3_TM_SERIES ≥ 1")

# EXACT binomial coefficients C(1/2, k) in rational arithmetic for the root
_BINOM_HALF = [Fraction(1)]
for _k in range(1, UNARY_SERIES_DEG + 1):
    _BINOM_HALF.append(_BINOM_HALF[-1]
                       * (Fraction(1, 2) - (_k - 1)) / _k)

# --- base monomiale ε^α, |α| ≤ N -------------------------------------------------
MONO = [m for m in itertools.product(range(TM_ORDER + 1), repeat=NG)
        if sum(m) <= TM_ORDER]
MONO.sort(key=lambda m: (sum(m), m))
NM = len(MONO)
MIDX = {m: i for i, m in enumerate(MONO)}
DEG = [sum(m) for m in MONO]
GRADE = [[i for i in range(NM) if DEG[i] == k]
         for k in range(TM_ORDER + 1)]
EVEN = [all(e % 2 == 0 for e in m) for m in MONO]
# splits α + β = γ (ordre fixe : α croissant)
SPLITS = []
for g, gm in enumerate(MONO):
    lst = []
    for am in itertools.product(*[range(e + 1) for e in gm]):
        bm = tuple(gm[k] - am[k] for k in range(NG))
        lst.append((MIDX[am], MIDX[bm]))
    SPLITS.append(lst)


def iv_absmax(x):
    lo, hi = mp.mpf(x.a), mp.mpf(x.b)
    return iv.mpf(max(abs(lo), abs(hi)))


def _dist0(x):
    lo, hi = mp.mpf(x.a), mp.mpf(x.b)
    if lo > 0:
        return iv.mpf(lo)
    if hi < 0:
        return iv.mpf(-hi)
    return IV0


def civ_absmax(c: CIV):
    return iv.sqrt(iv_absmax(c.re) ** 2 + iv_absmax(c.im) ** 2)


def civ_absmin(c: CIV):
    return iv.sqrt(_dist0(c.re) ** 2 + _dist0(c.im) ** 2)


def _hi(x):
    return mp.mpf(x.b)


# ===========================================================================
#  the range-aware guard — sérialisation des diagnostics de branche
#  A guard refusal must be reconstructible from the artefact ALONE.
#  Floats serve the JSON, 25-digit strings serve the audit: a
#  float rounded TO NEAREST can lie about an inclusion, and
#  it is precisely the inclusion "the guard disc avoids the cut" that
#  est en cause ici.
# ===========================================================================
def _iv_diag(x):
    lo, hi = mp.mpf(x.a), mp.mpf(x.b)
    return {"lo": float(lo), "hi": float(hi),
            "lo_exact": mp.nstr(lo, 25), "hi_exact": mp.nstr(hi, 25)}


def _civ_diag(c: CIV):
    return {"re": _iv_diag(c.re), "im": _iv_diag(c.im),
            "absmin": _iv_diag(civ_absmin(c)),
            "absmax": _iv_diag(civ_absmax(c))}


# C85 : instrumentation SÉPARÉE des queues (produits vs séries unaires)
# + `UNARY_TAIL_SCALE` : diviseur de test des SEULES queues unaires —
# the targeted negative control (the truth MUST come out if the disc is
# uniquement cette borne-là). Jamais ≠ 1 en production.
STATS = {}
UNARY_TAIL_SCALE = 1.0

# the cut guard : quelle branche de la garde √ a autorisé chaque appel. Compté
# separately, so that "the range-aware branch helps or does not help" is
# a measurement and not an impression.
GUARD_STATS = {"sqrt_disc": 0, "sqrt_range": 0}


def reset_guard_stats():
    GUARD_STATS.update({"sqrt_disc": 0, "sqrt_range": 0})


def reset_stats():
    STATS.clear()
    STATS.update({"n_inv": 0, "n_sqrt": 0, "n_mul": 0,
                  "max_q_inv": 0.0, "max_q_sqrt": 0.0,
                  "sum_inv_series_tails": 0.0,
                  "sum_sqrt_series_tails": 0.0,
                  "sum_product_tails": 0.0})


reset_stats()


def _mid_iv(x):
    """POINT (intervalle dégénéré) au milieu de x — ancrage the pointwise anchoring."""
    return iv.mpf(mp.mpf((mp.mpf(x.a) + mp.mpf(x.b)) / 2))


def _mid_civ(c: CIV) -> CIV:
    return CIV(_mid_iv(c.re), _mid_iv(c.im))


# ===========================================================================
#  Taylor-modèle réel
# ===========================================================================
class TMR:
    """P(ε) + I : p[NM] intervalles, rem ≥ 0 (|I| ≤ rem)."""
    __slots__ = ("p", "rem", "_gr")

    def __init__(self, p, rem=None):
        self.p = p
        self.rem = rem if rem is not None else IV0
        self._gr = None

    @staticmethod
    def const(v):
        p = [IV0] * NM
        p[0] = v
        return TMR(p)

    def grades(self):
        """normes par degré A_k = Σ_{|α|=k} |p_α| (cachées)."""
        if self._gr is None:
            self._gr = [sum((iv_absmax(self.p[i]) for i in GRADE[k]),
                            IV0) for k in range(TM_ORDER + 1)]
        return self._gr

    def norm(self):
        return sum(self.grades(), IV0) + self.rem

    def dev(self):
        """borne de |x − p₀| (déviation autour du terme constant)."""
        g = self.grades()
        return sum(g[1:], IV0) + self.rem

    def to_iv(self):
        acc = self.p[0]
        for i in range(1, NM):
            acc = acc + self.p[i] * (IV01 if EVEN[i] else IVPM)
        return acc + IVPM * self.rem

    def __add__(a, b):
        if isinstance(b, TMR):
            return TMR([x + y for x, y in zip(a.p, b.p)],
                       a.rem + b.rem)
        p = list(a.p)
        p[0] = p[0] + b
        return TMR(p, a.rem)

    def __sub__(a, b):
        if isinstance(b, TMR):
            return TMR([x - y for x, y in zip(a.p, b.p)],
                       a.rem + b.rem)
        p = list(a.p)
        p[0] = p[0] - b
        return TMR(p, a.rem)

    def __neg__(a):
        return TMR([-x for x in a.p], a.rem)

    def __mul__(a, b):
        if not isinstance(b, TMR):
            return TMR([x * b for x in a.p], a.rem * iv_absmax(b))
        p = []
        for g in range(NM):
            acc = IV0
            for ia, ib in SPLITS[g]:
                acc = acc + a.p[ia] * b.p[ib]
            p.append(acc)
        A, B = a.grades(), b.grades()
        tail = IV0
        for i in range(TM_ORDER + 1):
            for j in range(TM_ORDER + 1 - i, TM_ORDER + 1):
                tail = tail + A[i] * B[j]
        nA = sum(A, IV0)
        nB = sum(B, IV0)
        STATS["n_mul"] += 1
        STATS["sum_product_tails"] += float(mp.mpf(tail.b))
        rem = tail + nA * b.rem + nB * a.rem + a.rem * b.rem
        return TMR(p, rem)

    def inv(a):
        p0 = a.p[0]
        m = _dist0(p0)
        # C83 : borne INFÉRIEURE certifiée de la distance à 0 (.a),
        # not the upper one: an enclosure [0, eps] is not "nonzero"
        if not (mp.mpf(m.a) > 0):
            raise BranchCutError("TM inv : terme constant contenant 0")
        # POINTWISE ANCHORING. p_0 is an INTERVAL, so p_0 - p_0 is not {0}
        # and q = |1/p_0|.dev(a) did NOT bound the z actually built.
        # We anchor on a POINT c (a = c(1+z) is then an explicit
        # analytic identity) and bound q by the norm of the built z.
        c = _mid_iv(p0)
        u = IV1 / c
        U = iv_absmax(u)
        z = (a - c) * u
        q = z.norm()                          # ≥ sup |z| par construction
        if not (_hi(q) < 1):
            raise BranchCutError("TM inv : déviation ≥ rayon de "
                                 "convergence (q ≥ 1)")
        # série géométrique tronquée à z^K (C74 : K = UNARY_SERIES_DEG,
        # plus z³ en dur), queue |1/p₀|·q^{K+1}/(1−q)
        acc = TMR.const(IV1)
        zk = TMR.const(IV1)
        for _ in range(UNARY_SERIES_DEG):
            zk = zk * z
            acc = acc - zk if _ % 2 == 0 else acc + zk
        s = acc * u
        qh = iv.mpf(_hi(q))
        tail = U * qh ** (UNARY_SERIES_DEG + 1) / (IV1 - qh)
        if UNARY_TAIL_SCALE != 1.0:
            tail = tail / iv.mpf(UNARY_TAIL_SCALE)
        STATS["n_inv"] += 1
        STATS["max_q_inv"] = max(STATS["max_q_inv"], float(qh))
        STATS["sum_inv_series_tails"] += float(mp.mpf(tail.b))
        s.rem = s.rem + tail
        s._gr = None
        return s


# ===========================================================================
#  Taylor-modèle complexe
# ===========================================================================
class TMC:
    __slots__ = ("p", "rem", "_gr")

    def __init__(self, p, rem=None):
        self.p = p
        self.rem = rem if rem is not None else IV0
        self._gr = None

    @staticmethod
    def const(v: CIV):
        p = [CZERO] * NM
        p[0] = v
        return TMC(p)

    def grades(self):
        if self._gr is None:
            self._gr = [sum((civ_absmax(self.p[i]) for i in GRADE[k]),
                            IV0) for k in range(TM_ORDER + 1)]
        return self._gr

    def dev(self):
        g = self.grades()
        return sum(g[1:], IV0) + self.rem

    def __add__(a, b):
        return TMC([x + y for x, y in zip(a.p, b.p)], a.rem + b.rem)

    def __sub__(a, b):
        return TMC([x - y for x, y in zip(a.p, b.p)], a.rem + b.rem)

    def __neg__(a):
        return TMC([-x for x in a.p], a.rem)

    def __mul__(a, b):
        p = []
        for g in range(NM):
            acc = CZERO
            for ia, ib in SPLITS[g]:
                acc = acc + a.p[ia] * b.p[ib]
            p.append(acc)
        A, B = a.grades(), b.grades()
        tail = IV0
        for i in range(TM_ORDER + 1):
            for j in range(TM_ORDER + 1 - i, TM_ORDER + 1):
                tail = tail + A[i] * B[j]
        nA = sum(A, IV0)
        nB = sum(B, IV0)
        STATS["n_mul"] += 1
        STATS["sum_product_tails"] += float(mp.mpf(tail.b))
        rem = tail + nA * b.rem + nB * a.rem + a.rem * b.rem
        return TMC(p, rem)

    def mul_real(a, r):
        return TMC([x.mul_real(r) for x in a.p], a.rem * iv_absmax(r))

    def mul_rtm(a, r: TMR):
        """produit par un TM RÉEL (miroir de mul_rt2)."""
        p = []
        for g in range(NM):
            acc = CZERO
            for ia, ib in SPLITS[g]:
                acc = acc + a.p[ia].mul_real(r.p[ib])
            p.append(acc)
        A, B = a.grades(), r.grades()
        tail = IV0
        for i in range(TM_ORDER + 1):
            for j in range(TM_ORDER + 1 - i, TM_ORDER + 1):
                tail = tail + A[i] * B[j]
        STATS["n_mul"] += 1
        STATS["sum_product_tails"] += float(mp.mpf(tail.b))
        rem = (tail + sum(A, IV0) * r.rem + sum(B, IV0) * a.rem
               + a.rem * r.rem)
        return TMC(p, rem)

    def conj(a):
        return TMC([x.conj() for x in a.p], a.rem)

    def re_tm(a) -> TMR:
        return TMR([x.re for x in a.p], a.rem)

    def im_tm(a) -> TMR:
        return TMR([x.im for x in a.p], a.rem)

    def abs2_tm(a) -> TMR:
        return (a * a.conj()).re_tm()

    def to_iv_pair(a):
        return a.re_tm().to_iv(), a.im_tm().to_iv()

    def norm(self):
        return sum(self.grades(), IV0) + self.rem

    def inv(a):
        p0 = a.p[0]
        if not (mp.mpf(civ_absmin(p0).a) > 0):        # C83 : .a
            raise BranchCutError(
                "TM inv (C) : p₀ peut contenir 0",
                {"guard": "inv_p0_contains_zero", "kernel": "mpmath",
                 "p0": _civ_diag(p0)})
        c = _mid_civ(p0)                              # the pointwise anchoring : ancrage POINT
        if not (mp.mpf(civ_absmin(c).a) > 0):
            raise BranchCutError(
                "TM inv (C) : ancre nulle",
                {"guard": "inv_anchor_zero", "kernel": "mpmath",
                 "p0": _civ_diag(p0), "anchor": _civ_diag(c)})
        u = CONE.div(c)
        U = civ_absmax(u)
        z = (a - TMC.const(c)).mul_civ(u)
        q = z.norm()                                  # ≥ sup |z|
        if not (_hi(q) < 1):
            raise BranchCutError(
                "TM inv (C) : q ≥ 1",
                {"guard": "inv_q_ge_1", "kernel": "mpmath",
                 "p0": _civ_diag(p0), "anchor": _civ_diag(c),
                 "q": _iv_diag(q)})
        acc = TMC.const(CONE)
        zk = TMC.const(CONE)
        for _ in range(UNARY_SERIES_DEG):
            zk = zk * z
            acc = acc - zk if _ % 2 == 0 else acc + zk
        s = acc.mul_civ(u)
        qh = iv.mpf(_hi(q))
        tail = U * qh ** (UNARY_SERIES_DEG + 1) / (IV1 - qh)
        if UNARY_TAIL_SCALE != 1.0:
            tail = tail / iv.mpf(UNARY_TAIL_SCALE)
        STATS["n_inv"] += 1
        STATS["max_q_inv"] = max(STATS["max_q_inv"], float(qh))
        STATS["sum_inv_series_tails"] += float(mp.mpf(tail.b))
        s.rem = s.rem + tail
        s._gr = None
        return s

    def mul_civ(a, c: CIV):
        return TMC([x * c for x in a.p], a.rem * civ_absmax(c))

    def div(a, b):
        return a * b.inv()

    def sqrt_principal(a):
        """Principal root: a containment guard on p_0 plus a RANGE guard (the range
        entière évite (−∞, 0]) ; série binomiale tronquée + queue."""
        p0 = a.p[0]
        # the pointwise anchoring : ancrage POINT — a = c(1+z), q = norm(z) borne |z|
        c = _mid_civ(p0)
        if not (mp.mpf(civ_absmin(c).a) > 0):
            raise BranchCutError(
                "TM √ : ancre nulle",
                {"guard": "sqrt_anchor_zero", "kernel": "mpmath",
                 "p0": _civ_diag(p0), "anchor": _civ_diag(c)})
        u = CONE.div(c)
        U = civ_absmax(u)
        z = (a - TMC.const(c)).mul_civ(u)
        q = z.norm()
        if not (_hi(q) < 1):
            raise BranchCutError(
                "TM √ : q ≥ 1",
                {"guard": "sqrt_q_ge_1", "kernel": "mpmath",
                 "p0": _civ_diag(p0), "anchor": _civ_diag(c),
                 "q": _iv_diag(q)})
        # cut guard on the WHOLE RANGE: |a - c| <= |c|.q, so
        # it suffices that the disc of radius |c|q around c avoid
        # (−∞, 0] — Re(c) > |c|q OU |Im(c)| > |c|q
        rad = civ_absmax(c) * q
        rh = _hi(rad)
        re_ok = mp.mpf(c.re.a) > rh
        im_ok = (mp.mpf(c.im.a) > rh) or (mp.mpf(c.im.b) < -rh)
        disc_ok = re_ok or im_ok
        # the cut guard — branche RANGE-AWARE. Le lemme de convexité (note
        # `k3_cap_r12b_c103_range_guard_proof_2026_07_28.md` §2) autorise
        # any CONVEX K containing the anchor, avoiding the cut and inside
        # the closed disc of radius rho|c| with rho < 1. The isotropic disc is only ONE choice of
        # K; here K = R intersected with that disc, where R is the rectangular enclosure
        # of the model. The two branches do not imply each other (measured:
        # 240 of 400 shell cells pass by range and not by
        # disc), hence the DISJUNCTION, and therefore monotone safety:
        # aucune cellule qui passait ne peut se mettre à échouer, et
        # aucune borne en aval ne change.
        range_diag = None
        if not disc_ok:
            Rre, Rim = a.re_tm().to_iv(), a.im_tm().to_iv()
            r_lo, r_hi = mp.mpf(Rre.a), mp.mpf(Rre.b)
            i_lo, i_hi = mp.mpf(Rim.a), mp.mpf(Rim.b)
            cr, ci = mp.mpf(c.re.a), mp.mpf(c.im.a)
            # (G2b) le rectangle rencontre (−∞,0] ssi im ∋ 0 ET re_lo ≤ 0
            cut_free = (i_lo > 0) or (i_hi < 0) or (r_lo > 0)
            # (G1) the anchor is in R, exactly, since c is a POINT. It carries
            # hypothesis (1) of the lemma AND, with (G2b), forbids an anchor
            # on the cut (which the disc test excluded as a side
            # effect and which this branch would no longer receive).
            anchor_in = (r_lo <= cr <= r_hi) and (i_lo <= ci <= i_hi)
            # (G1b) hypothèse (3) : ρ = rad/|c| < 1. NON impliquée par
            # q < 1, because rad is rounded UPWARDS.
            rho_ok = rh < mp.mpf(civ_absmin(c).a)
            range_ok = bool(cut_free and anchor_in and rho_ok)
            range_diag = {"cut_free": bool(cut_free),
                          "anchor_in_R": bool(anchor_in),
                          "rho_lt_1": bool(rho_ok),
                          "accepted": range_ok,
                          "R_re": _iv_diag(Rre), "R_im": _iv_diag(Rim)}
        else:
            range_ok = False
        if not (disc_ok or range_ok):
            # This is THE guard that produces the branch shell.
            # Everything that makes it reproducible starts here, including
            # the three signed slacks that say BY HOW MUCH the disc
            # générique déborde, et désormais POURQUOI la branche
            # range-aware form could not save it either.
            raise BranchCutError(
                "TM √ : la plage touche la coupure (−∞, 0]",
                {"guard": "sqrt_disc_touches_cut", "kernel": "mpmath",
                 "p0": _civ_diag(p0), "anchor": _civ_diag(c),
                 "q": _iv_diag(q), "guard_radius": _iv_diag(rad),
                 "re_ok": bool(re_ok), "im_ok": bool(im_ok),
                 "slack_re": _iv_diag(iv.mpf(mp.mpf(c.re.a) - rh)),
                 "slack_im_pos": _iv_diag(
                     iv.mpf(mp.mpf(c.im.a) - rh)),
                 "slack_im_neg": _iv_diag(
                     iv.mpf(-rh - mp.mpf(c.im.b))),
                 "range_guard": range_diag})
        GUARD_STATS["sqrt_disc"] += int(disc_ok)
        GUARD_STATS["sqrt_range"] += int(range_ok)
        w0 = civ_sqrt_principal(c)
        # binôme tronqué à z^K, coefficients C(1/2,k) EXACTS ;
        # queue ≤ |√p₀|·q^{K+1}/(8(1−q)) car |C(1/2,k)| ≤ 1/8 (k ≥ 2)
        acc = TMC.const(CONE)
        zk = TMC.const(CONE)
        for k in range(1, UNARY_SERIES_DEG + 1):
            zk = zk * z
            acc = acc + zk.mul_real(riv(_BINOM_HALF[k]))
        s = acc.mul_civ(w0)
        qh = iv.mpf(_hi(q))
        tail = civ_absmax(w0) * qh ** (
            UNARY_SERIES_DEG + 1) / (iv.mpf(8) * (IV1 - qh))
        if UNARY_TAIL_SCALE != 1.0:
            tail = tail / iv.mpf(UNARY_TAIL_SCALE)
        STATS["n_sqrt"] += 1
        STATS["max_q_sqrt"] = max(STATS["max_q_sqrt"], float(qh))
        STATS["sum_sqrt_series_tails"] += float(mp.mpf(tail.b))
        s.rem = s.rem + tail
        s._gr = None
        return s


def _pow(z: TMC, k: int) -> TMC:
    r = z
    for _ in range(k - 1):
        r = r * z
    return r


# ===========================================================================
#  Section et métrique — miroir strict de t2_chart_metric
# ===========================================================================
def rotated_sigma_from_coeffs(a1, a2, ur, ui, vr, vi):
    """The COMPONENT of `Im R`, determined from the SIGNS
    of the coefficients and coordinates, not read off an enclosure which,
    on a face-aligned cell, contains 0 by construction.

    `Im R = 2a₁·Re(u)·Im(u) + 2a₂·Re(v)·Im(v)`.

    **The identical zero.** An interval factor `[0, 0]` makes
    its term IDENTICALLY ZERO; classifying it as "positive" (which
    `lo >= 0` did) is wrong as a sign algebra. A term identically
    zero does not contribute to the sum: it must be SET ASIDE, not counted.
    And if BOTH terms are identically zero then `Im R` vanishes: the
    radicand is real on the whole box and the component does not exist,
    so **refusal**, not a default sign.

    Retourne +1 (composante supérieure), −1 (inférieure), ou 0 =
    INDÉTERMINÉ — auquel cas la continuation doit être REFUSÉE.
    """
    def factor_sign(lo, hi):
        """-1 / +1 / 0 (signe non déterminé) / None (IDENTIQUEMENT nul)."""
        if lo == 0 and hi == 0:
            return None
        if lo >= 0:
            return 1
        if hi <= 0:
            return -1
        return 0

    terms = []
    for a, (rlo, rhi), (ilo, ihi) in ((a1, ur, ui), (a2, vr, vi)):
        if a == 0:
            continue                      # coefficient nul : terme absent
        sr, si = factor_sign(rlo, rhi), factor_sign(ilo, ihi)
        if sr is None or si is None:
            continue                      # terme IDENTIQUEMENT nul : écarté
        if sr == 0 or si == 0:
            return 0                      # facteur de signe indéterminé
        terms.append(sr * si * (1 if a > 0 else -1))
    if not terms:
        return 0                          # Im R ≡ 0 : pas de composante
    first = terms[0]
    return first if all(x == first for x in terms) else 0


def tm_sqrt_rotated(a: TMC, sigma: int) -> TMC:
    """the rotated continuation — la DÉTERMINATION TOURNÉE : `√_rot(R) = σ·i·√_principal(−R)`.

    `R` avoids the non-negative reals exactly when `-R` avoids the non-positive ones, so **the guard
    existante s'applique verbatim à `−R`** : il n'y a aucune nouvelle
    guard to prove, and the convexity lemma carries over without
    change (it does not depend on the orientation of the ray).

    The sign is NOT free: the argument computation gives
    `i·√_p(−R) = √_p(R)` si `Im R > 0` et `= −√_p(R)` si `Im R < 0`,
    so **it is the sign of `Im R`**, the component. That is why
    `rotated_sigma_from_coeffs` refuse de deviner quand elle est
    indéterminée.

    On a cell whose interior lies in the lower component and
    whose boundary meets the slice, the rotated root is CONTINUOUS on the closed
    cell and agrees with the principal root on the interior, which is
    exactement la continuation analytique cherchée, là où la
    détermination principale saute.
    """
    if sigma not in (-1, 1):
        raise BranchCutError(
            "TM √_rot : composante INDÉTERMINÉE (σ ∉ {−1, +1}) — la "
            "continuation ne doit pas être devinée",
            {"guard": "rotated_component_undetermined",
             "kernel": "mpmath", "sigma": sigma})
    w = a.mul_real(riv(-1.0)).sqrt_principal()
    return w.mul_civ(CIV(IV0, riv(float(sigma))))       # × (σ·i)


def section_radicands(S, g_col, u0: complex, v0: complex, h: float):
    """The enclosures of the THREE section radicands, computed
    WITHOUT going through the root, hence available even when the guard
    accepte et ne lève aucun diagnostic.

    This is what makes the range-aware prediction **recomputable** instead
    of being read back from an artefact produced under an earlier guard: the
    predicate "R avoids the cut" is a function of the cell alone.
    A strict mirror of the four rows of `tm_chart_cell_section` that
    construisent `R` (mêmes A, mêmes u², v²).
    """
    T = tuple(j for j in range(6) if j not in S)
    others = [c for c in T if c != g_col]
    A_exact = minor_inv_times_T_exact(S, T)
    perm = [list(T).index(g_col), list(T).index(others[0]),
            list(T).index(others[1])]
    A = [[riv(A_exact[r][c]) for c in perm] for r in range(3)]
    hr = riv(h)
    e = [MIDX[tuple(1 if k == j else 0 for k in range(NG))]
         for j in range(NG)]
    pu = [CZERO] * NM
    pu[0] = CIV.from_complex(u0)
    pu[e[0]] = CIV(hr, IV0)
    pu[e[1]] = CIV(IV0, hr)
    pv = [CZERO] * NM
    pv[0] = CIV.from_complex(v0)
    pv[e[2]] = CIV(hr, IV0)
    pv[e[3]] = CIV(IV0, hr)
    u, v = TMC(pu), TMC(pv)
    u2, v2 = u * u, v * v
    out = []
    for r in range(3):
        R = TMC.const(CIV(A[r][0])) + u2.mul_real(A[r][1]) \
            + v2.mul_real(A[r][2])
        re_iv, im_iv = R.re_tm().to_iv(), R.im_tm().to_iv()
        re_lo, re_hi = mp.mpf(re_iv.a), mp.mpf(re_iv.b)
        im_lo, im_hi = mp.mpf(im_iv.a), mp.mpf(im_iv.b)
        out.append({
            "row": r,
            "re": [float(re_lo), float(re_hi)],
            "im": [float(im_lo), float(im_hi)],
            # le prédicat (G2b) : le rectangle évite (−∞, 0]
            "cut_free": bool(im_lo > 0 or im_hi < 0 or re_lo > 0),
            "quadratic_coeffs_exact": [
                str(Fraction(A_exact[r][perm[j]])) for j in range(3)]})
    return out


def tm_chart_cell_section(S, g_col, eps, u0: complex, v0: complex,
                          h: float):
    T = tuple(j for j in range(6) if j not in S)
    others = [c for c in T if c != g_col]
    o1, o2 = others
    A_exact = minor_inv_times_T_exact(S, T)
    perm = [list(T).index(g_col), list(T).index(o1), list(T).index(o2)]
    A = [[riv(A_exact[r][c]) for c in perm] for r in range(3)]

    hr = riv(h)
    e = [MIDX[tuple(1 if k == j else 0 for k in range(NG))]
         for j in range(NG)]
    pu = [CZERO] * NM
    pu[0] = CIV.from_complex(u0)
    pu[e[0]] = CIV(hr, IV0)
    pu[e[1]] = CIV(IV0, hr)
    pv = [CZERO] * NM
    pv[0] = CIV.from_complex(v0)
    pv[e[2]] = CIV(hr, IV0)
    pv[e[3]] = CIV(IV0, hr)
    u, v = TMC(pu), TMC(pv)
    u2, v2 = u * u, v * v

    ZT = [TMC.const(CZERO) for _ in range(6)]
    WT = [[TMC.const(CZERO), TMC.const(CZERO)] for _ in range(6)]
    ZT[g_col] = TMC.const(CONE)
    ZT[o1], ZT[o2] = u, v
    WT[o1][0] = TMC.const(CONE)
    WT[o2][1] = TMC.const(CONE)
    for r, s_coord in enumerate(S):
        R = TMC.const(CIV(A[r][0])) + u2.mul_real(A[r][1]) \
            + v2.mul_real(A[r][2])
        # the root guard does not know WHICH section root it
        # refuses. The radicand is the EXPLICIT quadratic
        # R_r = a₀ + a₁·u² + a₂·v² (coefficients rationnels exacts) —
        # it is attached to the diagnostic, which makes the refusal reconstructible
        # depuis l'artefact seul, et donne à the cut guard son objet de travail.
        try:
            Zs = R.sqrt_principal().mul_real(riv(int(eps[r])))
        except BranchCutError as exc:
            exc.diag.update({
                "site": "section_root", "row": r, "s_coord": s_coord,
                "eps_r": int(eps[r]),
                "quadratic_coeffs_exact": [
                    str(Fraction(A_exact[r][perm[j]])) for j in range(3)],
                "quadratic_coeffs_float": [
                    float(A_exact[r][perm[j]]) for j in range(3)],
                "cell": {"u0": [u0.real, u0.imag],
                         "v0": [v0.real, v0.imag], "h": float(h),
                         "u0_hex": [u0.real.hex(), u0.imag.hex()],
                         "v0_hex": [v0.real.hex(), v0.imag.hex()],
                         "h_hex": float(h).hex()},
                "radicand_enclosure": {
                    "re": _iv_diag(R.re_tm().to_iv()),
                    "im": _iv_diag(R.im_tm().to_iv())}})
            raise
        ZT[s_coord] = Zs
        WT[s_coord][0] = u.mul_real(A[r][1]).div(Zs)
        WT[s_coord][1] = v.mul_real(A[r][2]).div(Zs)

    i_, j_, k_ = S
    VS = ((MU_INT[j_] - MU_INT[i_]) * (MU_INT[k_] - MU_INT[i_])
          * (MU_INT[k_] - MU_INT[j_]))
    det_MS = (ZT[i_] * ZT[j_] * ZT[k_]).mul_real(riv(8 * VS))
    return ZT, WT, det_MS


def tm_chart_metric(Z, W, M_civ, coeffs218, basis=B3, midx=B3_IDX,
                    multis=B3_MULTIS, rho_weight=None):
    s = TMR.const(IV0)
    for a in range(6):
        s = s + Z[a].abs2_tm()
    zW = []
    for al in range(2):
        acc = TMC.const(CZERO)
        for a in range(6):
            acc = acc + Z[a].conj() * W[a][al]
        zW.append(acc)
    WtWb = [[None, None], [None, None]]
    for A_ in range(2):
        for B_ in range(2):
            acc = TMC.const(CZERO)
            for a in range(6):
                acc = acc + W[a][A_] * W[a][B_].conj()
            WtWb[A_][B_] = acc

    m, p = {}, {}
    one = TMC.const(CONE)
    for Iu in multis:
        cnt = Counter(Iu)
        val = one
        for o, mo in cnt.items():
            val = val * _pow(Z[o], mo)
        m[Iu] = val
        pI = [TMC.const(CZERO), TMC.const(CZERO)]
        for a, ma in cnt.items():
            gv = TMC.const(CIV(riv(ma)))
            for o, mo in cnt.items():
                ee = mo - 1 if o == a else mo
                if ee:
                    gv = gv * _pow(Z[o], ee)
            pI[0] = pI[0] + gv * W[a][0]
            pI[1] = pI[1] + gv * W[a][1]
        p[Iu] = pI

    Md = [[TMC.const(M_civ[i][j]) for j in range(6)] for i in range(6)]
    rho = TMR.const(IV0)
    r_vec = []
    for a_ in range(6):
        acc = TMC.const(CZERO)
        for i_ in range(6):
            acc = acc + Z[i_].conj() * Md[i_][a_]
        r_vec.append(acc)
        rho = rho + (acc * Z[a_]).re_tm()
    if not (mp.mpf(rho.to_iv().a) > 0):
        raise BranchCutError(
            "rho non strictement positif (TM)",
            {"guard": "rho_not_positive", "kernel": "mpmath",
             "site": "chart_metric", "rho": _iv_diag(rho.to_iv())})

    v_vec = []
    for A_ in range(2):
        acc = TMC.const(CZERO)
        for a in range(6):
            acc = acc + W[a][A_] * r_vec[a]
        v_vec.append(acc)
    rho_inv = rho.inv()
    rho2_inv = rho_inv * rho_inv
    if rho_weight is not None:
        w_r = TMR.const(riv(rho_weight))
        rho_inv = rho_inv * w_r
        rho2_inv = rho2_inv * w_r
    G = [[None, None], [None, None]]
    for A_ in range(2):
        for B_ in range(2):
            T1 = TMC.const(CZERO)
            for a in range(6):
                for b in range(6):
                    T1 = T1 + W[a][A_] * Md[b][a] * W[b][B_].conj()
            G[A_][B_] = (T1.mul_rtm(rho_inv)
                         - (v_vec[A_] * v_vec[B_].conj())
                         .mul_rtm(rho2_inv))

    s_inv = s.inv()
    s3_inv = s_inv * s_inv * s_inv
    sd1 = s3_inv * s_inv
    sd2 = sd1 * s_inv
    cross = [TMC.const(CZERO), TMC.const(CZERO)]
    id1 = TMR.const(IV0)
    id2 = TMR.const(IV0)
    coeffs = np.asarray(coeffs218, float)
    for e_, be in enumerate(basis):
        c = float(coeffs[e_])
        if c == 0.0:
            continue
        Iu, Kk, typ = be["ij"], be["kl"], be["type"]
        d = len(Iu)
        cr = riv(c)
        mI, mK = m[Iu], m[Kk]
        pI, pK = p[Iu], p[Kk]
        if typ == "self":
            lead = [[pI[A_] * pI[B_].conj() for B_ in range(2)]
                    for A_ in range(2)]
            g1 = [pI[A_] * mK.conj() for A_ in range(2)]
            phi = (mI * mK.conj()).re_tm()
        elif typ == "real_pair":
            lead = [[pI[A_] * pK[B_].conj() + pK[A_] * pI[B_].conj()
                     for B_ in range(2)] for A_ in range(2)]
            g1 = [pI[A_] * mK.conj() + pK[A_] * mI.conj()
                  for A_ in range(2)]
            phi = (mI * mK.conj()).re_tm() * iv.mpf(2)
        else:
            Jc = TMC.const(CIV(IV0, IV1))
            lead = [[Jc * (pI[A_] * pK[B_].conj()
                           - pK[A_] * pI[B_].conj())
                     for B_ in range(2)] for A_ in range(2)]
            g1 = [Jc * (pI[A_] * mK.conj() - pK[A_] * mI.conj())
                  for A_ in range(2)]
            phi = (mI * mK.conj()).im_tm() * iv.mpf(-2)
        c_sd = s3_inv * cr
        c_d_sd1 = sd1 * (cr * d)
        for A_ in range(2):
            for B_ in range(2):
                G[A_][B_] = G[A_][B_] + lead[A_][B_].mul_rtm(c_sd)
            cross[A_] = cross[A_] + g1[A_].mul_rtm(c_d_sd1)
        id1 = id1 + phi * sd1 * (cr * d)
        id2 = id2 + phi * sd2 * (cr * (d * (d + 1)))

    for A_ in range(2):
        for B_ in range(2):
            G[A_][B_] = (G[A_][B_]
                         - cross[A_] * zW[B_].conj()
                         - zW[A_] * cross[B_].conj()
                         + (zW[A_] * zW[B_].conj()).mul_rtm(id2)
                         - WtWb[A_][B_].mul_rtm(id1))

    half = iv.mpf(1) / 2
    g00 = ((G[0][0] + G[0][0].conj()).re_tm()) * half
    g11 = ((G[1][1] + G[1][1].conj()).re_tm()) * half
    g01 = (G[0][1] + G[1][0].conj()).mul_real(half)
    return [g00, g11, g01.re_tm(), g01.im_tm()]


def det_packed_tm(g):
    """Determinant INSIDE the model algebra (correlations between components AND
    les annulations jusqu'au degré N survivent)."""
    return g[0] * g[1] - g[2] * g[2] - g[3] * g[3]


def tm_qfield_certificate(S, g_col, eps, u0, v0, h, M_civ, coeffs,
                          rho_w):
    import time as _time
    t1 = _time.time()
    try:
        Z, W, _ = tm_chart_cell_section(S, g_col, eps, u0, v0, h)
        q = tm_chart_metric(Z, W, M_civ, coeffs, rho_weight=rho_w)
    except BranchCutError as exc:
        return {"h": h, "status": "BRANCH", "error": str(exc)[:120],
                "branch_diag": exc.diag,          # the range-aware guard
                "t_call_s": _time.time() - t1}
    det = det_packed_tm(q)
    q00_lo, q00_hi = iv_bounds(q[0].to_iv())
    det_lo, det_hi = iv_bounds(det.to_iv())
    ok = q00_lo > 0 and det_lo > 0
    gr = det.grades()
    return {"h": h, "order": TM_ORDER,
            "poly_deg": TM_ORDER, "unary_series_deg": UNARY_SERIES_DEG,
            "n_monomials": NM,
            "status": "PASS" if ok
            else ("FAIL_DET" if q00_lo > 0 else "FAIL_PIVOT"),
            "q00": [q00_lo, q00_hi], "det": [det_lo, det_hi],
            "w_det_final": det_hi - det_lo,
            "det_grade_norms": [float(mp.mpf(x.b)) for x in gr],
            "det_remainder": float(mp.mpf(det.rem.b)),
            "w_components": [float(mp.mpf(x.to_iv().delta))
                             for x in q],
            "t_call_s": _time.time() - t1}


# ===========================================================================
#  Self-test
# ===========================================================================
def _selftest():
    import json
    import time
    fails = []
    RES = Path(os.environ.get(
        "K3_RES_DIR", Path(__file__).resolve().parent / "data"))
    print(f"      TM_ORDER = {TM_ORDER} ({NM} monômes)")

    # --- M1 : exactitude polynomiale + négatif ---------------------------------------
    x = TMR.const(IV0)
    x.p[MIDX[tuple(1 if k == 0 else 0 for k in range(NG))]] = IV1
    y = TMR.const(IV0)
    y.p[MIDX[tuple(1 if k == 1 else 0 for k in range(NG))]] = IV1
    pol = x * y * x + TMR.const(iv.mpf(2))       # degré 3
    ok_exact = mp.mpf(pol.rem.b) == 0
    over = pol * x                                # degré 4 > N=3
    ok_neg = mp.mpf(over.rem.b) > 0 if TM_ORDER < 4 else True
    t1 = ok_exact and ok_neg
    fails.append(not t1)
    print(f"[{'PASS' if t1 else 'FAIL'}] M1 polynôme degré ≤ N exact "
          f"(rem = 0) ; NÉGATIF degré N+1 : rem = "
          f"{float(mp.mpf(over.rem.b)):.2e} > 0")

    # --- M2 : identités -------------------------------------------------------------------
    fz = TMC.const(CIV(iv.mpf(2), iv.mpf(1)))
    fz.p[MIDX[tuple(1 if k == 0 else 0 for k in range(NG))]] = \
        CIV(riv(0.05), IV0)
    fz.p[MIDX[tuple(1 if k == 1 else 0 for k in range(NG))]] = \
        CIV(IV0, riv(0.03))
    one = fz * fz.inv()
    b_re, b_im = one.to_iv_pair()
    okinv = (iv_bounds(b_re)[0] <= 1.0 <= iv_bounds(b_re)[1]
             and iv_bounds(b_im)[0] <= 0.0 <= iv_bounds(b_im)[1])
    w = fz.sqrt_principal()
    d = w * w - fz
    d_re, d_im = d.to_iv_pair()
    oksq = (iv_bounds(d_re)[0] <= 0 <= iv_bounds(d_re)[1]
            and iv_bounds(d_im)[0] <= 0 <= iv_bounds(d_im)[1])
    t2 = okinv and oksq
    fails.append(not t2)
    print(f"[{'PASS' if t2 else 'FAIL'}] M2 identités : f·inv(f) ∋ 1 ; "
          f"(√f)² − f ∋ 0 (largeur re "
          f"{float(mp.mpf(d_re.delta)):.1e})")

    # --- setup réel ------------------------------------------------------------------------
    from .witness_registry import load_canonical_MH
    from .interval_arithmetic import build_M_civ
    from .width_attribution import (GAMMA, float_G_pair,
                                            load_boxes_from_direct)
    reg = load_canonical_MH()
    M_H = reg["M_H_canonical"]
    c218 = reg["coeffs218"]
    M_civ = build_M_civ(M_H)
    direct = json.loads(
        (RES / "k3_cap_b1e2iii_p0a2_direct.json").read_text(
            encoding="utf-8"))
    boxes = load_boxes_from_direct(direct)
    bB, bC = boxes[0], boxes[2]
    S, g_col, eps = bB["S"], bB["g"], bB["eps"]
    u0, v0 = bB["u0"], bB["v0"]
    rw = 1.0 - GAMMA

    # --- M3 : soundness réelle -----------------------------------------------------------
    rng = np.random.default_rng(23)
    inside = True
    reps = {}
    for hh in (8e-3, 1.7e-2):
        t0 = time.time()
        rep = tm_qfield_certificate(S, g_col, eps, u0, v0, hh, M_civ,
                                    c218, rw)
        rep["t_meas"] = time.time() - t0
        reps[hh] = rep
        Z, W, _ = tm_chart_cell_section(S, g_col, eps, u0, v0, hh)
        qa = tm_chart_metric(Z, W, M_civ, c218, rho_weight=rw)
        comp = [iv_bounds(x.to_iv()) for x in qa]
        for _ in range(200):
            d = rng.uniform(-hh, hh, 4)
            Gf, Gr = float_G_pair(S, g_col, eps,
                                  u0 + complex(d[0], d[1]),
                                  v0 + complex(d[2], d[3]), M_H, c218)
            Q = Gf - GAMMA * Gr
            vals = [Q[0, 0].real, Q[1, 1].real, Q[0, 1].real,
                    Q[0, 1].imag]
            for c in range(4):
                inside = (inside and comp[c][0] <= vals[c]
                          <= comp[c][1])
            det = float((Q[0, 0] * Q[1, 1]).real - abs(Q[0, 1]) ** 2)
            inside = inside and rep["det"][0] <= det <= rep["det"][1]
    fails.append(not inside)
    print(f"[{'PASS' if inside else 'FAIL'}] M3 soundness : 400 pts ⊆ "
          f"(4 comp + det) @8e-3 et @1.7e-2")
    for hh, rep in reps.items():
        print(f"      h={hh:g} : {rep['status']} det ∈ "
              f"[{rep['det'][0]:.4e}, {rep['det'][1]:.4e}] "
              f"(rem {rep['det_remainder']:.2e}) "
              f"{rep['t_meas']:.1f}s")

    # --- M4 : dégénéré -------------------------------------------------------------------
    Z, W, _ = tm_chart_cell_section(S, g_col, eps, u0, v0, 0.0)
    q0 = tm_chart_metric(Z, W, M_civ, c218, rho_weight=rw)
    Gf, Gr = float_G_pair(S, g_col, eps, u0, v0, M_H, c218)
    Qf = Gf - GAMMA * Gr
    ref = [Qf[0, 0].real, Qf[1, 1].real, Qf[0, 1].real, Qf[0, 1].imag]
    mids = [sum(iv_bounds(x.to_iv())) / 2 for x in q0]
    rel = max(abs(m - r) / max(abs(r), 1e-30)
              for m, r in zip(mids, ref))
    t4 = rel < 5e-12
    fails.append(not t4)
    print(f"[{'PASS' if t4 else 'FAIL'}] M4 dégénéré h=0 ≡ float "
          f"(rel {rel:.2e})")

    # --- M5 : négatifs ---------------------------------------------------------------------
    rep_neg = tm_qfield_certificate(S, g_col, eps, u0, v0, 1e-3,
                                    M_civ, c218, -1.0)
    probe = json.loads(
        (RES / "k3_cap_b1e2iii_p0a_probe.json").read_text(
            encoding="utf-8"))
    from .interval_arithmetic import leaf_of_float_point
    b = probe["b_r_sampled"]["argmin"]
    S_w, g_w = tuple(b["S"]), int(b["g"])
    Z_w = np.array([re + 1j * im
                    for re, im in zip(b["Z"], b["Z_imag"])])
    eps_w = leaf_of_float_point(S_w, g_w, Z_w)
    others = [c for c in range(6) if c not in S_w and c != g_w]
    t_bad = 1.25 * probe["d_discriminant_target"]["t_crit_sampled"]
    rep_bad = tm_qfield_certificate(
        S_w, g_w, eps_w, complex(Z_w[others[0]]),
        complex(Z_w[others[1]]), 0.0, M_civ, c218 * t_bad, 1.0)
    t5 = rep_neg["status"] != "PASS" and rep_bad["det"][1] < 0
    fails.append(not t5)
    print(f"[{'PASS' if t5 else 'FAIL'}] M5 négatifs : w=−1 → "
          f"{rep_neg['status']} ; t_bad det.hi = "
          f"{rep_bad['det'][1]:.3e} < 0")

    # --- M6 : garde de branche --------------------------------------------------------------
    t4o = json.loads(
        (RES / "k3_cap_b1e2iii_t4o_run.json").read_text(
            encoding="utf-8"))
    c63 = t4o["r4"]["c63"]
    u_star = complex(c63["u_star"][0], c63["u_star"][1])
    side = c63["side"]

    def st(u_c, v_c):
        try:
            tm_chart_cell_section(bC["S"], bC["g"], bC["eps"], u_c,
                                  v_c, 1e-3)
            return "OK"
        except BranchCutError:
            return "BRANCH"

    cross = st(u_star, 0j)
    okside = st(complex(side["u"][0], side["u"][1]),
                complex(side["v"][0], side["v"][1]))
    t6 = cross == "BRANCH" and okside == "OK"
    fails.append(not t6)
    print(f"[{'PASS' if t6 else 'FAIL'}] M6 branche : traversée "
          f"{cross}, côté {okside}")

    # --- M7: order monotonicity (diagnostic on the determinant at hard h) ---------------
    r17 = reps[1.7e-2]
    print(f"      M7 (info) : det @1.7e-2 largeur "
          f"{r17['w_det_final']:.3e}, normes par degré "
          f"{['%.2e' % x for x in r17['det_grade_norms']]}, reste "
          f"{r17['det_remainder']:.2e}")

    print("-" * 78)
    print("SELF-TEST:", "FAIL" if any(fails) else "ALL PASS")
    return 1 if any(fails) else 0


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--selftest":
        raise SystemExit(_selftest())
    raise SystemExit("usage: k3_cap_tm_kernel.py --selftest")
