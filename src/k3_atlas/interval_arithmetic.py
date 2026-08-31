#!/usr/bin/env python3
"""
k3_cap_r4_interval_kernel.py — R4 : noyau d'arithmétique d'intervalle pour le
witness v2 (datum R3, coeffs218 natifs, moteur Kähler convention holomorphe).

Transcription intervalle FIDÈLE de `k3_cap_kahler_engine.chart_metric_kahler`
(bloc CONVENTION LOCK — chain rule holomorphe Wᵀ, potentiel
K̃ = log(Z†MZ) + Σ coeffs218·q̃_e, q̃_e ∈ B₃), sur les charts radicaux
explicites : Z_S = ε·√(A₀ + A₁u² + A₂v²), A = −V_S⁻¹V_T EXACTE en rationnels
(Vandermonde entière, μ = 1,2,3,5,7,11).

Garanties d'enclosure (chaque étape est une extension d'intervalle valide) :
 - entrées float64 convergées EXACTEMENT en mpf (53 bits ⊂ prec bits) ;
 - A rationnelle exacte (pivot de Gauss sur Fractions) ;
 - √ complexe rectangulaire par module/argument (iv.atan2), avec GARDE DE
   BRANCHE : le rectangle du radicande doit éviter la coupure (−∞, 0]
   (sinon BranchCutError — la feuille ε n'est pas continue sur la cellule) ;
 - division complexe par conj/|·|² (dénominateur > 0 exigé) ;
 - hermitisation (G+G†)/2 : enclosure valide (les deux enclosent le même G).

Le selftest (--selftest) contient les tests que la mauvaise réponse échoue :
 K2 box dégénérée (h=0) ≡ moteur float à 5e-12 relatif (multi-charts, incl.
    charts boostés des loci radicaux) — toute erreur de transcription ou de
    convention casse cette égalité ;
 K3 containment Monte-Carlo : G_float(point) ∈ G_intervalle(boîte) pour des
    points tirés DANS la boîte (toute enclosure invalide casse l'inclusion).

Consommé par : k3_cap_r4a_cell_probe.py (cellules témoins, h*, scaling h→h/2).
Witness chargé UNIQUEMENT via k3_cap_witness_registry.load_active_witness().
"""
from __future__ import annotations

import io
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np
from mpmath import iv, mp

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

from .spectral_basis import TRIPLES, sample_chart          # noqa: E402
from .witness_parametrisation import B3, B3_MULTIS, B3_IDX              # noqa: E402
from .witness_registry import load_active_witness           # noqa: E402

iv.prec = 200
mp.prec = 200

MU_INT = (1, 2, 3, 5, 7, 11)          # entiers EXACTS (LAMBDA float ailleurs)


class BranchCutError(RuntimeError):
    """Radicande dont le rectangle touche la coupure (−∞, 0] ou 0.

    C101 : l'exception porte un diagnostic STRUCTURÉ optionnel (`diag`,
    toujours un dict, vide par défaut) — site de la garde, ancre, q,
    rayon de garde, enclosure du radicande, coefficients du quadratique
    de section. Sans lui, un refus de branche n'est reconstructible que
    par un probe externe non versionné : le message seul ne dit ni QUELLE
    garde a levé ni SUR QUELS NOMBRES. Les appelants historiques
    (`raise BranchCutError("...")`, `except BranchCutError`) sont
    inchangés ; `str(exc)` reste le message.
    """

    def __init__(self, message, diag=None):
        super().__init__(message)
        self.diag = dict(diag) if diag else {}


# ===========================================================================
#  Intervalle complexe (CIV) — repris du pilote NS-1c, + sqrt/div rigoureux
# ===========================================================================
def riv(x) -> "iv.mpf":
    """Intervalle réel dégénéré EXACT depuis float64/int/Fraction."""
    if isinstance(x, Fraction):
        return iv.mpf(x.numerator) / iv.mpf(x.denominator)
    return iv.mpf(x)


class CIV:
    __slots__ = ("re", "im")

    def __init__(self, re, im=None):
        self.re = re
        self.im = im if im is not None else iv.mpf(0)

    @staticmethod
    def from_complex(c: complex):
        """Dégénéré exact (composantes float64 ⊂ mpf prec)."""
        return CIV(iv.mpf(c.real), iv.mpf(c.imag))

    @staticmethod
    def box(c: complex, h: float):
        """Boîte c ± h : bornes par addition d'intervalle (arrondi extérieur
        garanti — JAMAIS c.real ± h en float, qui peut rétrécir la boîte)."""
        radius = iv.mpf([-h, h])
        return CIV(iv.mpf(c.real) + radius, iv.mpf(c.imag) + radius)

    def __add__(a, b): return CIV(a.re + b.re, a.im + b.im)
    def __sub__(a, b): return CIV(a.re - b.re, a.im - b.im)
    def __neg__(a):    return CIV(-a.re, -a.im)

    def __mul__(a, b):
        return CIV(a.re * b.re - a.im * b.im,
                   a.re * b.im + a.im * b.re)

    def mul_real(a, r): return CIV(a.re * r, a.im * r)
    def div_real(a, r): return CIV(a.re / r, a.im / r)

    def conj(a):  return CIV(a.re, -a.im)
    def abs2(a):  return a.re ** 2 + a.im ** 2          # ** : puissance paire

    def div(a, b):
        """a/b = a·conj(b)/|b|² ; exige 0 < |b|² (borne inf > 0)."""
        d = b.abs2()
        if not (mp.mpf(d.a) > 0):
            raise BranchCutError("division complexe par un intervalle "
                                 "contenant potentiellement 0")
        return (a * b.conj()).div_real(d)

    def contains(a, c: complex) -> bool:
        return (c.real in a.re) and (c.imag in a.im)

    def mid(a) -> complex:
        return complex((mp.mpf(a.re.a) + mp.mpf(a.re.b)) / 2,
                       (mp.mpf(a.im.a) + mp.mpf(a.im.b)) / 2)

    def max_width(a) -> float:
        return float(max(mp.mpf(a.re.delta), mp.mpf(a.im.delta)))


CZERO = CIV(iv.mpf(0), iv.mpf(0))
CONE = CIV(iv.mpf(1), iv.mpf(0))


def civ_sqrt_principal(R: CIV) -> CIV:
    """√ principale d'un rectangle complexe ÉVITANT la coupure (−∞, 0].

    Garde de branche : le rectangle doit satisfaire re > 0 OU im > 0 OU
    im < 0 (strictement, sur toute la boîte). L'argument y est continu et
    iv.atan2 en donne une enclosure ; √ = |R|^{1/4}·(cos θ/2 + i sin θ/2)."""
    re_pos = mp.mpf(R.re.a) > 0
    im_pos = mp.mpf(R.im.a) > 0
    im_neg = mp.mpf(R.im.b) < 0
    if not (re_pos or im_pos or im_neg):
        raise BranchCutError(
            f"radicande sur la coupure : re=[{mp.mpf(R.re.a)}, "
            f"{mp.mpf(R.re.b)}], im=[{mp.mpf(R.im.a)}, {mp.mpf(R.im.b)}]")
    modulus4 = iv.sqrt(iv.sqrt(R.abs2()))          # |R|^{1/2} via (|R|²)^{1/4}
    half_arg = iv.atan2(R.im, R.re) / 2
    return CIV(modulus4 * iv.cos(half_arg), modulus4 * iv.sin(half_arg))


# ===========================================================================
#  Chart radical exact : A = −V_S⁻¹·V_T (rationnel), section (Z, W) intervalle
# ===========================================================================
def minor_inv_times_T_exact(S, T) -> list[list[Fraction]]:
    """A[s, t] rationnelle exacte : solve V_S·A = −V_T (Vandermonde entière)."""
    VS = [[Fraction(MU_INT[a]) ** m for a in S] for m in range(3)]
    VT = [[-Fraction(MU_INT[a]) ** m for a in T] for m in range(3)]
    # Gauss-Jordan exact sur la matrice augmentée [VS | VT]
    aug = [VS[r] + VT[r] for r in range(3)]
    for col in range(3):
        piv = next(r for r in range(col, 3) if aug[r][col] != 0)
        aug[col], aug[piv] = aug[piv], aug[col]
        p = aug[col][col]
        aug[col] = [x / p for x in aug[col]]
        for r in range(3):
            if r != col and aug[r][col] != 0:
                f = aug[r][col]
                aug[r] = [x - f * y for x, y in zip(aug[r], aug[col])]
    return [row[3:] for row in aug]                  # (3 lignes S) × (3 cols T)


def chart_cell_section(S, g_col, eps, u0: complex, v0: complex, h: float):
    """Section holomorphe intervalle du chart (S, g_col), feuille ε, sur la
    cellule (u, v) ∈ (u0 ± h) × (v0 ± h) (boîte sur les 4 dims réelles).

    Retourne Z (6 CIV), W (6×2 CIV), det_MS (CIV). Convention IDENTIQUE à
    sample_chart : Z_g = 1, Z_{o1} = u, Z_{o2} = v, Z_S = ε·√(A₀+A₁u²+A₂v²),
    W[o1,0] = W[o2,1] = 1, W[S,0] = A₁·u/Z_S, W[S,1] = A₂·v/Z_S."""
    T = tuple(j for j in range(6) if j not in S)
    others = [c for c in T if c != g_col]
    o1, o2 = others
    A_exact = minor_inv_times_T_exact(S, T)
    perm = [list(T).index(g_col), list(T).index(o1), list(T).index(o2)]
    A = [[riv(A_exact[r][c]) for c in perm] for r in range(3)]  # (S,3) réels

    u = CIV.box(u0, h)
    v = CIV.box(v0, h)
    u2, v2 = u * u, v * v

    Z = [CZERO] * 6
    W = [[CZERO, CZERO] for _ in range(6)]
    Z[g_col] = CONE
    Z[o1], Z[o2] = u, v
    W[o1][0] = CONE
    W[o2][1] = CONE
    for r, s_coord in enumerate(S):
        R = CIV(A[r][0]) + u2.mul_real(A[r][1]) + v2.mul_real(A[r][2])
        w0 = civ_sqrt_principal(R)
        Zs = w0.mul_real(riv(int(eps[r])))
        Z[s_coord] = Zs
        W[s_coord][0] = u.mul_real(A[r][1]).div(Zs)
        W[s_coord][1] = v.mul_real(A[r][2]).div(Zs)

    i_, j_, k_ = S
    VS = ((MU_INT[j_] - MU_INT[i_]) * (MU_INT[k_] - MU_INT[i_])
          * (MU_INT[k_] - MU_INT[j_]))
    det_MS = (Z[i_] * Z[j_] * Z[k_]).mul_real(riv(8 * VS))
    return Z, W, det_MS


# ===========================================================================
#  Monômes + gradients holomorphes projetés (transcription holomorphic_grads)
# ===========================================================================
def interval_monomials(Z, W, multis):
    """m[I] (CIV), p[I][α] (CIV) : p = Wᵀ·∇z^I (chain rule holomorphe)."""
    m, p = {}, {}
    for I in multis:
        cnt = Counter(I)
        val = CONE
        for o, mo in cnt.items():
            val = val * _civ_pow(Z[o], mo)
        m[I] = val
        pI = [CZERO, CZERO]
        for a, ma in cnt.items():
            gv = CIV(riv(ma))
            for o, mo in cnt.items():
                e = mo - 1 if o == a else mo
                if e:
                    gv = gv * _civ_pow(Z[o], e)
            pI[0] = pI[0] + gv * W[a][0]
            pI[1] = pI[1] + gv * W[a][1]
        p[I] = pI
    return m, p


def _civ_pow(z: CIV, k: int) -> CIV:
    r = z
    for _ in range(k - 1):
        r = r * z
    return r


# ===========================================================================
#  Métrique de chart intervalle (transcription chart_metric_kahler, K = 1)
# ===========================================================================
def interval_chart_metric(Z, W, M_civ, coeffs218, basis=B3, midx=B3_IDX,
                          multis=B3_MULTIS):
    """G packé intervalle [g00, g11, Re g01, Im g01] (4 iv réels) + s, zW.

    Miroir terme à terme de chart_metric_kahler (bloc ρ, blocs φ, termes
    croisés, hermitisation) — voir CONVENTION LOCK du moteur."""
    s = iv.mpf(0)
    for a in range(6):
        s = s + Z[a].abs2()
    zW = [CZERO, CZERO]
    for al in range(2):
        acc = CZERO
        for a in range(6):
            acc = acc + Z[a].conj() * W[a][al]
        zW[al] = acc
    WtWb = [[CZERO, CZERO], [CZERO, CZERO]]
    for A_ in range(2):
        for B_ in range(2):
            acc = CZERO
            for a in range(6):
                acc = acc + W[a][A_] * W[a][B_].conj()
            WtWb[A_][B_] = acc

    m, p = interval_monomials(Z, W, multis)

    # --- bloc rho ---------------------------------------------------------
    rho = iv.mpf(0)
    r_vec = [CZERO] * 6                                # r_a = (Z†M)_a
    for a_ in range(6):
        acc = CZERO
        for i_ in range(6):
            acc = acc + Z[i_].conj() * M_civ[i_][a_]
        r_vec[a_] = acc
        rho = rho + (acc * Z[a_]).re
    if not (mp.mpf(rho.a) > 0):
        raise BranchCutError("rho = Z†MZ non strictement positif sur la boîte")

    G = [[CZERO, CZERO], [CZERO, CZERO]]
    v_vec = [CZERO, CZERO]
    for A_ in range(2):
        acc = CZERO
        for a in range(6):
            acc = acc + W[a][A_] * r_vec[a]
        v_vec[A_] = acc
    rho2 = rho ** 2
    for A_ in range(2):
        for B_ in range(2):
            T1 = CZERO                                  # (WᵀMᵀW̄)_{AB}
            for a in range(6):
                for b in range(6):
                    T1 = T1 + W[a][A_] * M_civ[b][a] * W[b][B_].conj()
            G[A_][B_] = (T1.div_real(rho)
                         - (v_vec[A_] * v_vec[B_].conj()).div_real(rho2))

    # --- blocs phi ----------------------------------------------------------
    s3 = s * s * s
    sd = 1 / s3
    sd1 = sd / s
    sd2 = sd1 / s
    cross = [CZERO, CZERO]
    id1 = iv.mpf(0)
    id2 = iv.mpf(0)
    coeffs = np.asarray(coeffs218, float)
    for e, be in enumerate(basis):
        c = float(coeffs[e])
        if c == 0.0:
            continue
        I, Kk, typ = be["ij"], be["kl"], be["type"]
        d = len(I)                                     # = 3 sur B3
        cr = riv(c)
        mI, mK = m[I], m[Kk]
        pI, pK = p[I], p[Kk]
        if typ == "self":
            lead = [[pI[A_] * pI[B_].conj() for B_ in range(2)]
                    for A_ in range(2)] if I == Kk else None
            # self ⟹ I == Kk toujours (enumerate_sector) ; garde générale :
            if lead is None:
                lead = [[pI[A_] * pK[B_].conj() for B_ in range(2)]
                        for A_ in range(2)]
            g1 = [pI[A_] * mK.conj() for A_ in range(2)]
            phi = (mI * mK.conj()).re
        elif typ == "real_pair":
            lead = [[pI[A_] * pK[B_].conj() + pK[A_] * pI[B_].conj()
                     for B_ in range(2)] for A_ in range(2)]
            g1 = [pI[A_] * mK.conj() + pK[A_] * mI.conj() for A_ in range(2)]
            phi = 2 * (mI * mK.conj()).re
        else:                                          # imag_pair
            J = CIV(iv.mpf(0), iv.mpf(1))
            lead = [[J * (pI[A_] * pK[B_].conj() - pK[A_] * pI[B_].conj())
                     for B_ in range(2)] for A_ in range(2)]
            g1 = [J * (pI[A_] * mK.conj() - pK[A_] * mI.conj())
                  for A_ in range(2)]
            phi = -2 * (mI * mK.conj()).im
        c_sd = cr * sd                                 # c·s^{-d} (iv réel)
        c_d_sd1 = cr * sd1 * d
        for A_ in range(2):
            for B_ in range(2):
                G[A_][B_] = G[A_][B_] + lead[A_][B_].mul_real(c_sd)
            cross[A_] = cross[A_] + g1[A_].mul_real(c_d_sd1)
        id1 = id1 + c * d * phi * sd1
        id2 = id2 + c * d * (d + 1) * phi * sd2

    for A_ in range(2):
        for B_ in range(2):
            G[A_][B_] = (G[A_][B_]
                         - cross[A_] * zW[B_].conj()
                         - zW[A_] * cross[B_].conj()
                         + (zW[A_] * zW[B_].conj()).mul_real(id2)
                         - WtWb[A_][B_].mul_real(id1))

    # hermitisation (G + G†)/2 — enclosure valide du G hermitien vrai
    g00 = ((G[0][0] + G[0][0].conj()).re) / 2
    g11 = ((G[1][1] + G[1][1].conj()).re) / 2
    g01 = (G[0][1] + G[1][0].conj()).div_real(iv.mpf(2))
    return [g00, g11, g01.re, g01.im], s, zW, WtWb


def interval_fs_metric(Z, W, s, zW, WtWb):
    """G_FS packé intervalle : (WᵀW̄)/s − zW⊗z̄W/s² (contrôle + q)."""
    s2 = s ** 2
    G = [[WtWb[A_][B_].div_real(s)
          - (zW[A_] * zW[B_].conj()).div_real(s2)
          for B_ in range(2)] for A_ in range(2)]
    g00 = ((G[0][0] + G[0][0].conj()).re) / 2
    g11 = ((G[1][1] + G[1][1].conj()).re) / 2
    g01 = (G[0][1] + G[1][0].conj()).div_real(iv.mpf(2))
    return [g00, g11, g01.re, g01.im]


def det_packed_iv(g):
    return g[0] * g[1] - g[2] ** 2 - g[3] ** 2


def residual_iv(g_packed, det_MS):
    """r = log det G + 2 log|det M_S| = log det G + log |det M_S|².
    Exige det G > 0 (borne inf) — sinon None (positivité non certifiée)."""
    detG = det_packed_iv(g_packed)
    if not (mp.mpf(detG.a) > 0):
        return None, detG
    ms2 = det_MS.abs2()
    if not (mp.mpf(ms2.a) > 0):
        return None, detG
    return iv.log(detG) + iv.log(ms2), detG


def sylvester_positive(g_packed) -> bool:
    """PD certifiée sur la boîte : g00 > 0 ET det G > 0 (bornes inf)."""
    detG = det_packed_iv(g_packed)
    return (mp.mpf(g_packed[0].a) > 0) and (mp.mpf(detG.a) > 0)


def build_M_civ(M: np.ndarray) -> list[list[CIV]]:
    """M complex128 → CIV dégénérés exacts, hermitisée (M+M†)/2 en mpf."""
    Mc = [[CIV.from_complex(complex(M[i, j])) for j in range(6)]
          for i in range(6)]
    return [[(Mc[i][j] + Mc[j][i].conj()).div_real(iv.mpf(2))
             for j in range(6)] for i in range(6)]


def iv_bounds(x) -> tuple[float, float]:
    return float(mp.mpf(x.a)), float(mp.mpf(x.b))


def iv_width(x) -> float:
    return float(mp.mpf(x.delta))


# ===========================================================================
#  Couche duale (AD forward, 4 directions réelles u_re, u_im, v_re, v_im)
#  — transposition du certificateur valeur-moyenne DualCIV du mini-cover
#  (k3_cap_dualciv, C+) au moteur Kähler v2. Forme valeur-moyenne :
#     F(t) ∈ F(t₀) + Σ_a ∂F/∂t_a(boîte)·[−h, h]   (F réel, boîte convexe)
#  ⟹ le terme linéaire coûte |∇F|·h (vrai gradient O(1)) au lieu de
#  ‖coeffs218‖·h ≈ 3.6e4·h en évaluation naïve.
# ===========================================================================
NG = 4                                     # directions de dérivation


class DCIV:
    """DualCIV du kernel : (val CIV, grads[4] CIV), ε_a·ε_b = 0."""
    __slots__ = ("val", "grads")

    def __init__(self, val, grads):
        self.val = val
        self.grads = grads

    @staticmethod
    def const(val: CIV):
        return DCIV(val, [CZERO] * NG)

    def __add__(a, b): return DCIV(a.val + b.val,
                                   [x + y for x, y in zip(a.grads, b.grads)])
    def __sub__(a, b): return DCIV(a.val - b.val,
                                   [x - y for x, y in zip(a.grads, b.grads)])
    def __neg__(a):    return DCIV(-a.val, [-x for x in a.grads])

    def __mul__(a, b):
        return DCIV(a.val * b.val,
                    [ga * b.val + a.val * gb
                     for ga, gb in zip(a.grads, b.grads)])

    def mul_real(a, r):                     # r : iv const
        return DCIV(a.val.mul_real(r), [g.mul_real(r) for g in a.grads])

    def conj(a):
        return DCIV(a.val.conj(), [g.conj() for g in a.grads])

    def inv(a):
        """1/a : val par conj/|·|², grads = −da·(1/a)² (enclosures valides)."""
        iv_val = CONE.div(a.val)
        m2 = iv_val * iv_val
        return DCIV(iv_val, [-(g * m2) for g in a.grads])

    def div(a, b):
        return a * b.inv()

    def mul_rdual(a, r: "DIV"):             # r : DualIV réel
        return DCIV(a.val.mul_real(r.val),
                    [g.mul_real(r.val) + a.val.mul_real(gr)
                     for g, gr in zip(a.grads, r.grads)])

    def div_rdual(a, r: "DIV"):
        return a.mul_rdual(r.inv())

    def re_dual(a) -> "DIV":
        return DIV(a.val.re, [g.re for g in a.grads])

    def im_dual(a) -> "DIV":
        return DIV(a.val.im, [g.im for g in a.grads])

    def abs2_dual(a) -> "DIV":
        """|a|² : d = 2 Re(ā·da)."""
        return DIV(a.val.abs2(),
                   [2 * (a.val.conj() * g).re for g in a.grads])


class DIV:
    """DualIV du kernel : (val iv réel, grads[4] iv réels)."""
    __slots__ = ("val", "grads")

    def __init__(self, val, grads):
        self.val = val
        self.grads = grads

    @staticmethod
    def const(val):
        return DIV(val, [iv.mpf(0)] * NG)

    def __add__(a, b):
        if isinstance(b, DIV):
            return DIV(a.val + b.val,
                       [x + y for x, y in zip(a.grads, b.grads)])
        return DIV(a.val + b, list(a.grads))

    def __sub__(a, b):
        if isinstance(b, DIV):
            return DIV(a.val - b.val,
                       [x - y for x, y in zip(a.grads, b.grads)])
        return DIV(a.val - b, list(a.grads))

    def __mul__(a, b):
        if isinstance(b, DIV):
            return DIV(a.val * b.val,
                       [x * b.val + a.val * y
                        for x, y in zip(a.grads, b.grads)])
        return DIV(a.val * b, [x * b for x in a.grads])   # b : iv/float const

    __rmul__ = __mul__

    def __neg__(a):
        return DIV(-a.val, [-x for x in a.grads])

    def inv(a):
        iv_inv = 1 / a.val
        m2 = iv_inv * iv_inv
        return DIV(iv_inv, [-(g * m2) for g in a.grads])


def _dciv_pow(z: DCIV, k: int) -> DCIV:
    r = z
    for _ in range(k - 1):
        r = r * z
    return r


def dual_chart_cell_section(S, g_col, eps, u0: complex, v0: complex, h: float):
    """Section duale : Z (6 DCIV), W (6×2 DCIV), det_MS (DCIV), seeds
    ∂/∂(u_re, u_im, v_re, v_im). Miroir strict de chart_cell_section."""
    T = tuple(j for j in range(6) if j not in S)
    others = [c for c in T if c != g_col]
    o1, o2 = others
    A_exact = minor_inv_times_T_exact(S, T)
    perm = [list(T).index(g_col), list(T).index(o1), list(T).index(o2)]
    A = [[riv(A_exact[r][c]) for c in perm] for r in range(3)]

    J = CIV(iv.mpf(0), iv.mpf(1))
    u = DCIV(CIV.box(u0, h), [CONE, J, CZERO, CZERO])
    v = DCIV(CIV.box(v0, h), [CZERO, CZERO, CONE, J])
    u2, v2 = u * u, v * v

    ZD = [DCIV.const(CZERO)] * 6
    WD = [[DCIV.const(CZERO), DCIV.const(CZERO)] for _ in range(6)]
    ZD = list(ZD)
    ZD[g_col] = DCIV.const(CONE)
    ZD[o1], ZD[o2] = u, v
    WD[o1][0] = DCIV.const(CONE)
    WD[o2][1] = DCIV.const(CONE)
    for r, s_coord in enumerate(S):
        R = DCIV.const(CIV(A[r][0])) + u2.mul_real(A[r][1]) \
            + v2.mul_real(A[r][2])
        w_val = civ_sqrt_principal(R.val)
        half_inv = CONE.div(w_val.mul_real(iv.mpf(2)))   # 1/(2√R) enclosure
        w = DCIV(w_val, [g * half_inv for g in R.grads])
        Zs = w.mul_real(riv(int(eps[r])))
        ZD[s_coord] = Zs
        WD[s_coord][0] = u.mul_real(A[r][1]).div(Zs)
        WD[s_coord][1] = v.mul_real(A[r][2]).div(Zs)

    i_, j_, k_ = S
    VS = ((MU_INT[j_] - MU_INT[i_]) * (MU_INT[k_] - MU_INT[i_])
          * (MU_INT[k_] - MU_INT[j_]))
    det_MS = (ZD[i_] * ZD[j_] * ZD[k_]).mul_real(riv(8 * VS))
    return ZD, WD, det_MS


def dual_interval_monomials(Z, W, multis):
    m, p = {}, {}
    one = DCIV.const(CONE)
    for I in multis:
        cnt = Counter(I)
        val = one
        for o, mo in cnt.items():
            val = val * _dciv_pow(Z[o], mo)
        m[I] = val
        pI = [DCIV.const(CZERO), DCIV.const(CZERO)]
        for a, ma in cnt.items():
            gv = DCIV.const(CIV(riv(ma)))
            for o, mo in cnt.items():
                e = mo - 1 if o == a else mo
                if e:
                    gv = gv * _dciv_pow(Z[o], e)
            pI[0] = pI[0] + gv * W[a][0]
            pI[1] = pI[1] + gv * W[a][1]
        p[I] = pI
    return m, p


def dual_chart_metric(Z, W, M_civ, coeffs218, basis=B3, midx=B3_IDX,
                      multis=B3_MULTIS):
    """G packé DUAL [g00, g11, Re g01, Im g01] (4 DIV) — miroir strict de
    interval_chart_metric avec propagation des 4 dérivées."""
    s = DIV.const(iv.mpf(0))
    for a in range(6):
        s = s + Z[a].abs2_dual()
    zW = []
    for al in range(2):
        acc = DCIV.const(CZERO)
        for a in range(6):
            acc = acc + Z[a].conj() * W[a][al]
        zW.append(acc)
    WtWb = [[None, None], [None, None]]
    for A_ in range(2):
        for B_ in range(2):
            acc = DCIV.const(CZERO)
            for a in range(6):
                acc = acc + W[a][A_] * W[a][B_].conj()
            WtWb[A_][B_] = acc

    m, p = dual_interval_monomials(Z, W, multis)

    Md = [[DCIV.const(M_civ[i][j]) for j in range(6)] for i in range(6)]
    rho = DIV.const(iv.mpf(0))
    r_vec = []
    for a_ in range(6):
        acc = DCIV.const(CZERO)
        for i_ in range(6):
            acc = acc + Z[i_].conj() * Md[i_][a_]
        r_vec.append(acc)
        rho = rho + (acc * Z[a_]).re_dual()
    if not (mp.mpf(rho.val.a) > 0):
        raise BranchCutError("rho non strictement positif (dual)")

    v_vec = []
    for A_ in range(2):
        acc = DCIV.const(CZERO)
        for a in range(6):
            acc = acc + W[a][A_] * r_vec[a]
        v_vec.append(acc)
    rho_inv = rho.inv()
    rho2_inv = rho_inv * rho_inv
    G = [[None, None], [None, None]]
    for A_ in range(2):
        for B_ in range(2):
            T1 = DCIV.const(CZERO)
            for a in range(6):
                for b in range(6):
                    T1 = T1 + W[a][A_] * Md[b][a] * W[b][B_].conj()
            G[A_][B_] = (T1.mul_rdual(rho_inv)
                         - (v_vec[A_] * v_vec[B_].conj()).mul_rdual(rho2_inv))

    s_inv = s.inv()
    s3_inv = s_inv * s_inv * s_inv
    sd1 = s3_inv * s_inv
    sd2 = sd1 * s_inv
    cross = [DCIV.const(CZERO), DCIV.const(CZERO)]
    id1 = DIV.const(iv.mpf(0))
    id2 = DIV.const(iv.mpf(0))
    coeffs = np.asarray(coeffs218, float)
    for e, be in enumerate(basis):
        c = float(coeffs[e])
        if c == 0.0:
            continue
        I, Kk, typ = be["ij"], be["kl"], be["type"]
        d = len(I)
        cr = riv(c)
        mI, mK = m[I], m[Kk]
        pI, pK = p[I], p[Kk]
        if typ == "self":
            lead = [[pI[A_] * pI[B_].conj() for B_ in range(2)]
                    for A_ in range(2)]
            g1 = [pI[A_] * mK.conj() for A_ in range(2)]
            phi = (mI * mK.conj()).re_dual()
        elif typ == "real_pair":
            lead = [[pI[A_] * pK[B_].conj() + pK[A_] * pI[B_].conj()
                     for B_ in range(2)] for A_ in range(2)]
            g1 = [pI[A_] * mK.conj() + pK[A_] * mI.conj() for A_ in range(2)]
            phi = (mI * mK.conj()).re_dual() * 2
        else:
            Jc = DCIV.const(CIV(iv.mpf(0), iv.mpf(1)))
            lead = [[Jc * (pI[A_] * pK[B_].conj() - pK[A_] * pI[B_].conj())
                     for B_ in range(2)] for A_ in range(2)]
            g1 = [Jc * (pI[A_] * mK.conj() - pK[A_] * mI.conj())
                  for A_ in range(2)]
            phi = (mI * mK.conj()).im_dual() * (-2)
        c_sd = s3_inv * cr                          # DIV
        c_d_sd1 = sd1 * (cr * d)
        for A_ in range(2):
            for B_ in range(2):
                G[A_][B_] = G[A_][B_] + lead[A_][B_].mul_rdual(c_sd)
            cross[A_] = cross[A_] + g1[A_].mul_rdual(c_d_sd1)
        id1 = id1 + phi * sd1 * (cr * d)
        id2 = id2 + phi * sd2 * (cr * (d * (d + 1)))

    for A_ in range(2):
        for B_ in range(2):
            G[A_][B_] = (G[A_][B_]
                         - cross[A_] * zW[B_].conj()
                         - zW[A_] * cross[B_].conj()
                         + (zW[A_] * zW[B_].conj()).mul_rdual(id2)
                         - WtWb[A_][B_].mul_rdual(id1))

    half = iv.mpf(1) / 2
    g00 = ((G[0][0] + G[0][0].conj()).re_dual()) * half
    g11 = ((G[1][1] + G[1][1].conj()).re_dual()) * half
    g01 = (G[0][1] + G[1][0].conj()).mul_real(half)
    return [g00, g11, g01.re_dual(), g01.im_dual()]


def mean_value_enclose(center_iv, dual: DIV, h: float):
    """Forme valeur-moyenne : F(boîte) ⊆ F(centre) + Σ_a ∂F(boîte)·[−h,h].
    center_iv : enclosure de F au CENTRE (boîte dégénérée). Retourne iv."""
    rad = iv.mpf([-h, h])
    acc = center_iv
    for g in dual.grads:
        acc = acc + g * rad
    return acc


def mean_value_metric(S, g_col, eps, u0, v0, h, M_civ, coeffs218):
    """G packé (4 iv réels) + det_MS (CIV) par forme valeur-moyenne,
    INTERSECTÉE avec l'évaluation naïve (les deux enclosent G ⟹ ∩ valide).

    Retourne (g_packed, det_MS, diag) — diag contient les largeurs des
    deux formes pour le probe."""
    # centre : boîte dégénérée (naïf, largeur ~ arrondi seulement)
    Zc, Wc, dMSc = chart_cell_section(S, g_col, eps, u0, v0, 0.0)
    g_center, *_ = interval_chart_metric(Zc, Wc, M_civ, coeffs218)
    # duals sur la boîte pleine
    ZD, WD, dMS_dual = dual_chart_cell_section(S, g_col, eps, u0, v0, h)
    g_dual = dual_chart_metric(ZD, WD, M_civ, coeffs218)
    # naïf sur la boîte pleine (pour intersection + comparaison)
    Zb, Wb, dMS_box = chart_cell_section(S, g_col, eps, u0, v0, h)
    g_naive, *_ = interval_chart_metric(Zb, Wb, M_civ, coeffs218)

    g_mv, w_mv, w_naive = [], [], []
    for c in range(4):
        mv = mean_value_enclose(g_center[c], g_dual[c], h)
        w_mv.append(float(mp.mpf(mv.delta)))
        w_naive.append(iv_width(g_naive[c]))
        inter = _iv_intersect(mv, g_naive[c])
        g_mv.append(inter)
    # det_MS : valeur-moyenne composante par composante, ∩ naïf
    dMS_center_re, dMS_center_im = dMSc.re, dMSc.im
    ms_re = _iv_intersect(
        mean_value_enclose(dMS_center_re, dMS_dual.re_dual(), h), dMS_box.re)
    ms_im = _iv_intersect(
        mean_value_enclose(dMS_center_im, dMS_dual.im_dual(), h), dMS_box.im)
    det_MS = CIV(ms_re, ms_im)
    diag = {"w_mv": w_mv, "w_naive": w_naive}
    return g_mv, det_MS, diag


def _iv_intersect(a, b):
    """Intersection de deux enclosures du même réel (toujours non vide si
    les deux sont valides ; sinon on le SAURA — erreur du kernel)."""
    lo = max(mp.mpf(a.a), mp.mpf(b.a))
    hi = min(mp.mpf(a.b), mp.mpf(b.b))
    if lo > hi:
        raise RuntimeError("intersection vide : une des deux enclosures "
                           "est invalide (bug kernel)")
    return iv.mpf([lo, hi])


# ===========================================================================
#  Couche ordre 2 (jets tronqués : val + 4 gradients + 10 hessiennes)
#  — R4-A passe 2. Forme de Taylor ordre 2 :
#     F(t₀+δ) = F(t₀) + ∇F(t₀)·δ + ½ δᵀ·H(ξ)·δ,  ξ ∈ boîte (Lagrange)
#  Le terme linéaire prend le gradient EXACT au centre (largeur ~arrondi) ;
#  seule la Hessienne est enclosée naïvement sur la boîte ⟹ la cancellation
#  ‖coeffs218‖ est reléguée au terme h² (avec constante vraie |∇²F|) + h³.
# ===========================================================================
HPAIRS = [(0, 0), (0, 1), (0, 2), (0, 3), (1, 1),
          (1, 2), (1, 3), (2, 2), (2, 3), (3, 3)]
NH = len(HPAIRS)


class T2CIV:
    """Jet ordre 2 complexe : (val CIV, g[4] CIV, h[10] CIV, paires a≤b)."""
    __slots__ = ("val", "g", "h")

    def __init__(self, val, g, h):
        self.val = val
        self.g = g
        self.h = h

    @staticmethod
    def const(val: CIV):
        return T2CIV(val, [CZERO] * NG, [CZERO] * NH)

    def __add__(a, b): return T2CIV(a.val + b.val,
                                    [x + y for x, y in zip(a.g, b.g)],
                                    [x + y for x, y in zip(a.h, b.h)])
    def __sub__(a, b): return T2CIV(a.val - b.val,
                                    [x - y for x, y in zip(a.g, b.g)],
                                    [x - y for x, y in zip(a.h, b.h)])
    def __neg__(a):    return T2CIV(-a.val, [-x for x in a.g],
                                    [-x for x in a.h])

    def __mul__(f, g):
        val = f.val * g.val
        gr = [f.g[a] * g.val + f.val * g.g[a] for a in range(NG)]
        hs = [f.h[i] * g.val + f.val * g.h[i]
              + f.g[a] * g.g[b] + f.g[b] * g.g[a]
              for i, (a, b) in enumerate(HPAIRS)]
        return T2CIV(val, gr, hs)

    def mul_real(f, r):
        return T2CIV(f.val.mul_real(r), [x.mul_real(r) for x in f.g],
                     [x.mul_real(r) for x in f.h])

    def conj(f):
        return T2CIV(f.val.conj(), [x.conj() for x in f.g],
                     [x.conj() for x in f.h])

    def inv(f):
        """1/f : u_a = −f_a u², u_ab = −f_ab u² + 2 f_a f_b u³."""
        u = CONE.div(f.val)
        u2 = u * u
        u3 = u2 * u
        gr = [-(f.g[a] * u2) for a in range(NG)]
        hs = [-(f.h[i] * u2) + (f.g[a] * f.g[b] * u3).mul_real(iv.mpf(2))
              for i, (a, b) in enumerate(HPAIRS)]
        return T2CIV(u, gr, hs)

    def div(f, g):
        return f * g.inv()

    def mul_rt2(f, r: "T2IV"):
        """f · r, r jet réel (règle produit complète)."""
        val = f.val.mul_real(r.val)
        gr = [f.g[a].mul_real(r.val) + f.val.mul_real(r.g[a])
              for a in range(NG)]
        hs = [f.h[i].mul_real(r.val) + f.val.mul_real(r.h[i])
              + f.g[a].mul_real(r.g[b]) + f.g[b].mul_real(r.g[a])
              for i, (a, b) in enumerate(HPAIRS)]
        return T2CIV(val, gr, hs)

    def re_t2(f) -> "T2IV":
        return T2IV(f.val.re, [x.re for x in f.g], [x.re for x in f.h])

    def im_t2(f) -> "T2IV":
        return T2IV(f.val.im, [x.im for x in f.g], [x.im for x in f.h])

    def abs2_t2(f) -> "T2IV":
        return (f * f.conj()).re_t2()

    def sqrt_principal(f) -> "T2CIV":
        """w = √f (branche principale, garde de branche sur f.val) :
        w_a = f_a/(2w), w_ab = (f_ab − 2 w_a w_b)/(2w)."""
        w = civ_sqrt_principal(f.val)
        inv2w = CONE.div(w.mul_real(iv.mpf(2)))
        gr = [f.g[a] * inv2w for a in range(NG)]
        hs = [(f.h[i] - (gr[a] * gr[b]).mul_real(iv.mpf(2))) * inv2w
              for i, (a, b) in enumerate(HPAIRS)]
        return T2CIV(w, gr, hs)


class T2IV:
    """Jet ordre 2 réel : (val iv, g[4] iv, h[10] iv)."""
    __slots__ = ("val", "g", "h")

    def __init__(self, val, g, h):
        self.val = val
        self.g = g
        self.h = h

    @staticmethod
    def const(val):
        return T2IV(val, [iv.mpf(0)] * NG, [iv.mpf(0)] * NH)

    def __add__(a, b):
        if isinstance(b, T2IV):
            return T2IV(a.val + b.val, [x + y for x, y in zip(a.g, b.g)],
                        [x + y for x, y in zip(a.h, b.h)])
        return T2IV(a.val + b, list(a.g), list(a.h))

    def __sub__(a, b):
        if isinstance(b, T2IV):
            return T2IV(a.val - b.val, [x - y for x, y in zip(a.g, b.g)],
                        [x - y for x, y in zip(a.h, b.h)])
        return T2IV(a.val - b, list(a.g), list(a.h))

    def __mul__(f, g):
        if isinstance(g, T2IV):
            val = f.val * g.val
            gr = [f.g[a] * g.val + f.val * g.g[a] for a in range(NG)]
            hs = [f.h[i] * g.val + f.val * g.h[i]
                  + f.g[a] * g.g[b] + f.g[b] * g.g[a]
                  for i, (a, b) in enumerate(HPAIRS)]
            return T2IV(val, gr, hs)
        return T2IV(f.val * g, [x * g for x in f.g], [x * g for x in f.h])

    __rmul__ = __mul__

    def __neg__(f):
        return T2IV(-f.val, [-x for x in f.g], [-x for x in f.h])

    def inv(f):
        u = 1 / f.val
        u2 = u * u
        u3 = u2 * u
        gr = [-(f.g[a] * u2) for a in range(NG)]
        hs = [-(f.h[i] * u2) + 2 * f.g[a] * f.g[b] * u3
              for i, (a, b) in enumerate(HPAIRS)]
        return T2IV(u, gr, hs)


def t2_chart_cell_section(S, g_col, eps, u0: complex, v0: complex, h: float):
    """Section jet ordre 2 : Z (6 T2CIV), W (6×2), det_MS (T2CIV).
    u, v sont LINÉAIRES dans les 4 params ⟹ leurs hessiennes sont nulles."""
    T = tuple(j for j in range(6) if j not in S)
    others = [c for c in T if c != g_col]
    o1, o2 = others
    A_exact = minor_inv_times_T_exact(S, T)
    perm = [list(T).index(g_col), list(T).index(o1), list(T).index(o2)]
    A = [[riv(A_exact[r][c]) for c in perm] for r in range(3)]

    J = CIV(iv.mpf(0), iv.mpf(1))
    u = T2CIV(CIV.box(u0, h), [CONE, J, CZERO, CZERO], [CZERO] * NH)
    v = T2CIV(CIV.box(v0, h), [CZERO, CZERO, CONE, J], [CZERO] * NH)
    u2, v2 = u * u, v * v

    ZT = [T2CIV.const(CZERO) for _ in range(6)]
    WT = [[T2CIV.const(CZERO), T2CIV.const(CZERO)] for _ in range(6)]
    ZT[g_col] = T2CIV.const(CONE)
    ZT[o1], ZT[o2] = u, v
    WT[o1][0] = T2CIV.const(CONE)
    WT[o2][1] = T2CIV.const(CONE)
    for r, s_coord in enumerate(S):
        R = T2CIV.const(CIV(A[r][0])) + u2.mul_real(A[r][1]) \
            + v2.mul_real(A[r][2])
        Zs = R.sqrt_principal().mul_real(riv(int(eps[r])))
        ZT[s_coord] = Zs
        WT[s_coord][0] = u.mul_real(A[r][1]).div(Zs)
        WT[s_coord][1] = v.mul_real(A[r][2]).div(Zs)

    i_, j_, k_ = S
    VS = ((MU_INT[j_] - MU_INT[i_]) * (MU_INT[k_] - MU_INT[i_])
          * (MU_INT[k_] - MU_INT[j_]))
    det_MS = (ZT[i_] * ZT[j_] * ZT[k_]).mul_real(riv(8 * VS))
    return ZT, WT, det_MS


def _t2_pow(z: T2CIV, k: int) -> T2CIV:
    r = z
    for _ in range(k - 1):
        r = r * z
    return r


def t2_chart_metric(Z, W, M_civ, coeffs218, basis=B3, midx=B3_IDX,
                    multis=B3_MULTIS, want_fs=False, rho_weight=None):
    """G packé jet ordre 2 [g00, g11, Re g01, Im g01] (4 T2IV) — miroir
    strict de dual_chart_metric. want_fs=True : retourne (G, G_FS) packés
    jets — G_FS réutilise s/zW/WᵀW̄ des MÊMES jets, donc le pinceau
    H_α = G − α·G_FS (R4-B0) garde les corrélations entre les deux.

    rho_weight (P0a-2, review GPT p0a_probe §5) : si non-None, le bloc ρ
    (pullback G_ρ = ∂∂̄ log ρ) est multiplié par ce scalaire AVANT l'ajout
    du bloc φ — l'assemblage retourne alors le champ COMBINÉ
    Q = rho_weight·G_ρ + H_Φ en un seul jet (cancellations préservées ;
    avec rho_weight = 1−γ c'est le certificat Loewner direct
    G ⪰ γ·G_ρ ⟺ Q ≻ 0, sans inverse ni racine matricielle).
    None (défaut) : comportement STRICTEMENT identique à avant."""
    s = T2IV.const(iv.mpf(0))
    for a in range(6):
        s = s + Z[a].abs2_t2()
    zW = []
    for al in range(2):
        acc = T2CIV.const(CZERO)
        for a in range(6):
            acc = acc + Z[a].conj() * W[a][al]
        zW.append(acc)
    WtWb = [[None, None], [None, None]]
    for A_ in range(2):
        for B_ in range(2):
            acc = T2CIV.const(CZERO)
            for a in range(6):
                acc = acc + W[a][A_] * W[a][B_].conj()
            WtWb[A_][B_] = acc

    m, p = {}, {}
    one = T2CIV.const(CONE)
    for I in multis:
        cnt = Counter(I)
        val = one
        for o, mo in cnt.items():
            val = val * _t2_pow(Z[o], mo)
        m[I] = val
        pI = [T2CIV.const(CZERO), T2CIV.const(CZERO)]
        for a, ma in cnt.items():
            gv = T2CIV.const(CIV(riv(ma)))
            for o, mo in cnt.items():
                e = mo - 1 if o == a else mo
                if e:
                    gv = gv * _t2_pow(Z[o], e)
            pI[0] = pI[0] + gv * W[a][0]
            pI[1] = pI[1] + gv * W[a][1]
        p[I] = pI

    Md = [[T2CIV.const(M_civ[i][j]) for j in range(6)] for i in range(6)]
    rho = T2IV.const(iv.mpf(0))
    r_vec = []
    for a_ in range(6):
        acc = T2CIV.const(CZERO)
        for i_ in range(6):
            acc = acc + Z[i_].conj() * Md[i_][a_]
        r_vec.append(acc)
        rho = rho + (acc * Z[a_]).re_t2()
    if not (mp.mpf(rho.val.a) > 0):
        raise BranchCutError("rho non strictement positif (t2)")

    v_vec = []
    for A_ in range(2):
        acc = T2CIV.const(CZERO)
        for a in range(6):
            acc = acc + W[a][A_] * r_vec[a]
        v_vec.append(acc)
    rho_inv = rho.inv()
    rho2_inv = rho_inv * rho_inv
    if rho_weight is not None:
        w_rt2 = T2IV.const(riv(rho_weight))
        rho_inv = rho_inv * w_rt2          # poids appliqué UNE fois
        rho2_inv = rho2_inv * w_rt2        # (pas au carré)
    G = [[None, None], [None, None]]
    for A_ in range(2):
        for B_ in range(2):
            T1 = T2CIV.const(CZERO)
            for a in range(6):
                for b in range(6):
                    T1 = T1 + W[a][A_] * Md[b][a] * W[b][B_].conj()
            G[A_][B_] = (T1.mul_rt2(rho_inv)
                         - (v_vec[A_] * v_vec[B_].conj()).mul_rt2(rho2_inv))

    s_inv = s.inv()
    s3_inv = s_inv * s_inv * s_inv
    sd1 = s3_inv * s_inv
    sd2 = sd1 * s_inv
    cross = [T2CIV.const(CZERO), T2CIV.const(CZERO)]
    id1 = T2IV.const(iv.mpf(0))
    id2 = T2IV.const(iv.mpf(0))
    coeffs = np.asarray(coeffs218, float)
    for e, be in enumerate(basis):
        c = float(coeffs[e])
        if c == 0.0:
            continue
        I, Kk, typ = be["ij"], be["kl"], be["type"]
        d = len(I)
        cr = riv(c)
        mI, mK = m[I], m[Kk]
        pI, pK = p[I], p[Kk]
        if typ == "self":
            lead = [[pI[A_] * pI[B_].conj() for B_ in range(2)]
                    for A_ in range(2)]
            g1 = [pI[A_] * mK.conj() for A_ in range(2)]
            phi = (mI * mK.conj()).re_t2()
        elif typ == "real_pair":
            lead = [[pI[A_] * pK[B_].conj() + pK[A_] * pI[B_].conj()
                     for B_ in range(2)] for A_ in range(2)]
            g1 = [pI[A_] * mK.conj() + pK[A_] * mI.conj() for A_ in range(2)]
            phi = (mI * mK.conj()).re_t2() * 2
        else:
            Jc = T2CIV.const(CIV(iv.mpf(0), iv.mpf(1)))
            lead = [[Jc * (pI[A_] * pK[B_].conj() - pK[A_] * pI[B_].conj())
                     for B_ in range(2)] for A_ in range(2)]
            g1 = [Jc * (pI[A_] * mK.conj() - pK[A_] * mI.conj())
                  for A_ in range(2)]
            phi = (mI * mK.conj()).im_t2() * (-2)
        c_sd = s3_inv * cr
        c_d_sd1 = sd1 * (cr * d)
        for A_ in range(2):
            for B_ in range(2):
                G[A_][B_] = G[A_][B_] + lead[A_][B_].mul_rt2(c_sd)
            cross[A_] = cross[A_] + g1[A_].mul_rt2(c_d_sd1)
        id1 = id1 + phi * sd1 * (cr * d)
        id2 = id2 + phi * sd2 * (cr * (d * (d + 1)))

    for A_ in range(2):
        for B_ in range(2):
            G[A_][B_] = (G[A_][B_]
                         - cross[A_] * zW[B_].conj()
                         - zW[A_] * cross[B_].conj()
                         + (zW[A_] * zW[B_].conj()).mul_rt2(id2)
                         - WtWb[A_][B_].mul_rt2(id1))

    half = iv.mpf(1) / 2
    g00 = ((G[0][0] + G[0][0].conj()).re_t2()) * half
    g11 = ((G[1][1] + G[1][1].conj()).re_t2()) * half
    g01 = (G[0][1] + G[1][0].conj()).mul_real(half)
    g_packed = [g00, g11, g01.re_t2(), g01.im_t2()]
    if not want_fs:
        return g_packed
    s2_inv = s_inv * s_inv
    Gfs = [[WtWb[A_][B_].mul_rt2(s_inv)
            - (zW[A_] * zW[B_].conj()).mul_rt2(s2_inv)
            for B_ in range(2)] for A_ in range(2)]
    f00 = ((Gfs[0][0] + Gfs[0][0].conj()).re_t2()) * half
    f11 = ((Gfs[1][1] + Gfs[1][1].conj()).re_t2()) * half
    f01 = (Gfs[0][1] + Gfs[1][0].conj()).mul_real(half)
    return g_packed, [f00, f11, f01.re_t2(), f01.im_t2()]


def det_packed_t2(g):
    return g[0] * g[1] - g[2] * g[2] - g[3] * g[3]


def taylor2_enclose(center: T2IV, box: T2IV, h: float):
    """Enclosure ordre 2 : val(t₀) + ∇(t₀)·δ + ½δᵀH(boîte)δ, ∩ MV, ∩ naïf.

    center : jet évalué en boîte DÉGÉNÉRÉE (val/grads serrés) ;
    box    : jet évalué sur la boîte pleine (val naïve, grads MV, hess).
    Les trois formes enclosent F(boîte) ⟹ l'intersection est valide."""
    rad = iv.mpf([-h, h])
    sq_diag = iv.mpf([0, h * h])
    sq_off = iv.mpf([-h * h, h * h])
    half = iv.mpf(1) / 2
    # forme de Taylor ordre 2
    t2 = center.val
    for a in range(NG):
        t2 = t2 + center.g[a] * rad
    for i, (a, b) in enumerate(HPAIRS):
        if a == b:
            t2 = t2 + box.h[i] * sq_diag * half
        else:
            t2 = t2 + box.h[i] * sq_off
    # forme valeur-moyenne (grads sur la boîte)
    mv = center.val
    for a in range(NG):
        mv = mv + box.g[a] * rad
    return _iv_intersect(_iv_intersect(t2, mv), box.val)


def taylor2_metric(S, g_col, eps, u0, v0, h, M_civ, coeffs218):
    """G packé (4 iv), det G (iv), det_MS (CIV), |det_MS|² (iv) par forme de
    Taylor ordre 2 (∩ MV ∩ naïf, composante par composante ET sur det G
    directement — le jet de det G garde les corrélations entre composantes).

    diag : largeurs des 3 formes pour g00 et det G."""
    Zc, Wc, dMSc = t2_chart_cell_section(S, g_col, eps, u0, v0, 0.0)
    g_c = t2_chart_metric(Zc, Wc, M_civ, coeffs218)
    det_c = det_packed_t2(g_c)
    ms2_c = dMSc.abs2_t2()

    Zb, Wb, dMSb = t2_chart_cell_section(S, g_col, eps, u0, v0, h)
    g_b = t2_chart_metric(Zb, Wb, M_civ, coeffs218)
    det_b = det_packed_t2(g_b)
    ms2_b = dMSb.abs2_t2()

    g_packed = [taylor2_enclose(g_c[c], g_b[c], h) for c in range(4)]
    detG = _iv_intersect(taylor2_enclose(det_c, det_b, h),
                         det_packed_iv(g_packed))
    ms2 = taylor2_enclose(ms2_c, ms2_b, h)
    det_MS = CIV(taylor2_enclose(dMSc.re_t2(), dMSb.re_t2(), h),
                 taylor2_enclose(dMSc.im_t2(), dMSb.im_t2(), h))
    diag = {
        "w_t2_g00": float(mp.mpf(g_packed[0].delta)),
        "w_naive_g00": iv_width(g_b[0].val),
        "w_mv_g00": None,   # inclus dans l'intersection ; suivi via passes 1
        "w_t2_detG": float(mp.mpf(detG.delta)),
        "w_naive_detG": iv_width(det_b.val),
    }
    return g_packed, detG, det_MS, ms2, diag


def residual_t2(detG, ms2):
    """r = log det G + log |det M_S|² depuis les enclosures ordre 2."""
    if not (mp.mpf(detG.a) > 0 and mp.mpf(ms2.a) > 0):
        return None
    return iv.log(detG) + iv.log(ms2)


# ===========================================================================
#  Récupération de la feuille ε d'un point float du sampler
# ===========================================================================
def leaf_of_float_point(S, g_col, Z_float: np.ndarray) -> tuple[int, ...]:
    """Retrouve ε ∈ {±1}³ : Z_S = ε·√(R) branche principale (float)."""
    T = tuple(j for j in range(6) if j not in S)
    others = [c for c in T if c != g_col]
    o1, o2 = others
    A_exact = minor_inv_times_T_exact(S, T)
    perm = [list(T).index(g_col), list(T).index(o1), list(T).index(o2)]
    u, v = complex(Z_float[o1]), complex(Z_float[o2])
    eps = []
    for r, s_coord in enumerate(S):
        a0, a1, a2 = (float(A_exact[r][c]) for c in perm)
        w0 = np.sqrt(complex(a0 + a1 * u * u + a2 * v * v))
        ratio = complex(Z_float[s_coord]) / w0
        eps.append(1 if ratio.real > 0 else -1)
    return tuple(eps)


# ===========================================================================
#  Selftest
# ===========================================================================
def _selftest() -> int:
    import time
    from .kahler_metric import chart_metric_kahler, fs_pullback
    from .ricci_functional import pack_herm
    from .spectral_basis import detMS_on_block

    witness = load_active_witness()
    coeffs218 = np.asarray(witness["coeffs218"], float)
    M = np.asarray(witness["M"], complex)
    M_civ = build_M_civ(M)
    failures = []

    def check(name, ok, detail=""):
        print(f"[{'PASS' if ok else 'FAIL'}] {name}" +
              (f" — {detail}" if detail else ""), flush=True)
        if not ok:
            failures.append(name)

    # K1 — A exacte vs A float du sampler
    from .spectral_basis import minor_inv_times_T_float
    worst = 0.0
    for S in TRIPLES:
        T = tuple(j for j in range(6) if j not in S)
        A_f = minor_inv_times_T_float(S, T)
        A_e = minor_inv_times_T_exact(S, T)
        diff = max(abs(A_f[r, c] - float(A_e[r][c]))
                   for r in range(3) for c in range(3))
        worst = max(worst, diff)
    check("K1_exact_A_vs_float", worst < 1e-11, f"max|ΔA| = {worst:.2e}")

    # charts témoins : 1 générique + les 2 boostés (loci radicaux)
    cases = [((0, 1, 2), 3), ((1, 3, 5), 2), ((0, 2, 5), 1)]
    rng = np.random.default_rng(2026)
    t0 = time.time()
    worst_rel, worst_r = 0.0, 0.0
    n_pts = 0
    samples = []
    for S, g_col in cases:
        res = sample_chart(rng, S, g_col, 400)
        if res is None:
            continue
        Z_f, W_f, UV = res
        idx = rng.choice(Z_f.shape[0], size=min(4, Z_f.shape[0]),
                         replace=False)
        for i in idx:
            eps = leaf_of_float_point(S, g_col, Z_f[i])
            u0, v0 = complex(UV[i, 0]), complex(UV[i, 1])
            samples.append((S, g_col, eps, u0, v0))
            # K2 : box dégénérée h=0 ≡ moteur float
            Z, W, dMS = chart_cell_section(S, g_col, eps, u0, v0, 0.0)
            g_iv, s, zW, WtWb = interval_chart_metric(Z, W, M_civ, coeffs218)
            G_ref = chart_metric_kahler(Z_f[i:i + 1], W_f[i:i + 1], M,
                                        coeffs218, B3, B3_MULTIS, B3_IDX)
            g_ref = pack_herm(G_ref)[0]
            mids = np.array([(a + b) / 2 for a, b in map(iv_bounds, g_iv)])
            rel = float(np.abs(mids - g_ref).max() / np.abs(g_ref).max())
            worst_rel = max(worst_rel, rel)
            # r float vs r intervalle (box dégénérée)
            blk = {"S": S, "Z": Z_f[i:i + 1]}
            dMS_f = detMS_on_block(blk)[0]
            r_ref = float(np.log(g_ref[0] * g_ref[1] - g_ref[2] ** 2
                                 - g_ref[3] ** 2) + 2 * np.log(abs(dMS_f)))
            r_iv_, _ = residual_iv(g_iv, dMS)
            lo, hi = iv_bounds(r_iv_)
            worst_r = max(worst_r, abs((lo + hi) / 2 - r_ref))
            n_pts += 1
    check("K2_degenerate_box_equals_float", worst_rel < 5e-12,
          f"max rel(G) = {worst_rel:.2e}, max |Δr| = {worst_r:.2e} "
          f"sur {n_pts} pts × {len(cases)} charts ({time.time() - t0:.0f}s)")

    # K3 — containment Monte-Carlo : float(point ∈ boîte) ∈ intervalle(boîte)
    t0 = time.time()
    h = 1e-3
    n_contained, n_tested, n_branch = 0, 0, 0
    for S, g_col, eps, u0, v0 in samples[:6]:
        try:
            Z, W, dMS = chart_cell_section(S, g_col, eps, u0, v0, h)
            g_iv, s, zW, WtWb = interval_chart_metric(Z, W, M_civ, coeffs218)
        except BranchCutError:
            n_branch += 1
            continue
        r_iv_, detG_iv = residual_iv(g_iv, dMS)
        for _ in range(8):
            du = rng.uniform(-h, h, 4)
            up = u0 + complex(du[0], du[1])
            vp = v0 + complex(du[2], du[3])
            Z_p, W_p = _float_section(S, g_col, eps, up, vp)
            if Z_p is None:
                continue
            G_p = chart_metric_kahler(Z_p, W_p, M, coeffs218,
                                      B3, B3_MULTIS, B3_IDX)
            g_p = pack_herm(G_p)[0]
            ok = all(g_p[c] in g_iv[c] for c in range(4))
            if r_iv_ is not None:
                i_, j_, k_ = S
                VS = ((MU_INT[j_] - MU_INT[i_]) * (MU_INT[k_] - MU_INT[i_])
                      * (MU_INT[k_] - MU_INT[j_]))
                dMS_p = 8.0 * Z_p[0, i_] * Z_p[0, j_] * Z_p[0, k_] * VS
                detg_p = g_p[0] * g_p[1] - g_p[2] ** 2 - g_p[3] ** 2
                r_p = float(np.log(detg_p) + 2 * np.log(abs(dMS_p)))
                ok = ok and (r_p in r_iv_)
            n_contained += ok
            n_tested += 1
    check("K3_containment_MC", n_tested > 0 and n_contained == n_tested,
          f"{n_contained}/{n_tested} points contenus (h = {h:g}, "
          f"{n_branch} cellules rejetées par garde de branche, "
          f"{time.time() - t0:.0f}s)")

    # K4 — largeur décroissante h → h/2
    S, g_col, eps, u0, v0 = samples[0]
    widths = []
    for hh in (1e-3, 5e-4):
        Z, W, dMS = chart_cell_section(S, g_col, eps, u0, v0, hh)
        g_iv, *_ = interval_chart_metric(Z, W, M_civ, coeffs218)
        widths.append(max(iv_width(g) for g in g_iv))
    check("K4_width_shrinks", widths[1] < widths[0],
          f"width(h)={widths[0]:.3e} → width(h/2)={widths[1]:.3e} "
          f"(ratio {widths[0] / widths[1]:.2f})")

    # K5 — adverse : coefficients altérés ⟹ K2 échoue (le test a des dents)
    bad = coeffs218.copy()
    bad[7] += 1e-6 * np.linalg.norm(coeffs218)
    S, g_col, eps, u0, v0 = samples[0]
    Z, W, dMS = chart_cell_section(S, g_col, eps, u0, v0, 0.0)
    g_bad, *_ = interval_chart_metric(Z, W, M_civ, bad)
    Z_p, W_p = _float_section(S, g_col, eps, u0, v0)
    g_ref = pack_herm(chart_metric_kahler(Z_p, W_p, M, coeffs218,
                                          B3, B3_MULTIS, B3_IDX))[0]
    mids = np.array([(a + b) / 2 for a, b in map(iv_bounds, g_bad)])
    rel_bad = float(np.abs(mids - g_ref).max() / np.abs(g_ref).max())
    check("K5_adverse_tamper_detected", rel_bad > 1e-10,
          f"coeffs altérés ⟹ rel = {rel_bad:.2e} ≫ tolérance K2")

    # D1 — dérivées duales (h=0) vs différences finies du moteur float
    S, g_col, eps, u0, v0 = samples[0]
    ZD, WD, _ = dual_chart_cell_section(S, g_col, eps, u0, v0, 0.0)
    g_dual = dual_chart_metric(ZD, WD, M_civ, coeffs218)
    fd_step = 1e-6
    worst_d1 = 0.0
    for a, (du, dv) in enumerate([(1, 0), (1j, 0), (0, 1), (0, 1j)]):
        gp = _float_packed(S, g_col, eps, u0 + fd_step * du,
                           v0 + fd_step * dv, M, coeffs218)
        gm = _float_packed(S, g_col, eps, u0 - fd_step * du,
                           v0 - fd_step * dv, M, coeffs218)
        fd = (gp - gm) / (2 * fd_step)
        for c in range(4):
            lo, hi = iv_bounds(g_dual[c].grads[a])
            mid = (lo + hi) / 2
            scale = max(1.0, abs(fd[c]))
            worst_d1 = max(worst_d1, abs(mid - fd[c]) / scale)
    check("D1_dual_grads_vs_FD", worst_d1 < 1e-5,
          f"max rel(∂G dual − FD) = {worst_d1:.2e} (pas FD {fd_step:g})")

    # D2 — containment MC de la forme valeur-moyenne
    h = 2.5e-4
    n_ok, n_tot = 0, 0
    for S, g_col, eps, u0, v0 in samples[:4]:
        try:
            g_mv, dMS_mv, diag = mean_value_metric(
                S, g_col, eps, u0, v0, h, M_civ, coeffs218)
        except (BranchCutError, RuntimeError):
            continue
        for _ in range(10):
            du = rng.uniform(-h, h, 4)
            up = u0 + complex(du[0], du[1])
            vp = v0 + complex(du[2], du[3])
            g_p = _float_packed(S, g_col, eps, up, vp, M, coeffs218)
            if g_p is None:
                continue
            n_ok += all(g_p[c] in g_mv[c] for c in range(4))
            n_tot += 1
    check("D2_mean_value_containment_MC", n_tot > 0 and n_ok == n_tot,
          f"{n_ok}/{n_tot} points contenus dans l'enclosure valeur-moyenne "
          f"(h = {h:g})")

    # D3 — la forme valeur-moyenne gagne sur la forme naïve
    S, g_col, eps, u0, v0 = samples[0]
    _, _, diag = mean_value_metric(S, g_col, eps, u0, v0, h, M_civ, coeffs218)
    gain = min(n / m for n, m in zip(diag["w_naive"], diag["w_mv"]))
    check("D3_mean_value_gain", gain > 3.0,
          f"largeur naïve / valeur-moyenne ≥ {gain:.1f}× (h = {h:g}, "
          f"w_mv g00 = {diag['w_mv'][0]:.2e} vs naïf {diag['w_naive'][0]:.2e})")

    # T1 — jets ordre 2 (h=0) : val/grads ≡ dual, hessiennes vs FD 2nd float
    import time as _time
    t0 = _time.time()
    S, g_col, eps, u0, v0 = samples[0]
    ZT, WT, dMS_t2 = t2_chart_cell_section(S, g_col, eps, u0, v0, 0.0)
    g_t2 = t2_chart_metric(ZT, WT, M_civ, coeffs218)
    t_eval = _time.time() - t0
    worst_g = 0.0
    for c in range(4):
        for a in range(NG):
            lo_d, hi_d = iv_bounds(g_dual[c].grads[a])
            lo_t, hi_t = iv_bounds(g_t2[c].g[a])
            worst_g = max(worst_g, abs((lo_t + hi_t) / 2 - (lo_d + hi_d) / 2))
    # Hessiennes vs FD des GRADIENTS DUAUX (prec 200 — la référence float
    # est plafonnée par le bruit de cancellation ‖c‖·ε/h² ≈ 4e-4)
    worst_h = 0.0
    deltas = {0: (1, 0), 1: (1j, 0), 2: (0, 1), 3: (0, 1j)}

    def dual_grads_at(du, dv):
        ZD_, WD_, _ = dual_chart_cell_section(S, g_col, eps, u0 + du,
                                              v0 + dv, 0.0)
        gd = dual_chart_metric(ZD_, WD_, M_civ, coeffs218)
        return np.array([[sum(iv_bounds(gd[c].grads[a])) / 2
                          for a in range(NG)] for c in range(4)])

    step = 1e-4
    H_fd = {}
    for b in range(NG):
        db, dvb = deltas[b]
        d1 = (dual_grads_at(step * db, step * dvb)
              - dual_grads_at(-step * db, -step * dvb)) / (2 * step)
        d2 = (dual_grads_at(step / 2 * db, step / 2 * dvb)
              - dual_grads_at(-step / 2 * db, -step / 2 * dvb)) / step
        H_fd[b] = (4.0 * d2 - d1) / 3.0          # Richardson O(step⁴)
    for i, (a, b) in enumerate(HPAIRS):
        for c in range(4):
            lo, hi = iv_bounds(g_t2[c].h[i])
            ref = H_fd[b][c][a]
            scale = max(10.0, abs(ref))
            worst_h = max(worst_h, abs((lo + hi) / 2 - ref) / scale)
    check("T1_t2_jets_vs_dual_and_FD", worst_g < 1e-12 and worst_h < 1e-6,
          f"grads t2 ≡ dual à {worst_g:.2e} ; hess vs FD(grads duaux) "
          f"Richardson rel {worst_h:.2e} ({t_eval:.1f}s/éval t2)")

    # T2 — containment MC de la forme de Taylor ordre 2
    h = 1e-3
    n_ok, n_tot = 0, 0
    for S, g_col, eps, u0, v0 in samples[:4]:
        try:
            g_pk, detG_t, dMS_t, ms2_t, diag = taylor2_metric(
                S, g_col, eps, u0, v0, h, M_civ, coeffs218)
        except (BranchCutError, RuntimeError):
            continue
        for _ in range(10):
            du = rng.uniform(-h, h, 4)
            up = u0 + complex(du[0], du[1])
            vp = v0 + complex(du[2], du[3])
            g_p = _float_packed(S, g_col, eps, up, vp, M, coeffs218)
            if g_p is None:
                continue
            detg_p = g_p[0] * g_p[1] - g_p[2] ** 2 - g_p[3] ** 2
            ok = all(g_p[c] in g_pk[c] for c in range(4))
            ok = ok and (detg_p in detG_t)
            n_ok += ok
            n_tot += 1
    check("T2_taylor2_containment_MC", n_tot > 0 and n_ok == n_tot,
          f"{n_ok}/{n_tot} points contenus (G packé + det G, h = {h:g})")

    # T3 — l'ordre 2 gagne sur la valeur-moyenne
    S, g_col, eps, u0, v0 = samples[0]
    _, _, diag_mv = mean_value_metric(S, g_col, eps, u0, v0, h,
                                      M_civ, coeffs218)
    _, _, _, _, diag_t2 = taylor2_metric(S, g_col, eps, u0, v0, h,
                                         M_civ, coeffs218)
    gain = diag_mv["w_mv"][0] / diag_t2["w_t2_g00"]
    check("T3_taylor2_gain", gain > 3.0,
          f"w_mv(g00) / w_t2(g00) = {gain:.1f}× (h = {h:g} ; "
          f"w_t2 = {diag_t2['w_t2_g00']:.2e})")

    print(f"\nselftest kernel : {11 - len(failures)}/11 PASS")
    return 1 if failures else 0


def _float_packed(S, g_col, eps, u, v, M, coeffs218):
    from .kahler_metric import chart_metric_kahler
    from .ricci_functional import pack_herm
    Z_p, W_p = _float_section(S, g_col, eps, u, v)
    if Z_p is None:
        return None
    return pack_herm(chart_metric_kahler(Z_p, W_p, M, coeffs218,
                                         B3, B3_MULTIS, B3_IDX))[0]


def _float_section(S, g_col, eps, u: complex, v: complex):
    """Section float exacte du sampler pour (u, v) donnés (feuille ε)."""
    T = tuple(j for j in range(6) if j not in S)
    others = [c for c in T if c != g_col]
    o1, o2 = others
    A_exact = minor_inv_times_T_exact(S, T)
    perm = [list(T).index(g_col), list(T).index(o1), list(T).index(o2)]
    A = np.array([[float(A_exact[r][c]) for c in perm] for r in range(3)])
    R = A[:, 0] + A[:, 1] * u * u + A[:, 2] * v * v
    if np.abs(R).min() < 1e-15:
        return None, None
    w0 = np.sqrt(R.astype(complex))
    Z = np.zeros((1, 6), dtype=complex)
    Z[0, g_col] = 1.0
    Z[0, o1], Z[0, o2] = u, v
    Z[0, list(S)] = np.array(eps, float) * w0
    W = np.zeros((1, 6, 2), dtype=complex)
    W[0, o1, 0] = 1.0
    W[0, o2, 1] = 1.0
    W[0, list(S), 0] = A[:, 1] * u / Z[0, list(S)]
    W[0, list(S), 1] = A[:, 2] * v / Z[0, list(S)]
    return Z, W


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--selftest":
        raise SystemExit(_selftest())
    raise SystemExit("usage: k3_cap_r4_interval_kernel.py --selftest")
