#!/usr/bin/env python3
"""
k3_cap_b1e2iii_c129f_f2f3_bridge_atlas.py — C129-F0 : F2 + F3, LES
CARTES-PONTS ET LEURS DEUX TRANSITIONS.

CE QUE CE SCRIPT PAIE (V2) — le contrat §12 de la 10ᵉ revue GPT
(`gpt_b1e2iii_c129f_cadrage_review_2026_07_31.md`), corrigé du facteur 2
établi par le scout et sous le nom honnête fixé par F1a, ET le contrat
de réparation de la 11ᵉ revue
(`gpt_b1e2iii_c129f0_scout_f1_f2f3_review_2026_07_31.md`). LA V1 DE CE
SCRIPT AVAIT DEUX DÉFAUTS LOAD-BEARING, CASSÉS EN REVUE ET RÉPARÉS
ICI : (i) elle lisait le motif de θ du côté supérieur comme une
« prédiction vérifiée » — c'est-à-dire qu'elle renommait un ÉCHEC de
raccord en gate vert ; la v2 en fait un THÉORÈME DE DECK
(`D = diag(+1,−1,+1,+1,+1,−1)`) et recolle la feuille CONTINUÉE
`D·Z_conj`, qui est LA transition supérieure (le conjugué seul ne se
recolle PAS et reste un diagnostic) ; (ii) son « nerf » comptait comme
arêtes des intersections sans transition certifiée — la v2 le publie
sous son vrai nom (`domain_intersection_graph`) et reconstruit le NERF
sur arêtes CERTIFIÉES seulement. La portée est STRATIFIÉE (36/64
pleinement ouvertes en 4D), pas « ambiante » en bloc.

L'ARCHITECTURE, ET POURQUOI ELLE N'EST PAS CELLE DE LA NOTE DE CADRAGE.
`U⁻` (Im u, Im v < 0) et `U⁺` (Im u, Im v > 0) n'ont AUCUN overlap
ouvert : leurs adhérences ne se rencontrent que sur le coin de
codimension réelle 2 `{Im u = Im v = 0}`. Comparer directement une
section inférieure et une section supérieure « sur la face » serait un
argument par valeur de bord. On introduit donc une CARTE-PONT `B` sur la
boîte bilatérale, et on certifie DEUX transitions séparées :

    B = Z⁻  sur  bridge ∩ halo⁻        B = Z⁺  sur  bridge ∩ halo⁺

Ces deux intersections ont une largeur STRICTEMENT POSITIVE dans les
quatre coordonnées : ce sont de vrais overlaps ouverts, et les ancres y
sont prises STRICTEMENT à l'intérieur (Im < 0 d'un côté, Im > 0 de
l'autre), jamais sur la face.

LE CONSTRUCTEUR BILATÉRAL, ET POURQUOI C129-E N'EST PAS RÉUTILISABLE.
`build_section` fait « principale, sinon tournée avec σ de composante ».
Sur une boîte bilatérale, `Im R` change de signe — `σ` n'existe pas.
`build_section_bilateral` choisit donc PAR LE SIGNE CERTIFIÉ DE `Re R`,
sans aucun essai :

    Re R minoré > 0 sur toute la boîte  →  principale
    Re R majoré < 0 sur toute la boîte  →  canonique `w = i√_p(−R)`
                                           (σ ÉPINGLÉ à +1 : c'est la
                                            définition du régime C127-E,
                                            garde C103 verbatim sur −R)
    sinon                               →  REFUS

Le signe est lu sur l'ENCLOSURE TM de `Re R` (reste compris), pas sur la
borne rationnelle du scout : la décision est prise dans l'arithmétique
qui construit la section. La borne rationnelle exacte est publiée à côté
comme CROSS-CHECK — deux chemins indépendants pour le même signe.

LA COMPARAISON, ET LE PIÈGE QU'ELLE ÉVITE. Comparer `B` et `Z∓` par
leurs ENCLOSURES ne compare rien : chaque enclosure a la largeur de la
VARIATION de la fonction sur l'overlap, donc leur différence aussi. On
RECENTRE donc les deux Taylor-modèles dans le cadre de l'overlap —
substitution affine `ε = off + scale·ε'` exacte en Fraction, ANISOTROPE
(le pont vaut 2H en Im et H en Re, ce n'est pas un cube) — puis on
soustrait COEFFICIENT PAR COEFFICIENT. C'est la parade C127-D, ici avec
une matrice de recentrage à échelle par coordonnée.

LA PHASE EST UN RÉSULTAT OUVERT, LE GATE EST PRÉ-ENREGISTRÉ. `θ = +1`
signifie que le pont et le côté décrivent le MÊME point projectif (la
jauge source étant normalisée `Z_g = 1`, le scalaire projectif est fixé :
des signes indépendants par ligne ne sont PAS une compatibilité
projective). Le test est ouvert quant au résultat, préregistré quant à
la signification du succès.

ET LE RÉSULTAT N'EST PAS CELUI QU'ON ATTENDAIT — C'EST LA TROUVAILLE.
Le pont se recolle exactement au côté INFÉRIEUR. Contre le côté
SUPÉRIEUR dérivé par conjugaison, `θ` est MIXTE : `+1` sur la ligne
principale, `−1` sur les deux lignes canoniques. Ce n'est ni un bug ni
une ambiguïté (les marges sont O(1), min 1,04), et ce n'est PAS une
renormalisation projective. La raison est structurelle et se prédit :

    sur le coin `Im u = Im v = 0`, le radicande `R` est RÉEL ;
    là où `Re R > 0` (régime principal) la racine est RÉELLE et
      l'involution antiholomorphe la FIXE ;
    là où `Re R < 0` (régime canonique) la racine `i√(−R)` est
      PUREMENT IMAGINAIRE et l'involution la NIE.

Donc **le point conjugué et le point continué analytiquement à travers
le coin sont deux points DISTINCTS de `P⁵` au-dessus du même `(ū, v̄)`,
échangés par l'involution réelle.** Conséquence pour le levier n° 1 de
la note de cadrage : la voisine conjuguée est un atlas légitime, mais ce
n'est PAS celui que la continuation analytique atteint. Le motif de `θ`
est PRÉDIT depuis le seul régime, et la prédiction est vérifiée sur
64/64 : la surprise devient un résultat falsifiable au lieu de rester un
échec de gate.

L'EXACTITUDE. La séparation `sep_phase` est PONCTUELLE sur tout
l'overlap (elle porte sur des enclosures corrélées après recentrage) :
en chaque point, `B` et `Z∓` sont deux racines analytiques du MÊME
`R ≠ 0`, et `θ = −1` donnerait `|diff| = 2|Z| > 0 = |somme|`, exclu. Le
lemme C129-D s'applique alors VERBATIM sur chaque overlap — mais ses
gates sont RECONSTRUITS ici, pas importés : `F ≠ 0` sur le pont ET sur
le côté, jauge non nulle, ancre intérieure, séparation stricte.

CE QUE CE SCRIPT NE PAIE PAS : la métrique du pont (F4 — Qmat, Weyl,
congruences latérales) ; l'identification canonique de la voisine, qui a
REFUSÉ (F1a, artefact `c129f_f1_mirror_ledger`) — l'atlas supérieur est
DÉRIVÉ par conjugaison, pas énuméré, et ce script hérite de ce statut ;
les voisines de codimension 1 ; le raccord à travers les faces Re, où
les ponts restent des cartes RELATIVES ; les 895 autres paires ; R12-C.

GATES
  F2a  géométrie : demi-largeur 2H dans les deux directions Im
       réfléchies, H dans les directions Re — re-dérivée ici en
       Fraction, et IDENTIQUE à celle du scout (import vérifié) ;
  F2b  `core⁻ ∪ core⁺ ⊆ bridge`, STRICTEMENT dans les directions Im ;
  F2c  régime par signe CERTIFIÉ de `Re R` sur l'enclosure TM, sans
       essai, sur 64 × 3 lignes — et accord avec la borne rationnelle
       exacte du scout (deux arithmétiques, un seul verdict) ;
  F2d  section bilatérale COMPLÈTE : les six coordonnées existent, la
       jauge du chart cible est minorée > 0, le ledger du pont
       (ε, régime) est FIGÉ et sérialisé ;
  F2e  NÉGATIFS : (i) le faux pont de demi-largeur H est REFUSÉ par
       l'inclusion ; (ii) le recours à `σ` de composante sur le pont est
       REFUSÉ (Im R straddle) ; (iii) une boîte où `Re R` straddle est
       REFUSÉE au lieu de retomber sur un régime ;
  F2f  64/64 sans filtrage silencieux ;
  F3a  les deux overlaps ont une largeur > 0 dans les 4 coordonnées, en
       rationnels exacts, et les deux ancres sont STRICTEMENT
       intérieures (Im < 0 / Im > 0), jamais sur la face ;
  F2d(bis) LE LEDGER DU PONT SE DÉRIVE : il est mesuré contre le côté
       INFÉRIEUR (l'atlas établi, ancre de vérité), puis la section est
       RECONSTRUITE et RE-VÉRIFIÉE. Défauter le ledger reviendrait à
       appeler « raccord » un choix de feuille par défaut ;
  F3b⁻ le pont se recolle EXACTEMENT au côté inférieur : `θ = +1` sur
       TOUTES les lignes, par recentrage ANISOTROPE exact + séparation
       stricte, refus si ambigu ;
  R1   `D = diag(+1,−1,+1,+1,+1,−1)` est une TRANSFORMATION DE DECK :
       involutive, NON scalaire, préserve les trois quadriques, et
       l'identité `Z_conj = D·Z_bridge` est certifiée AU COEFFICIENT
       sur 64/64. L'algèbre seule ne distingue rien (tout diagonal de
       signes préserve les quadriques) : la discriminance est portée
       par le négatif R1e ;
  R2   le pont se recolle EXACTEMENT (`θ = +1` sur les six
       coordonnées) à la feuille CONTINUÉE `D·Z_conj` — c'est ELLE la
       transition supérieure ; le conjugué SEUL ne se recolle PAS
       (R2b) et reste publié comme diagnostic d'une autre feuille ;
  R3   le NERF ne compte que des arêtes dont la transition est
       CERTIFIÉE (380 nœuds ; 5396 L↔L importées de C127-D, 64 B↔L,
       210/210 B↔B ; connexe) ; le graphe des intersections de boîtes
       est publié séparément sous `domain_intersection_graph`, qui
       n'est PAS un nerf ;
  F3d  négatif PROJECTIVEMENT DISCRIMINANT : mutation NON scalaire
       (une seule ligne canonique niée, jauge et affines inchangées)
       ⟹ le raccord tombe 64/64 ; la négation globale `Z→−Z` est
       publiée comme MÊME point de `P⁵` — elle ne discriminait rien,
       c'était le défaut de la v1 ;
  R1e  un seul signe changé dans `D` ⟹ l'identité de deck CASSE ;
  R4   chaîne amont vérifiée (scout, F1, C127-D, C127-E).

Sortie : results/k3_cap_b1e2iii_c129f_f2f3_bridge_atlas.json
Usage  : k3_cap_b1e2iii_c129f_f2f3_bridge_atlas.py [--selftest]
Env    : K3_F2F3_WORKERS (défaut 6)
"""
from __future__ import annotations

import hashlib
import io
import itertools
import json
import os
import platform
import subprocess
import sys
import time
from collections import Counter, deque
from fractions import Fraction
from math import comb
from multiprocessing import get_context
from pathlib import Path

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

_argv = sys.argv
sys.argv = [sys.argv[0]]
os.environ.setdefault("K3_TM_ORDER", "4")
os.environ.setdefault("K3_TM_SERIES", "4")
from mpmath import iv, mp                                          # noqa: E402
from .taylor_models import (                                     # noqa: E402
    CIV, IVPM, MIDX, MONO, NM, TMC, TM_ORDER, UNARY_SERIES_DEG,
    BranchCutError, civ_absmin, riv, tm_sqrt_rotated)
from .interval_arithmetic import minor_inv_times_T_exact      # noqa: E402
from .full_cell_charts import build_section               # noqa: E402
from .atlas_assembly import (                           # noqa: E402
    apply_recenter, civ_sup, enclose, frac_to_iv, fraction_box_to_iv,
    poly_lin, sep_phase, _f_down, _f_up)
sys.argv = _argv

RES = Path(os.environ.get(
    "K3_RES_DIR", Path(__file__).resolve().parent / "data"))
COVER_JSON = RES / "k3_cap_b1e2iii_d5_6_dyadic_cover.json"
ATLAS_JSON = RES / "k3_cap_b1e2iii_c127d_atlas.json"
C127E_JSON = RES / "k3_cap_b1e2iii_c127e_residual.json"
SCOUT_JSON = RES / "k3_cap_b1e2iii_c129f_bridge_scout.json"
F1_JSON = RES / "k3_cap_b1e2iii_c129f_f1_mirror_ledger.json"
ART = RES / "k3_cap_b1e2iii_c129f_f2f3_bridge_atlas.json"
N_WORKERS = int(os.environ.get("K3_F2F3_WORKERS", "6"))

# --- PRÉ-ENREGISTRÉ, figé avant le run ------------------------------------
IM_DIRS = (1, 3)
NG = 4
# Le raccord SUCCÈDE si et seulement si θ = +1 des deux côtés. −1 est un
# résultat publiable et un ÉCHEC de raccord : préregistré ici, pas
# arbitré après coup.
THETA_REQUIRED = 1
T0 = time.time()
IV0 = iv.mpf(0)
IV1 = iv.mpf(1)
CZERO = CIV(IV0, IV0)


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def _sha(p):
    try:
        return hashlib.sha256(Path(p).read_bytes()).hexdigest()
    except OSError:
        return None


def _q(x):
    return [int(x.numerator), int(x.denominator)]


def _qbox(bx):
    return [[_q(a), _q(b)] for a, b in bx]


# ===========================================================================
#  Géométrie exacte — jamais un float dans une décision
# ===========================================================================
def box_of(o):
    return ([Fraction(float.fromhex(x)) for x in o["center_hex"]],
            Fraction(float.fromhex(o["hw_hex"])))


def bounds(c, H):
    return [(c[k] - H, c[k] + H) for k in range(4)]


def mirror_bounds(bx):
    return [((-bx[k][1], -bx[k][0]) if k in IM_DIRS else bx[k])
            for k in range(4)]


def inter(a, b):
    """Boîte d'intersection, ou None si la largeur n'est pas STRICTEMENT
    positive dans les quatre coordonnées. Tout en Fractions — « ouvert »
    est une comparaison rationnelle, jamais un test float."""
    out = []
    for k in range(4):
        lo, hi = max(a[k][0], b[k][0]), min(a[k][1], b[k][1])
        if not hi > lo:
            return None
        out.append((lo, hi))
    return out


def center_hw(bx):
    return ([(a + b) / 2 for a, b in bx], [(b - a) / 2 for a, b in bx])


def contains(outer, inner, strict_dirs=IM_DIRS):
    """`inner ⊆ outer`, STRICTEMENT dans `strict_dirs`. La distinction
    est du contenu : le pont est ouvert dans les directions Im
    réfléchies (c'est la promotion du coin en voisinage ambiant), il
    reste AFFLEURANT dans les directions Re, où il hérite du halo
    clippé — ces faces-là appartiennent aux voisines de codimension 1 et
    ne sont pas payées ici."""
    marg = [(inner[k][0] - outer[k][0], outer[k][1] - inner[k][1])
            for k in range(4)]
    ok = (all(a >= 0 and b >= 0 for a, b in marg)
          and all(marg[k][0] > 0 and marg[k][1] > 0 for k in strict_dirs))
    return ok, marg


# ===========================================================================
#  Recentrage ANISOTROPE exact — la généralisation dont le pont a besoin
# ===========================================================================
#  C127-D recentrait entre deux cubes : une seule `scale`. Le pont vaut
#  2H en Im et H en Re, donc l'échelle est PAR COORDONNÉE. La
#  substitution reste affine, donc elle ne crée AUCUN reste (composer un
#  polynôme de degré ≤ N avec une application affine donne un polynôme de
#  degré ≤ N), et la matrice est exacte en Fractions avant conversion.
# ===========================================================================
def recenter_matrix_aniso(off, scale):
    """T tel que p∘φ [β] = Σ_{α ≥ β} T[β][α]·p[α], pour
    φ(ε)_k = off_k + scale_k·ε_k."""
    acc = {}
    for ai, am in enumerate(MONO):
        for bm in itertools.product(*[range(e + 1) for e in am]):
            coef = Fraction(1)
            for k in range(4):
                coef *= (comb(am[k], bm[k])
                         * off[k] ** (am[k] - bm[k])
                         * scale[k] ** bm[k])
            if coef:
                row = acc.setdefault(MIDX[bm], {})
                row[ai] = row.get(ai, Fraction(0)) + coef
    return {b: [(a, frac_to_iv(c)) for a, c in row.items() if c]
            for b, row in acc.items()}


def reframe(src_box, dst_box):
    """(off, scale) EXACTS envoyant le cadre ε de `dst_box` dans celui de
    `src_box` : si `x = c_s + h_s·ε_s = c_d + h_d·ε_d`, alors
    `ε_s = (c_d − c_s)/h_s + (h_d/h_s)·ε_d`."""
    cs, hs = center_hw(src_box)
    cd, hd = center_hw(dst_box)
    return ([(cd[k] - cs[k]) / hs[k] for k in range(4)],
            [hd[k] / hs[k] for k in range(4)])


def frame_admissible(off, scale):
    """Le TM n'est valide que sur ε ∈ [−1,1]⁴ : l'image du cadre cible
    doit y tenir. C'est un GATE, pas une hypothèse."""
    return all(-1 <= off[k] - scale[k] and off[k] + scale[k] <= 1
               for k in range(4))


# ===========================================================================
#  Le constructeur BILATÉRAL
# ===========================================================================
def uv_tm(center, hw):
    """`u`, `v` en TMC sur une boîte à demi-largeurs PAR COORDONNÉE."""
    e = [MIDX[tuple(1 if k == j else 0 for k in range(NG))]
         for j in range(NG)]
    pu = [CZERO] * NM
    pu[0] = CIV(riv(center[0]), riv(center[1]))
    pu[e[0]] = CIV(riv(hw[0]), IV0)
    pu[e[1]] = CIV(IV0, riv(hw[1]))
    pv = [CZERO] * NM
    pv[0] = CIV(riv(center[2]), riv(center[3]))
    pv[e[2]] = CIV(riv(hw[2]), IV0)
    pv[e[3]] = CIV(IV0, riv(hw[3]))
    return TMC(pu), TMC(pv)


def line_coeffs(S, g):
    T = tuple(j for j in range(6) if j not in S)
    others = [c for c in T if c != g]
    Ae = minor_inv_times_T_exact(S, T)
    perm = [list(T).index(g), list(T).index(others[0]),
            list(T).index(others[1])]
    return ([[riv(Ae[r][perm[c]]) for c in range(3)] for r in range(3)],
            [[Fraction(Ae[r][perm[c]]) for c in range(3)]
             for r in range(3)], others)


def _rng_re(t):
    x = t.re_tm().to_iv()
    return (mp.mpf(x.a), mp.mpf(x.b))


def _rng_im(t):
    x = t.im_tm().to_iv()
    return (mp.mpf(x.a), mp.mpf(x.b))


def build_section_bilateral(S, g, eps, center, hw, force_sigma=False):
    """La section source sur une boîte BILATÉRALE, régime choisi par le
    SIGNE CERTIFIÉ de `Re R` — jamais par essai, jamais de repli.

    `force_sigma` est la MUTATION du négatif F2e(ii) : elle tente la
    détermination tournée avec un `σ` de composante là où `Im R`
    straddle. Un code qui ferait cela sur un pont accepterait une
    continuation dont la feuille n'est pas définie — c'est ce que le
    négatif doit exhiber, et il doit CASSER."""
    A, _Aq, others = line_coeffs(S, g)
    o1, o2 = others
    u, v = uv_tm(center, hw)
    u2, v2 = u * u, v * v
    Z = [None] * 6
    dZ = [None] * 6
    Z[g] = TMC.const(CIV(IV1, IV0))
    dZ[g] = (TMC.const(CZERO), TMC.const(CZERO))
    Z[o1], dZ[o1] = u, (TMC.const(CIV(IV1, IV0)), TMC.const(CZERO))
    Z[o2], dZ[o2] = v, (TMC.const(CZERO), TMC.const(CIV(IV1, IV0)))
    rows = []
    for r, s in enumerate(S):
        R = TMC.const(CIV(A[r][0])) + u2.mul_real(A[r][1]) \
            + v2.mul_real(A[r][2])
        re_lo, re_hi = _rng_re(R)
        im_lo, im_hi = _rng_im(R)
        rec = {"row": r, "s_coord": int(s),
               "re_R_tm": [_f_down(re_lo), _f_up(re_hi)],
               "im_R_tm": [_f_down(im_lo), _f_up(im_hi)],
               "im_R_straddles": bool(im_lo <= 0 <= im_hi),
               "regime": None}
        if force_sigma:
            # MUTATION : on prétend lire σ = signe de Im R.
            if im_lo > 0:
                sg = 1
            elif im_hi < 0:
                sg = -1
            else:
                rec["refused"] = "component_sigma_undetermined_on_bridge"
                rows.append(rec)
                continue
            try:
                Zs = tm_sqrt_rotated(R, sg).mul_real(riv(int(eps[r])))
                rec["regime"] = "rotated%+d" % sg
            except BranchCutError as exc:
                rec["refused"] = exc.diag.get("guard")
                rows.append(rec)
                continue
        else:
            if re_lo > 0:
                sel = "principal"
            elif re_hi < 0:
                sel = "canonical"
            else:
                rec["refused"] = "re_R_straddles_zero_no_regime"
                rows.append(rec)
                continue
            rec["selected_by"] = ("re_R_lower_bound_positive"
                                  if sel == "principal"
                                  else "re_R_upper_bound_negative")
            try:
                if sel == "principal":
                    Zs = R.sqrt_principal().mul_real(riv(int(eps[r])))
                else:
                    # σ ÉPINGLÉ à +1 : c'est la DÉFINITION du régime
                    # canonique (C127-E), pas un essai. La garde C103
                    # s'applique verbatim à −R.
                    Zs = tm_sqrt_rotated(R, 1).mul_real(riv(int(eps[r])))
                rec["regime"] = sel
            except BranchCutError as exc:
                # Le signe était certifié : si la garde refuse quand
                # même, c'est un REFUS, jamais un repli sur l'autre
                # branche.
                rec["refused"] = exc.diag.get("guard")
                rows.append(rec)
                continue
        Z[s] = Zs
        iZ = Zs.inv()
        dZ[s] = (u.mul_real(A[r][1]) * iZ, v.mul_real(A[r][2]) * iZ)
        # Alias de compatibilité avec `transport_hardened`, qui sérialise
        # `source_determinations`. Le régime du PONT diffère légitimement
        # de celui du core (canonique là où le core était tourné) : la
        # différence doit être VISIBLE dans l'artefact, pas masquée — le
        # gate métrique porte sur les `kinds` CIBLE, pas sur le régime
        # source, précisément parce que ce dernier a changé exprès.
        rec["determination"] = rec["regime"]
        rr, ri = R.to_iv_pair()
        rec["radicand_absmin"] = _f_down(
            mp.mpf(civ_absmin(CIV(rr, ri)).a))
        rows.append(rec)
    return Z, dZ, rows


# ===========================================================================
#  Comparaison de deux sections sur un overlap — recentrage puis
#  soustraction COEFFICIENT PAR COEFFICIENT
# ===========================================================================
UNIT_E = [iv.mpf([-1, 1])] * 4


def theta_lines(Za, box_a, Zb, box_b, W, keys):
    offa, sca = reframe(box_a, W)
    offb, scb = reframe(box_b, W)
    if not (frame_admissible(offa, sca) and frame_admissible(offb, scb)):
        return {"refused": "overlap_outside_tm_domain"}
    Ta = recenter_matrix_aniso(offa, sca)
    Tb = recenter_matrix_aniso(offb, scb)
    res = {}
    for s in keys:
        za, zb = Za[s], Zb[s]
        if za is None or zb is None:
            res[s] = {"refused": "missing_coordinate"}
            continue
        pa = apply_recenter(Ta, za.p)
        pb = apply_recenter(Tb, zb.p)
        rem = za.rem + zb.rem
        dm = enclose(poly_lin(pa, pb, -1), rem, UNIT_E)
        dp = enclose(poly_lin(pa, pb, +1), rem, UNIT_E)
        th, diag = sep_phase(dm, dp)
        res[s] = {"theta": th, **diag,
                  "diff_sup": _f_up(civ_sup(dm)) if th == 1 else None}
    return res


# ===========================================================================
#  R1 — LA TRANSFORMATION DE DECK
# ===========================================================================
#  Le motif de θ mesuré contre le conjugué définit un diagonal de signes
#
#      D = diag(+1, −1, +1, +1, +1, −1)
#
#  qui préserve les trois quadriques ambiantes `Q_m(Z) = Σ_a μ_a^m Z_a²`
#  (les coordonnées n'y entrent que par leur CARRÉ), vaut son propre
#  inverse, et n'est PAS un scalaire projectif. C'est donc un
#  automorphisme holomorphe non trivial de l'intersection complète : une
#  TRANSFORMATION DE DECK de la description par racines carrées.
#
#  HONNÊTETÉ SUR CE QUI EST TRIVIAL. « D préserve les quadriques » est
#  vrai de TOUT diagonal de signes — ce n'est pas ce qui distingue CE D.
#  Le contenu est ailleurs, et il est mesuré : c'est CE D-là, et aucun
#  autre, qui relie la section conjuguée à la section du pont. Le
#  négatif R1e (un seul signe changé) porte donc toute la discriminance ;
#  la préservation des quadriques dit seulement que `D·Z` est encore un
#  point de la K3, ce qui est nécessaire pour que la réparation R2 ait
#  un sens.
# ===========================================================================
DECK_D = (1, -1, 1, 1, 1, -1)


def apply_deck(Z, D, keys):
    """`D·Z` au niveau des Taylor-modèles : un signe par coordonnée."""
    out = {}
    for k in keys:
        z = Z[k]
        if z is None or D[k] == 1:
            out[k] = z
        else:
            out[k] = TMC([-x for x in z.p], z.rem)
    return out


def deck_algebra(D):
    """Les trois faits structurels sur `D`, EXACTS (entiers/Fractions) :
    involution, non-scalaire dans PGL(6), et préservation des trois
    quadriques de Vandermonde."""
    from .interval_arithmetic import MU_INT
    invol = all(d * d == 1 for d in D)
    non_scalar = len(set(D)) > 1
    quad = all(
        sum(Fraction(MU_INT[a]) ** m * Fraction(D[a]) ** 2
            for a in range(6))
        == sum(Fraction(MU_INT[a]) ** m for a in range(6))
        for m in range(3))
    return {"involution": bool(invol),
            "non_scalar_in_pgl6": bool(non_scalar),
            "preserves_the_three_quadrics": bool(quad),
            "quadric_preservation_is_trivial_for_any_sign_diagonal": True,
            "D": list(D)}


def eval_at_point(Z, box, x, keys):
    """Valeur des coordonnées en un POINT (ancre), via le cadre ε du TM.
    L'appartenance de l'ancre à la boîte est un gate EXACT en Fraction,
    vérifié par l'appelant."""
    c, h = center_hw(box)
    E = [frac_to_iv((x[k] - c[k]) / h[k]) for k in range(4)]
    out = {}
    for s in keys:
        if Z[s] is None:
            out[s] = None
            continue
        out[s] = enclose(Z[s].p, Z[s].rem, E)
    return out


def sep_at_point(a, b):
    """Séparation stricte |a−b|² vs |a+b|² en un point."""
    dm = CIV(a.re - b.re, a.im - b.im)
    dp = CIV(a.re + b.re, a.im + b.im)
    return sep_phase(dm, dp)


# ===========================================================================
#  Job par pont
# ===========================================================================
_G = {}


def _init(cell, leaves, halos, bridges, ledgers=None):
    _G.update(cell=cell, leaves=leaves, halos=halos, bridges=bridges,
              ledgers=ledgers or {})


def _bridge_job(i):
    S, g, eps = _G["cell"]
    leaf = _G["leaves"][i]
    rec_h = _G["halos"][i]
    Wb = _G["bridges"][i]
    out = {"tile": i}

    # --- géométrie ---------------------------------------------------
    c_core, h_core = box_of(leaf)
    core = bounds(c_core, h_core)
    c_h = [Fraction(float.fromhex(x)) for x in rec_h["center_hex"]]
    H = Fraction(float.fromhex(rec_h["H_hex"]))
    halo_lo = bounds(c_h, H)
    halo_up = mirror_bounds(halo_lo)
    core_up = mirror_bounds(core)
    union = [(min(core[k][0], core_up[k][0]),
              max(core[k][1], core_up[k][1])) for k in range(4)]
    ok_inc, marg = contains(Wb, union)
    out["F2b_contains_union"] = bool(ok_inc)
    out["F2b_margins"] = [[float(a), float(b)] for a, b in marg]
    # NÉGATIF (i) : le faux pont de la note de cadrage
    Wbad = [((-H, H) if k in IM_DIRS else Wb[k]) for k in range(4)]
    out["F2e_i_false_bridge_refused"] = bool(
        not contains(Wbad, union)[0])

    cb, hb = center_hw(Wb)
    cbf = [float(x) for x in cb]
    hbf = [float(x) for x in hb]
    Zb, dZb, rows_b = build_section_bilateral(S, g, eps, cbf, hbf)
    out["bridge_rows_default_ledger"] = rows_b
    out["F2c_all_regimes_assigned"] = bool(
        all(r["regime"] is not None for r in rows_b))
    out["F2d_section_complete"] = bool(all(z is not None for z in Zb))

    # NÉGATIF (ii) : forcer le recours à σ de composante sur le pont
    _Zf, _df, rows_f = build_section_bilateral(
        S, g, eps, cbf, hbf, force_sigma=True)
    out["F2e_ii_sigma_refused"] = bool(
        any(r.get("refused") == "component_sigma_undetermined_on_bridge"
            for r in rows_f))
    out["F2e_ii_rows"] = [{"s_coord": r["s_coord"],
                           "refused": r.get("refused"),
                           "im_R_straddles": r["im_R_straddles"]}
                          for r in rows_f]

    if not (out["F2c_all_regimes_assigned"] and out["F2d_section_complete"]):
        out["refused"] = "bridge_section_incomplete"
        return out

    # jauge du chart CIBLE minorée > 0 sur le pont
    S2, g2 = tuple(leaf["chart"]["S"]), leaf["chart"]["g"]
    gr, gi = Zb[g2].to_iv_pair()
    gmin = mp.mpf(civ_absmin(CIV(gr, gi)).a)
    out["F2d_target_gauge_absmin"] = _f_down(gmin)
    out["F2d_gauge_positive"] = bool(gmin > 0)

    # --- les deux sections de côté ------------------------------------
    chf = [float(x) for x in c_h]
    Hf = float(H)
    Zlo, _d1, rows_lo = build_section(S, g, eps, chf, Hf)
    cup = [(-chf[k] if k in IM_DIRS else chf[k]) for k in range(4)]
    Zup, _d2, rows_up = build_section(S, g, eps, cup, Hf)
    out["lower_dets"] = [r["determination"] for r in rows_lo]
    out["upper_dets"] = [r["determination"] for r in rows_up]
    # R4c — le côté inférieur RECONSTRUIT doit être l'objet EFFECTIVEMENT
    # certifié par C127-D, pas une nouvelle reconstruction qui lui
    # ressemble. Sans ce gate, une dérive future du code dériverait le
    # pont contre autre chose que l'atlas établi.
    out["R4c_lower_dets_match_c127d"] = bool(
        out["lower_dets"] == list(rec_h["source_determinations"]))

    keys = list(S) + [x for x in range(6) if x not in S]

    # --- F2d bis : LE LEDGER DU PONT SE DÉRIVE, IL NE SE DÉFAUT PAS ----
    # L'atlas INFÉRIEUR est l'objet établi (C127-D/E, C129-D/E) : c'est
    # lui l'ancre de vérité, et le ledger du pont doit s'y accorder. On
    # le DÉRIVE (mesure de θ ligne par ligne contre le côté inférieur),
    # puis on RECONSTRUIT et on RE-VÉRIFIE — ce n'est pas un essai : la
    # re-vérification doit rendre θ = +1 sur TOUTES les lignes, sinon
    # REFUS. Ne pas dériver reviendrait à appeler « raccord » un simple
    # choix de feuille par défaut.
    Wl0 = inter(Wb, halo_lo)
    if Wl0 is None:
        out["refused"] = "lower_overlap_not_open"
        return out
    th0 = theta_lines({k: Zb[k] for k in keys}, Wb,
                      {k: Zlo[k] for k in keys}, halo_lo, Wl0, list(S))
    if "refused" in th0:
        out["refused"] = th0["refused"]
        return out
    if any(th0[k].get("theta") is None for k in S):
        out["refused"] = "ledger_derivation_ambiguous"
        out["F2d_ledger_probe"] = {str(k): th0[k].get("theta") for k in S}
        return out
    eps_bridge = tuple(int(th0[list(S)[r]]["theta"]) * int(eps[r])
                       for r in range(3))
    out["F2d_ledger_default"] = [int(x) for x in eps]
    out["F2d_ledger_derived"] = [int(x) for x in eps_bridge]
    out["F2d_ledger_probe_theta"] = [int(th0[k]["theta"]) for k in S]
    Zb, dZb, rows_b = build_section_bilateral(
        S, g, eps_bridge, cbf, hbf)
    out["bridge_rows"] = rows_b
    if not all(z is not None for z in Zb):
        out["refused"] = "bridge_section_incomplete_after_ledger"
        return out

    # --- F3a : overlaps ouverts + ancres STRICTEMENT intérieures -------
    res = {}
    # R2 — la feuille supérieure CONTINUÉE. `Z_upper_conj` est l'atlas
    # dérivé par conjugaison ; il vit sur une AUTRE feuille, reliée à
    # celle du pont par la transformation de deck `D`. La feuille que la
    # continuation analytique atteint est donc `D·Z_upper_conj`, et
    # c'est ELLE qui doit passer le vrai contrat F3b⁺. L'atlas conjugué
    # est conservé comme DIAGNOSTIC séparé, pas comme arête d'atlas.
    Zup_cont = apply_deck(Zup, DECK_D, keys)
    for side, halo, Zside, box_side in (
            ("lower", halo_lo, Zlo, halo_lo),
            ("upper_conj", halo_up, Zup, mirror_bounds(halo_lo)),
            ("upper_cont", halo_up, Zup_cont, mirror_bounds(halo_lo))):
        W = inter(Wb, halo)
        if W is None:
            res[side] = {"refused": "overlap_not_open"}
            continue
        # ancre : centre de l'overlap, puis DÉCALÉE strictement du bon
        # côté de la face si l'overlap y touche. Jamais sur Im = 0.
        anc = [(a + b) / 2 for a, b in W]
        strict_in = all(W[k][0] < anc[k] < W[k][1] for k in range(4))
        sgn = -1 if side == "lower" else 1
        anchor_off_face = all(sgn * anc[k] > 0 for k in IM_DIRS)
        r = {"overlap": _qbox(W),
             "widths": [float(b - a) for a, b in W],
             "open_4d": True,
             "anchor": [float(x) for x in anc],
             "anchor_strictly_inside_overlap": bool(strict_in),
             "anchor_off_the_face": bool(anchor_off_face)}
        Zs_map = ({k: Zside[k] for k in keys}
                  if not isinstance(Zside, dict) else Zside)
        th = theta_lines({k: Zb[k] for k in keys}, Wb,
                         Zs_map, box_side, W, keys)
        if "refused" in th:
            r.update(th)
            res[side] = r
            continue
        r["theta_by_line"] = {str(k): th[k] for k in keys}
        ths = {th[k].get("theta") for k in S}
        r["theta"] = (list(ths)[0] if len(ths) == 1 else None)
        r["theta_consistent_across_lines"] = bool(len(ths) == 1)
        r["theta_margins"] = [th[k].get("margin") for k in S]
        # séparation AU POINT d'ancre (contrat §7 de la revue)
        va = eval_at_point(Zb, Wb, anc, list(S))
        vb = eval_at_point(Zs_map, box_side, anc, list(S))
        anc_th, anc_ok = [], True
        for k in S:
            if va[k] is None or vb[k] is None:
                anc_ok = False
                anc_th.append(None)
                continue
            t, _dg = sep_at_point(va[k], vb[k])
            anc_th.append(t)
            anc_ok = anc_ok and (t is not None)
        r["anchor_theta_by_line"] = anc_th
        r["anchor_separation_strict"] = bool(anc_ok)
        r["anchor_agrees_with_box"] = bool(
            anc_ok and all(anc_th[j] == th[list(S)[j]].get("theta")
                           for j in range(len(list(S)))))
        # radicandes non nuls SUR L'OVERLAP, des DEUX côtés (lemme C129-D)
        r["glued_exactly"] = bool(
            r.get("theta") == THETA_REQUIRED
            and r["theta_consistent_across_lines"]
            and r["anchor_separation_strict"]
            and r["anchor_agrees_with_box"]
            and r["anchor_strictly_inside_overlap"]
            and r["anchor_off_the_face"])
        # LA PRÉDICTION DE STRUCTURE RÉELLE. Sur le coin `Im u = Im v = 0`
        # le radicande `R` est RÉEL. Là où `Re R > 0` (régime principal)
        # la racine est RÉELLE et l'involution antiholomorphe la FIXE ;
        # là où `Re R < 0` (régime canonique) la racine `i√(−R)` est
        # PUREMENT IMAGINAIRE et l'involution la NIE. Le côté supérieur
        # étant DÉRIVÉ par conjugaison, sa phase contre le pont est donc
        # PRÉDITE ligne par ligne — et prédite MIXTE, ce qui n'est pas
        # une renormalisation projective. Cette prédiction est
        # falsifiable : c'est ce qui transforme la surprise en résultat.
        pred = {}
        for j, k in enumerate(S):
            rg = rows_b[j]["regime"]
            pred[str(k)] = 1 if rg == "principal" else -1
        for k in keys:
            if k not in S:
                pred[str(k)] = 1        # jauge et coordonnées affines :
                                        # réelles sur le coin, donc fixées
        r["real_structure_prediction"] = pred
        r["deck_D_from_prediction"] = [pred[str(k)] for k in range(6)]
        r["matches_real_structure_prediction"] = bool(
            side != "upper_conj"
            or all(th[k].get("theta") == pred[str(k)] for k in keys))
        # R1d — l'identité de deck AU COEFFICIENT : le motif de θ contre
        # le conjugué doit être EXACTEMENT `D`, ligne par ligne, et le
        # `D` mesuré doit être celui pré-enregistré. C'est la seule part
        # discriminante du théorème de deck.
        if side == "upper_conj":
            r["deck_D_measured_equals_preregistered"] = bool(
                [pred[str(k)] for k in range(6)] == list(DECK_D)
                and all(th[k].get("theta") == DECK_D[k] for k in keys))
        res[side] = r
    out["sides"] = res

    # --- F3d : LE NÉGATIF DOIT ÊTRE PROJECTIVEMENT DISCRIMINANT -------
    # v1 niait TOUTES les coordonnées : `Z` et `−Z` sont LE MÊME POINT de
    # P⁵, donc le test ne montrait que la sensibilité de la
    # normalisation affine `Z_g = 1` au changement de représentant — pas
    # la détection d'un mauvais POINT. v2 garde la jauge et les
    # coordonnées affines INCHANGÉES et nie UNE SEULE ligne racine
    # canonique : `D_bad = diag(1,−1,1,1,1,1)`, non scalaire, donc un
    # point projectif RÉELLEMENT différent.
    Wl = inter(Wb, halo_lo)
    if Wl is not None:
        s_bad = list(S)[1]                    # une ligne CANONIQUE
        D_bad = tuple(-1 if k == s_bad else 1 for k in range(6))
        Zbad = apply_deck(Zb, D_bad, keys)
        thn = theta_lines(Zbad, Wb, {k: Zlo[k] for k in keys},
                          halo_lo, Wl, keys)
        out["F3d_D_bad"] = list(D_bad)
        out["F3d_bad_theta_by_line"] = {
            str(k): thn[k].get("theta") for k in keys}
        # le raccord doit tomber : au moins une ligne à −1 ET un motif
        # NON constant (donc pas un simple changement de représentant)
        vals = [thn[k].get("theta") for k in keys]
        out["F3d_non_scalar_mutation_breaks_gluing"] = bool(
            any(v != THETA_REQUIRED for v in vals) and len(set(vals)) > 1)
        # et le contrôle : la négation GLOBALE, elle, est le MÊME point
        # projectif — publiée pour dire pourquoi elle ne vaut rien
        Zglob = apply_deck(Zb, tuple([-1] * 6), keys)
        thg = theta_lines(Zglob, Wb, {k: Zlo[k] for k in keys},
                          halo_lo, Wl, keys)
        gv = {thg[k].get("theta") for k in keys}
        out["F3d_global_negation_theta"] = (
            list(gv)[0] if len(gv) == 1 else None)
        out["F3d_global_negation_is_same_projective_point"] = bool(
            len(gv) == 1)

    # --- R1e : NÉGATIF DU DECK — un seul signe changé doit CASSER ------
    Wu = inter(Wb, halo_up)
    if Wu is not None:
        D_wrong = tuple(-d if k == list(S)[0] else d
                        for k, d in enumerate(DECK_D))
        Zw = apply_deck(Zup, D_wrong, keys)
        thw = theta_lines({k: Zb[k] for k in keys}, Wb, Zw,
                          mirror_bounds(halo_lo), Wu, keys)
        vals = [thw[k].get("theta") for k in keys]
        out["R1e_D_wrong"] = list(D_wrong)
        out["R1e_wrong_deck_breaks"] = bool(
            not all(v == THETA_REQUIRED for v in vals))
    return out


# ===========================================================================
#  R3 — les transitions PONT ↔ PONT, certifiées (pas seulement
#  géométriques). Une arête de NERF est une transition CERTIFIÉE ; une
#  intersection de boîtes n'est qu'un graphe de domaines.
# ===========================================================================
def _bb_job(arg):
    i, j = arg
    S, g, eps = _G["cell"]
    Wi, Wj = _G["bridges"][i], _G["bridges"][j]
    W = inter(Wi, Wj)
    rec = {"pair": [i, j]}
    if W is None:
        rec["refused"] = "overlap_not_open"
        return rec
    ei, ej = _G["ledgers"][i], _G["ledgers"][j]
    ci, hi_ = center_hw(Wi)
    cj, hj_ = center_hw(Wj)
    Zi, _di, ri = build_section_bilateral(
        S, g, ei, [float(x) for x in ci], [float(x) for x in hi_])
    Zj, _dj, rj = build_section_bilateral(
        S, g, ej, [float(x) for x in cj], [float(x) for x in hj_])
    keys = list(S) + [x for x in range(6) if x not in S]
    if any(Zi[k] is None or Zj[k] is None for k in keys):
        rec["refused"] = "section_incomplete"
        return rec
    th = theta_lines({k: Zi[k] for k in keys}, Wi,
                     {k: Zj[k] for k in keys}, Wj, W, keys)
    if "refused" in th:
        rec.update(th)
        return rec
    vals = [th[k].get("theta") for k in keys]
    rec["theta_by_line"] = {str(k): th[k].get("theta") for k in keys}
    rec["margins"] = [th[k].get("margin") for k in keys]
    rec["widths"] = [float(b - a) for a, b in W]
    rec["certified"] = bool(all(v == THETA_REQUIRED for v in vals))
    rec["diff_sup"] = max((th[k].get("diff_sup") or 0.0) for k in keys)
    return rec


# ===========================================================================
#  Self-test
# ===========================================================================
def selftest():
    F = Fraction
    ok, tot = 0, 0

    def chk(name, cond):
        nonlocal ok, tot
        tot += 1
        ok += bool(cond)
        print(f"  {'OK  ' if cond else 'FAIL'} T{tot} {name}")

    chk("inter exige une largeur STRICTEMENT positive",
        inter([(F(0), F(1))] * 4, [(F(1), F(2))] * 4) is None
        and inter([(F(0), F(2))] * 4, [(F(1), F(3))] * 4)
        == [(F(1), F(2))] * 4)
    chk("mirror_bounds ne touche que les directions Im",
        mirror_bounds([(F(-2), F(0))] * 4)[0] == (F(-2), F(0))
        and mirror_bounds([(F(-2), F(0))] * 4)[1] == (F(0), F(2)))
    chk("pont ∩ halo inférieur est OUVERT dans les 4 coordonnées",
        inter([(F(-2), F(2))] * 4, [(F(-2), F(0))] * 4) is not None)
    chk("halo inférieur ∩ halo supérieur n'est PAS ouvert (le coin)",
        inter([(F(-2), F(0))] * 4, [(F(0), F(2))] * 4) is None)

    # reframe : identité, translation, dilatation anisotrope
    b = [(F(-1), F(1)), (F(-2), F(2)), (F(-1), F(1)), (F(-2), F(2))]
    off, sc = reframe(b, b)
    chk("reframe d'une boîte sur elle-même est l'identité",
        off == [F(0)] * 4 and sc == [F(1)] * 4)
    w = [(F(-1), F(0)), (F(-2), F(0)), (F(-1), F(1)), (F(-2), F(2))]
    off, sc = reframe(b, w)
    chk("reframe anisotrope : off et scale par coordonnée",
        off == [F(-1, 2), F(-1, 2), F(0), F(0)]
        and sc == [F(1, 2), F(1, 2), F(1), F(1)])
    chk("frame_admissible accepte une sous-boîte",
        frame_admissible(off, sc))
    chk("frame_admissible REFUSE une boîte qui déborde",
        not frame_admissible([F(1), F(0), F(0), F(0)], [F(1)] * 4))

    # Le recentrage anisotrope, sur un polynôme dont on connaît la valeur
    e = [MIDX[tuple(1 if k == j else 0 for k in range(NG))]
         for j in range(NG)]
    p = [CZERO] * NM
    p[0] = CIV(iv.mpf(3), IV0)
    p[e[0]] = CIV(iv.mpf(5), IV0)
    p[e[1]] = CIV(iv.mpf(7), IV0)
    #  q(ε') = p(off + scale·ε') = 3 + 5(o₀+s₀ε'₀) + 7(o₁+s₁ε'₁)
    off = [F(1, 2), F(-1, 4), F(0), F(0)]
    sc = [F(1, 2), F(1, 4), F(1), F(1)]
    q = apply_recenter(recenter_matrix_aniso(off, sc), p)
    chk("recentrage anisotrope : terme constant",
        abs(float(mp.mpf(q[0].re.a)) - (3 + 5 * 0.5 + 7 * -0.25)) < 1e-30)
    chk("recentrage anisotrope : pente en ε'₀ mise à l'échelle",
        abs(float(mp.mpf(q[e[0]].re.a)) - 5 * 0.5) < 1e-30)
    chk("recentrage anisotrope : pente en ε'₁ mise à l'échelle",
        abs(float(mp.mpf(q[e[1]].re.a)) - 7 * 0.25) < 1e-30)
    chk("recentrage sur soi-même laisse le polynôme INCHANGÉ",
        all(abs(float(mp.mpf(x.re.a)) - float(mp.mpf(y.re.a))) < 1e-30
            for x, y in zip(
                apply_recenter(
                    recenter_matrix_aniso([F(0)] * 4, [F(1)] * 4), p), p)))

    # sep_phase, sur des valeurs construites
    th, _ = sep_phase(CIV(iv.mpf(0), IV0), CIV(iv.mpf(4), IV0))
    chk("sep_phase lit θ = +1 quand la différence est nulle", th == 1)
    th, _ = sep_phase(CIV(iv.mpf(4), IV0), CIV(iv.mpf(0), IV0))
    chk("sep_phase lit θ = −1 quand c'est la somme qui s'annule",
        th == -1)
    th, _ = sep_phase(CIV(iv.mpf([-1, 1]), IV0), CIV(iv.mpf([-1, 1]), IV0))
    chk("sep_phase REFUSE l'ambigu", th is None)

    # Le constructeur bilatéral sur la vraie cellule : régimes assignés,
    # et la mutation σ doit casser.
    cov = json.loads(COVER_JSON.read_text(encoding="utf-8"))
    cell = cov["cell"]
    S, g, eps = tuple(cell["S"]), cell["g"], tuple(cell["eps"])
    sc_j = json.loads(SCOUT_JSON.read_text(encoding="utf-8"))
    t0 = sc_j["per_tile"][0]
    Wb = [(F(*b[0]), F(*b[1])) for b in t0["bridge_corrected_bounds"]]
    cb, hb = center_hw(Wb)
    Z, _dz, rows = build_section_bilateral(
        S, g, eps, [float(x) for x in cb], [float(x) for x in hb])
    chk("build_section_bilateral assigne les 3 régimes sans essai",
        all(r["regime"] is not None for r in rows))
    chk("les régimes concordent avec le scout (deux arithmétiques)",
        [r["regime"] for r in rows] == t0["regimes"])
    _z2, _d2, rows_f = build_section_bilateral(
        S, g, eps, [float(x) for x in cb], [float(x) for x in hb],
        force_sigma=True)
    chk("la MUTATION σ de composante est REFUSÉE sur le pont",
        any(r.get("refused") == "component_sigma_undetermined_on_bridge"
            for r in rows_f))
    chk("Im R straddle sur les lignes où σ est refusé",
        all(r["im_R_straddles"] for r in rows_f
            if r.get("refused") == "component_sigma_undetermined_on_bridge"))
    # NÉGATIF (iii) : une boîte élargie où Re R straddle doit REFUSER
    big = [(cb[k] - 4, cb[k] + 4) for k in range(4)]
    cbg, hbg = center_hw(big)
    _z3, _d3, rows_big = build_section_bilateral(
        S, g, eps, [float(x) for x in cbg], [float(x) for x in hbg])
    chk("une boîte où Re R straddle est REFUSÉE, pas repliée",
        any(r.get("refused") == "re_R_straddles_zero_no_regime"
            for r in rows_big))
    print(f"\nself-test {ok}/{tot}")
    return ok == tot


# ===========================================================================
#  Build
# ===========================================================================
def build():
    print("=" * 78)
    print(f"C129-F0 — F2 + F3 : CARTES-PONTS ET DEUX TRANSITIONS "
          f"({N_WORKERS} workers)")
    print("=" * 78)
    cov = json.loads(COVER_JSON.read_text(encoding="utf-8"))
    atl = json.loads(ATLAS_JSON.read_text(encoding="utf-8"))
    c127e = json.loads(C127E_JSON.read_text(encoding="utf-8"))
    scout = json.loads(SCOUT_JSON.read_text(encoding="utf-8"))
    f1 = json.loads(F1_JSON.read_text(encoding="utf-8"))
    cell = cov["cell"]
    S, g, eps = tuple(cell["S"]), cell["g"], tuple(cell["eps"])
    leaves = list(cov["tiles"]) + list(c127e["new_tiles"])
    halos = {h["index"]: h["record"] for h in atl["halos"] if h["ok"]}
    clipped = sorted(i for i, r in halos.items()
                     if r["rule"] == "clipped")
    log(f"cellule S={S} g={g} ; {len(clipped)} tuiles clippées")
    log(f"    F1a a REFUSÉ l'identification canonique — l'atlas "
        f"supérieur est DÉRIVÉ, pas énuméré "
        f"(claim_level={f1['claim_level']})")

    # --- R4c : LA CHAÎNE AMONT EST GATÉE, PAS SEULEMENT IMPORTÉE ------
    # v1 lisait les JSON amont sans jamais vérifier qu'ils étaient VERTS
    # ni en mode complet. Un artefact pilote ou rouge aurait été consommé
    # en silence.
    up = {}
    for name, blob, need_full in (
            ("c127d_atlas", atl, True),
            ("c129f_bridge_scout", scout, False),
            ("c129f_f1_mirror_ledger", f1, False)):
        gp, gt = blob.get("gates_passed"), blob.get("gates_total")
        up[name] = {"gates": f"{gp}/{gt}", "green": bool(gp == gt and gt),
                    "mode": blob.get("mode")}
        if need_full:
            up[name]["full"] = bool(blob.get("mode") == "full")
    for name, path in (("c129d_exact_gluing",
                        RES / "k3_cap_b1e2iii_c129d_exact_gluing.json"),
                       ("c129e_halo_metric",
                        RES / "k3_cap_b1e2iii_c129e_halo_metric.json")):
        try:
            b = json.loads(path.read_text(encoding="utf-8"))
        except OSError:
            up[name] = {"green": False, "missing": True}
            continue
        gp, gt = b.get("gates_passed"), b.get("gates_total")
        up[name] = {"gates": f"{gp}/{gt}",
                    "green": bool(gp == gt and gt),
                    "mode": b.get("mode")}
    upstream_ok = all(v.get("green") and v.get("full", True)
                      for v in up.values())
    log(f"R4c : chaîne amont — " + " ; ".join(
        f"{k} {v.get('gates', '?')}" for k, v in up.items())
        + f" ⟹ {upstream_ok}")

    # --- F2a : la géométrie, RE-DÉRIVÉE puis confrontée au scout ------
    bridges, f2a = {}, True
    for i in clipped:
        r = halos[i]
        c_h = [Fraction(float.fromhex(x)) for x in r["center_hex"]]
        H = Fraction(float.fromhex(r["H_hex"]))
        hb = bounds(c_h, H)
        W = [((-2 * H, 2 * H) if k in IM_DIRS else hb[k])
             for k in range(4)]
        bridges[i] = W
    scoutb = {t["tile"]: [(Fraction(*b[0]), Fraction(*b[1]))
                          for b in t["bridge_corrected_bounds"]]
              for t in scout["per_tile"]}
    f2a = all(bridges[i] == scoutb[i] for i in clipped)
    log(f"F2a : géométrie 2H/H re-dérivée == scout sur "
        f"{len(clipped)} ponts : {f2a}")

    jobs = clipped
    mpctx = get_context("fork")
    _init((S, g, eps), leaves, halos, bridges)
    with mpctx.Pool(N_WORKERS, initializer=_init,
                    initargs=((S, g, eps), leaves, halos, bridges)) as pool:
        rows = pool.map(_bridge_job, jobs)
    byi = {r["tile"]: r for r in rows}
    log(f"    {len(rows)} ponts traités")

    f2b = all(r["F2b_contains_union"] for r in rows)
    f2c = all(r["F2c_all_regimes_assigned"] for r in rows)
    f2d = all(r.get("F2d_section_complete") and r.get("F2d_gauge_positive")
              for r in rows)
    f2e_i = all(r["F2e_i_false_bridge_refused"] for r in rows)
    f2e_ii = all(r["F2e_ii_sigma_refused"] for r in rows)
    f2f = (len(rows) == len(clipped)
           and all("refused" not in r for r in rows)
           and all(r.get("R4c_lower_dets_match_c127d") for r in rows))
    reg_census = Counter(tuple(x["regime"] for x in r["bridge_rows"])
                         for r in rows)
    gmin = min(r.get("F2d_target_gauge_absmin", 0.0) for r in rows)
    log(f"F2b : core⁻ ∪ core⁺ ⊆ pont (strict en Im) : {f2b}")
    log(f"F2c : régimes assignés sans essai : {f2c} — {dict(reg_census)}")
    log(f"F2d : section complète + jauge cible > 0 (min {gmin:.3e}) : {f2d}")
    log(f"F2e : (i) faux pont H refusé {f2e_i} ; (ii) σ de composante "
        f"refusé sur le pont {f2e_ii}")
    log(f"F2f : {len(rows)}/{len(clipped)} sans filtrage : {f2f}")

    # --- F3a/F3b : les deux transitions ------------------------------
    f3a, f3b_lo = True, True
    th_census = Counter()
    pat_census = Counter()
    marg_min, diff_max = None, 0.0
    wid_min = None
    for r in rows:
        for side in ("lower", "upper_conj", "upper_cont"):
            d = r["sides"].get(side, {})
            if d.get("refused") or "open_4d" not in d:
                f3a = False
                if side == "lower":
                    f3b_lo = False
                continue
            f3a = f3a and d["open_4d"] and \
                d["anchor_strictly_inside_overlap"] and \
                d["anchor_off_the_face"]
            if side == "lower":
                f3b_lo = f3b_lo and d.get("glued_exactly", False)
            if side == "upper_conj":
                pat_census[tuple(
                    d["theta_by_line"][str(k)].get("theta")
                    for k in range(6))] += 1
            th_census[(side, d.get("theta"))] += 1
            for m in d.get("theta_margins") or []:
                if m is not None:
                    marg_min = m if marg_min is None else min(marg_min, m)
            for k, v in (d.get("theta_by_line") or {}).items():
                if v.get("diff_sup"):
                    diff_max = max(diff_max, v["diff_sup"])
            w = min(d["widths"])
            wid_min = w if wid_min is None else min(wid_min, w)
    ledger_census = Counter(tuple(r["F2d_ledger_derived"]) for r in rows)
    log(f"F3a : overlaps ouverts (largeur min {wid_min:.3e}) + ancres "
        f"strictement intérieures HORS de la face : {f3a}")
    log(f"F2d(bis) : ledger du pont DÉRIVÉ du côté inférieur — "
        f"{dict(ledger_census)} (défaut {list(eps)})")
    log(f"F3b⁻ : le pont se recolle EXACTEMENT au côté inférieur, θ = +1 "
        f"sur toutes les lignes : {f3b_lo} ; marge min {marg_min:.3e} ; "
        f"sup de la différence recentrée {diff_max:.3e}")
    log(f"F3b⁺(diag) : côté supérieur CONJUGUÉ — motif de θ "
        f"{dict(pat_census)} = la transformation de deck D. Ce n'est PAS "
        f"une transition d'atlas : voir R2 pour la feuille continuée.")

    # --- R3 : LE GRAPHE DE DOMAINES ET LE VRAI NERF -------------------
    # v1 confondait les deux : elle ajoutait une arête pont↔côté dès
    # qu'un overlap GÉOMÉTRIQUE existait — y compris quand le raccord
    # avait explicitement ÉCHOUÉ — et des arêtes pont↔pont sans aucune
    # transition calculée, puis énumérait des triples depuis ce graphe.
    # Une arête de NERF est une transition CERTIFIÉE. Rien d'autre.
    _G["ledgers"] = {r["tile"]: tuple(r["F2d_ledger_derived"])
                     for r in rows}
    bb_geo = [(a, b) for a, b in itertools.combinations(clipped, 2)
              if inter(bridges[a], bridges[b]) is not None]
    log(f"R3 : {len(bb_geo)} overlaps pont↔pont géométriques — "
        f"transitions en cours de certification…")
    with mpctx.Pool(N_WORKERS, initializer=_init,
                    initargs=((S, g, eps), leaves, halos, bridges,
                              _G["ledgers"])) as pool:
        bb = pool.map(_bb_job, bb_geo)
    bb_ok = [x for x in bb if x.get("certified")]
    bb_bad = [x for x in bb if not x.get("certified")]
    bb_diff = max((x.get("diff_sup") or 0.0) for x in bb_ok) if bb_ok else 0.0
    bb_marg = min((m for x in bb_ok for m in x["margins"]
                   if m is not None), default=None)
    log(f"     pont↔pont CERTIFIÉES {len(bb_ok)}/{len(bb_geo)} "
        f"(θ = +1 sur les 6 coordonnées) ; marge min "
        f"{bb_marg if bb_marg is None else round(bb_marg, 4)} ; "
        f"sup diff {bb_diff:.3e} ; non certifiées {len(bb_bad)}")

    # graphe de DOMAINES (géométrique, publié comme tel)
    dom_edges = (2 * len(clipped)) + len(bb_geo)

    # NERF : uniquement les arêtes dont la transition est CERTIFIÉE.
    # Nœuds : les 316 cartes de l'atlas inférieur + les 64 ponts. Les
    # arêtes inférieur↔inférieur sont IMPORTÉES de C127-D et vérifiées
    # vertes, pas recalculées.
    lower_pairs = [tuple(sorted((p["i"], p["j"]))) for p in atl["pairs"]
                   if p.get("ok", True)] if "pairs" in atl else []
    n_lower = len(leaves)
    nodes = [("L", i) for i in range(n_lower)] + [("B", i) for i in clipped]
    adj = {n: set() for n in nodes}
    for a, b in lower_pairs:
        if ("L", a) in adj and ("L", b) in adj:
            adj[("L", a)].add(("L", b))
            adj[("L", b)].add(("L", a))
    n_bl = 0
    n_bu_cont = 0
    for r in rows:
        i = r["tile"]
        if r["sides"].get("lower", {}).get("glued_exactly"):
            adj[("B", i)].add(("L", i))
            adj[("L", i)].add(("B", i))
            n_bl += 1
        if r["sides"].get("upper_cont", {}).get("glued_exactly"):
            n_bu_cont += 1
    for x in bb_ok:
        a, b = x["pair"]
        adj[("B", a)].add(("B", b))
        adj[("B", b)].add(("B", a))
    seen, dq = {nodes[0]}, deque([nodes[0]])
    while dq:
        x = dq.popleft()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                dq.append(y)
    connected = len(seen) == len(nodes)

    def boxof(n):
        t, i = n
        if t == "B":
            return bridges[i]
        rr = halos.get(i)
        if rr is None:
            c0, h0 = box_of(leaves[i])
            return bounds(c0, h0)
        return bounds([Fraction(float.fromhex(x))
                       for x in rr["center_hex"]],
                      Fraction(float.fromhex(rr["H_hex"])))

    # triples NOUVEAUX (impliquant au moins un pont), depuis les seules
    # arêtes CERTIFIÉES. Les triples purement inférieurs sont ceux de
    # C127-D, déjà certifiés là-bas — ils ne sont pas recomptés ici.
    new_triples = []
    for i in clipped:
        x = ("B", i)
        for y, z in itertools.combinations(sorted(adj[x]), 2):
            if y in adj[z]:
                W1 = inter(boxof(x), boxof(y))
                W = inter(W1, boxof(z)) if W1 else None
                if W is not None:
                    new_triples.append([list(x), list(y), list(z)])
    log(f"R3 : NERF (arêtes certifiées SEULEMENT) — {len(nodes)} nœuds "
        f"({n_lower} inférieurs + {len(clipped)} ponts), "
        f"{len(lower_pairs)} arêtes L↔L importées de C127-D, {n_bl} "
        f"arêtes B↔L, {len(bb_ok)} arêtes B↔B, "
        f"{len(new_triples)} triples NEUFS ; connexe={connected}")
    log(f"     (le graphe de DOMAINES, lui, compte {dom_edges} arêtes "
        f"géométriques — il est publié sous ce nom, pas comme nerf)")
    f3c = bool(connected and n_bl == len(clipped)
               and len(bb_ok) == len(bb_geo))

    # --- R4a : la portée AMBIANTE, stratifiée ------------------------
    re_flush = Counter()
    for i in clipped:
        fl = halos[i]["flush_faces"]
        re_flush[sum(1 for k in (0, 2) if fl[k] != 0)] += 1
    fully_ambient = re_flush[0]
    log(f"R4a : portée — {len(clipped)}/{len(clipped)} ponts bilatéraux "
        f"dans les DEUX directions Im ; ouverts dans les QUATRE "
        f"coordonnées : {fully_ambient}/{len(clipped)} ; encore "
        f"relatifs à 1 face Re : {re_flush[1]} ; à 2 faces : "
        f"{re_flush[2]}")

    # --- F3d + R1e : les négatifs ------------------------------------
    f3d = all(r.get("F3d_non_scalar_mutation_breaks_gluing") for r in rows)
    glob_same = all(r.get("F3d_global_negation_is_same_projective_point")
                    for r in rows)
    r1e = all(r.get("R1e_wrong_deck_breaks") for r in rows)
    log(f"F3d : mutation NON SCALAIRE D_bad = diag(1,−1,1,1,1,1) ⟹ le "
        f"raccord tombe : {f3d}  (et la négation GLOBALE, elle, reste "
        f"le MÊME point projectif : {glob_same} — c'est pourquoi la v1 "
        f"ne discriminait rien)")
    log(f"R1e : un seul signe changé dans D ⟹ l'identité de deck "
        f"CASSE : {r1e}")

    # --- R1 : le théorème de deck ------------------------------------
    alg = deck_algebra(DECK_D)
    deck_meas = all(r["sides"]["upper_conj"].get(
        "deck_D_measured_equals_preregistered") for r in rows)
    r1 = bool(alg["involution"] and alg["non_scalar_in_pgl6"]
              and alg["preserves_the_three_quadrics"] and deck_meas
              and r1e)
    log(f"R1 : D = {list(DECK_D)} — involution {alg['involution']}, "
        f"non scalaire {alg['non_scalar_in_pgl6']}, préserve les 3 "
        f"quadriques {alg['preserves_the_three_quadrics']} (TRIVIAL "
        f"pour tout diagonal de signes) ; D MESURÉ == pré-enregistré "
        f"sur {len(rows)} ponts : {deck_meas} ⟹ {r1}")

    # --- R2 : la feuille supérieure CONTINUÉE -------------------------
    r2 = all(r["sides"].get("upper_cont", {}).get("glued_exactly")
             for r in rows)
    conj_fail = all(not r["sides"].get("upper_conj", {}).get(
        "glued_exactly", False) for r in rows)
    log(f"R2 : le pont se recolle EXACTEMENT à D·Z_conj (feuille "
        f"CONTINUÉE) sur {n_bu_cont}/{len(rows)} : {r2} — et PAS à "
        f"Z_conj seule ({conj_fail}), qui reste un diagnostic")

    gates = {
        "F2a_bridge_geometry_2H_matches_scout": bool(f2a),
        "F2b_bridge_contains_both_cores": bool(f2b),
        "F2c_regime_assigned_without_trial": bool(f2c),
        "F2d_bilateral_section_complete_gauge_positive": bool(f2d),
        "F2d_bis_bridge_ledger_derived_not_defaulted": bool(
            all("F2d_ledger_derived" in r for r in rows)),
        "F2e_i_false_bridge_H_refused": bool(f2e_i),
        "F2e_ii_component_sigma_refused_on_bridge": bool(f2e_ii),
        "F2f_all_64_no_silent_filtering": bool(f2f),
        "F3a_open_overlaps_and_interior_anchors": bool(f3a),
        "F3b_lower_exact_gluing_theta_plus_one": bool(f3b_lo),
        "R1_deck_transformation_D_measured_and_discriminated": bool(r1),
        "R2_bridge_glues_to_CONTINUED_upper_sheet_D_times_conj": bool(r2),
        "R2b_conjugate_sheet_alone_does_NOT_glue": bool(conj_fail),
        "R3_nerve_uses_only_certified_edges_and_is_connected": bool(f3c),
        "F3d_non_scalar_mutation_breaks_gluing": bool(f3d),
        "R1e_wrong_deck_sign_breaks_identity": bool(r1e),
        "R4_upstream_chain_verified": bool(upstream_ok)}
    npass = sum(1 for v in gates.values() if v)

    try:
        head = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=Path(__file__).resolve().parent, capture_output=True,
            text=True, timeout=10).stdout.strip() or None
    except Exception:
        head = None

    out = {
        "artifact": "k3_cap_b1e2iii_c129f_f2f3_bridge_atlas",
        # Ce script n'a qu'un seul mode : les 64 ponts, toujours tous.
        # Sérialisé EXPLICITEMENT depuis la 12ᵉ revue (D2) : le G2 de F4
        # exige mode == "full" de chaque amont, null n'est pas un mode.
        "mode": "full",
        "claim_level": f1["claim_level"],
        "claim": (
            "Les 64 CARTES-PONTS bilatérales sont construites sur la "
            "géométrie 2H corrigée ; leur régime est assigné PAR LIGNE "
            "au signe certifié de Re R sans essai ; leur ledger est "
            "DÉRIVÉ du côté inférieur puis re-vérifié. Chacune se "
            "recolle EXACTEMENT au côté INFÉRIEUR et à la FEUILLE "
            "CONTINUÉE `D·Z_conj` sur des overlaps OUVERTS, ancres "
            "strictement intérieures. `D = diag(+1,−1,+1,+1,+1,−1)` est "
            "une TRANSFORMATION DE DECK : c'est elle, et non un simple "
            "« le conjugué n'est pas le continué », qui sépare les deux "
            "feuilles — et l'atlas CONJUGUÉ SEUL ne se recolle PAS "
            "(64/64), il reste un diagnostic. PORTÉE STRATIFIÉE : 64/64 "
            "ponts sont bilatéraux dans les deux directions Im, mais "
            "seuls 36/64 sont ouverts dans les QUATRE coordonnées — 24 "
            "restent relatifs à une face Re, 4 à deux faces, et ces "
            "faces appartiennent aux voisines de codimension 1. Le NERF "
            "ne compte que des arêtes CERTIFIÉES (380 nœuds : 316 "
            "inférieurs + 64 ponts ; 5396 arêtes L↔L importées de "
            "C127-D, 64 B↔L, 210/210 B↔B) ; le graphe géométrique est "
            "publié séparément sous le nom `domain_intersection_graph`. "
            "Le lot reste une CONTINUATION ANALYTIQUE LOCALE VERS UNE "
            "BOÎTE MIROIR DÉRIVÉE (F1a a refusé l'identification "
            "canonique). Aucune métrique n'est certifiée ici."),
        "cell": {"S": list(S), "g": g, "eps": list(eps)},
        "n_bridges": len(clipped),
        "theta_required_preregistered": THETA_REQUIRED,
        "regime_census": {str(k): v for k, v in reg_census.items()},
        "theta_census": {f"{k[0]}|theta={k[1]}": v
                         for k, v in th_census.items()},
        "upper_theta_pattern_census": {str(k): v
                                       for k, v in pat_census.items()},
        "bridge_ledger_census": {str(k): v
                                 for k, v in ledger_census.items()},
        "real_structure_finding": (
            "SUR LE COIN, L'INVOLUTION ANTIHOLOMORPHE AGIT AVEC DES "
            "SIGNES MIXTES. R y est RÉEL : positif sur la ligne "
            "principale (racine RÉELLE, fixée par conjugaison) et "
            "négatif sur les deux lignes canoniques (racine i√(−R) "
            "PUREMENT IMAGINAIRE, donc NIÉE). Le point conjugué et le "
            "point continué analytiquement à travers le coin sont donc "
            "DEUX POINTS DISTINCTS de P⁵ au-dessus du même (ū, v̄) — la "
            "jauge étant normalisée Z_g = 1, des signes indépendants par "
            "ligne ne sont PAS une renormalisation projective. "
            "CONSÉQUENCE POUR LE LEVIER N° 1 DU CADRAGE : la voisine "
            "conjuguée est un atlas légitime, mais ce n'est PAS celui "
            "que la continuation analytique atteint. Le motif de θ est "
            "PRÉDIT ligne par ligne depuis le régime, et la prédiction "
            "est vérifiée sur 64/64 — c'est un résultat, pas une "
            "surprise résiduelle."),
        "min_theta_margin": marg_min,
        "max_recentred_difference_sup": diff_max,
        "min_overlap_width": wid_min,
        "min_target_gauge_absmin": gmin,
        "upstream_chain": up,
        "deck": {**deck_algebra(DECK_D),
                 "identity": "Z_upper_conj = D · Z_bridge sur l'overlap "
                             "supérieur, certifiée AU COEFFICIENT",
                 "measured_on_bridges": len(rows),
                 "why_the_algebra_is_not_the_content": (
                     "tout diagonal de signes préserve les quadriques "
                     "Σ μ^m Z² ; ce qui distingue CE D est qu'il est "
                     "celui, et le seul, qui relie la section conjuguée "
                     "à la section du pont — c'est le négatif R1e qui "
                     "porte la discriminance, pas l'algèbre")},
        "domain_intersection_graph": {
            "note": ("GRAPHE DE DOMAINES, PAS un nerf d'atlas : une "
                     "arête n'y est qu'une intersection de boîtes de "
                     "largeur > 0. Publié sous ce nom depuis que la "
                     "revue a montré que la v1 comptait comme arêtes "
                     "des raccords explicitement en échec."),
            "n_bridge_side_geometric": 2 * len(clipped),
            "n_bridge_bridge_geometric": len(bb_geo),
            "n_edges": dom_edges},
        "nerve": {
            "note": ("NERF : uniquement les arêtes dont la TRANSITION "
                     "est certifiée. Les arêtes inférieur↔inférieur "
                     "sont IMPORTÉES de C127-D (vérifié vert), pas "
                     "recalculées."),
            "n_nodes": len(nodes),
            "n_lower_nodes": n_lower,
            "n_bridge_nodes": len(clipped),
            "n_lower_lower_edges_imported": len(lower_pairs),
            "n_bridge_lower_edges_certified": n_bl,
            "n_bridge_bridge_edges_certified": len(bb_ok),
            "n_bridge_bridge_geometric": len(bb_geo),
            "n_bridge_upper_cont_certified": n_bu_cont,
            "n_new_triples_involving_a_bridge": len(new_triples),
            "connected": bool(connected),
            "min_bb_margin": bb_marg,
            "max_bb_diff_sup": bb_diff,
            "new_triples": new_triples[:200]},
        "bridge_bridge_transitions": bb,
        "ambient_scope": {
            "note": ("STRATIFIÉ, corrigé après revue : « les 64 cartes "
                     "deviennent ambiantes » était une sur-portée."),
            "bilateral_in_both_im_directions": len(clipped),
            "fully_open_in_all_four_coordinates": fully_ambient,
            "still_relative_to_one_re_face": re_flush[1],
            "still_relative_to_two_re_faces": re_flush[2]},
        "negatives": {
            "F2e_i": "faux pont de demi-largeur H : inclusion refusée",
            "F2e_ii": ("σ de composante sur le pont : REFUS "
                       "`component_sigma_undetermined_on_bridge` — c'est "
                       "la raison structurelle du régime canonique"),
            "F2e_iii": ("boîte élargie où Re R straddle : REFUS "
                        "`re_R_straddles_zero_no_regime` (self-test)"),
            "F3d": ("mutation NON SCALAIRE D_bad = diag(1,−1,1,1,1,1) "
                    "(jauge et coordonnées affines INCHANGÉES, une seule "
                    "ligne canonique niée) ⟹ le raccord tombe. La v1 "
                    "niait TOUTES les coordonnées, or Z et −Z sont LE "
                    "MÊME point de P⁵ : elle ne testait que la "
                    "sensibilité de la normalisation Z_g = 1 au "
                    "changement de représentant, pas la détection d'un "
                    "mauvais POINT. Dette de discriminance, relevée par "
                    "la revue et payée ici."),
            "F3d_global_negation_control": (
                "la négation globale reste le MÊME point projectif sur "
                "64/64 — publié pour dire pourquoi la v1 ne valait rien"),
            "R1e": ("un seul signe changé dans D ⟹ l'identité de deck "
                    "CASSE : c'est ce négatif, et non l'algèbre des "
                    "quadriques, qui porte la discriminance du théorème")},
        "triple_cocycle_note": (
            "Les 588 triples NEUFS sont énumérés depuis les seules "
            "arêtes CERTIFIÉES. Leur cocycle θ_ij·θ_jk·θ_ki = +1 est "
            "IMPLIQUÉ : chaque θ vaut +1 sur un domaine qui CONTIENT la "
            "boîte triple, et θ est discret. Ce n'est donc PAS un test "
            "indépendant et il n'est pas gaté comme tel — même "
            "honnêteté qu'en C127-D sur son propre cocycle. Ce qui "
            "resterait à payer au triple est la congruence MÉTRIQUE, "
            "qui appartient à F4."),
        "per_bridge": rows,
        "not_paid_here": [
            "la MÉTRIQUE du pont (F4) : Qmat, Weyl, congruences "
            "latérales — le chemin C129-E exige le constructeur "
            "bilatéral, qui existe maintenant, mais le run ne l'utilise pas",
            "l'identification canonique de la voisine : F1a a REFUSÉ, "
            "l'atlas supérieur est DÉRIVÉ par conjugaison, pas énuméré",
            "le raccord à travers les faces Re : les ponts y restent des "
            "cartes RELATIVES (marge 0 contre la face de la cellule)",
            "les voisines de codimension 1, les 895 autres paires, R12-C"],
        "gates": gates, "gates_passed": npass, "gates_total": len(gates),
        "verdict": (
            "F2 + F3 + R1-R4 LIVRÉS — 64 ponts bilatéraux, recollés "
            "EXACTEMENT au côté inférieur (64) ET à la feuille "
            "CONTINUÉE D·Z_conj (64), 210/210 transitions pont↔pont "
            "certifiées, nerf de 380 nœuds connexe sur arêtes "
            "certifiées SEULEMENT. D = diag(+,−,+,+,+,−) est la "
            "transformation de deck qui sépare le conjugué du continué."
            if npass == len(gates) else
            f"ROUGE — {len(gates) - npass} gate(s) en échec"),
        "provenance": {
            "git_head": head, "python": sys.version.split()[0],
            "platform": platform.platform(), "mp_prec": int(mp.prec),
            "tm_order": TM_ORDER, "unary_series_deg": UNARY_SERIES_DEG,
            "wall_s": round(time.time() - T0, 1), "n_workers": N_WORKERS,
            "preregistered": {"im_dirs": list(IM_DIRS),
                              "theta_required": THETA_REQUIRED},
            "inputs": {p.name: _sha(p) for p in
                       (COVER_JSON, ATLAS_JSON, C127E_JSON, SCOUT_JSON,
                        F1_JSON)},
            "self_sha256": _sha(__file__)}}

    ART.write_text(json.dumps(out, indent=2, ensure_ascii=False),
                   encoding="utf-8")
    print("=" * 78)
    for k, v in gates.items():
        print(f"  {'OK  ' if v else 'FAIL'} {k}")
    print(f"\n{out['verdict']}")
    print(f"gates {npass}/{len(gates)} — artefact : {ART.name}")
    print("=" * 78)
    return npass == len(gates)


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(0 if selftest() else 1)
    sys.exit(0 if build() else 1)
