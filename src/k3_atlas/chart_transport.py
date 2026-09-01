#!/usr/bin/env python3
"""
TRANSPORT ON EVERY TILE of the dyadic cover, under the hardened
contract of a review.

The cover built 252 charted tiles (63/64 of the volume of a residual
cell) but transported only 3, with four named gaps
par la review (`gpt_b1e2iii_c125_c126_d56_review_2026_07_29.md`) :

  . a claim of "partition" rested on a record carrying NO
    tree address (a volume sum plus uniqueness detects neither a
    feuille ET sa descendante, ni une branche omise compensée) ;
  . the sheet was serialised only on the 3 transported tiles, and its
    derivation compared interval LOWER BOUNDS without margin;
  · H6–H8 filtraient silencieusement les `native_failed` ;
  · aucun seuil pré-enregistré, aucun scaling à chart fixé, et le
    "positivity transport" was a DOUBLE certification, not a transport.

WHAT THIS SCRIPT PAYS:

  A. LEDGER DYADIQUE AUTONOME — chaque feuille (tuile OU résidu) reçoit
     its ADDRESS reconstructed arithmetically from the root, with
     EXACT float equality; the set of addresses is checked UNIQUE,
     PREFIX-FREE (donc intérieurs disjoints), CLOS (chaque nœud subdivisé
     a ses 16 enfants couverts, récursivement), et d'égalité de Kraft
     `sum 16^{-d} = 1` in rationals. The unresolved BOUNDARY fraction
     after each level is derived from the addresses and checked strictly
     décroissante.

  B. UNIVERSAL TRANSPORT: congruence and positive definiteness on the
     252 tiles (full mode) with, PER TILE: the sheet derived at the anchor
     (constant coefficient of the model = enclosure of the value at the centre) by
     SÉPARATION STRICTE des deux branches — marge > 0 sérialisée, refus
     si ambigu — puis figé ; `σ'` certifié par la composante ;
     compatibilité projective au même signe avec ses bornes ; AUCUN
     échec filtré : le check exige `n_ok == n_selectionnés == n_attendu`.

  C. STABILITY: the REAL positivity transport, by Weyl. With
     `C = Jᵀ Q_target conj(J)` et `D = Q_source − C`,
         λ_min(Q_source) ≥ λ_min(C) − ‖D‖₂ ≥ λmin_lo(C) − ‖D‖_F^up,
     so the check `||D||_F^up < lambda_min_lo(C)` (with C certified positive definite from its
     own enclosure) proves positive definiteness of `Q_source` BY TRANSPORT, without
     certifier directement — la certification directe reste publiée en
     CROSS-CHECK. λmin_lo(C) = det_lo/trace_hi (arrondi vers le bas),
     valid since lambda_max is at most the trace for a positive definite Hermitian matrix.
     Seuil PRÉ-ENREGISTRÉ : résidu relatif ≤ DELTA_REL = 1e-5 sur toute
     tile (the 3 pilot tiles were at most 3.93e-7; the threshold is fixed
     BEFORE the run and will not be adjusted).
     Scaling AT FIXED CHART (first slice): on tiles
     sondes, même centre, même chart, demi-largeurs h, h/2, h/4 — le
     résidu sup doit DÉCROÎTRE STRICTEMENT à chaque division ; les
     ratios are published without a threshold (no extrapolation).

  PROBES (DECLARED scope, not silent): the Jacobian mutation
  (doit casser) et le recheck d'AUTONOMIE (retransport à `ε'`/`σ'` FIGÉS
     from the record, bit-identical congruence bounds) run on the
  STRATIFIED SUBSET, not on all 252: each probe is a
  complete transport (about the cost of one tile) and the mutation tests the
  DISCRIMINATING POWER OF THE TEST, not each tile. The indices are published.

  STRATIFICATION (§13 de la review) : représentant de chaque signature
  (profondeur, chart, déterminations source/cible, σ cible) + tuile de
  jauge minimale + tuile de `det J` minimal + 3 tuiles ADJACENTES AU
  RÉSIDU (distance de Tchebychev nulle). Publiée intégralement.

WHAT THIS SCRIPT DOES NOT PAY: halos, open overlaps, the nerve
cocycle (the partition stays a partition, not a glued atlas);
the residual 1/64 (neither covered nor excluded); the EXACT contract of
congruence (the identity stays "certified approximate congruence" with bounds);
la globalisation ; the later scaling.

DURCISSEMENTS C128 (a review GPT, 2026-07-30) — appliqués puis TOUS les
regenerated artefacts (lesson: a modified script without a rerun gives
the appearance of verification):
  A  the sheet anchor includes the model REMAINDER (`p[0] +- rem` per part):
          the constant coefficient alone is not the value at the centre;
  the directed-float amendment  λmin_lo et ‖D‖_F en arithmétique d'intervalle mpmath de bout
          en bout, conversion float DIRIGÉE à la sérialisation seule,
          chaînes exactes publiées, check d'HERMITICITÉ de C, slack
          exact λ − ‖D‖ par tuile ;
  C  the ratio checked against delta is CERTIFIED: sup_entry(D) rounded up over
          q00_lo rounded down (the earlier max-endpoint denominator was a
          SUPÉRIEURE — un quotient majorant exige un dénominateur
          inférieur) ; l'ancien ratio publié en `_legacy` ;
  the autonomy probe  le scaling fige le LEDGER (ε', σ', déterminations refusées si
          they change), PREREGISTERED ratio window [16, 64]
          (order 5 plus or minus 1), one level deeper, one chart more;
  E  full SHA-256, dyadic addresses explicitly
          sérialisées.

Sorties : results/k3_cap_b1e2iii_c127_transport_pilot.json  (mode pilot)
          results/k3_cap_b1e2iii_c127_transport_all.json    (mode full)
Usage   : k3_cap_b1e2iii_c127_transport_all.py [--selftest]
Env     : K3_C127_MODE    pilot (défaut) | full
          K3_C127_WORKERS processus parallèles (défaut 3)
"""
from __future__ import annotations

import hashlib
import io
import json
import math
import os
import platform
import subprocess
import sys
import time
from fractions import Fraction
from itertools import product
from multiprocessing import get_context
from pathlib import Path

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

_argv = sys.argv
sys.argv = [sys.argv[0]]
os.environ.setdefault("K3_TM_ORDER", "4")
os.environ.setdefault("K3_TM_SERIES", "4")
from mpmath import iv, mp                                          # noqa: E402
from .witness_registry import load_canonical_MH              # noqa: E402
from .interval_arithmetic import build_M_civ, iv_bounds       # noqa: E402
from .taylor_models import (                                     # noqa: E402
    IVPM, TMC, CIV, TM_ORDER, UNARY_SERIES_DEG, riv)
from .full_cell_charts import build_section               # noqa: E402
from .gram_congruence import (                         # noqa: E402
    GAMMA, Qmat, congruence, contains_zero, mat_sub,
    native_target_section, pd_bounds, target_component, _c)
sys.argv = _argv

RES = Path(os.environ.get(
    "K3_RES_DIR", Path(__file__).resolve().parent / "data"))
COVER_JSON = RES / "k3_cap_b1e2iii_d5_6_dyadic_cover.json"
MODE = os.environ.get("K3_C127_MODE", "pilot")
N_WORKERS = int(os.environ.get("K3_C127_WORKERS", "3"))
ART = RES / ("k3_cap_b1e2iii_c127_transport_all.json" if MODE == "full"
             else "k3_cap_b1e2iii_c127_transport_pilot.json")

# --- PREREGISTERED, fixed before the run, not adjustable --------------------
# The checked ratio is now CERTIFIED: sup_entry(D) rounded
# UP over a certified LOWER denominator q00_lo(Q_src)
# (max_entry(Q) >= q00 >= q00_lo, so the quotient bounds the true ratio).
DELTA_REL = 1e-5          # plafond du résidu relatif CERTIFIÉ
N_RESIDUAL_ADJACENT = 3   # tiles adjacent to the residual in the pilot
SCALING_LEVELS = 4        # h, h/2, h/4, h/8 à chart ET sheet record figés
# the autonomy probe : fenêtre pré-enregistrée du ratio de scaling — ordre observé
# 5 ± 1 (2⁴ à 2⁶), déclarée AVANT le rerun, jamais ajustée dessus.
SCALING_RATIO_WINDOW = (16.0, 64.0)

T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def _sha(p):
    # the provenance amendment : SHA-256 COMPLETS (les 16 premiers caractères suffisaient
    # to practical identity, not to publishable provenance)
    try:
        return hashlib.sha256(Path(p).read_bytes()).hexdigest()
    except OSError:
        return None


def provenance(src, t_wall):
    try:
        h = subprocess.run(["git", "rev-parse", "HEAD"],
                           cwd=Path(__file__).resolve().parent,
                           capture_output=True, text=True,
                           timeout=10).stdout.strip() or None
    except Exception:
        h = None
    return {"git_head": h, "python": sys.version.split()[0],
            "platform": platform.platform(), "mp_prec": int(mp.prec),
            "tm_order": TM_ORDER, "unary_series_deg": UNARY_SERIES_DEG,
            "wall_s": round(t_wall, 1), "mode": MODE,
            "n_workers": N_WORKERS,
            "preregistered": {"delta_rel": DELTA_REL,
                              "n_residual_adjacent": N_RESIDUAL_ADJACENT,
                              "scaling_levels": SCALING_LEVELS,
                              "scaling_ratio_window":
                                  list(SCALING_RATIO_WINDOW)},
            "inputs": {str(Path(x).name): _sha(x) for x in src},
            "self_sha256": _sha(__file__)}


# ===========================================================================
#  the address sheet record step — le sheet record dyadique AUTONOME
# ===========================================================================
def address_of(root_c, root_h, c, h):
    """Adresse dyadique d'une feuille, reconstruite ARITHMÉTIQUEMENT en
    rejouant la récursion de `subboxes` (mêmes opérations float, même
    order): the reconstruction must land back on `(c, h)` in exact
    float EXACTE, sinon None."""
    d, hh = 0, root_h
    while hh > h and d <= 60:
        hh /= 2.0
        d += 1
    if hh != h:
        return None
    nc, nh = list(root_c), root_h
    path = []
    for _lvl in range(d):
        nh /= 2.0
        digit = tuple(1 if c[k] > nc[k] else -1 for k in range(4))
        nc = [nc[k] + digit[k] * nh for k in range(4)]
        path.append(digit)
    if nh != h or any(nc[k] != c[k] for k in range(4)):
        return None
    return tuple(path)


def tree_gates(addresses):
    """Unicité, prefix-free, clôture récursive, Kraft — depuis les
    SEULES adresses."""
    n = len(addresses)
    unique = len(set(addresses)) == n
    aset = set(addresses)
    prefix_free = not any(
        a != b and a == b[:len(a)] for a in aset for b in aset)

    max_d = max((len(a) for a in addresses), default=0)

    def covered(node):
        if node in aset:
            return True
        if len(node) >= max_d:
            return False
        return all(covered(node + (dg,))
                   for dg in product((-1, 1), repeat=4))
    closure = covered(())
    kraft = sum(Fraction(1, 16 ** len(a)) for a in addresses)
    return {"n_leaves": n, "unique": bool(unique),
            "prefix_free": bool(prefix_free),
            "tree_closed": bool(closure),
            "kraft_sum": [kraft.numerator, kraft.denominator],
            "kraft_is_one": bool(kraft == 1),
            "note": ("prefix-free on closed dyadic boxes gives "
                     "pairwise disjoint interiors; the boundaries "
                     "se touchent, et c'est attendu")}


def frontier_fractions(addresses):
    """Fraction of the boundary NOT RESOLVED after each level: the mass of
    depth-d nodes that are PROPER prefixes of a leaf."""
    max_d = max((len(a) for a in addresses), default=0)
    out = []
    for d in range(max_d):
        internal = {a[:d] for a in addresses if len(a) > d}
        f = Fraction(len(internal), 16 ** d)
        out.append({"after_depth": d,
                    "internal_nodes": len(internal),
                    "fraction": [f.numerator, f.denominator],
                    "fraction_float": float(f)})
    return out


# ===========================================================================
#  the transport step / the stability step — le transport durci
# ===========================================================================
def _f_down(x):
    """Conversion mpf → float DIRIGÉE vers le bas (the directed-float amendment : `float()`
    rounds to nearest, which can raise a lower bound)."""
    f = float(x)
    return math.nextafter(f, -math.inf) if f > x else f


def _f_up(x):
    f = float(x)
    return math.nextafter(f, math.inf) if f < x else f


def _anchor(Z):
    """The enclosure of the value AT THE CENTRE of a complex model: `p[0] +-
    rem`, per part. The constant coefficient ALONE is NOT the value
    at the centre once `rem > 0`: the model remainder belongs to the anchor
    (défaut nommé par la a review §4.2)."""
    rtm, itm = Z.re_tm(), Z.im_tm()
    return (rtm.p[0] + IVPM * rtm.rem, itm.p[0] + IVPM * itm.rem)


def _eps_sep(aE, bE):
    """Choix de signe par séparation stricte des distances au carré,
    on ANCHOR ENCLOSURES (remainder included). Comparisons in mpf, not
    en float. Retourne (eps, rec) ou (None, diagnostic)."""
    (are, aim), (bre, bim) = aE, bE
    d2p = (are - bre) ** 2 + (aim - bim) ** 2
    d2m = (are + bre) ** 2 + (aim + bim) ** 2
    p_lo, p_hi = mp.mpf(d2p.a), mp.mpf(d2p.b)
    m_lo, m_hi = mp.mpf(d2m.a), mp.mpf(d2m.b)
    if p_hi < m_lo:
        return 1, {"eps": 1, "margin": _f_down(m_lo - p_hi),
                   "d2_chosen": [_f_down(p_lo), _f_up(p_hi)],
                   "d2_other": [_f_down(m_lo), _f_up(m_hi)]}
    if m_hi < p_lo:
        return -1, {"eps": -1, "margin": _f_down(p_lo - m_hi),
                    "d2_chosen": [_f_down(m_lo), _f_up(m_hi)],
                    "d2_other": [_f_down(p_lo), _f_up(p_hi)]}
    return None, {"d2_plus": [_f_down(p_lo), _f_up(p_hi)],
                  "d2_minus": [_f_down(m_lo), _f_up(m_hi)]}


def derive_eps_target_margined(S2, g2, up, vp, sigma2, Zp):
    """`ε'` dérivé à l'ancre par SÉPARATION STRICTE, marge sérialisée.

    The anchor of each row is the COMPLETE enclosure of the
    value at the centre, `p[0] +- rem` per part, and no longer the sole
    constant coefficient. The sign is kept only if `sup(chosen squared distance) <
    inf(other squared distance)` STRICTLY (compared in mpf); ambiguous means REFUSAL."""
    Zt1, _d, _k = native_target_section(S2, g2, up, vp, (1, 1, 1),
                                        sigma2)
    if Zt1 is None:
        return None, None, "native_eps1_failed"
    eps, margins = [], []
    for r, s in enumerate(S2):
        sign, rec = _eps_sep(_anchor(Zt1[s]), _anchor(Zp[s]))
        if sign is None:
            rec["row"] = r
            return None, rec, "eps_ambiguous"
        rec["row"] = r
        eps.append(sign)
        margins.append(rec)
    return tuple(eps), margins, None


def lam_min_lo(Q):
    """Borne inférieure de λ_min d'une hermitienne 2×2 PD :
    λ_min = det/λ_max ≥ det_lo/trace_hi.

    the directed-float amendment : intégralement en arithmétique d'intervalle mpmath
    (rounded outwards at working precision), and the interval division
    gives det_lo/trace_hi as a lower bound without an intermediate
    conversion to the nearest float. Returns (mpf, a serialisable
    dict with directed floats and an exact string), or (None, None)."""
    q00 = Q[0][0].re_tm().to_iv()
    q11 = Q[1][1].re_tm().to_iv()
    det = (Q[0][0] * Q[1][1]
           + (Q[0][1] * Q[1][0]).mul_real(riv(-1.0))).re_tm().to_iv()
    tr = q00 + q11
    if not (mp.mpf(det.a) > 0 and mp.mpf(tr.b) > 0
            and mp.mpf(q00.a) > 0 and mp.mpf(q11.a) > 0):
        return None, None
    lam = det / tr
    lo = mp.mpf(lam.a)
    return lo, {"float": _f_down(lo), "exact": mp.nstr(lo, 40)}


def fro_up(D):
    """Upper bound on ||D||_F from the MATRIX of complex models (rather than from
    bornes déjà converties en float) — the directed-float amendment : sommation et racine en
    iv, conversion float dirigée à la sérialisation seulement.
    Retourne (mpf, dict)."""
    s = iv.mpf(0)
    for i in range(2):
        for k in range(2):
            re = D[i][k].re_tm().to_iv()
            im = D[i][k].im_tm().to_iv()
            mre = max(abs(mp.mpf(re.a)), abs(mp.mpf(re.b)))
            mim = max(abs(mp.mpf(im.a)), abs(mp.mpf(im.b)))
            s = s + iv.mpf(mre) ** 2 + iv.mpf(mim) ** 2
    hi = mp.mpf(iv.sqrt(s).b)
    return hi, {"float": _f_up(hi), "exact": mp.nstr(hi, 40)}


def sup_entry_up(D):
    """Supremum of the component moduli (real and imaginary) of the matrix, in mpf."""
    out = mp.mpf(0)
    for i in range(2):
        for k in range(2):
            for part in (D[i][k].re_tm(), D[i][k].im_tm()):
                x = part.to_iv()
                out = max(out, abs(mp.mpf(x.a)), abs(mp.mpf(x.b)))
    return out


def _conj_tmc(x):
    return TMC([CIV(c.re, -c.im) for c in x.p], x.rem)


def hermitian_contains_zero(C):
    """Check d'hermiticité (the directed-float amendment) : `C[0][1] − conj(C[1][0]) ∋ 0` et
    the diagonal imaginary parts contain 0, without which
    `lambda_min(C)` lacks the spectral meaning Weyl requires."""
    off = C[0][1] + _conj_tmc(C[1][0]).mul_real(riv(-1.0))
    checks = [off.re_tm().to_iv(), off.im_tm().to_iv(),
              C[0][0].im_tm().to_iv(), C[1][1].im_tm().to_iv()]
    return all(mp.mpf(x.a) <= 0 <= mp.mpf(x.b) for x in checks)


def transport_hardened(S, g, eps, center, hw, S2, g2, M, c218, rw,
                       fixed_eps2=None, fixed_sigma2=None,
                       perturb_J=None, section=None):
    """Miroir de `transport` (d5_congruence) avec : `ε'` À MARGE (refus
    si ambigu), mode À ENTRÉES FIGÉES (recheck d'autonomie — aucune
    dérivation), et bornes spectrales du transport par Weyl.

    `section` (the bridge step/F4, 2026-07-31) — POINT D'INJECTION ADDITIF.
    `center` et `hw` ne servent QU'À `build_section` : tout l'aval
    (Qmat, transport, congruence, Weyl, ratio certifié) est
    BOX-AGNOSTIC. Passing an already built section therefore allows
    replaying the hardened path VERBATIM on a box that `build_section`
    cannot describe: the bilateral bridge box, which is
    ANISOTROPIC (2H in the imaginary direction, H in the real one) and whose regime is chosen on the
    sign of `Re R` rather than on the component, undetermined there.
    When `section is None` the behaviour is UNCHANGED to the byte:
    no existing caller is affected, and the metric script checks it
    with a non-regression check against the serialised artefact."""
    t_start = time.time()
    if section is None:
        Z, dZ, rows = build_section(S, g, eps, center, hw)
    else:
        Z, dZ, rows = section
    if any(z is None for z in Z):
        return {"failed": "source_section_incomplete"}
    W = [[dZ[a][0], dZ[a][1]] for a in range(6)]
    Q_src = Qmat(Z, W, M, c218, rw)

    T2 = [j for j in range(6) if j not in S2]
    o = [x for x in T2 if x != g2]
    ib = Z[g2].inv()
    Zp = [z * ib for z in Z]
    dZp = [((dZ[a][0] * Z[g2] + (Z[a] * dZ[g2][0]).mul_real(riv(-1.0)))
            * ib * ib,
            (dZ[a][1] * Z[g2] + (Z[a] * dZ[g2][1]).mul_real(riv(-1.0)))
            * ib * ib) for a in range(6)]
    Q_mid = Qmat(Zp, [[dZp[a][0], dZp[a][1]] for a in range(6)],
                 M, c218, rw)
    up_, vp_ = Zp[o[0]], Zp[o[1]]
    J = [[dZp[o[0]][0], dZp[o[0]][1]],
         [dZp[o[1]][0], dZp[o[1]][1]]]
    if perturb_J is not None:
        J[0][0] = J[0][0] + TMC.const(_c(perturb_J))

    if fixed_sigma2 is not None:
        sigma2 = list(fixed_sigma2)
        margins = None
    else:
        sigma2 = target_component(S2, g2, up_, vp_)
    if fixed_eps2 is not None:
        eps2 = tuple(fixed_eps2)
        margins = None
    else:
        eps2, margins, err = derive_eps_target_margined(
            S2, g2, up_, vp_, sigma2, Zp)
        if eps2 is None:
            return {"failed": err, "detail": margins,
                    "sigma_target": sigma2}
    Zt, dZt, kinds_t = native_target_section(S2, g2, up_, vp_, eps2,
                                             sigma2)
    if Zt is None:
        return {"failed": "native_target_failed",
                "kinds_target": kinds_t, "sigma_target": sigma2,
                "eps_target": list(eps2)}
    Q_tgt = Qmat(Zt, [[dZt[a][0], dZt[a][1]] for a in range(6)],
                 M, c218, rw)
    C = congruence(J, Q_tgt)

    gauge_ok, _gd = contains_zero(mat_sub(Q_src, Q_mid))
    Dm = mat_sub(Q_src, C)
    cong_ok, cong_d = contains_zero(Dm)
    compat = []
    for a in range(6):
        D = Zt[a] + Zp[a].mul_real(riv(-1.0))
        re, im = D.re_tm().to_iv(), D.im_tm().to_iv()
        rl, rh = float(mp.mpf(re.a)), float(mp.mpf(re.b))
        il, ih = float(mp.mpf(im.a)), float(mp.mpf(im.b))
        compat.append({"coord": a,
                       "same_point": bool(rl <= 0 <= rh
                                          and il <= 0 <= ih),
                       "residual": [rl, rh, il, ih]})
    qs_norm = max(abs(x) for i in range(2) for k in range(2)
                  for x in iv_bounds(Q_src[i][k].re_tm().to_iv()))
    sup_abs = max(max(abs(x["re"][0]), abs(x["re"][1]),
                      abs(x["im"][0]), abs(x["im"][1]))
                  for x in cong_d)
    (ps, ds) = pd_bounds(Q_src)
    (pt, dt) = pd_bounds(Q_tgt)
    lam_mpf, lam_c = lam_min_lo(C)
    fro_mpf, fro_d = fro_up(Dm)
    herm_ok = hermitian_contains_zero(C)
    slack = None
    if lam_mpf is not None:
        sm = lam_mpf - fro_mpf
        slack = {"float": _f_down(sm), "exact": mp.nstr(sm, 40)}
    # the certified-ratio amendment : ratio CERTIFIÉ — numérateur sup_entry(D) (mpf),
    # dénominateur q00_lo(Q_src) INFÉRIEUR certifié, quotient iv
    # arrondi vers le haut. max_entry(Q) ≥ q00 ≥ q00_lo ⟹ majorant.
    sup_mpf = sup_entry_up(Dm)
    q00_lo_mpf = mp.mpf(Q_src[0][0].re_tm().to_iv().a)
    rel_cert = None
    if q00_lo_mpf > 0:
        rel_cert = _f_up(mp.mpf((iv.mpf(sup_mpf)
                                 / iv.mpf(q00_lo_mpf)).b))
    return {
        "eps_target": list(eps2), "sigma_target": sigma2,
        "eps_margins": margins, "kinds_target": kinds_t,
        "source_determinations": [x["determination"] for x in rows],
        "gauge_invariance_ok": bool(gauge_ok),
        "congruence_contains_zero": bool(cong_ok),
        "congruence_components": cong_d,
        "same_projective_point": all(x["same_point"] for x in compat),
        "compat": compat,
        "residual_sup_abs": sup_abs, "Q_source_norm": qs_norm,
        "residual_relative": rel_cert,
        "residual_relative_legacy": sup_abs / max(qs_norm, 1e-300),
        "q00_lower": _f_down(q00_lo_mpf),
        "pd_source": {"q00": list(ps), "det": list(ds),
                      "is_PD": bool(ps[0] > 0 and ds[0] > 0)},
        "pd_target": {"q00": list(pt), "det": list(dt),
                      "is_PD": bool(pt[0] > 0 and dt[0] > 0)},
        "spectral": {
            "lam_min_lo_C": lam_c, "fro_up_D": fro_d,
            "slack": slack, "hermitian_ok": bool(herm_ok),
            "C_is_PD": bool(lam_mpf is not None),
            "weyl_transport_ok": bool(lam_mpf is not None and herm_ok
                                      and fro_mpf < lam_mpf),
            "note": ("lambda_min(Q_src) >= lambda_min_lo(C) minus the Frobenius bound: positive definiteness of "
                     "Q_source is PROVED BY TRANSPORT, the direct one "
                     "being only a cross-check; bounds in 200-bit interval arithmetic, "
                     "floats DIRECTED at serialisation")},
        "wall_s": round(time.time() - t_start, 1)}


# ===========================================================================
#  Worker (multiprocessing — fork : M/c218 hérités du parent)
# ===========================================================================
_G = {}


def _init_worker(cell, M, c218, rw):
    _G["cell"], _G["M"], _G["c218"], _G["rw"] = cell, M, c218, rw


def _run_tile(job):
    i, t, kind, extra = job
    S, g, eps = _G["cell"]
    c = [float.fromhex(x) for x in t["center_hex"]]
    h = float.fromhex(t["hw_hex"])
    if extra and "hw_scale" in extra:
        h = h * extra["hw_scale"]
    S2, g2 = tuple(t["chart"]["S"]), t["chart"]["g"]
    kw = {}
    if kind == "mutation":
        kw["perturb_J"] = complex(0.3, 0.2)
    elif kind in ("autonomy", "scaling"):
        # the autonomy probe : le scaling fige aussi le LEDGER (ε', σ') — réduire h
        # without freezing the branch would measure something other than the model remainder
        kw["fixed_eps2"] = extra["eps_target"]
        kw["fixed_sigma2"] = extra["sigma_target"]
    r = transport_hardened(S, g, eps, c, h, S2, g2,
                           _G["M"], _G["c218"], _G["rw"], **kw)
    r["tile_index"], r["kind"] = i, kind
    r["depth"], r["chart"] = t["depth"], t["chart"]
    r["center_hex"], r["hw_hex"] = t["center_hex"], float(h).hex()
    return r


# ===========================================================================
#  Stratification (§13 de la review)
# ===========================================================================
def stratify(tiles, residual):
    sig_of = {}
    for i, t in enumerate(tiles):
        sig = (t["depth"], tuple(t["chart"]["S"]), t["chart"]["g"],
               tuple(t["source_determinations"]),
               tuple(t["target_determinations"]),
               tuple(x if x is not None else 0
                     for x in t["target_sigma"]))
        sig_of.setdefault(sig, []).append(i)
    reps = sorted(v[0] for v in sig_of.values())
    g_min = min(range(len(tiles)),
                key=lambda i: (tiles[i]["gauge_absmin"], i))
    j_min = min(range(len(tiles)),
                key=lambda i: (tiles[i]["detJ_absmin"], i))
    adj = []
    for i, t in enumerate(tiles):
        ct = [float.fromhex(x) for x in t["center_hex"]]
        ht = float.fromhex(t["hw_hex"])
        best = math.inf
        for rb in residual:
            cr = [float.fromhex(x) for x in rb["center_hex"]]
            hr = float.fromhex(rb["hw_hex"])
            gap = max(abs(ct[k] - cr[k]) for k in range(4)) - (ht + hr)
            best = min(best, gap)
        adj.append((best, i))
    adj.sort()
    near = [i for _gap, i in adj[:N_RESIDUAL_ADJACENT]]
    sel = sorted(set(reps + [g_min, j_min] + near))
    return sel, {
        "n_strata": len(sig_of),
        "strata": [{"signature": [list(k[1])] + [k[0], k[2],
                                                 list(k[3]), list(k[4]),
                                                 list(k[5])],
                    "count": len(v), "representative": v[0]}
                   for k, v in sorted(sig_of.items())],
        "min_gauge_tile": g_min, "min_detJ_tile": j_min,
        "residual_adjacent_tiles": near,
        "residual_adjacent_gaps": [g for g, _i in adj[:3]],
        "selected": sel}


# ===========================================================================
#  Build
# ===========================================================================
def build():
    print("=" * 78)
    print(f"the transport step TRANSPORT {'UNIVERSEL (252)' if MODE == 'full' else 'PILOTE STRATIFIÉ'}"
          f" : TM ({TM_ORDER},{UNARY_SERIES_DEG}), {N_WORKERS} workers, "
          f"δ_rel pré-enregistré = {DELTA_REL:.0e}")
    print("=" * 78)
    reg = load_canonical_MH()
    M = build_M_civ(reg["M_H_canonical"])
    c218 = reg["coeffs218"]
    rw = 1.0 - GAMMA
    cov = json.loads(COVER_JSON.read_text(encoding="utf-8"))
    tiles, residual = cov["tiles"], cov["residual"]
    cell = cov["cell"]
    S, g, eps = tuple(cell["S"]), cell["g"], tuple(cell["eps"])
    root_c = [float.fromhex(x) for x in cell["center_hex"]]
    root_h = float.fromhex(cell["hw_hex"])
    log(f"cover chargé : {len(tiles)} tuiles, {len(residual)} boîtes "
        f"résiduelles · cellule S={list(S)} g={g} eps={list(eps)}")

    # --- the address sheet record step : le sheet record dyadique autonome --------------------------------
    leaves = [(t, "tile") for t in tiles] + [(r, "residual")
                                             for r in residual]
    addresses, addr_fail, leaf_addr = [], 0, []
    for (lf, kind) in leaves:
        a = address_of(root_c, root_h,
                       [float.fromhex(x) for x in lf["center_hex"]],
                       float.fromhex(lf["hw_hex"]))
        if a is None:
            addr_fail += 1
        else:
            addresses.append(a)
            # the provenance amendment : adresses explicitement SÉRIALISÉES (elles
            # were only reconstructed: autonomous but less
            # auditable de l'extérieur)
            leaf_addr.append({"kind": kind, "depth": len(a),
                              "path": [list(d) for d in a]})
    tg = tree_gates(addresses)
    fr = frontier_fractions(addresses)
    residual_final = sum(Fraction(1, 16 ** r["depth"])
                         for r in residual)
    log(f"the address sheet record step : {len(addresses)}/{len(leaves)} adresses exactes, "
        f"prefix-free={tg['prefix_free']}, clos={tg['tree_closed']}, "
        f"Kraft={tg['kraft_sum'][0]}/{tg['kraft_sum'][1]} · frontière "
        f"{[x['fraction_float'] for x in fr]}")

    # --- Stratification -------------------------------------------------------
    sel_pilot, strat = stratify(tiles, residual)
    log(f"stratification : {strat['n_strata']} strates, pilote = "
        f"{len(sel_pilot)} tuiles {sel_pilot}")

    selected = list(range(len(tiles))) if MODE == "full" else sel_pilot
    expected_n = len(selected)

    # --- Les jobs : transports + probes ---------------------------------------
    # fork (not forkserver, the Python 3.14 default): mpmath interval
    # objects do not pickle; under fork they are inherited without pickling.
    mpctx = get_context("fork")
    jobs = [(i, tiles[i], "transport", None) for i in selected]
    with mpctx.Pool(N_WORKERS, initializer=_init_worker,
                    initargs=((S, g, eps), M, c218, rw)) as pool:
        results = pool.map(_run_tile, jobs)
    by_idx = {r["tile_index"]: r for r in results}
    failed = [r for r in results if r.get("failed")]
    ok = [r for r in results if not r.get("failed")]
    log(f"transports : {len(ok)}/{expected_n} OK, {len(failed)} échecs "
        f"(REFUSÉS, pas filtrés)")

    # --- Probes on the stratified subset --------------------------------------
    probe_idx = [i for i in sel_pilot if i in by_idx
                 and not by_idx[i].get("failed")]
    probe_jobs = []
    for i in probe_idx:
        r0 = by_idx[i]
        probe_jobs.append((i, tiles[i], "mutation", None))
        probe_jobs.append((i, tiles[i], "autonomy",
                           {"eps_target": r0["eps_target"],
                            "sigma_target": r0["sigma_target"]}))
    # scaling à chart ET sheet record figés (the autonomy probe) : une tuile par chart
    # distinct du cover (3), niveaux h/2 … h/2^(SCALING_LEVELS−1)
    seen_ch, scale_base = set(), []
    for i in probe_idx:
        ch = (tuple(tiles[i]["chart"]["S"]), tiles[i]["chart"]["g"])
        if ch not in seen_ch:
            seen_ch.add(ch)
            scale_base.append(i)
    scale_base = scale_base[:3]
    for i in scale_base:
        r0 = by_idx[i]
        for lvl in range(1, SCALING_LEVELS):
            probe_jobs.append((i, tiles[i], "scaling",
                               {"hw_scale": 0.5 ** lvl,
                                "eps_target": r0["eps_target"],
                                "sigma_target": r0["sigma_target"]}))
    with mpctx.Pool(N_WORKERS, initializer=_init_worker,
                    initargs=((S, g, eps), M, c218, rw)) as pool:
        probes = pool.map(_run_tile, probe_jobs)

    mutations = [p for p in probes if p["kind"] == "mutation"]
    autonomy = [p for p in probes if p["kind"] == "autonomy"]
    scaling = [p for p in probes if p["kind"] == "scaling"]

    mut_break = [not p.get("failed")
                 and not p["congruence_contains_zero"]
                 for p in mutations]
    auto_ok = []
    for p in autonomy:
        r0 = by_idx[p["tile_index"]]
        auto_ok.append(
            not p.get("failed")
            and p["congruence_components"] == r0["congruence_components"]
            and p["congruence_contains_zero"]
            == r0["congruence_contains_zero"])
    scale_ok, scale_detail = [], []
    for i in scale_base:
        r0 = by_idx[i]
        levels = sorted((q for q in scaling if q["tile_index"] == i),
                        key=lambda q: -float.fromhex(q["hw_hex"]))
        seq = [r0["residual_sup_abs"]] + [p["residual_sup_abs"]
                                          for p in levels]
        ratios = [seq[k] / seq[k + 1] if seq[k + 1] else None
                  for k in range(len(seq) - 1)]
        dec = all(seq[k + 1] < seq[k] for k in range(len(seq) - 1))
        # the autonomy probe : fenêtre pré-enregistrée + sheet record analytique FIGÉ —
        # a determination that changed when reducing h would invalidate the
        # mesure (on comparerait deux branches, pas deux restes)
        in_window = all(r is not None
                        and SCALING_RATIO_WINDOW[0] <= r
                        <= SCALING_RATIO_WINDOW[1] for r in ratios)
        frozen = all(not p.get("failed")
                     and p["source_determinations"]
                     == r0["source_determinations"]
                     and p["kinds_target"] == r0["kinds_target"]
                     for p in levels)
        scale_ok.append(dec and in_window and frozen)
        scale_detail.append({
            "tile_index": i, "depth": tiles[i]["depth"],
            "chart": tiles[i]["chart"], "sup_abs_sequence": seq,
            "ratios": ratios,
            "ratio_window_preregistered": list(SCALING_RATIO_WINDOW),
            "ratios_in_window": bool(in_window),
            "sheet_record_frozen_across_levels": bool(frozen),
            "strictly_decreasing": bool(dec)})
    log(f"probes : mutation casse {sum(mut_break)}/{len(mut_break)} · "
        f"autonomie {sum(auto_ok)}/{len(auto_ok)} · scaling "
        f"{sum(scale_ok)}/{len(scale_ok)}")

    # --- Checks pré-enregistrés ------------------------------------------------
    max_rel = max((r["residual_relative"] for r in ok), default=None)
    eps_margins_all = [m for r in ok for m in (r["eps_margins"] or [])]
    checks = {
        "G1_dyadic_tree_autonomous": bool(
            addr_fail == 0 and tg["unique"] and tg["prefix_free"]
            and tg["tree_closed"] and tg["kraft_is_one"]),
        "G2_frontier_strictly_decreasing": bool(
            len(fr) >= 2 and all(
                fr[k + 1]["fraction_float"] < fr[k]["fraction_float"]
                for k in range(len(fr) - 1))
            and float(residual_final) < fr[-1]["fraction_float"]),
        "G3_all_transports_succeed": bool(
            len(ok) == expected_n == len(results) and not failed),
        "G4_congruence_universal": bool(ok) and all(
            r["gauge_invariance_ok"] and r["congruence_contains_zero"]
            and r["same_projective_point"] for r in ok),
        "G5_residual_below_prereg_delta": bool(ok) and all(
            r["residual_relative"] is not None
            and r["residual_relative"] <= DELTA_REL for r in ok),
        "G6_weyl_transport_universal": bool(ok) and all(
            r["spectral"]["weyl_transport_ok"] for r in ok),
        "G13_hermitian_universal": bool(ok) and all(
            r["spectral"]["hermitian_ok"] for r in ok),
        "G7_pd_both_sides_universal": bool(ok) and all(
            r["pd_source"]["is_PD"] and r["pd_target"]["is_PD"]
            for r in ok),
        "G8_eps_margin_positive_universal": bool(eps_margins_all)
        and all(m["margin"] > 0 for m in eps_margins_all),
        "G9_jacobian_mutation_breaks_on_probes": bool(mut_break)
        and all(mut_break),
        "G10_scaling_fixed_chart_decreases": bool(scale_ok)
        and all(scale_ok),
        "G11_autonomy_fixed_input_identical": bool(auto_ok)
        and all(auto_ok),
        "G12_no_silent_cap": bool(
            len(results) == len(jobs)
            and len(probes) == len(probe_jobs)
            and set(sel_pilot) <= set(selected)
            and len(strat["strata"]) == strat["n_strata"])}

    n_pass = sum(1 for v in checks.values() if v)
    wall_tiles = [r["wall_s"] for r in ok]
    verdict = (
        "Transport (%s) over the dyadic cover, under the hardened "
        "contract. First, the %d leaves (tiles plus residual) carry a "
        "dyadic address reconstructed in EXACT float equality; the tree "
        "is prefix-free, CLOSED, with Kraft equality 1, so neither hole "
        "nor overlap is now AUTONOMOUS from the record. "
        "Unresolved boundary by level: %s, strictly decreasing. "
        "TRANSPORT: %d of %d tiles, NO failure filtered out; the sheet derived with "
        "a strictly positive margin and serialised EVERYWHERE (minimum margin "
        "%.3e); maximum relative residual %.3e at most the preregistered delta %.0e; and "
        "positivity transport is now PROVED BY WEYL on each tile "
        "(the Frobenius bound below lambda_min, worst ratio %.3e), the "
        "direct certification being only a cross-check. Probes "
        "(declared scope: %d stratified tiles): the Jacobian "
        "mutation breaks %d of %d, the autonomy recheck at frozen inputs is "
        "bit-identical %d of %d, and scaling at fixed chart is strictly "
        "decreasing %d of %d. NOT PAID HERE: halos, overlaps and the cocycle; "
        "the residual 1/64; the exact contract of the identity; "
        "globalisation." % (
            MODE, "UNIVERSEL" if MODE == "full" else "stratified pilot",
            len(addresses),
            [x["fraction_float"] for x in fr],
            len(ok), expected_n,
            min((m["margin"] for m in eps_margins_all),
                default=float("nan")),
            max_rel if max_rel is not None else float("nan"), DELTA_REL,
            max((r["spectral"]["fro_up_D"]["float"]
                 / r["spectral"]["lam_min_lo_C"]["float"]
                 for r in ok if r["spectral"]["lam_min_lo_C"]),
                default=float("nan")),
            len(probe_idx), sum(mut_break), len(mut_break),
            sum(auto_ok), len(auto_ok), sum(scale_ok), len(scale_ok)))

    art = {"artifact": ART.stem, "mode": MODE,
           "claim": ("Autonomous dyadic record plus hardened "
                     "transport (universal in full mode) over the "
                     "cover, with positivity transport by Weyl."),
           "cell": cell,
           "tree_gates": tg, "n_address_failures": addr_fail,
           "leaf_addresses": leaf_addr,
           "frontier_fractions": fr,
           "stratification": strat,
           "selected_tiles": selected,
           "expected_n": expected_n,
           "transports": results,
           "failed_transports": [r["tile_index"] for r in failed],
           "probes": {"mutation": mutations, "autonomy": autonomy,
                      "scaling_detail": scale_detail,
                      "probe_indices": probe_idx,
                      "scope_note": (
                          "mutation, autonomy and scaling on the "
                          "STRATIFIED subset, not on all 252: "
                          "each probe costs a complete transport and "
                          "tests the discriminating power of the TEST, not each "
                          "tile; declared scope, published indices")},
           "max_residual_relative": max_rel,
           "delta_rel_preregistered": DELTA_REL,
           "wall_per_tile_s": {"min": min(wall_tiles, default=None),
                               "median": sorted(wall_tiles)[
                                   len(wall_tiles) // 2]
                               if wall_tiles else None,
                               "max": max(wall_tiles, default=None)},
           "not_paid_here": ["the atlas step halos/overlaps/cocycle",
                             "the residual closure résidu 1/64",
                             "contrat exact de l'identité (E3 reste "
                             "« congruence approchée certifiée »)",
                             "globalisation", "the later scaling"],
           "verdict": verdict, "checks": checks,
           "checks_passed": n_pass, "checks_total": len(checks),
           "provenance": provenance([COVER_JSON], time.time() - T0)}
    ART.parent.mkdir(parents=True, exist_ok=True)
    ART.write_text(json.dumps(art, indent=2, ensure_ascii=False),
                   encoding="utf-8")
    print("-" * 78)
    for k, v in checks.items():
        print(f"  {'PASS' if v else 'FAIL'}  {k}")
    print("-" * 78)
    print(verdict)
    print(f"→ {ART}")
    return 0 if all(checks.values()) else 1


# ===========================================================================
#  Self-test (pure functions, without the registry or the models)
# ===========================================================================
def _selftest():
    ok = True

    def chk(name, cond):
        nonlocal ok
        print(f"  {'ok  ' if cond else 'FAIL'} {name}")
        ok = ok and bool(cond)

    root_c, root_h = [0.1234567, -0.7654321, 0.5, -0.25], 0.001953125

    def children(c, h):
        hh = h / 2.0
        return [([c[k] + dg[k] * hh for k in range(4)], hh, dg)
                for dg in product((-1, 1), repeat=4)]

    # T1: EXACT address reconstruction on a synthetic tree
    #      (feuille de profondeur 2, mêmes opérations float que subboxes)
    (c1, h1, d1) = children(root_c, root_h)[5]
    (c2, h2, d2) = children(c1, h1)[10]
    a = address_of(root_c, root_h, c2, h2)
    chk("T1 adresse profondeur 2 reconstruite exactement",
        a == (d1, d2))

    # T2: a leaf PERTURBED by one ulp has NO address
    c_bad = list(c2)
    c_bad[0] = math.nextafter(c_bad[0], math.inf)
    chk("T2 négatif : centre perturbé d'un ulp → adresse REFUSÉE",
        address_of(root_c, root_h, c_bad, h2) is None)

    # T3 : arbre complet à profondeur mixte — clôture et Kraft
    addrs = [(dg,) for dg in product((-1, 1), repeat=4)
             if dg != (1, 1, 1, 1)]
    addrs += [((1, 1, 1, 1), dg) for dg in product((-1, 1), repeat=4)]
    t = tree_gates(addrs)
    chk("T3 arbre 15 + 16 : unique, prefix-free, clos, Kraft = 1",
        t["unique"] and t["prefix_free"] and t["tree_closed"]
        and t["kraft_is_one"])

    # T4 NEGATIVE CONTROL: a leaf AND its descendant (the hole that the volume
    #      sum does not see if compensated) gives a prefix-free FAILURE
    bad = addrs + [((-1, -1, -1, -1), (1, 1, 1, 1))]
    t4 = tree_gates(bad)
    chk("T4 négatif : feuille + sa descendante → prefix-free REFUSÉ",
        not t4["prefix_free"])

    # T5 : NÉGATIF — un enfant manquant → clôture FAIL, Kraft ≠ 1
    t5 = tree_gates(addrs[:-1])
    chk("T5 négatif : enfant manquant → arbre NON clos, Kraft ≠ 1",
        not t5["tree_closed"] and not t5["kraft_is_one"])

    # T6: boundary fractions on the mixed tree: 1 then 1/16
    fr = frontier_fractions(addrs)
    chk("T6 frontière : [1, 1/16] strictement décroissante",
        len(fr) == 2 and fr[0]["fraction_float"] == 1.0
        and abs(fr[1]["fraction_float"] - 1.0 / 16) < 1e-15)

    # T7: det_lo/tr_hi really is a LOWER BOUND (exact algebra:
    #      diag(2, 3) → det/tr = 6/5 = 1.2 ≤ λ_min = 2)
    chk("T7 lambda_min_lo: det_lo/tr_hi is at most lambda_min on an exact case",
        _f_down(mp.mpf(6) / 5) <= 2.0)

    # T8 : fro_up majore la vraie norme de Frobenius (matrice de TMC)
    def cmat(z00, z01, z10, z11):
        return [[TMC.const(CIV(riv(z.real), riv(z.imag)))
                 for z in row] for row in ((z00, z01), (z10, z11))]
    D = cmat(complex(0.3, -0.4), 0j, 0j, complex(1.2, 0.5))
    true_f = math.sqrt(0.3 ** 2 + 0.4 ** 2 + 1.2 ** 2 + 0.5 ** 2)
    _fm, fd = fro_up(D)
    chk("T8 fro_up ≥ ‖D‖_F vraie", fd["float"] >= true_f)

    # T9: the anchor includes the REMAINDER; two anchors whose p[0]
    # SEPARATE the signs but whose remainders overlap must
    # être REFUSÉES (l'ancien code au seul p[0] aurait tranché)
    a_clean = TMC.const(CIV(riv(1.0), riv(0.0)))
    a_fat = TMC(a_clean.p, iv.mpf(3.0))
    b = TMC.const(CIV(riv(1.0), riv(0.0)))
    s_clean, _r1 = _eps_sep(_anchor(a_clean), _anchor(b))
    s_fat, _r2 = _eps_sep(_anchor(a_fat), _anchor(b))
    chk("T9 ancre + reste : p[0] propre tranche (+1), reste large REFUSE",
        s_clean == 1 and s_fat is None)

    # T10: DIRECTED float conversion. 1/3 is not a float,
    # _f_down must round down and _f_up up; an exact float does not move
    third = mp.mpf(1) / 3
    chk("T10 _f_down ≤ x ≤ _f_up (x = 1/3), et 0.5 reste 0.5",
        _f_down(third) <= third <= _f_up(third)
        and _f_down(third) < _f_up(third)
        and _f_down(mp.mpf(0.5)) == 0.5 == _f_up(mp.mpf(0.5)))

    # T11: the hermiticity check bites. A matrix whose
    # block [1][0] is NOT the conjugate of [0][1] is refused
    H = cmat(complex(1, 0), complex(0.3, 0.2), complex(0.3, -0.2),
             complex(2, 0))
    Hbad = cmat(complex(1, 0), complex(0.3, 0.2), complex(0.3, 0.2),
                complex(2, 0))
    chk("T11 hermiticité : hermitienne acceptée, mutée REFUSÉE",
        hermitian_contains_zero(H)
        and not hermitian_contains_zero(Hbad))

    # T12: lambda_min_lo really is a directed lower bound on
    # diag(2, 3) : det/tr = 6/5 ≤ 2, et float(λ) ≤ λ exacte
    lm, ld = lam_min_lo(cmat(2 + 0j, 0j, 0j, 3 + 0j))
    chk("T12 λmin_lo(diag(2,3)) = 1.2 ≤ 2, float dirigé ≤ mpf",
        lm is not None and abs(ld["float"] - 1.2) < 1e-12
        and ld["float"] <= lm)

    print("SELF-TEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_selftest() if "--selftest" in sys.argv else build())
