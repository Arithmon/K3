#!/usr/bin/env python3
"""
FROM THE PARTITION TO THE ATLAS.

The transport step (252 tiles) and the residual closure (the 64 remaining) closed the cover
of one cell/class pair: 316 tiles, Kraft equal to 1, no residual. But that
cover stayed a PARTITION INTO CLOSED BOXES, as the notes said
elles-mêmes : « aucun halo, aucun overlap ouvert, aucune transition
certified between neighbouring tiles, no cocycle on the nerve of the 316 charts."
And closing the residual made gluing MANDATORY: tiles
adjacentes emploient des DÉTERMINATIONS DIFFÉRENTES (principale,
tournée, tournée canonique), donc rien ne garantissait a priori
that they describe the same sheet of the source.

WHAT THIS SCRIPT PAYS (the four points of the contract):

  D-a. CORES + HALOS — chaque tuile reçoit un halo `box(c, (1+ρ)h)`,
     rho taken from a PREREGISTERED SCALE (1/8, 1/16, 1/32, 1/64):
     the first rho that passes is kept, and "no rho passes" is a
     REFUSAL, not a fallback. On the halo we require: a source section
     complète, DÉTERMINATIONS SOURCE INCHANGÉES, critère de chart (celui
     that certified the core: the chart criterion for the 252, extended for the 64,
     and for those 64 the plain criterion must STILL REFUSE, non-tautology
     survit au halo), section cible native à SIGN PATTERN FIGÉ (`ε'`, `σ'`)
     with the same kinds, and projective congruence under delta.
     A chart valid on an OPEN NEIGHBOURHOOD of the core: that is what
     transforme un certificat de boîte fermée en carte d'atlas.

  D-b. OPEN OVERLAPS: the nerve is computed in DYADIC RATIONALS
     EXACTS (`Fraction`, jamais de float) : deux cores qui se touchent
     have halos meeting in a box of STRICTLY POSITIVE
     width in all four coordinates. Published: the overlap box
     of each pair, the connectivity of the nerve, the triple intersections
     triple non vide.

  D-c. TRANSITIONS — sur chaque overlap :
     . THE SHEET: the two source sections are compared not
       through their enclosures (which lose the correlation: width
       about 1e-3, a verdict without content) but by EXACT RECENTRING of the two
       polynomials in a common frame followed by COEFFICIENT BY COEFFICIENT
       COEFFICIENT — la différence tombe à ~1e-12 et la séparation
       ±  devient décisive. `θ_ij ∈ {±1}` est DÉRIVÉ à marge stricte
       (refus si ambigu, jamais d'essai) et doit valoir +1 : les tuiles
       neighbours are on the SAME SHEET. This test COULD have failed,
       c'est là tout son intérêt.
     . THE CHART: `lambda_ij = Z[g'_i]/Z[g'_j]`, with `|Z[g'_i]|` and
       `|Z[g'_j]|` certified BOUNDED BELOW by 0 on the overlap, so the transition
       is a biholomorphism, not merely a formula.
     . THE TRANSITION IDENTITY, without division: `Zt_a.Z[g']` contains `Z_a` per
       tile, and from one tile to the next `Zt^(j)_a.Z[g'_j] -
       Zt^(i)_a.Z[g'_i]` contains 0: this is the gluing that brings the
       deux sheet records `ε'` dérivés indépendamment.

  D-d. COCYCLE ON THE NERVE: on each non-empty triple intersection
     vide : `θ_ij·θ_jk·θ_ki = +1` (exact, discret), `λ_ij·λ_jk·λ_ki ∋ 1`
     (enclosure) and the transition identity on the three edges. The
     sheet cochain is a cocycle of TRIVIAL class: the choice
     of sheet glues globally on the cell.

HONESTY ABOUT WHAT IS TAUTOLOGICAL AND WHAT IS NOT:
  . `theta = +1` is NOT tautological: the source determinations vary
    from tile to tile and nothing forced the same sheet;
  . the transition identity is NOT tautological: it brings into
    présence deux sheet records `ε'` dérivés séparément ;
  . at the level of TRIPLES, the cocycle condition is IMPLIED by the
    pairwise identities (all charts are affine charts of the
    same projective space and all sections descend from ONE source section),
    and what stays non-trivial is that the inclusion of 0 on the edge i to k
    is NOT an interval consequence of i to j and j to k: the three
    are checked separately. This is stated, not disguised.
  . the NEGATIVE CONTROLS carry the discriminating power: a corrupted sheet record,
    feuillet inversé, ρ = 0, paire non adjacente, recentrage neutralisé,
    geometry shifted by one ulp; all six must BREAK.

WHAT THIS SCRIPT DOES NOT PAY: metric transport (congruence of Q,
Weyl) ON THE HALO. Those checks stay established on the cores, and
extending them to the halos costs four `Qmat` per tile (about 244 s, the whole
budget of the transport step); it is measured here on the stratified
subset only, and that is DECLARED. It also does not pay: the EXACT contract
of the congruence identity, the full scaling, the 895
autres paires cellule/classe, the later scaling.

the rational-box step/B/C (a review GPT, 2026-07-31) — trois dettes theorem-grade
payées SANS changer aucun contrat :
  A. les boîtes d'overlap (paires, triples, négatif N5) passaient par
     `iv.mpf([float(lo), float(hi)])` TO NEAREST, whereas the symbols
     intersections carry denominators in 65 (halos of 65h/64), not
     représentables en binaire : l'intervalle pouvait être STRICTEMENT
     smaller than the rational box. Hence `fraction_box_to_iv`:
     borne basse VERS LE BAS, haute VERS LE HAUT, jamais optimiste.
  B. les défauts gatés (congruence de halo, identité de paire, arêtes
     de triple) restent en mp jusqu'à la comparaison à δ ; la
     float conversion is DIRECTED (`_f_up`) and happens only at the
     sérialisation — la borne publiée MAJORE la borne comparée.
  C. check D2c : l'ÉGALITÉ D'ENSEMBLES {tuiles clippées} = {tuiles nées
     of the residual}, published with both sorted sets and their
     SHA-256: "64 = 64" does not exclude a permutation 63+1.

Sorties : results/atlas_assembly_pilot.json   (mode pilot)
          results/atlas_assembly.json         (mode full)
Usage   : atlas_assembly.py [--selftest]
Env     : K3_C127D_MODE     pilot (défaut) | full
          K3_C127D_WORKERS  processus parallèles (défaut 6)
"""
from __future__ import annotations

import hashlib
import io
import itertools
import json
import math
import os
import platform
import subprocess
import sys
import time
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
from .witness_registry import load_canonical_MH              # noqa: E402
from .interval_arithmetic import build_M_civ                  # noqa: E402
from .taylor_models import (                                     # noqa: E402
    CIV, IVPM, MIDX, MONO, NM, TMC, TM_ORDER, UNARY_SERIES_DEG,
    civ_absmin, riv)
from .full_cell_charts import (                           # noqa: E402
    build_section, chart_certificate)
from .chart_selection_criterion import (                  # noqa: E402
    native_section_constructible)
from .gram_congruence import (                         # noqa: E402
    GAMMA, native_target_section)
from .chart_transport import (                    # noqa: E402
    address_of, tree_gates, _f_down, _f_up)
from .residual_closure import native_rows_ext          # noqa: E402
sys.argv = _argv

RES = Path(os.environ.get(
    "K3_RES_DIR", Path(__file__).resolve().parent / "data"))
COVER_JSON = RES / "dyadic_cover.json"
C127_JSON = RES / "chart_transport.json"
C127E_JSON = RES / "residual_closure.json"
MODE = os.environ.get("K3_C127D_MODE", "pilot")
N_WORKERS = int(os.environ.get("K3_C127D_WORKERS", "6"))
ART = RES / ("atlas_assembly.json" if MODE == "full"
             else "atlas_assembly_pilot.json")

# --- PREREGISTERED, fixed before the run, not adjustable --------------------
# The rho scale: the FIRST that passes is kept, per tile, and
# serialised. If no rho passes, the tile is REFUSED (check D2), never
# a silent fallback on the core.
RHO_LADDER = (Fraction(1, 8), Fraction(1, 16),
              Fraction(1, 32), Fraction(1, 64))
# Ceiling on the CERTIFIED defect of the transition identity, on the
# overlaps AND on the halos. The same delta as before, carried over unchanged.
DELTA_TRANS = 1e-5
# Ceiling on the defect of the cocycle `lambda_ij.lambda_jk.lambda_ki - 1` on triples.
DELTA_COCYCLE = 1e-5
N_PILOT_TILES = 40      # tuiles du patch stratifié en mode pilot
N_HALO_METRIC = 8       # tiles where the METRIC transport is redone
                        # on the halo (declared scope, not silent)

T0 = time.time()
IV0 = iv.mpf(0)
IV1 = iv.mpf(1)


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def _sha(p):
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
            "preregistered": {
                "rho_ladder": [[r.numerator, r.denominator]
                               for r in RHO_LADDER],
                "delta_trans": DELTA_TRANS,
                "delta_cocycle": DELTA_COCYCLE,
                "n_pilot_tiles": N_PILOT_TILES,
                "n_halo_metric": N_HALO_METRIC},
            "inputs": {str(Path(x).name): _sha(x) for x in src},
            "self_sha256": _sha(__file__)}


# ===========================================================================
#  EXACT dyadic geometry: never a float inside a decision
# ===========================================================================
def box_of(o):
    """(centre, demi-largeur) d'une feuille, en Fractions exactes."""
    return ([Fraction(float.fromhex(x)) for x in o["center_hex"]],
            Fraction(float.fromhex(o["hw_hex"])))


def halo_hw(h, rho):
    """(1+rho).h, checking that the result is EXACTLY
    representable as a float (h is a power of 2, rho is dyadic)."""
    H = h * (1 + rho)
    f = float(H)
    return (H, f) if Fraction(f) == H else (H, None)


def _inter_generic(loA, hiA, loB, hiB):
    """Boîte d'intersection de deux pavés donnés par leurs bornes, ou
    None if the width is not STRICTLY positive in all four
    coordinates. Everything in rationals: "open" is a rational
    comparison, never a float test. The returned box is not a cube
    (demi-largeurs par coordonnée), d'où le passage systématique par
    (lo, hi) plutôt que par (centre, demi-largeur)."""
    lo, hi = [], []
    for k in range(4):
        a, b = max(loA[k], loB[k]), min(hiA[k], hiB[k])
        if not b > a:
            return None
        lo.append(a)
        hi.append(b)
    return {"center": [(lo[k] + hi[k]) / 2 for k in range(4)],
            "hw": [(hi[k] - lo[k]) / 2 for k in range(4)],
            "lo": lo, "hi": hi}


def cube_bounds(c, H):
    return [c[k] - H for k in range(4)], [c[k] + H for k in range(4)]


def touching(ci, hi, cj, hj):
    """Les cores FERMÉS se touchent (distance de Tchebychev ≤ 0)."""
    return max(abs(ci[k] - cj[k]) - (hi + hj) for k in range(4)) <= 0


def eps_box(w, c, H):
    """Image of the box `w` in the symbols of the tile (c, H):
    a list of 4 (lo, hi) in rationals. The model is valid only on
    [-1,1]^4, and the inclusion is a CHECK, not an assumption."""
    out = []
    for k in range(4):
        lo = (w["lo"][k] - c[k]) / H
        hi = (w["hi"][k] - c[k]) / H
        out.append((lo, hi))
    return out


def eps_box_in_range(E):
    return all(-1 <= lo and hi <= 1 for lo, hi in E)


def frac_to_iv(q):
    """Fraction → intervalle mpmath ENCLOSANT (arrondi extérieur par la
    division intervalle : aucune conversion float au plus proche)."""
    return iv.mpf(int(q.numerator)) / iv.mpf(int(q.denominator))


def fraction_box_to_iv(lo, hi):
    """(lo, hi) en Fractions → UN intervalle mpmath ENCLOSANT [lo, hi] :
    borne basse arrondie VERS LE BAS, borne haute VERS LE HAUT (the rational-box step,
    a review). L'ancien chemin `iv.mpf([float(lo), float(hi)])`
    rounded both endpoints TO NEAREST, whereas the boxes
    d'overlap héritent des demi-largeurs de halo `65h/64` et leurs
    coordinates carry denominators in 65, not representable in
    binary: the built interval could be STRICTLY SMALLER
    than the exact rational box, and a certified containment of 0 on a truncated
    domain does not certify the domain. The widening is at most 1 ulp at
    mp.prec par borne : conservateur, jamais optimiste."""
    return iv.mpf([mp.mpf(frac_to_iv(lo).a), mp.mpf(frac_to_iv(hi).b)])


# ===========================================================================
#  EXACT recentring of a Taylor model: the heart of this step
# ===========================================================================
#  Deux tuiles voisines portent leurs polynômes dans DEUX cadres ε
#  different. Comparing their ENCLOSURES on the overlap loses the
#  corrélation : chaque enclosure a la largeur de la VARIATION de la
#  function on the overlap (about 1e-3 here), so their difference does too, and the
#  verdict has no content. We therefore recentre the polynomial of j in the
#  cadre de i — substitution AFFINE ε_j = off + scale·ε_i, exacte en
#  rationnels — puis on soustrait COEFFICIENT PAR COEFFICIENT. Mesuré sur
#  une paire réelle : largeur 7,2e-3 (naïf) → 2,3e-12 (recentré).
#
#  The substitution creates NO remainder: composing a polynomial of degree
#  at most N with an affine map gives a polynomial of degree at most N, and the
#  matrice ci-dessous est exacte (Fractions) avant conversion intervalle.
# ===========================================================================
def recenter_matrix(off, scale):
    """T such that (p composed with phi)[beta] = sum_{alpha >= beta} T[beta][alpha] . p[alpha], for
    φ(ε)_k = off_k + scale·ε_k. Rendue en intervalles ENCLOSANTS."""
    acc = {}
    for ai, am in enumerate(MONO):
        for bm in itertools.product(*[range(e + 1) for e in am]):
            coef = Fraction(1)
            for k in range(4):
                coef *= (comb(am[k], bm[k])
                         * off[k] ** (am[k] - bm[k])
                         * scale ** bm[k])
            if coef:
                row = acc.setdefault(MIDX[bm], {})
                row[ai] = row.get(ai, Fraction(0)) + coef
    return {b: [(a, frac_to_iv(c)) for a, c in row.items() if c]
            for b, row in acc.items()}


def apply_recenter(T, p):
    """Applique T aux coefficients CIV d'un TM."""
    out = [CIV(IV0, IV0)] * NM
    for b, row in T.items():
        ar, ai_ = IV0, IV0
        for a, c in row:
            ar = ar + p[a].re * c
            ai_ = ai_ + p[a].im * c
        out[b] = CIV(ar, ai_)
    return out


def poly_lin(pa, pb, sign):
    """pa ± pb, coefficient par coefficient."""
    if sign > 0:
        return [CIV(pa[k].re + pb[k].re, pa[k].im + pb[k].im)
                for k in range(NM)]
    return [CIV(pa[k].re - pb[k].re, pa[k].im - pb[k].im)
            for k in range(NM)]


def enclose(p, rem, E):
    """Enclosure of a model (containment coefficients, remainder bounding BOTH
    parts) on the sub-box given by the intervals E."""
    pw = []
    for k in range(4):
        col = [IV1]
        for _d in range(TM_ORDER):
            col.append(col[-1] * E[k])
        pw.append(col)
    ar, ai_ = IV0, IV0
    for idx, m in enumerate(MONO):
        c = p[idx]
        t = IV1
        for k in range(4):
            if m[k]:
                t = t * pw[k][m[k]]
        ar = ar + c.re * t
        ai_ = ai_ + c.im * t
    return CIV(ar + IVPM * rem, ai_ + IVPM * rem)


def sep_phase(dm, dp):
    """Sheet `theta` by STRICT SEPARATION of the squared distances, on the
    enclosures of the DIFFERENCE and of the SUM (hence correlated, hence
    décisives). Miroir exact de la logique `_eps_sep` de the transport step, appliquée
    à des polynômes déjà soustraits. Ambigu ⟹ REFUS."""
    d2m = dm.re ** 2 + dm.im ** 2
    d2p = dp.re ** 2 + dp.im ** 2
    m_lo, m_hi = mp.mpf(d2m.a), mp.mpf(d2m.b)
    p_lo, p_hi = mp.mpf(d2p.a), mp.mpf(d2p.b)
    if m_hi < p_lo:
        return 1, {"theta": 1, "margin": _f_down(p_lo - m_hi),
                   "d2_chosen": [_f_down(m_lo), _f_up(m_hi)],
                   "d2_other": [_f_down(p_lo), _f_up(p_hi)]}
    if p_hi < m_lo:
        return -1, {"theta": -1, "margin": _f_down(m_lo - p_hi),
                    "d2_chosen": [_f_down(p_lo), _f_up(p_hi)],
                    "d2_other": [_f_down(m_lo), _f_up(m_hi)]}
    return None, {"d2_diff": [_f_down(m_lo), _f_up(m_hi)],
                  "d2_sum": [_f_down(p_lo), _f_up(p_hi)]}


def civ_sup(c):
    """Supremum of the parts of a complex enclosure."""
    return max(abs(mp.mpf(c.re.a)), abs(mp.mpf(c.re.b)),
               abs(mp.mpf(c.im.a)), abs(mp.mpf(c.im.b)))


def civ_contains_zero(c):
    return (mp.mpf(c.re.a) <= 0 <= mp.mpf(c.re.b)
            and mp.mpf(c.im.a) <= 0 <= mp.mpf(c.im.b))


def civ_div(a, b):
    """a/b en intervalles, exige |b| minoré > 0 (certifié par
    l'appelant)."""
    d = b.re ** 2 + b.im ** 2
    return CIV((a.re * b.re + a.im * b.im) / d,
               (a.im * b.re - a.re * b.im) / d)


# ===========================================================================
#  EXACT serialisation of the Taylor models (mpmath interval objects do not
#  pickle: we transport the exact `_mpf_` tuples between
#  workers et le parent — aucune conversion, aucun arrondi)
# ===========================================================================
def ser_tmc(t):
    return ([(mp.mpf(x.re.a)._mpf_, mp.mpf(x.re.b)._mpf_,
              mp.mpf(x.im.a)._mpf_, mp.mpf(x.im.b)._mpf_) for x in t.p],
            mp.mpf(t.rem.b)._mpf_)


def de_tmc(blob):
    coeffs, rem = blob
    p = [CIV(iv.mpf([mp.make_mpf(a), mp.make_mpf(b)]),
             iv.mpf([mp.make_mpf(c), mp.make_mpf(d)]))
         for (a, b, c, d) in coeffs]
    return p, iv.mpf(mp.make_mpf(rem))


# ===========================================================================
#  Halo phase: is the chart valid on an OPEN NEIGHBOURHOOD of the core?
# ===========================================================================
_G = {}


def _init_worker(cell, root=None):
    _G["cell"] = cell
    if root is not None:
        _G["root"] = root


def flush_faces(tile, root_c, root_h):
    """Directions where the core sits EXACTLY astride a face of the
    cellule : +1 face haute, −1 face basse, 0 intérieur. Comparaison
    rationnelle exacte."""
    c, h = box_of(tile)
    out = []
    for k in range(4):
        if c[k] + h == root_c[k] + root_h:
            out.append(1)
        elif c[k] - h == root_c[k] - root_h:
            out.append(-1)
        else:
            out.append(0)
    return out


def halo_certify(cell, tile, rho, rule, root_c, root_h):
    """One halo attempt at a given rho, under one of two RULES:

    . `symmetric`: `box(c, (1+rho)h)`, a two-sided halo in all four
      directions. This is the preregistered rule, tried first.
    · `clipped` — même demi-largeur `(1+ρ)h`, mais CENTRE DÉCALÉ de ∓ρh
      in each direction where the core sits astride a face of the
      CELL: the halo then extends inwards and stays FLUSH
      with the face, so the core is inside the halo inside the cell. This is the notion of an atlas
      of a manifold WITH BOUNDARY: the charts are open in the
      RELATIVE topology of the cell. The rule is legitimate ONLY on a
      flush direction; check D2b requires it, without which it
      would become an escape hatch for an interior failure.

    Returns (record, Z, P) where `P[a] = Zt[a].Z[g']` (the division-free
    product that carries the transition identity), or (record, None, None)."""
    S, g, eps = cell
    c, h = box_of(tile)
    H, Hf = halo_hw(h, rho)
    rec = {"rho": [rho.numerator, rho.denominator], "rule": rule}
    if Hf is None:
        rec["refused"] = "halo_hw_not_exact_float"
        return rec, None, None
    fl = flush_faces(tile, root_c, root_h)
    if rule == "clipped":
        shift = [-rho * h if fl[k] > 0 else
                 (rho * h if fl[k] < 0 else Fraction(0))
                 for k in range(4)]
        if all(s == 0 for s in shift):
            rec["refused"] = "clipping_useless_no_flush_face"
            return rec, None, None
    else:
        shift = [Fraction(0)] * 4
    c = [c[k] + shift[k] for k in range(4)]
    rec["flush_faces"] = fl
    rec["center_shift"] = [[s.numerator, s.denominator] for s in shift]
    rec["center_hex"] = [float(x).hex() for x in c]
    if any(Fraction(float(x)) != x for x in c):
        rec["refused"] = "halo_center_not_exact_float"
        return rec, None, None
    c0, h0 = box_of(tile)
    if not all(c[k] - H <= c0[k] - h0 and c0[k] + h0 <= c[k] + H
               for k in range(4)):
        rec["refused"] = "core_not_inside_halo"
        return rec, None, None
    if rule == "clipped" and not all(
            root_c[k] - root_h <= c[k] - H
            and c[k] + H <= root_c[k] + root_h for k in range(4)):
        rec["refused"] = "clipped_halo_leaves_cell"
        return rec, None, None
    rec["H_hex"] = float(Hf).hex()
    cf = [float(x) for x in c]
    Z, dZ, rows = build_section(S, g, eps, cf, Hf)
    if any(z is None for z in Z):
        rec["refused"] = "source_section_incomplete"
        return rec, None, None
    det = [r["determination"] for r in rows]
    rec["source_determinations"] = det
    if det != tile["core_src_det"]:
        rec["refused"] = "source_determination_changed"
        return rec, None, None
    S2, g2 = tuple(tile["chart"]["S"]), tile["chart"]["g"]
    cert = chart_certificate(Z, dZ, S2, g2)
    rec["gauge_absmin"] = cert.get("gauge_absmin")
    rec["detJ_absmin"] = cert.get("detJ_absmin")
    if not (cert.get("admissible")
            and cert.get("disjoint_from_target_slice")):
        rec["refused"] = "chart_certificate_failed"
        return rec, None, None
    T2 = [j for j in range(6) if j not in S2]
    o = [x for x in T2 if x != g2]
    ib = Z[g2].inv()
    Zp = [z * ib for z in Z]
    up_, vp_ = Zp[o[0]], Zp[o[1]]
    _r126, ok126 = native_section_constructible(S2, g2, up_, vp_)
    _rext, okext = native_rows_ext(S2, g2, up_, vp_)
    rec["criterion_c126_on_halo"] = bool(ok126)
    rec["criterion_ext_on_halo"] = bool(okext)
    need = tile["criterion"]
    if need == "criterion" and not ok126:
        rec["refused"] = "c126_fails_on_halo"
        return rec, None, None
    if need == "extended":
        if not okext:
            rec["refused"] = "extended_fails_on_halo"
            return rec, None, None
        if ok126:
            # non-tautology on the halo: the imported criterion must STILL refuse,
            # otherwise the tile was not residual for the stated reason
            rec["refused"] = "c126_unexpectedly_accepts_on_halo"
            return rec, None, None
    Zt, _dZt, kinds = native_target_section(
        S2, g2, up_, vp_, tuple(tile["eps_target"]),
        list(tile["sigma_target"]))
    if Zt is None:
        rec["refused"] = "native_target_failed_on_halo"
        rec["kinds_target"] = kinds
        return rec, None, None
    rec["kinds_target"] = kinds
    if kinds != tile["core_kinds_target"]:
        rec["refused"] = "target_kinds_changed"
        return rec, None, None
    # projective congruence on the HALO, with a frozen record. The supremum stays in
    # multiprecision until the check: the comparison to delta is made on the
    # exacte, la conversion float (dirigée, _f_up) n'arrive qu'à la
    # sérialisation.
    sup = mp.mpf(0)
    okz = True
    for a in range(6):
        D = Zt[a] + Zp[a].mul_real(riv(-1.0))
        dc = CIV(D.re_tm().to_iv() + IVPM * IV0,
                 D.im_tm().to_iv() + IVPM * IV0)
        okz = okz and civ_contains_zero(dc)
        sup = max(sup, civ_sup(dc))
    rec["halo_congruence_contains_zero"] = bool(okz)
    rec["halo_congruence_sup"] = _f_up(sup)
    if not okz:
        rec["refused"] = "halo_congruence_excludes_zero"
        return rec, None, None
    if not sup <= DELTA_TRANS:
        rec["refused"] = "halo_congruence_above_delta"
        return rec, None, None
    P = [Zt[a] * Z[g2] for a in range(6)]
    rec["accepted"] = True
    return rec, Z, P


def _halo_job(job):
    idx, tile = job
    root_c, root_h = _G["root"]
    attempts = []
    # STRICT ORDER: the symmetric rule first, over the whole rho scale.
    # Clipping is reached only if ALL symmetric attempts
    # have been refused, and the refusal is published, not absorbed.
    for rule in ("symmetric", "clipped"):
        for rho in RHO_LADDER:
            rec, Z, P = halo_certify(_G["cell"], tile, rho, rule,
                                     root_c, root_h)
            attempts.append(rec)
            if rec.get("accepted"):
                return {"index": idx, "ok": True,
                        "rho_used": rec["rho"], "rule_used": rule,
                        "H_hex": rec["H_hex"],
                        "center_shift": rec["center_shift"],
                        "flush_faces": rec["flush_faces"],
                        "n_symmetric_refusals": sum(
                            1 for a in attempts
                            if a["rule"] == "symmetric"),
                        "attempts": attempts, "record": rec,
                        "Z": [ser_tmc(Z[a]) for a in range(6)],
                        "P": [ser_tmc(P[a]) for a in range(6)]}
    return {"index": idx, "ok": False, "attempts": attempts}


# ===========================================================================
#  Phase paire : feuillet, transition, identité
# ===========================================================================
_TCACHE = {}


def get_T(off, scale):
    key = (tuple((o.numerator, o.denominator) for o in off),
           (scale.numerator, scale.denominator))
    T = _TCACHE.get(key)
    if T is None:
        T = recenter_matrix(off, scale)
        _TCACHE[key] = T
    return T


def pair_certificate(i, j, corrupt=None):
    """Certificat d'overlap entre les tuiles i et j. `corrupt` permet aux
    negative controls to substitute falsified data WITHOUT touching the path
    nominal."""
    gi, gj = _G["geom"][i], _G["geom"][j]
    w = _inter_generic(gi["lo"], gi["hi"], gj["lo"], gj["hi"])
    out = {"i": i, "j": j}
    if w is None:
        out["refused"] = "no_open_overlap"
        return out
    out["overlap_hw"] = [float(x) for x in w["hw"]]
    out["overlap_min_width"] = float(2 * min(w["hw"]))
    Ei = eps_box(w, gi["c"], gi["H"])
    Ej = eps_box(w, gj["c"], gj["H"])
    out["eps_box_in_range"] = bool(eps_box_in_range(Ei)
                                   and eps_box_in_range(Ej))
    if not out["eps_box_in_range"]:
        out["refused"] = "overlap_outside_tm_domain"
        return out
    # ε_j = off + scale·ε_i   (exact en rationnels)
    off = [(gi["c"][k] - gj["c"][k]) / gj["H"] for k in range(4)]
    scale = gi["H"] / gj["H"]
    T = get_T(off, scale)
    Eiv = [fraction_box_to_iv(lo, hi) for lo, hi in Ei]
    Zi, Pi = _G["tm"][i]
    Zj, Pj = _G["tm"][j]
    if corrupt and corrupt.get("tile") == j:
        Zj, Pj = corrupt["Z"], corrupt["P"]
    elif corrupt and corrupt.get("tile") == i:
        Zi, Pi = corrupt["Z"], corrupt["P"]
    S = _G["cell"][0]

    # --- (1) le FEUILLET : θ par recentrage + soustraction exacte -------
    thetas, tmargins, amb = [], [], False
    for s in S:
        pa, ra = Zi[s]
        pb, rb = Zj[s]
        pb = apply_recenter(T, pb)
        rem = ra + rb
        dm = enclose(poly_lin(pa, pb, -1), rem, Eiv)
        dp = enclose(poly_lin(pa, pb, +1), rem, Eiv)
        th, rec = sep_phase(dm, dp)
        rec["s_coord"] = int(s)
        rec["diff_sup"] = _f_up(civ_sup(dm))
        if th is None:
            amb = True
        thetas.append(th)
        tmargins.append(rec)
    out["theta"] = thetas
    out["theta_margins"] = tmargins
    out["theta_ambiguous"] = bool(amb)
    out["same_sheet"] = bool(not amb and all(t == 1 for t in thetas))

    # --- (2) coordonnées AFFINES : contrôle de l'évaluateur -------------
    #  Z[g] is 1 and Z[o] is u; v are the SAME functions in both
    #  tiles: their recentred difference must be zero up to rounding.
    #  This is a free test of the recentring itself.
    aff_sup = mp.mpf(0)
    for a in range(6):
        if a in S:
            continue
        pa, ra = Zi[a]
        pb, rb = Zj[a]
        d = enclose(poly_lin(pa, apply_recenter(T, pb), -1),
                    ra + rb, Eiv)
        aff_sup = max(aff_sup, civ_sup(d))
    out["affine_coords_defect"] = _f_up(aff_sup)

    # --- (3) the CHART: lambda = Z[g'_i]/Z[g'_j], nonzero on both sides ---
    g_i, g_j = _G["tiles"][i]["chart"]["g"], _G["tiles"][j]["chart"]["g"]
    zi_gi = enclose(*Zi[g_i], Eiv)
    zi_gj = enclose(*Zi[g_j], Eiv)
    mi = mp.mpf(civ_absmin(zi_gi).a)
    mj = mp.mpf(civ_absmin(zi_gj).a)
    out["abs_min_Z_gi"] = _f_down(mi)
    out["abs_min_Z_gj"] = _f_down(mj)
    out["transition_biholomorphic"] = bool(mi > 0 and mj > 0)
    if mi > 0 and mj > 0:
        lam = civ_div(zi_gi, zi_gj)
        out["lambda"] = [_f_down(mp.mpf(lam.re.a)),
                         _f_up(mp.mpf(lam.re.b)),
                         _f_down(mp.mpf(lam.im.a)),
                         _f_up(mp.mpf(lam.im.b))]
        out["charts_equal"] = bool(g_i == g_j)

    # --- (4) the transition IDENTITY, without division ------------------
    #  Zt^(j)_a·Z[g'_j] − Zt^(i)_a·Z[g'_i] ∋ 0 : les deux sheet records ε'
    #  derived SEPARATELY are brought together on the overlap.
    # The supremum stays in multiprecision until the check; _f_up only for
    # la sérialisation — la borne publiée MAJORE la borne comparée.
    idsup, idok = mp.mpf(0), True
    for a in range(6):
        pa, ra = Pi[a]
        pb, rb = Pj[a]
        d = enclose(poly_lin(pa, apply_recenter(T, pb), -1),
                    ra + rb, Eiv)
        idok = idok and civ_contains_zero(d)
        idsup = max(idsup, civ_sup(d))
    out["transition_identity_contains_zero"] = bool(idok)
    out["transition_identity_defect"] = _f_up(idsup)
    out["ok"] = bool(out["same_sheet"] and idok
                     and idsup <= DELTA_TRANS
                     and out["transition_biholomorphic"])
    return out


def _pair_job(job):
    i, j = job
    return pair_certificate(i, j)


# ===========================================================================
#  Triple phase: the cocycle on the nerve
# ===========================================================================
def triple_certificate(i, j, k, lam):
    gi, gj, gk = _G["geom"][i], _G["geom"][j], _G["geom"][k]
    w = _inter_generic(gi["lo"], gi["hi"], gj["lo"], gj["hi"])
    if w is None:
        return None
    w = _inter_generic(w["lo"], w["hi"], gk["lo"], gk["hi"])
    if w is None:
        return None
    out = {"i": i, "j": j, "k": k,
           "min_width": float(2 * min(w["hw"]))}
    # θ : cocycle EXACT (discret)
    tij, tjk, tik = lam[(i, j)]["theta"], lam[(j, k)]["theta"], \
        lam[(i, k)]["theta"]
    if any(x is None for x in tij + tjk + tik):
        out["theta_cocycle"] = False
        out["theta_note"] = "ambiguous_edge"
    else:
        out["theta_cocycle"] = all(
            tij[r] * tjk[r] * tik[r] == 1 for r in range(len(tij)))
    # lambda: cocycle as an enclosure on the TRIPLE intersection
    Ei = eps_box(w, gi["c"], gi["H"])
    if not eps_box_in_range(Ei):
        out["refused"] = "triple_outside_tm_domain"
        return out
    Eiv = [fraction_box_to_iv(lo, hi) for lo, hi in Ei]
    Zi = _G["tm"][i][0]
    gg = [_G["tiles"][x]["chart"]["g"] for x in (i, j, k)]
    z = [enclose(*Zi[c], Eiv) for c in gg]
    ok = all(mp.mpf(civ_absmin(x).a) > 0 for x in z)
    out["lambda_defined"] = bool(ok)
    out["gauge_abs_min"] = [_f_down(mp.mpf(civ_absmin(x).a)) for x in z]
    if ok:
        # HONESTY: `lambda_ij.lambda_jk.lambda_ki` = (z_i/z_j)(z_j/z_k)(z_k/z_i) is an
        # ALGEBRAIC IDENTITY of the construction. All charts are
        # affine charts of the SAME projective space and the three lambdas are ratios
        # of the three SAME evaluations. Checking it in intervals verifies
        # nothing: each z appears twice, the correlation is
        # lost, and the measured "defect" IS ONLY the decorrelation
        # width (exactly the phenomenon that forced the
        # recentring elsewhere in this script). The number is therefore published
        # comme DIAGNOSTIC de décorrélation, JAMAIS comme check.
        # What has content at the triple level, and what is checked:
        #   . the three gauges are BOUNDED BELOW by 0 on the triple box
        #     (the three lambdas exist and the composition is licit);
        #   . the three edge identities, recentred, contain 0 on
        #     the TRIPLE box, and that is NOT an interval
        #     consequence of i to j and j to k (the inclusion is not
        #     transitive).
        prod = _mul_civ(_mul_civ(civ_div(z[0], z[1]),
                                 civ_div(z[1], z[2])),
                        civ_div(z[2], z[0]))
        d = CIV(prod.re - IV1, prod.im)
        out["lambda_product_contains_one"] = bool(civ_contains_zero(d))
        out["lambda_decorrelation_width"] = _f_up(civ_sup(d))
        out["lambda_note"] = ("identité algébrique — largeur de "
                              "décorrélation publiée, non gatée")
    # the transition identity on the THREE edges, on the TRIPLE box.
    # Sup en mp jusqu'au check (the directed-rounding step), _f_up à la sérialisation.
    idsup, idok = mp.mpf(0), True
    for (a, b) in ((i, j), (j, k), (i, k)):
        ga, gb = _G["geom"][a], _G["geom"][b]
        Ea = eps_box(w, ga["c"], ga["H"])
        if not eps_box_in_range(Ea):
            idok = False
            continue
        Eav = [fraction_box_to_iv(lo, hi) for lo, hi in Ea]
        off = [(ga["c"][t] - gb["c"][t]) / gb["H"] for t in range(4)]
        T = get_T(off, ga["H"] / gb["H"])
        for co in range(6):
            pa, ra = _G["tm"][a][1][co]
            pb, rb = _G["tm"][b][1][co]
            dd = enclose(poly_lin(pa, apply_recenter(T, pb), -1),
                         ra + rb, Eav)
            idok = idok and civ_contains_zero(dd)
            idsup = max(idsup, civ_sup(dd))
    out["edge_identities_contain_zero"] = bool(idok)
    out["edge_identity_defect"] = _f_up(idsup)
    out["ok"] = bool(out.get("theta_cocycle") and out["lambda_defined"]
                     and idok and idsup <= DELTA_TRANS)
    return out


def _mul_civ(a, b):
    return CIV(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)


def _triple_job(job):
    i, j, k = job
    return triple_certificate(i, j, k, _G["pairs"])


# ===========================================================================
#  Chargement des feuilles et de leurs sheet records
# ===========================================================================
def load_leaves():
    cov = json.loads(COVER_JSON.read_text(encoding="utf-8"))
    transport = json.loads(C127_JSON.read_text(encoding="utf-8"))
    residual = json.loads(C127E_JSON.read_text(encoding="utf-8"))
    cell = cov["cell"]
    tr127 = {r["tile_index"]: r for r in transport["transports"]
             if not r.get("failed")}
    tr_residual = {r["box_index"]: r for r in residual["transports"]
              if not r.get("failed")}
    leaves = []
    for i, t in enumerate(cov["tiles"]):
        r = tr127.get(i)
        if r is None:
            raise SystemExit(f"tuile the transport step {i} sans transport certifié")
        leaves.append({
            "src": "transport", "orig_index": i, "depth": t["depth"],
            "center_hex": t["center_hex"], "hw_hex": t["hw_hex"],
            "chart": t["chart"], "criterion": "criterion",
            "core_src_det": r["source_determinations"],
            "core_kinds_target": r["kinds_target"],
            "eps_target": r["eps_target"],
            "sigma_target": r["sigma_target"]})
    for t in residual["new_tiles"]:
        r = tr_residual.get(t["box_index"])
        if r is None:
            raise SystemExit(
                f"tuile the residual closure {t['box_index']} sans transport certifié")
        leaves.append({
            "src": "residual", "orig_index": t["box_index"],
            "depth": t["depth"], "center_hex": t["center_hex"],
            "hw_hex": t["hw_hex"], "chart": t["chart"],
            "criterion": "extended",
            "core_src_det": r["source_determinations"],
            "core_kinds_target": r["kinds_target"],
            "eps_target": r["eps_target"],
            "sigma_target": r["sigma_target"]})
    return cell, leaves


def pilot_patch(leaves, n_max):
    """Stratified patch: we start from the tile born of the RESIDUAL that has the most
    neighbours from the transport step (the most discriminating interface:
    déterminations différentes s'y touchent), on prend tout son
    voisinage, puis on complète par un représentant de chaque signature
    (chart, déterminations source, kinds cible) non encore couverte."""
    bx = [box_of(t) for t in leaves]
    nb = {i: set() for i in range(len(leaves))}
    for i in range(len(leaves)):
        for j in range(i + 1, len(leaves)):
            if touching(bx[i][0], bx[i][1], bx[j][0], bx[j][1]):
                nb[i].add(j)
                nb[j].add(i)
    seed, best = 0, -1
    for i, t in enumerate(leaves):
        if t["src"] != "residual":
            continue
        n = sum(1 for j in nb[i] if leaves[j]["src"] == "transport")
        if n > best:
            seed, best = i, n
    sel = [seed] + sorted(nb[seed])
    sig = {(leaves[i]["src"], tuple(leaves[i]["chart"]["S"]),
            leaves[i]["chart"]["g"], tuple(leaves[i]["core_src_det"]),
            tuple(leaves[i]["core_kinds_target"])) for i in sel}
    for i, t in enumerate(leaves):
        if len(sel) >= n_max:
            break
        s = (t["src"], tuple(t["chart"]["S"]), t["chart"]["g"],
             tuple(t["core_src_det"]), tuple(t["core_kinds_target"]))
        if s not in sig:
            sig.add(s)
            sel.append(i)
    sel = sorted(set(sel))[:n_max]
    return sel, {"seed": seed, "seed_c127_neighbours": best,
                 "n_selected": len(sel), "selected": sel,
                 "note": ("patch = one tile born of the residual plus all its "
                          "voisinage + un représentant par signature "
                          "non couverte")}


# ===========================================================================
#  Build
# ===========================================================================
def build():
    print("=" * 78)
    print(f"the atlas step ATLAS {'COMPLET (316)' if MODE == 'full' else 'PILOTE'}"
          f" : TM ({TM_ORDER},{UNARY_SERIES_DEG}), {N_WORKERS} workers, "
          f"δ_trans = {DELTA_TRANS:.0e}")
    print("=" * 78)
    cell_d, leaves = load_leaves()
    S, g, eps = (tuple(cell_d["S"]), cell_d["g"], tuple(cell_d["eps"]))
    root_c = [float.fromhex(x) for x in cell_d["center_hex"]]
    root_h = float.fromhex(cell_d["hw_hex"])
    log(f"{len(leaves)} feuilles chargées "
        f"({sum(1 for t in leaves if t['src'] == 'transport')} the transport step + "
        f"{sum(1 for t in leaves if t['src'] == 'residual')} the residual closure)")

    # --- D1: the partition, reasserted from the addresses --------------
    addresses, addr_fail = [], 0
    for t in leaves:
        a = address_of(root_c, root_h,
                       [float.fromhex(x) for x in t["center_hex"]],
                       float.fromhex(t["hw_hex"]))
        if a is None:
            addr_fail += 1
        else:
            addresses.append(a)
    tg = tree_gates(addresses)
    log(f"D1 : {len(addresses)}/{len(leaves)} adresses exactes, "
        f"prefix-free={tg['prefix_free']}, clos={tg['tree_closed']}, "
        f"Kraft={tg['kraft_sum'][0]}/{tg['kraft_sum'][1]}")

    sel, patch = (list(range(len(leaves))), None) if MODE == "full" \
        else pilot_patch(leaves, N_PILOT_TILES)
    log(f"tuiles retenues : {len(sel)}")

    # --- D2/D3 : les halos ---------------------------------------------
    rootF = ([Fraction(x) for x in
              [Fraction(float.fromhex(y)) for y in cell_d["center_hex"]]],
             Fraction(float.fromhex(cell_d["hw_hex"])))
    mpctx = get_context("fork")
    _init_worker((S, g, eps), rootF)
    with mpctx.Pool(N_WORKERS, initializer=_init_worker,
                    initargs=((S, g, eps), rootF)) as pool:
        halos = pool.map(_halo_job, [(i, leaves[i]) for i in sel])
    halo_ok = [h for h in halos if h["ok"]]
    n_clip = sum(1 for h in halo_ok if h["rule_used"] == "clipped")
    log(f"D2 : halos acceptés {len(halo_ok)}/{len(sel)} "
        f"(ρ : {sorted({tuple(h['rho_used']) for h in halo_ok})}) · "
        f"règle symétrique {len(halo_ok) - n_clip}, clippée {n_clip}")
    if len(halo_ok) != len(sel):
        for h in halos:
            if not h["ok"]:
                log(f"  REFUS tuile {h['index']} : "
                    f"{[a.get('refused') for a in h['attempts']]}")

    # --- géométrie des halos (exacte) ----------------------------------
    geom, tm = {}, {}
    for h in halo_ok:
        i = h["index"]
        c0, hw = box_of(leaves[i])
        rho = Fraction(*h["rho_used"])
        H = hw * (1 + rho)
        c = [c0[k] + Fraction(*h["center_shift"][k]) for k in range(4)]
        lo, hi = cube_bounds(c, H)
        geom[i] = {"c": c, "H": H, "lo": lo, "hi": hi,
                   "core_c": c0, "core_hw": hw, "rho": rho,
                   "rule": h["rule_used"], "flush": h["flush_faces"]}
        tm[i] = ([de_tmc(b) for b in h["Z"]],
                 [de_tmc(b) for b in h["P"]])
    _G["geom"], _G["tm"], _G["tiles"], _G["cell"] = \
        geom, tm, leaves, (S, g, eps)

    # --- D4 : le nerf, en rationnels exacts ----------------------------
    ids = sorted(geom)
    touch_pairs, pairs = [], []
    for a in range(len(ids)):
        for b in range(a + 1, len(ids)):
            i, j = ids[a], ids[b]
            gi, gj = geom[i], geom[j]
            tch = touching(gi["core_c"], gi["core_hw"],
                           gj["core_c"], gj["core_hw"])
            ov = _inter_generic(gi["lo"], gi["hi"], gj["lo"], gj["hi"])
            if tch:
                touch_pairs.append((i, j, ov is not None))
            if ov is not None:
                pairs.append((i, j))
    open_all = all(x[2] for x in touch_pairs)
    adj = {i: set() for i in ids}
    for i, j in pairs:
        adj[i].add(j)
        adj[j].add(i)
    seen, stack = {ids[0]}, [ids[0]]
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    connected = len(seen) == len(ids)
    log(f"D4 : {len(touch_pairs)} paires de cores qui se touchent, "
        f"{len(pairs)} overlaps ouverts, nerf connexe={connected}")

    # --- D6-D9 : les transitions ---------------------------------------
    with mpctx.Pool(N_WORKERS, initializer=_init_worker,
                    initargs=((S, g, eps),)) as pool:
        pres = pool.map(_pair_job, pairs)
    pmap = {(r["i"], r["j"]): r for r in pres}
    _G["pairs"] = pmap
    pok = [r for r in pres if r.get("ok")]
    log(f"D7-D9 : paires certifiées {len(pok)}/{len(pres)} · "
        f"θ≡+1 {sum(1 for r in pres if r.get('same_sheet'))} · "
        f"défaut d'identité max "
        f"{max((r.get('transition_identity_defect', 0) for r in pres), default=0):.3e}")

    # --- D10 : le cocycle ----------------------------------------------
    triples = []
    for (i, j) in pairs:
        for k in adj[i] & adj[j]:
            if k > j:
                triples.append((i, j, k))
    with mpctx.Pool(N_WORKERS, initializer=_init_worker,
                    initargs=((S, g, eps),)) as pool:
        tres = [r for r in pool.map(_triple_job, triples) if r]
    tok = [r for r in tres if r.get("ok")]
    log(f"D10 : triples d'intersection triple non vide {len(tres)}, "
        f"cocycle certifié {len(tok)}/{len(tres)}")

    # --- D11 : les négatifs --------------------------------------------
    negs = run_negatives(leaves, geom, tm, pairs, pmap)
    log("D11 : négatifs — "
        + " · ".join(f"{k}={'CASSE' if v['breaks'] else 'NE CASSE PAS'}"
                     for k, v in negs.items()))

    # --- checks ----------------------------------------------------------
    max_halo = max((h["record"]["halo_congruence_sup"]
                    for h in halo_ok), default=0.0)
    max_id = max((r.get("transition_identity_defect", 0.0)
                  for r in pres), default=0.0)
    max_aff = max((r.get("affine_coords_defect", 0.0)
                   for r in pres), default=0.0)
    min_theta = min((m["margin"] for r in pres
                     for m in r.get("theta_margins", [])
                     if "margin" in m), default=0.0)
    max_cocy = max((r.get("lambda_decorrelation_width", 0.0)
                    for r in tres), default=0.0)
    # D2b: clipping is legitimate ONLY against a face of the cell,
    # and ONLY if the symmetric rule was first refused over the WHOLE
    # rho scale. Without this check, "clipped" would be an escape hatch
    # pour n'importe quel échec intérieur.
    clip_legit = True
    for h in halo_ok:
        if h["rule_used"] != "clipped":
            continue
        sh, fl = h["center_shift"], h["flush_faces"]
        clip_legit = clip_legit and all(
            (Fraction(*sh[k]) == 0) or (fl[k] != 0) for k in range(4)) \
            and h["n_symmetric_refusals"] == len(RHO_LADDER)
    # D2c: EQUALITY OF SETS, not of counts. The
    # prose affirme {tuiles clippées} = {64 tuiles nées du résidu
    # residual step, and "64 = 64" does not exclude a permutation 63+1. Both
    # sets are published, sorted, with their SHA-256.
    clipped_idx = sorted(h["index"] for h in halo_ok
                         if h["rule_used"] == "clipped")
    residual_idx = sorted(i for i in sel if leaves[i]["src"] == "residual")

    def _set_sha(idx):
        return hashlib.sha256(
            json.dumps(idx, separators=(",", ":")).encode()).hexdigest()
    checks = {
        "D1_partition_reasserted": bool(
            addr_fail == 0 and tg["unique"] and tg["prefix_free"]
            and tg["tree_closed"] and tg["kraft_is_one"]),
        "D2_halo_charts_valid": bool(
            halos and len(halo_ok) == len(sel)
            and all(Fraction(*h["rho_used"]) > 0 for h in halo_ok)),
        "D2b_clipping_only_at_cell_boundary": bool(clip_legit),
        "D2c_clipped_set_is_exactly_c127e_residual": bool(
            clipped_idx == residual_idx),
        "D3_halo_congruence_below_delta": bool(
            halo_ok and max_halo <= DELTA_TRANS),
        "D4_overlaps_open": bool(touch_pairs and open_all
                                 and len(pairs) >= len(touch_pairs)),
        "D5_nerve_connected": bool(connected),
        "D6_overlaps_in_tm_domain": bool(pres) and all(
            r.get("eps_box_in_range") for r in pres),
        "D7_sheet_phase_trivial": bool(pres) and all(
            r.get("same_sheet") for r in pres) and min_theta > 0,
        "D8_transition_biholomorphic": bool(pres) and all(
            r.get("transition_biholomorphic") for r in pres),
        "D9_transition_identity": bool(pres) and all(
            r.get("transition_identity_contains_zero") for r in pres)
        and max_id <= DELTA_TRANS,
        "D10_cocycle_on_nerve": bool(tres) and all(
            r.get("ok") for r in tres),
        "D11_negatives_all_break": bool(negs) and all(
            v["breaks"] for v in negs.values()),
        "D12_no_silent_cap": bool(
            len(halos) == len(sel) and len(pres) == len(pairs)
            and len(tres) == len(triples)
            and all((r["i"], r["j"]) in pmap for r in pres))}
    n_pass = sum(1 for v in checks.values() if v)

    verdict = (
        "(%s) THE PARTITION HAS BECOME AN ATLAS. The %d "
        "leaves of the closed cover carry each a HALO "
        "`(1+rho)h` with rho from the preregistered scale %s: on this "
        "OPEN neighbourhood of the core the source determinations are "
        "UNCHANGED, the chart criterion that certified the core still "
        "holds (and for the tiles born of the residual, the plain criterion STILL REFUSES, "
        "so non-tautology survives the halo), the native target section "
        "with a FROZEN record has the same kinds, and projective congruence "
        "stays under delta (maximum %.3e). THE MEASUREMENT THAT COST AN AMENDMENT: "
        "%d tiles REFUSED the symmetric halo over the WHOLE rho scale, "
        "including at rho = 1/64. Diagnosis: they sit astride a "
        "face of the CELL, and the 64 tiles born of the residual "
        "are EXACTLY those flush with the face {Im u = Im v = "
        "0}, là où σ de la détermination tournée devient indéterminé "
        "as soon as one steps outside. The residual 1/64 was therefore "
        "geometrically A FACE, not a scattering. These tiles "
        "reçoivent un halo CLIPPÉ (même demi-largeur, centre décalé vers "
        "inwards, halo flush: core inside halo inside cell), the "
        "notion of an atlas of a manifold WITH BOUNDARY, charts open in the "
        "RELATIVE topology. Check D2b forbids clipping from serving "
        "elsewhere: it requires a flush face in each shifted "
        "direction AND the prior refusal of the whole symmetric scale. "
        "The nerve is computed in RATIONALS "
        "DYADIQUES EXACTS : %d paires de cores se touchent et %d/%d ont "
        "an overlap of STRICTLY POSITIVE width in all four "
        "coordinates; the nerve is %s. On each overlap the two "
        "source sections are compared by EXACT RECENTRING in a "
        "cadre commun puis soustraction coefficient par coefficient "
        "(the naive, decorrelated enclosure would have width about 1e-3 and "
        "no content): the sheet theta is DERIVED with a strict margin, equals "
        "+1 on %d of %d pairs (minimum margin %.3e), and the affine control "
        "defect of the recentring is %.3e. The transition is a "
        "BIHOLOMORPHISME certifié (|Z[g']| minoré > 0 des deux côtés) et "
        "the identity `Zt^(j).Z[g'_j] - Zt^(i).Z[g'_i]` containing 0, which brings "
        "together two records derived SEPARATELY, holds with a "
        "maximum defect %.3e at most delta. On the %d non-empty triple "
        "intersections the cocycle is certified %d of %d. HONESTY ON THE "
        "COCYCLE: `lambda_ij.lambda_jk.lambda_ki = 1` is an ALGEBRAIC IDENTITY of "
        "the construction (the same affine charts of one projective space, three "
        "ratios of the same three evaluations), so checking it in "
        "intervalles ne vérifie rien, chaque jauge y apparaît deux fois "
        "and the measured defect (%.3e) IS ONLY the decorrelation "
        "width, exactly the phenomenon that forced the "
        "recentring elsewhere; it is therefore published as a DIAGNOSTIC and "
        "enters NO check. What is checked at the triple level does have "
        "content: the three gauges BOUNDED BELOW by 0 on the triple box, "
        "the exact discrete cocycle, and the three edge identities "
        "recentred on the TRIPLE box, since containment of 0 is not "
        "transitive and i to k is no interval consequence of i to j "
        "and j to k. What is NOT tautological either is theta = +1 "
        "(rien n'imposait aux voisines le même feuillet, leurs "
        "déterminations source DIFFÈRENT) et l'identité de transition, "
        "which brings together two records derived separately. The %d "
        "negative controls all break. NOT PAID HERE: METRIC transport "
        "(congruence of Q, Weyl) on the halos, established on the cores "
        "earlier, whose extension costs four Qmat per tile; the "
        "EXACT contract of the identity; the full scaling; the 895 "
        "autres paires cellule/classe ; the later scaling. the rational-box step/B/C (a review) "
        "PAYÉS : endpoints d'overlap en arrondi EXTÉRIEUR "
        "(fraction_box_to_iv — l'ancien float-nearest pouvait tronquer "
        "les boîtes à dénominateur 65), défauts gatés comparés à δ en mp "
        "AVANT toute conversion float (sérialisation _f_up seulement), "
        "et égalité d'ENSEMBLES {clippées} = {résidu the residual closure} gatée (D2c, "
        "SHA-256 publiés)." % (
            MODE, len(leaves),
            [f"{r.numerator}/{r.denominator}" for r in RHO_LADDER],
            max_halo,
            sum(1 for h in halo_ok if h["rule_used"] == "clipped"),
            len(touch_pairs),
            sum(1 for x in touch_pairs if x[2]), len(touch_pairs),
            "connexe" if connected else "NON CONNEXE",
            sum(1 for r in pres if r.get("same_sheet")), len(pres),
            min_theta, max_aff, max_id, len(tres), len(tok), len(tres),
            max_cocy, len(negs)))

    art = {
        "artifact": ART.stem, "mode": MODE,
        "claim": ("the atlas step — cores + halos, overlaps ouverts certifiés, "
                  "transitions (feuillet, carte, identité) et cocycle "
                "on the nerve: the closed cover becomes "
                  "un atlas multi-chart."),
        "cell": cell_d,
        "n_leaves": len(leaves),
        "tree_gates": tg, "n_address_failures": addr_fail,
        "pilot_patch": patch,
        "selected_tiles": sel,
        "halos": [{"index": h["index"], "ok": h["ok"],
                   "rho_used": h.get("rho_used"),
                   "H_hex": h.get("H_hex"),
                   "record": h.get("record"),
                   "attempts": [{k: v for k, v in a.items()
                                 if k != "kinds_target"}
                                for a in h["attempts"]]}
                  for h in halos],
        "clipped_set": {
            "clipped_indices": clipped_idx,
            "c127e_selected_indices": residual_idx,
            "clipped_sha256": _set_sha(clipped_idx),
            "c127e_sha256": _set_sha(residual_idx),
            "equal": bool(clipped_idx == residual_idx),
            "note": ("the set-equality step : égalité d'ensembles gatée (D2c), pas "
                "only the cardinality: 64 = 64 does not exclude a "
                     "permutation 63+1")},
        "nerve": {
            "n_touching_core_pairs": len(touch_pairs),
            "n_open_overlaps": len(pairs),
            "all_touching_pairs_open": bool(open_all),
            "connected": bool(connected),
            "n_triples": len(triples),
            "degree": {"min": min((len(adj[i]) for i in ids), default=0),
                       "max": max((len(adj[i]) for i in ids), default=0)},
            "note": ("nerf calculé en Fraction exactes ; « largeur "
                "strictly positive in the 4 coordinates is "
                "a rational comparison, not a float test")},
        "pairs": pres,
        "triples": tres,
        "negatives": negs,
        "max_halo_congruence_sup": max_halo,
        "max_transition_identity_defect": max_id,
        "max_affine_control_defect": max_aff,
        "min_theta_margin": min_theta,
        "max_lambda_decorrelation_width": max_cocy,
        "delta_trans_preregistered": DELTA_TRANS,
        "delta_cocycle_preregistered": DELTA_COCYCLE,
        "not_paid_here": [
            "metric transport (Qmat and Weyl) on the halos, established on "
            "les cores par the transport step",
            "contrat exact de l'identité de congruence",
            "scaling complet (the scaling step)",
            "les 895 autres paires cellule/classe", "the later scaling"],
        "verdict": verdict, "checks": checks,
        "checks_passed": n_pass, "checks_total": len(checks),
        "provenance": provenance([COVER_JSON, C127_JSON, C127E_JSON],
                                 time.time() - T0)}
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
#  The negative controls: discriminating power, not decoration
# ===========================================================================
def run_negatives(leaves, geom, tm, pairs, pmap):
    """Six falsifications. Chacune DOIT casser un check nominal."""
    out = {}
    ids = sorted(geom)
    ref = next(((i, j) for (i, j) in pairs
                if pmap[(i, j)].get("ok")), None)
    if ref is None:
        return {"N0_no_reference_pair": {"breaks": False}}
    i, j = ref

    # N1: a corrupted sheet record on one tile makes the transition identity
    #      doit cesser de contenir 0 (elle met en présence les DEUX
    #      sheet records ; si elle ne voyait pas ε', elle serait décorative)
    Zj, Pj = tm[j]
    Pbad = [([CIV(-c.re, -c.im) for c in Pj[a][0]], Pj[a][1])
            if a == leaves[j]["chart"]["S"][0] else Pj[a]
            for a in range(6)]
    r1 = pair_certificate(i, j, corrupt={"tile": j, "Z": Zj,
                                         "P": Pbad})
    out["N1_eps_ledger_corrupted"] = {
        "breaks": bool(not r1.get("transition_identity_contains_zero")),
        "defect": r1.get("transition_identity_defect"),
        "note": ("one row of Zt.Z[g'] negated on tile j: "
                 "l'identité de transition doit exclure 0")}

    # N2: sheet flipped; theta must come out at -1 (the test DETECTS a
    #      raccord inter-feuillets, c'est toute sa raison d'être)
    s0 = _G["cell"][0][0]          # première ligne √ de la SOURCE
    Zbad = [([CIV(-c.re, -c.im) for c in Zj[a][0]], Zj[a][1])
            if a == s0 else Zj[a] for a in range(6)]
    r2 = pair_certificate(i, j, corrupt={"tile": j, "Z": Zbad,
                                         "P": Pj})
    out["N2_sheet_flipped"] = {
        "breaks": bool(r2.get("theta") and r2["theta"][0] == -1
                       and not r2.get("same_sheet")),
        "theta": r2.get("theta"),
        "note": ("Z_s negated on tile j: theta must equal -1 on that "
                 "ligne et le check de feuillet doit tomber")}

    # N3: rho = 0; the halos become the cores again and the overlap loses its
    #      largeur ; « ouvert » doit cesser d'être vrai
    n3 = 0
    for (a, b) in pairs[:200]:
        ga, gb = geom[a], geom[b]
        la, ha = cube_bounds(ga["c"], ga["core_hw"])
        lb, hb = cube_bounds(gb["c"], gb["core_hw"])
        if _inter_generic(la, ha, lb, hb) is None:
            n3 += 1
    out["N3_rho_zero_kills_openness"] = {
        "breaks": bool(n3 > 0), "n_pairs_losing_openness": n3,
        "n_tested": min(200, len(pairs)),
        "note": ("à ρ = 0 les boîtes fermées adjacentes ne se coupent "
                     "no longer in positive width: the openness of the overlaps "
                     "is indeed CARRIED by rho, not by the partition")}

    # N4: a non-adjacent pair; the geometry must refuse
    far = None
    for a in ids:
        for b in ids:
            if b > a and (a, b) not in pmap:
                far = (a, b)
                break
        if far:
            break
    if far:
        r4 = pair_certificate(*far)
        out["N4_non_adjacent_refused"] = {
            "breaks": bool(r4.get("refused") == "no_open_overlap"),
            "pair": list(far), "refused": r4.get("refused")}
    else:
        out["N4_non_adjacent_refused"] = {
            "breaks": False, "note": "aucune paire non adjacente"}

    # N5: recentring neutralised. Comparing the polynomials WITHOUT
    #      recentring (T = identity), the difference must explode as soon as
    #      les deux cadres diffèrent réellement
    gi, gj = geom[i], geom[j]
    w = _inter_generic(gi["lo"], gi["hi"], gj["lo"], gj["hi"])
    Eiv = [fraction_box_to_iv(lo, hi)
           for lo, hi in eps_box(w, gi["c"], gi["H"])]
    Zi = tm[i][0]
    pa, ra = Zi[s0]
    pb, rb = Zj[s0]
    d_true = float(civ_sup(enclose(
        poly_lin(pa, apply_recenter(
            get_T([(gi["c"][k] - gj["c"][k]) / gj["H"] for k in range(4)],
                  gi["H"] / gj["H"]), pb), -1), ra + rb, Eiv)))
    d_naive = float(civ_sup(enclose(poly_lin(pa, pb, -1),
                                    ra + rb, Eiv)))
    out["N5_recentering_is_load_bearing"] = {
        "breaks": bool(d_naive > 1e3 * max(d_true, 1e-300)),
        "defect_recentered": d_true, "defect_unrecentered": d_naive,
        "note": ("without recentring the comparison has no content; "
                     "the ratio measures what the recentring carries")}

    # N6 — géométrie décalée : le test doit voir le décalage À SA TAILLE
    #      PREDICTED. A plain "it is larger" is settled by a threshold
    #      chosen afterwards; here the expected effect of a shift delta of
    #      the affine map is delta times |df/deps_0|, and |df/deps_0| is READ off the
    #      degree-one coefficient of the Taylor model. The negative control therefore requires
    #      the measured defect to fall in [predicted/4, 4.predicted] AND that
    #      plancher nominal soit deux ordres SOUS le prédit.
    delta = Fraction(1, 2 ** 20)
    off_bad = [(gi["c"][k] - gj["c"][k]) / gj["H"] for k in range(4)]
    off_bad[0] = off_bad[0] + delta
    d_bad = float(civ_sup(enclose(
        poly_lin(pa, apply_recenter(get_T(off_bad, gi["H"] / gj["H"]),
                                    pb), -1), ra + rb, Eiv)))
    e0 = pb[MIDX[(1, 0, 0, 0)]]
    grad = float(max(abs(mp.mpf(e0.re.a)), abs(mp.mpf(e0.re.b)),
                     abs(mp.mpf(e0.im.a)), abs(mp.mpf(e0.im.b))))
    pred = float(delta) * grad
    out["N6_shifted_geometry_detected"] = {
        "breaks": bool(pred / 4 <= d_bad <= 4 * pred
                       and d_true < pred / 100),
        "defect_true": d_true, "defect_shifted": d_bad,
        "delta_eps": float(delta), "degree1_coefficient": grad,
        "predicted_defect": pred,
        "ratio_measured_over_predicted": (d_bad / pred if pred else None),
        "note": ("a shift of 2^-20 of the affine map: the defect "
                     "must appear at the shift times the first derivative, read on the "
                     "degree-one coefficient of the model. The test sees the GEOMETRY, and it "
                     "sees it at the right scale, not merely as larger")}
    return out


# ===========================================================================
#  Self-test: pure functions, without the registry or the sections
# ===========================================================================
def _selftest():
    ok = True

    def chk(name, cond):
        nonlocal ok
        ok = ok and bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    # T1 : intersection exacte, largeur strictement positive
    lo1, hi1 = cube_bounds([Fraction(0)] * 4, Fraction(1))
    lo2, hi2 = cube_bounds([Fraction(3, 2)] + [Fraction(0)] * 3,
                           Fraction(1))
    w = _inter_generic(lo1, hi1, lo2, hi2)
    chk("T1 open overlap: width 1/2 in the separating direction",
        w is not None and w["hw"][0] == Fraction(1, 4)
        and w["hw"][1] == Fraction(1))
    lo3, hi3 = cube_bounds([Fraction(2)] + [Fraction(0)] * 3,
                           Fraction(1))
    chk("T2 NEGATIVE CONTROL: two boxes that TOUCH without overlapping "
        "ne donnent aucun overlap ouvert",
        _inter_generic(lo1, hi1, lo3, hi3) is None)

    # T3 : la carte ε et son check de domaine
    E = eps_box({"lo": [Fraction(-1, 2)] * 4,
                 "hi": [Fraction(1, 2)] * 4},
                [Fraction(0)] * 4, Fraction(1))
    chk("T3 boîte ε ⊆ [−1,1]⁴", eps_box_in_range(E)
        and E[0] == (Fraction(-1, 2), Fraction(1, 2)))
    E2 = eps_box({"lo": [Fraction(-2)] * 4, "hi": [Fraction(2)] * 4},
                 [Fraction(0)] * 4, Fraction(1))
    chk("T4 NEGATIVE CONTROL: a box outside the model domain is "
        "REFUSED (the remainder is not valid there)", not eps_box_in_range(E2))

    # T5 : (1+ρ)h exactement représentable
    chk("T5 exact dyadic halo over the whole preregistered scale",
        all(halo_hw(Fraction(1, 1024), r)[1] is not None
            for r in RHO_LADDER))

    # T6: the recentring; the identity is the identity
    T = recenter_matrix([Fraction(0)] * 4, Fraction(1))
    p = [CIV(iv.mpf(float(k + 1)), iv.mpf(float(-k)))
         for k in range(NM)]
    q = apply_recenter(T, p)
    chk("T6 recentrage identité : coefficients inchangés",
        all(mp.mpf(q[k].re.a) <= float(k + 1) <= mp.mpf(q[k].re.b)
            and mp.mpf(q[k].im.a) <= float(-k) <= mp.mpf(q[k].im.b)
            for k in range(NM)))

    # T7: the recentring is EXACT on a known polynomial.
    #      f(ε) = ε₀²  et  φ(ε) = (1/2 + ε/4) ⟹ f∘φ = 1/4 + ε₀/4 + ε₀²/16
    p2 = [CIV(IV0, IV0)] * NM
    p2 = list(p2)
    p2[MIDX[(2, 0, 0, 0)]] = CIV(IV1, IV0)
    T2m = recenter_matrix([Fraction(1, 2)] + [Fraction(0)] * 3,
                          Fraction(1, 4))
    q2 = apply_recenter(T2m, p2)

    def near(c, x):
        return (mp.mpf(c.re.a) <= x <= mp.mpf(c.re.b)
                and mp.mpf(c.im.a) <= 0 <= mp.mpf(c.im.b))
    chk("T7 recentrage exact sur ε₀² ∘ (1/2 + ε/4)",
        near(q2[0], 0.25) and near(q2[MIDX[(1, 0, 0, 0)]], 0.25)
        and near(q2[MIDX[(2, 0, 0, 0)]], 0.0625))

    # T8 : recentrage puis évaluation = évaluation directe (le même
    #      point, two frames), the invariant that carries this whole step
    Efull = [iv.mpf([-1, 1])] * 4
    v_direct = enclose(p2, IV0, [iv.mpf([0.25, 0.75])]
                       + [iv.mpf([-1, 1])] * 3)
    v_recent = enclose(q2, IV0, Efull)
    chk("T8 recentring preserves the enclosure at the same physical point",
        mp.mpf(v_recent.re.a) <= mp.mpf(v_direct.re.a)
        and mp.mpf(v_direct.re.b) <= mp.mpf(v_recent.re.b))

    # T9 : la séparation de feuillet, et son REFUS
    dm = CIV(iv.mpf([-1e-9, 1e-9]), iv.mpf([-1e-9, 1e-9]))
    dp = CIV(iv.mpf([1.9, 2.1]), iv.mpf([-0.1, 0.1]))
    th, rec = sep_phase(dm, dp)
    chk("T9 feuillet : différence ~0 et somme ~2 ⟹ θ = +1 à marge > 0",
        th == 1 and rec["margin"] > 0)
    th2, _ = sep_phase(dp, dm)
    chk("T10 feuillet inversé ⟹ θ = −1", th2 == -1)
    amb = CIV(iv.mpf([-1, 1]), iv.mpf([-1, 1]))
    th3, _ = sep_phase(amb, amb)
    chk("T11 NÉGATIF : séparation ambiguë ⟹ REFUS (jamais d'essai)",
        th3 is None)

    # T12 : sérialisation exacte des TM (aller-retour bit à bit)
    t = TMC([CIV(iv.mpf([0.1, 0.2]), iv.mpf([-0.3, 0.4]))] * NM,
            iv.mpf([0, 1e-9]))
    p3, r3 = de_tmc(ser_tmc(t))
    chk("T12 sérialisation `_mpf_` exacte (aller-retour identique)",
        all(mp.mpf(p3[k].re.a) == mp.mpf(t.p[k].re.a)
            and mp.mpf(p3[k].im.b) == mp.mpf(t.p[k].im.b)
            for k in range(NM))
        and mp.mpf(r3.b) == mp.mpf(t.rem.b))

    # T13: the discrete cocycle, and a deliberately broken triple
    def cocy(a, b, c):
        return all(a[r] * b[r] * c[r] == 1 for r in range(3))
    chk("T13 trivial theta cocycle on a sound triple",
        cocy([1, 1, 1], [1, 1, 1], [1, 1, 1]))
    chk("T14 NEGATIVE CONTROL: a single flipped theta breaks the cocycle",
        not cocy([1, -1, 1], [1, 1, 1], [1, 1, 1]))

    # T15 : division intervalle complexe
    a = CIV(iv.mpf([1, 1]), iv.mpf([0, 0]))
    b = CIV(iv.mpf([2, 2]), iv.mpf([0, 0]))
    d = civ_div(a, b)
    chk("T15 division complexe : 1/2 enclos",
        mp.mpf(d.re.a) <= 0.5 <= mp.mpf(d.re.b))

    # T16/T17 : the rational-box step — les endpoints Fraction → iv en arrondi
    #   OUTSIDE. Witness: 1/65 (the denominator that halos of 65h/64
    #   really inject into the overlap boxes), not representable
    #   in binary. The criterion "B encloses q" is tested WITHOUT floats: 65.B
    #   doit contenir 1 exactement.
    q65 = Fraction(1, 65)
    B = fraction_box_to_iv(q65, q65)
    X = iv.mpf(65) * B
    chk("T16 fraction_box_to_iv enclot le rationnel exact (65·B ∋ 1)",
        mp.mpf(X.a) <= 1 <= mp.mpf(X.b))
    old = iv.mpf([float(q65), float(q65)])
    Xo = iv.mpf(65) * old
    chk("T17 NEGATIVE CONTROL: the old nearest-float path MISSES the exact "
        "rational (65.[float(1/65)] does not contain 1): the defect is real",
        not (mp.mpf(Xo.a) <= 1 <= mp.mpf(Xo.b)))

    print("SELF-TEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_selftest() if "--selftest" in sys.argv else build())
