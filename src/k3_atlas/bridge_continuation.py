#!/usr/bin/env python3
"""
bridge_continuation.py — the bridge step: THE BRIDGE CHARTS AND THEIR TWO
TRANSITIONS.

WHAT THIS SCRIPT PAYS (v2): the contract of a review, corrected for the
factor 2 established by the preliminary computation and under the honest
name fixed by the mirror identification, AND the repair contract of a
second review. The v1 of this script HAD TWO LOAD-BEARING DEFECTS, BROKEN
IN REVIEW AND REPAIRED HERE: (i) it read the theta pattern on the upper
side as a "verified prediction", that is, it renamed a gluing FAILURE as a
green check; v2 turns it into a DECK THEOREM
(`D = diag(+1,-1,+1,+1,+1,-1)`) and glues the CONTINUED sheet `D.Z_conj`,
which IS the upper transition (the conjugate alone does not glue and stays
a diagnostic); (ii) its "nerve" counted as edges intersections without a
certified transition; v2 publishes that under its real name
(`domain_intersection_graph`) and rebuilds the NERVE on CERTIFIED edges
only. The scope is STRATIFIED (36 of 64 fully open in 4D), not "ambient"
as a block.

THE ARCHITECTURE, AND WHY IT IS NOT THE ONE OF THE SCOPING NOTE.
`U-` (Im u, Im v < 0) and `U+` (Im u, Im v > 0) have NO open overlap:
their closures meet only on the corner of real codimension 2
`{Im u = Im v = 0}`. Comparing a lower section and an upper section
directly "on the face" would be an argument by boundary value. We
therefore introduce a BRIDGE CHART `B` on the bilateral box, and certify
TWO separate transitions:

    B = Z-  on  bridge ∩ halo-        B = Z+  on  bridge ∩ halo+

These two intersections have STRICTLY POSITIVE width in all four
coordinates: they are genuine open overlaps, and the anchors there are
taken STRICTLY inside (Im < 0 on one side, Im > 0 on the other), never on
the face.

THE BILATERAL CONSTRUCTOR, AND WHY THE EARLIER ONE IS NOT REUSABLE.
`build_section` does "plain, otherwise rotated with the component sigma".
On a bilateral box `Im R` changes sign, so the component does not exist.
`build_section_bilateral` therefore chooses BY THE CERTIFIED SIGN OF
`Re R`, without any trial:

    Re R strictly positive on the whole box   ->  principal branch
    Re R strictly negative on the whole box   ->  canonical `w = i.sqrt_p(-R)`
                                        (component PINNED to +1: this is the
                                         definition of the regime of the
                                         residual closure, keeping the cut
                                         guard verbatim on -R)
    otherwise                                 ->  REFUSAL

The sign is read on the model ENCLOSURE of `Re R` (remainder included), not
on the rational bound of the scouting step: the decision is taken in the
arithmetic that builds the section. The exact rational bound is published
alongside as a CROSS-CHECK: two independent paths to the same sign.

THE COMPARISON, AND THE TRAP IT AVOIDS. Comparing `B` and `Z` through
their ENCLOSURES compares nothing: each enclosure has the width of the
VARIATION of the function on the overlap, so their difference does too. We
therefore RECENTRE both Taylor models in the frame of the overlap, an
affine substitution `eps = off + scale.eps'` exact in Fraction and
ANISOTROPIC (the bridge is 2H in the imaginary direction and H in the real
one, so not a cube), then subtract COEFFICIENT BY COEFFICIENT. This is the
same countermeasure as before, here with a recentring matrix whose scale
is per coordinate.

THE PHASE IS AN OPEN RESULT, THE CHECK IS PREREGISTERED. `theta = +1`
means that the bridge and the side describe the SAME projective point (the
source gauge being normalised to `Z_g = 1`, the projective scalar is
fixed: independent signs per row are NOT a projective compatibility). The
test is open as to its outcome, preregistered as to the meaning of
success.

AND THE RESULT IS NOT THE EXPECTED ONE, WHICH IS THE FINDING. The bridge
glues exactly to the LOWER side. Against the UPPER side derived by
conjugation, `theta` is MIXED: `+1` on the principal row, `-1` on the two
canonical rows. This is neither a bug nor an ambiguity (the margins are of
order 1, minimum 1.04), and it is NOT a projective renormalisation. The
reason is structural and can be predicted:

    on the corner `Im u = Im v = 0` the radicand `R` is REAL;
    where `Re R > 0` (principal regime) the root is REAL and
      the antiholomorphic involution FIXES it;
    where `Re R < 0` (canonical regime) the root `i.sqrt(-R)` is
      PURELY IMAGINARY and the involution NEGATES it.

So the conjugate point and the point continued analytically across the
corner are two DISTINCT points of projective space above the same base
point, exchanged by the real involution. The consequence for the first
lever of the scoping note: the conjugate neighbour is a legitimate atlas,
but it is NOT the one that analytic continuation reaches. The `theta`
pattern is PREDICTED from the regime alone, and the prediction is verified
on 64 of 64: the surprise becomes a falsifiable result instead of staying
a failed check.

EXACTNESS. The separation `sep_phase` is POINTWISE over the whole overlap
(it bears on correlated enclosures after recentring): at each point `B`
and `Z` are two analytic roots of the SAME `R != 0`, and `theta = -1`
would give `|diff| = 2|Z| > 0 = |sum|`, which is excluded. The lemma then
applies VERBATIM on each overlap, but its checks are REBUILT here, not
imported: `F` nonzero on the bridge AND on the side, nonzero gauge,
interior anchor, strict separation.

WHAT THIS SCRIPT DOES NOT PAY: the metric of the bridge (Qmat, Weyl,
lateral congruences); the canonical identification of the neighbour, which
REFUSED (see the mirror record artefact): the upper atlas is DERIVED by
conjugation, not enumerated, and this script inherits that status; the
codimension-1 neighbours; gluing across the Re faces, where the bridges
stay RELATIVE charts; the 895 other pairs.

CHECKS
  F2a  geometry: half-width 2H in the two reflected imaginary
       directions, H in the real ones, re-derived here in Fraction, and
       IDENTICAL to the one of the preliminary computation (import
       verified);
  F2b  the union of the two cores is inside the bridge, STRICTLY in the
       imaginary directions;
  F2c  regime by CERTIFIED sign of `Re R` on the model enclosure, without
       trial, on 64 x 3 rows, and agreement with the exact rational
       bound of the scouting step (two arithmetics, one verdict);
  F2d  COMPLETE bilateral section: the six coordinates exist, the gauge of
       the target chart is strictly bounded away from zero, and the bridge record
       (eps, regime) is FROZEN and serialised;
  F2e  NEGATIVE CONTROLS: (i) the false bridge of half-width H is REFUSED
       by the inclusion; (ii) resorting to the component rule on the
       bridge is REFUSED (Im R straddles); (iii) a box where `Re R`
       straddles is REFUSED instead of falling back on a regime;
  F2f  64/64 with no silent filtering;
  F3a  both overlaps have width > 0 in all 4 coordinates, in exact
       rationals, and both anchors are STRICTLY interior (Im < 0 and
       Im > 0), never on the face;
  F2d(bis) THE BRIDGE RECORD IS DERIVED: it is measured against the LOWER
       side (the established atlas, anchor of truth), then the section is
       REBUILT and RE-VERIFIED. Defaulting the sheet record would amount
       to calling a default choice of sheet a "gluing";
  F3b  the bridge glues EXACTLY to the lower side: `theta = +1` on ALL
       rows, by exact ANISOTROPIC recentring plus strict separation,
       with refusal if ambiguous;
  R1   `D = diag(+1,-1,+1,+1,+1,-1)` is a DECK TRANSFORMATION:
       involutive, NON scalar, preserving the three quadrics, and
       the identity `Z_conj = D.Z_bridge` is certified AT THE COEFFICIENT
       level on 64 of 64. Algebra alone distinguishes nothing (any
       diagonal of signs preserves the quadrics): the discriminating
       power is carried by the negative control R1e;
  R2   the bridge glues EXACTLY (`theta = +1` on all six coordinates) to
       the CONTINUED sheet `D.Z_conj`, and THAT is the upper transition;
       the conjugate ALONE does NOT glue (R2b) and stays published as the
       diagnostic of another sheet;
  R3   the NERVE counts only edges whose transition is CERTIFIED (380
       nodes; 5396 leaf-to-leaf imported from the atlas step, 64
       bridge-to-leaf, 210 of 210 bridge-to-bridge; connected). The box
       intersection graph is published separately as
       `domain_intersection_graph`, which is NOT a nerve;
  F3d  a PROJECTIVELY DISCRIMINATING negative control: a NON-scalar
       mutation (a single canonical row negated, gauge and affine
       coordinates unchanged) so the gluing falls 64 of 64; the global
       negation `Z -> -Z` is published as the SAME point of `P^5`, since
       it discriminated nothing, and that was the defect of v1;
  R1e  a single sign changed in `D` breaks the deck identity;
  R4   the upstream chain is verified (preliminary computation, the
       mirror record, the atlas step, the residual closure).

Output : results/bridge_continuation.json
Usage  : bridge_continuation.py [--selftest]
Env    : K3_F2F3_WORKERS (default 6)
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
COVER_JSON = RES / "dyadic_cover.json"
ATLAS_JSON = RES / "atlas_assembly.json"
C127E_JSON = RES / "residual_closure.json"
PRELIMINARY_JSON = RES / "bridge_preliminary.json"
F1_JSON = RES / "mirror_record.json"
ART = RES / "bridge_continuation.json"
N_WORKERS = int(os.environ.get("K3_F2F3_WORKERS", "6"))

# --- PREREGISTERED, frozen before the run ---------------------------------
IM_DIRS = (1, 3)
NG = 4
# The gluing SUCCEEDS if and only if theta = +1 on both sides. -1 is a
# publishable result and a gluing FAILURE: preregistered here, not
# adjudicated after the fact.
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
#  Exact geometry: never a float inside a decision
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
    """Intersection box, or None if the width is not STRICTLY
    positive in all four coordinates. Everything in rationals: "open"
    is a rational comparison, never a float test."""
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
    """`inner` inside `outer`, STRICTLY in `strict_dirs`. The distinction
    is content: the bridge is open in the reflected imaginary
    directions (this is the promotion of the corner to an ambient neighbourhood), it
    stays FLUSH in the real directions, where it inherits the halo
    clipped: those faces belong to the codimension-1 neighbours and
    are not paid here."""
    marg = [(inner[k][0] - outer[k][0], outer[k][1] - inner[k][1])
            for k in range(4)]
    ok = (all(a >= 0 and b >= 0 for a, b in marg)
          and all(marg[k][0] > 0 and marg[k][1] > 0 for k in strict_dirs))
    return ok, marg


# ===========================================================================
#  Exact ANISOTROPIC recentring: the generalisation the bridge needs
# ===========================================================================
#  The earlier step recentred between two cubes, with a single `scale`. The bridge is
#  2H in the imaginary directions and H in the real ones, so the scale is PER COORDINATE. The
#  substitution stays affine, so it creates NO remainder (composing a
#  polynomial of degree at most N with an affine map gives a polynomial of
#  degree at most N), and the matrix is exact in rationals before conversion.
# ===========================================================================
def recenter_matrix_aniso(off, scale):
    """T such that (p composed with phi)[beta] = sum_{alpha >= beta} T[beta][alpha].p[alpha], for
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
    """EXACT (off, scale) sending the frame of `dst_box` into that of
    `src_box`: if `x = c_s + h_s.eps_s = c_d + h_d.eps_d`, then
    `eps_s = (c_d - c_s)/h_s + (h_d/h_s).eps_d`."""
    cs, hs = center_hw(src_box)
    cd, hd = center_hw(dst_box)
    return ([(cd[k] - cs[k]) / hs[k] for k in range(4)],
            [hd[k] / hs[k] for k in range(4)])


def frame_admissible(off, scale):
    """The model is valid only on the unit cube of symbols: the image of the target frame
    must fit there. This is a CHECK, not an assumption."""
    return all(-1 <= off[k] - scale[k] and off[k] + scale[k] <= 1
               for k in range(4))


# ===========================================================================
#  The BILATERAL constructor
# ===========================================================================
def uv_tm(center, hw):
    """`u` and `v` as complex models on a box with PER-COORDINATE half-widths."""
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
    """The source section on a BILATERAL box, regime chosen by the
    CERTIFIED SIGN of `Re R`, never by trial and never by fallback.

    `force_sigma` is the MUTATION of negative control (ii): it attempts the
    rotated determination with a component rule where `Im R`
    straddles. Code doing that on a bridge would accept a
    continuation whose sheet is not defined, which is what the
    negative control must exhibit, and it must BREAK."""
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
            # MUTATION: we pretend to read sigma as the sign of Im R.
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
                    # component PINNED to +1: this is the DEFINITION of the
                    # canonical regime, not a trial. The guard
                    # applies verbatim to -R.
                    Zs = tm_sqrt_rotated(R, 1).mul_real(riv(int(eps[r])))
                rec["regime"] = sel
            except BranchCutError as exc:
                # The sign was certified: if the guard refuses
                # anyway, that is a REFUSAL, never a fallback on the other
                # branche.
                rec["refused"] = exc.diag.get("guard")
                rows.append(rec)
                continue
        Z[s] = Zs
        iZ = Zs.inv()
        dZ[s] = (u.mul_real(A[r][1]) * iZ, v.mul_real(A[r][2]) * iZ)
        # Compatibility alias with `transport_hardened`, which serialises
        # `source_determinations`. The regime of the BRIDGE legitimately differs
        # from that of the core (canonical where the core was rotated): the
        # difference must be VISIBLE in the artefact, not masked. The
        # metric check bears on the TARGET kinds, not on the source
        # regime, precisely because the latter changed on purpose.
        rec["determination"] = rec["regime"]
        rr, ri = R.to_iv_pair()
        rec["radicand_absmin"] = _f_down(
            mp.mpf(civ_absmin(CIV(rr, ri)).a))
        rows.append(rec)
    return Z, dZ, rows


# ===========================================================================
#  Comparing two sections on an overlap: recentring, then
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
#  R1 — THE DECK TRANSFORMATION
# ===========================================================================
#  The theta pattern measured against the conjugate defines a diagonal of signs
#
#      D = diag(+1, −1, +1, +1, +1, −1)
#
#  which preserves the three ambient quadrics `Q_m(Z) = sum_a mu_a^m Z_a^2`
#  (the coordinates enter only through their SQUARE), is its own
#  inverse, and is NOT a projective scalar. It is therefore a
#  non-trivial holomorphic automorphism of the complete intersection: a
#  DECK TRANSFORMATION of the description by square roots.
#
#  HONESTY ABOUT WHAT IS TRIVIAL. "D preserves the quadrics" is
#  true of ANY diagonal of signs; that is not what singles out THIS D.
#  The content lies elsewhere, and it is measured: it is THAT D, and no
#  other, that links the conjugate section to the bridge section. The
#  negative control (a single sign changed) therefore carries all the discriminating power;
#  preserving the quadrics says only that `D.Z` is still a
#  point of the surface, which is necessary for the repair to have
#  un sens.
# ===========================================================================
DECK_D = (1, -1, 1, 1, 1, -1)


def apply_deck(Z, D, keys):
    """`D.Z` at the level of Taylor models: one sign per coordinate."""
    out = {}
    for k in keys:
        z = Z[k]
        if z is None or D[k] == 1:
            out[k] = z
        else:
            out[k] = TMC([-x for x in z.p], z.rem)
    return out


def deck_algebra(D):
    """The three structural facts about `D`, EXACT (integers and rationals):
    involution, non-scalar in the projective group, and preservation of the three
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
    """Value of the coordinates at a POINT (anchor), through the symbol frame of the model.
    Membership of the anchor in the box is an EXACT rational check,
    checked by the caller."""
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
    """Strict separation of |a-b|^2 against |a+b|^2 at a point."""
    dm = CIV(a.re - b.re, a.im - b.im)
    dp = CIV(a.re + b.re, a.im + b.im)
    return sep_phase(dm, dp)


# ===========================================================================
#  Job par pont
# ===========================================================================
_G = {}


def _init(cell, leaves, halos, bridges, sheet_records=None):
    _G.update(cell=cell, leaves=leaves, halos=halos, bridges=bridges,
              sheet_records=sheet_records or {})


def _bridge_job(i):
    S, g, eps = _G["cell"]
    leaf = _G["leaves"][i]
    rec_h = _G["halos"][i]
    Wb = _G["bridges"][i]
    out = {"tile": i}

    # --- geometry -----------------------------------------------------
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
    # NEGATIVE CONTROL (i): the false bridge of the scoping note
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

    # NEGATIVE CONTROL (ii): forcing the component rule on the bridge
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

    # gauge of the TARGET chart strictly bounded away from zero on the bridge
    S2, g2 = tuple(leaf["chart"]["S"]), leaf["chart"]["g"]
    gr, gi = Zb[g2].to_iv_pair()
    gmin = mp.mpf(civ_absmin(CIV(gr, gi)).a)
    out["F2d_target_gauge_absmin"] = _f_down(gmin)
    out["F2d_gauge_positive"] = bool(gmin > 0)

    # --- the two side sections ----------------------------------------
    chf = [float(x) for x in c_h]
    Hf = float(H)
    Zlo, _d1, rows_lo = build_section(S, g, eps, chf, Hf)
    cup = [(-chf[k] if k in IM_DIRS else chf[k]) for k in range(4)]
    Zup, _d2, rows_up = build_section(S, g, eps, cup, Hf)
    out["lower_dets"] = [r["determination"] for r in rows_lo]
    out["upper_dets"] = [r["determination"] for r in rows_up]
    # R4c: the RECONSTRUCTED lower side must be the object ACTUALLY
    # certified by the atlas step, not a new reconstruction that
    # resembles it. Without this check, a future drift of the code would derive the
    # bridge against something other than the established atlas.
    out["R4c_lower_dets_match_c127d"] = bool(
        out["lower_dets"] == list(rec_h["source_determinations"]))

    keys = list(S) + [x for x in range(6) if x not in S]

    # --- THE BRIDGE RECORD IS DERIVED, IT IS NOT DEFAULTED ------------
    # The LOWER atlas is the established object: it is
    # the anchor of truth, and the bridge record must agree with it. We
    # DERIVE it (measuring theta row by row against the lower side),
    # then RECONSTRUCT and RE-CHECK. This is not a trial: the
    # re-check must give theta = +1 on ALL rows, otherwise
    # REFUSAL. Not deriving would amount to calling "gluing" a mere
    # default choice of sheet.
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
        out["refused"] = "sheet_record_derivation_ambiguous"
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

    # --- F3a: open overlaps plus STRICTLY interior anchors -------------
    res = {}
    # R2: the CONTINUED upper sheet. `Z_upper_conj` is the atlas
    # derived by conjugation; it lives on ANOTHER sheet, linked to
    # that of the bridge by the deck transformation `D`. The sheet that
    # analytic continuation reaches is therefore `D.Z_upper_conj`, and
    # THAT is what must pass the real upper contract. The conjugate atlas
    # is kept as a separate DIAGNOSTIC, not as an atlas edge.
    Zup_cont = apply_deck(Zup, DECK_D, keys)
    for side, halo, Zside, box_side in (
            ("lower", halo_lo, Zlo, halo_lo),
            ("upper_conj", halo_up, Zup, mirror_bounds(halo_lo)),
            ("upper_cont", halo_up, Zup_cont, mirror_bounds(halo_lo))):
        W = inter(Wb, halo)
        if W is None:
            res[side] = {"refused": "overlap_not_open"}
            continue
        # anchor: centre of the overlap, then SHIFTED strictly to the right
        # side of the face if the overlap touches it. Never on Im = 0.
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
        # separation AT THE ANCHOR POINT (as the review contract asks)
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
        # nonzero radicands ON THE OVERLAP, on BOTH sides (the lemma)
        r["glued_exactly"] = bool(
            r.get("theta") == THETA_REQUIRED
            and r["theta_consistent_across_lines"]
            and r["anchor_separation_strict"]
            and r["anchor_agrees_with_box"]
            and r["anchor_strictly_inside_overlap"]
            and r["anchor_off_the_face"])
        # THE REAL-STRUCTURE PREDICTION. On the corner `Im u = Im v = 0`
        # the radicand `R` is REAL. Where `Re R > 0` (principal regime)
        # the root is REAL and the antiholomorphic involution FIXES it;
        # where `Re R < 0` (canonical regime) the root `i.sqrt(-R)` is
        # PURELY IMAGINARY and the involution NEGATES it. The upper side
        # being DERIVED by conjugation, its phase against the bridge is therefore
        # PREDICTED row by row, and predicted MIXED, which is not
        # a projective renormalisation. This prediction is
        # falsifiable, and that is what turns the surprise into a result.
        pred = {}
        for j, k in enumerate(S):
            rg = rows_b[j]["regime"]
            pred[str(k)] = 1 if rg == "principal" else -1
        for k in keys:
            if k not in S:
                pred[str(k)] = 1        # gauge and affine coordinates:
                                        # real on the corner, hence fixed
        r["real_structure_prediction"] = pred
        r["deck_D_from_prediction"] = [pred[str(k)] for k in range(6)]
        r["matches_real_structure_prediction"] = bool(
            side != "upper_conj"
            or all(th[k].get("theta") == pred[str(k)] for k in keys))
        # R1d: the deck identity AT THE COEFFICIENT level. The theta pattern against
        # the conjugate must be EXACTLY `D`, row by row, and the measured
        # `D` must be the preregistered one. This is the only
        # discriminating part of the deck theorem.
        if side == "upper_conj":
            r["deck_D_measured_equals_preregistered"] = bool(
                [pred[str(k)] for k in range(6)] == list(DECK_D)
                and all(th[k].get("theta") == DECK_D[k] for k in keys))
        res[side] = r
    out["sides"] = res

    # --- F3d: THE NEGATIVE CONTROL MUST BE PROJECTIVELY DISCRIMINATING -
    # v1 negated ALL coordinates: `Z` and `-Z` are THE SAME POINT of
    # projective space, so the test showed only the sensitivity of the
    # affine normalisation `Z_g = 1` to a change of representative, not
    # detecting a wrong POINT. v2 keeps the gauge and the affine
    # coordinates UNCHANGED and negates A SINGLE canonical root row:
    # `D_bad = diag(1,-1,1,1,1,1)`, non scalar, hence a REALLY different
    # projective point.
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
        # the gluing must fall: at least one row at -1 AND a
        # NON constant pattern (so not a mere change of representative)
        vals = [thn[k].get("theta") for k in keys]
        out["F3d_non_scalar_mutation_breaks_gluing"] = bool(
            any(v != THETA_REQUIRED for v in vals) and len(set(vals)) > 1)
        # and the control: the GLOBAL negation is the SAME point
        # projective point, published to say why it is worthless
        Zglob = apply_deck(Zb, tuple([-1] * 6), keys)
        thg = theta_lines(Zglob, Wb, {k: Zlo[k] for k in keys},
                          halo_lo, Wl, keys)
        gv = {thg[k].get("theta") for k in keys}
        out["F3d_global_negation_theta"] = (
            list(gv)[0] if len(gv) == 1 else None)
        out["F3d_global_negation_is_same_projective_point"] = bool(
            len(gv) == 1)

    # --- R1e: DECK NEGATIVE CONTROL, a single sign changed must BREAK --
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
#  R3: the BRIDGE to BRIDGE transitions, certified (not merely
#  geometric). A NERVE edge is a CERTIFIED transition; an
#  intersection of boxes is only a graph of domains.
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
    ei, ej = _G["sheet_records"][i], _G["sheet_records"][j]
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
    chk("mirror_bounds touches only the imaginary directions",
        mirror_bounds([(F(-2), F(0))] * 4)[0] == (F(-2), F(0))
        and mirror_bounds([(F(-2), F(0))] * 4)[1] == (F(0), F(2)))
    chk("bridge intersected with the lower halo is OPEN in all 4 coordinates",
        inter([(F(-2), F(2))] * 4, [(F(-2), F(0))] * 4) is not None)
    chk("lower halo intersected with upper halo is NOT open (the corner)",
        inter([(F(-2), F(0))] * 4, [(F(0), F(2))] * 4) is None)

    # reframe: identity, translation, anisotropic dilation
    b = [(F(-1), F(1)), (F(-2), F(2)), (F(-1), F(1)), (F(-2), F(2))]
    off, sc = reframe(b, b)
    chk("reframing a box onto itself is the identity",
        off == [F(0)] * 4 and sc == [F(1)] * 4)
    w = [(F(-1), F(0)), (F(-2), F(0)), (F(-1), F(1)), (F(-2), F(2))]
    off, sc = reframe(b, w)
    chk("anisotropic reframe: off and scale per coordinate",
        off == [F(-1, 2), F(-1, 2), F(0), F(0)]
        and sc == [F(1, 2), F(1, 2), F(1), F(1)])
    chk("frame_admissible accepts a subbox",
        frame_admissible(off, sc))
    chk("frame_admissible REFUSES a box that overflows",
        not frame_admissible([F(1), F(0), F(0), F(0)], [F(1)] * 4))

    # The anisotropic recentring, on a polynomial whose value is known
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
    chk("anisotropic recentring: constant term",
        abs(float(mp.mpf(q[0].re.a)) - (3 + 5 * 0.5 + 7 * -0.25)) < 1e-30)
    chk("anisotropic recentring: the slope in the first symbol is scaled",
        abs(float(mp.mpf(q[e[0]].re.a)) - 5 * 0.5) < 1e-30)
    chk("anisotropic recentring: the slope in the second symbol is scaled",
        abs(float(mp.mpf(q[e[1]].re.a)) - 7 * 0.25) < 1e-30)
    chk("recentring onto itself leaves the polynomial UNCHANGED",
        all(abs(float(mp.mpf(x.re.a)) - float(mp.mpf(y.re.a))) < 1e-30
            for x, y in zip(
                apply_recenter(
                    recenter_matrix_aniso([F(0)] * 4, [F(1)] * 4), p), p)))

    # sep_phase, on constructed values
    th, _ = sep_phase(CIV(iv.mpf(0), IV0), CIV(iv.mpf(4), IV0))
    chk("sep_phase reads theta = +1 when the difference vanishes", th == 1)
    th, _ = sep_phase(CIV(iv.mpf(4), IV0), CIV(iv.mpf(0), IV0))
    chk("sep_phase reads theta = -1 when it is the sum that vanishes",
        th == -1)
    th, _ = sep_phase(CIV(iv.mpf([-1, 1]), IV0), CIV(iv.mpf([-1, 1]), IV0))
    chk("sep_phase REFUSE l'ambigu", th is None)

    # The bilateral constructor on the real cell: regimes assigned,
    # and the sigma mutation must break.
    cov = json.loads(COVER_JSON.read_text(encoding="utf-8"))
    cell = cov["cell"]
    S, g, eps = tuple(cell["S"]), cell["g"], tuple(cell["eps"])
    sc_j = json.loads(PRELIMINARY_JSON.read_text(encoding="utf-8"))
    t0 = sc_j["per_tile"][0]
    Wb = [(F(*b[0]), F(*b[1])) for b in t0["bridge_corrected_bounds"]]
    cb, hb = center_hw(Wb)
    Z, _dz, rows = build_section_bilateral(
        S, g, eps, [float(x) for x in cb], [float(x) for x in hb])
    chk("build_section_bilateral assigns the 3 regimes without trial",
        all(r["regime"] is not None for r in rows))
    chk("the regimes agree with the scouting bound (two arithmetics)",
        [r["regime"] for r in rows] == t0["regimes"])
    _z2, _d2, rows_f = build_section_bilateral(
        S, g, eps, [float(x) for x in cb], [float(x) for x in hb],
        force_sigma=True)
    chk("the component-rule MUTATION is REFUSED on the bridge",
        any(r.get("refused") == "component_sigma_undetermined_on_bridge"
            for r in rows_f))
    chk("Im R straddles on the rows where the component is refused",
        all(r["im_R_straddles"] for r in rows_f
            if r.get("refused") == "component_sigma_undetermined_on_bridge"))
    # NEGATIVE CONTROL (iii): an enlarged box where Re R straddles must REFUSE
    big = [(cb[k] - 4, cb[k] + 4) for k in range(4)]
    cbg, hbg = center_hw(big)
    _z3, _d3, rows_big = build_section_bilateral(
        S, g, eps, [float(x) for x in cbg], [float(x) for x in hbg])
    chk("a box where Re R straddles is REFUSED, not fallen back on",
        any(r.get("refused") == "re_R_straddles_zero_no_regime"
            for r in rows_big))
    print(f"\nself-test {ok}/{tot}")
    return ok == tot


# ===========================================================================
#  Build
# ===========================================================================
def build():
    print("=" * 78)
    print(f"BRIDGE CHARTS AND THEIR TWO TRANSITIONS "
          f"({N_WORKERS} workers)")
    print("=" * 78)
    cov = json.loads(COVER_JSON.read_text(encoding="utf-8"))
    atl = json.loads(ATLAS_JSON.read_text(encoding="utf-8"))
    residual = json.loads(C127E_JSON.read_text(encoding="utf-8"))
    preliminary = json.loads(PRELIMINARY_JSON.read_text(encoding="utf-8"))
    f1 = json.loads(F1_JSON.read_text(encoding="utf-8"))
    cell = cov["cell"]
    S, g, eps = tuple(cell["S"]), cell["g"], tuple(cell["eps"])
    leaves = list(cov["tiles"]) + list(residual["new_tiles"])
    halos = {h["index"]: h["record"] for h in atl["halos"] if h["ok"]}
    clipped = sorted(i for i, r in halos.items()
                     if r["rule"] == "clipped")
    log(f"cell S={S} g={g}; {len(clipped)} clipped tiles")
    log(f"    the canonical identification REFUSED; the upper atlas "
        f"is DERIVED, not enumerated "
        f"(claim_level={f1['claim_level']})")

    # --- R4c: THE UPSTREAM CHAIN IS CHECKED, NOT MERELY IMPORTED ------
    # v1 read the upstream JSON without ever checking that they were GREEN
    # nor in full mode. A pilot or red artefact would have been consumed
    # in silence.
    up = {}
    for name, blob, need_full in (
            ("c127d_atlas", atl, True),
            ("c129f_bridge_scout", preliminary, False),
            ("c129f_f1_mirror_ledger", f1, False)):
        gp, gt = blob.get("checks_passed"), blob.get("checks_total")
        up[name] = {"checks": f"{gp}/{gt}", "green": bool(gp == gt and gt),
                    "mode": blob.get("mode")}
        if need_full:
            up[name]["full"] = bool(blob.get("mode") == "full")
    for name, path in (("c129d_exact_gluing",
                        RES / "exact_gluing.json"),
                       ("c129e_halo_metric",
                        RES / "halo_metric.json")):
        try:
            b = json.loads(path.read_text(encoding="utf-8"))
        except OSError:
            up[name] = {"green": False, "missing": True}
            continue
        gp, gt = b.get("checks_passed"), b.get("checks_total")
        up[name] = {"checks": f"{gp}/{gt}",
                    "green": bool(gp == gt and gt),
                    "mode": b.get("mode")}
    upstream_ok = all(v.get("green") and v.get("full", True)
                      for v in up.values())
    log(f"R4c: upstream chain — " + " ; ".join(
        f"{k} {v.get('checks', '?')}" for k, v in up.items())
        + f" ⟹ {upstream_ok}")

    # --- F2a: the geometry, RE-DERIVED then confronted with the preliminary computation
    bridges, f2a = {}, True
    for i in clipped:
        r = halos[i]
        c_h = [Fraction(float.fromhex(x)) for x in r["center_hex"]]
        H = Fraction(float.fromhex(r["H_hex"]))
        hb = bounds(c_h, H)
        W = [((-2 * H, 2 * H) if k in IM_DIRS else hb[k])
             for k in range(4)]
        bridges[i] = W
    preliminary_boxes = {t["tile"]: [(Fraction(*b[0]), Fraction(*b[1]))
                          for b in t["bridge_corrected_bounds"]]
              for t in preliminary["per_tile"]}
    f2a = all(bridges[i] == preliminary_boxes[i] for i in clipped)
    log(f"F2a: 2H/H geometry re-derived, matching the preliminary bound on "
        f"{len(clipped)} ponts : {f2a}")

    jobs = clipped
    mpctx = get_context("fork")
    _init((S, g, eps), leaves, halos, bridges)
    with mpctx.Pool(N_WORKERS, initializer=_init,
                    initargs=((S, g, eps), leaves, halos, bridges)) as pool:
        rows = pool.map(_bridge_job, jobs)
    byi = {r["tile"]: r for r in rows}
    log(f"    {len(rows)} bridges processed")

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
    log(f"F2b: core- union core+ inside the bridge (strict in Im): {f2b}")
    log(f"F2c: regimes assigned without trial: {f2c} — {dict(reg_census)}")
    log(f"F2d: complete section, target gauge > 0 (min {gmin:.3e}): {f2d}")
    log(f"F2e: (i) false H bridge refused {f2e_i}; (ii) component sigma "
        f"refused on the bridge {f2e_ii}")
    log(f"F2f: {len(rows)}/{len(clipped)} with no filtering: {f2f}")

    # --- F3a/F3b: the two transitions --------------------------------
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
    sheet_record_census = Counter(tuple(r["F2d_ledger_derived"]) for r in rows)
    log(f"F3a: open overlaps (minimum width {wid_min:.3e}) plus anchors "
        f"strictly interior and OFF the face: {f3a}")
    log(f"F2d(bis): bridge sign pattern DERIVED from the lower side, "
        f"{dict(sheet_record_census)} (default {list(eps)})")
    log(f"F3b(lower): the bridge glues EXACTLY to the lower side, theta = +1 "
        f"on every row: {f3b_lo}; minimum margin {marg_min:.3e}; "
        f"supremum of the recentred difference {diff_max:.3e}")
    log(f"F3b(upper, diagnostic): CONJUGATE upper side, theta pattern "
        f"{dict(pat_census)} = the deck transformation D. This is NOT "
        f"an atlas transition: see R2 for the continued sheet.")

    # --- R3: THE DOMAIN GRAPH AND THE REAL NERVE ---------------------
    # v1 confused the two: it added a bridge-to-side edge as soon
    # as a GEOMETRIC overlap existed, including when the gluing
    # had explicitly FAILED, and bridge-to-bridge edges without any
    # computed transition, then enumerated triples from that graph.
    # A NERVE edge is a CERTIFIED transition. Nothing else.
    _G["sheet_records"] = {r["tile"]: tuple(r["F2d_ledger_derived"])
                     for r in rows}
    bb_geo = [(a, b) for a, b in itertools.combinations(clipped, 2)
              if inter(bridges[a], bridges[b]) is not None]
    log(f"R3: {len(bb_geo)} geometric bridge-to-bridge overlaps — "
        f"transitions en cours de certification…")
    with mpctx.Pool(N_WORKERS, initializer=_init,
                    initargs=((S, g, eps), leaves, halos, bridges,
                              _G["sheet_records"])) as pool:
        bb = pool.map(_bb_job, bb_geo)
    bb_ok = [x for x in bb if x.get("certified")]
    bb_bad = [x for x in bb if not x.get("certified")]
    bb_diff = max((x.get("diff_sup") or 0.0) for x in bb_ok) if bb_ok else 0.0
    bb_marg = min((m for x in bb_ok for m in x["margins"]
                   if m is not None), default=None)
    log(f"     bridge-to-bridge CERTIFIED {len(bb_ok)}/{len(bb_geo)} "
        f"(theta = +1 on all 6 coordinates); minimum margin "
        f"{bb_marg if bb_marg is None else round(bb_marg, 4)} ; "
        f"supremum difference {bb_diff:.3e}; uncertified {len(bb_bad)}")

    # the DOMAIN graph (geometric, published as such)
    dom_edges = (2 * len(clipped)) + len(bb_geo)

    # NERVE: only the edges whose transition is CERTIFIED.
    # Nodes: the 316 charts of the lower atlas plus the 64 bridges. The
    # lower-to-lower edges are IMPORTED from the atlas step and checked
    # green, not recomputed.
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

    # NEW triples (involving at least one bridge), from the certified
    # edges alone. The purely lower triples are those of the
    # atlas step, already certified there; they are not recounted here.
    new_triples = []
    for i in clipped:
        x = ("B", i)
        for y, z in itertools.combinations(sorted(adj[x]), 2):
            if y in adj[z]:
                W1 = inter(boxof(x), boxof(y))
                W = inter(W1, boxof(z)) if W1 else None
                if W is not None:
                    new_triples.append([list(x), list(y), list(z)])
    log(f"R3: NERVE (certified edges ONLY) — {len(nodes)} nodes "
        f"({n_lower} lower + {len(clipped)} bridges), "
        f"{len(lower_pairs)} lower-to-lower edges imported from the atlas step, {n_bl} "
        f"bridge-to-leaf edges, {len(bb_ok)} bridge-to-bridge edges, "
        f"{len(new_triples)} NEW triples; connected={connected}")
    log(f"     (the DOMAIN graph, for its part, counts {dom_edges} "
        f"geometric ones; it is published under that name, not as a nerve)")
    f3c = bool(connected and n_bl == len(clipped)
               and len(bb_ok) == len(bb_geo))

    # --- R4a: the AMBIENT scope, stratified --------------------------
    re_flush = Counter()
    for i in clipped:
        fl = halos[i]["flush_faces"]
        re_flush[sum(1 for k in (0, 2) if fl[k] != 0)] += 1
    fully_ambient = re_flush[0]
    log(f"R4a: scope — {len(clipped)}/{len(clipped)} bilateral bridges "
        f"in BOTH imaginary directions; open in ALL FOUR "
        f"coordinates: {fully_ambient}/{len(clipped)}; still "
        f"relative to 1 Re face: {re_flush[1]}; to 2 faces: "
        f"{re_flush[2]}")

    # --- F3d + R1e: the negative controls ----------------------------
    f3d = all(r.get("F3d_non_scalar_mutation_breaks_gluing") for r in rows)
    glob_same = all(r.get("F3d_global_negation_is_same_projective_point")
                    for r in rows)
    r1e = all(r.get("R1e_wrong_deck_breaks") for r in rows)
    log(f"F3d: NON SCALAR mutation D_bad = diag(1,-1,1,1,1,1) gives "
        f"the gluing to fall: {f3d}  (whereas the GLOBAL negation stays "
        f"the SAME projective point: {glob_same}, which is why v1 "
        f"discriminated nothing)")
    log(f"R1e: a single sign changed in D makes the deck identity "
        f"BREAK: {r1e}")

    # --- R1: the deck theorem ----------------------------------------
    alg = deck_algebra(DECK_D)
    deck_meas = all(r["sides"]["upper_conj"].get(
        "deck_D_measured_equals_preregistered") for r in rows)
    r1 = bool(alg["involution"] and alg["non_scalar_in_pgl6"]
              and alg["preserves_the_three_quadrics"] and deck_meas
              and r1e)
    log(f"R1 : D = {list(DECK_D)} — involution {alg['involution']}, "
        f"non scalar {alg['non_scalar_in_pgl6']}, preserves the 3 "
        f"quadrics {alg['preserves_the_three_quadrics']} (TRIVIAL "
        f"for any diagonal of signs); D MEASURED == preregistered "
        f"on {len(rows)} bridges: {deck_meas} => {r1}")

    # --- R2: the CONTINUED upper sheet --------------------------------
    r2 = all(r["sides"].get("upper_cont", {}).get("glued_exactly")
             for r in rows)
    conj_fail = all(not r["sides"].get("upper_conj", {}).get(
        "glued_exactly", False) for r in rows)
    log(f"R2: the bridge glues EXACTLY to D.Z_conj (the CONTINUED "
        f"sheet) on {n_bu_cont}/{len(rows)}: {r2}, and NOT to "
        f"Z_conj alone ({conj_fail}), which stays a diagnostic")

    checks = {
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
    npass = sum(1 for v in checks.values() if v)

    try:
        head = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=Path(__file__).resolve().parent, capture_output=True,
            text=True, timeout=10).stdout.strip() or None
    except Exception:
        head = None

    out = {
        "artifact": "bridge_continuation",
        # This script has only one mode: the 64 bridges, always all of them.
        # Serialised EXPLICITLY: the upstream check of the metric step
        # requires mode == "full" of every upstream, and null is not a mode.
        "mode": "full",
        "claim_level": f1["claim_level"],
        "claim": (
            "The 64 bilateral BRIDGE CHARTS are built on the "
            "corrected 2H geometry; their regime is assigned PER ROW "
            "from the certified sign of Re R without trial; their record is "
            "DERIVED from the lower side then re-checked. Each glues "
            "EXACTLY to the LOWER side and to the CONTINUED "
            "SHEET `D.Z_conj` on OPEN overlaps, with anchors "
            "strictly interior. `D = diag(+1,-1,+1,+1,+1,-1)` is "
            "a DECK TRANSFORMATION: it is what separates the two "
            "sheets, rather than a bare statement that the conjugate is not "
            "the continued sheet, and the CONJUGATE atlas ALONE does NOT glue "
            "(64 of 64); it stays a diagnostic. STRATIFIED SCOPE: 64 of 64 "
            "bridges are bilateral in both imaginary directions, but "
            "only 36 of 64 are open in ALL FOUR coordinates; 24 "
            "stay relative to one real face, 4 to two faces, and those "
            "faces belong to the codimension-1 neighbours. The NERVE "
            "counts only CERTIFIED edges (380 nodes: 316 "
            "lower plus 64 bridges; 5396 lower-to-lower edges imported "
            "from the atlas step, 64 bridge-to-lower, 210 of 210 bridge-to-bridge); the geometric graph is "
            "published separately as `domain_intersection_graph`. "
            "The result stays a LOCAL ANALYTIC CONTINUATION TOWARDS A "
            "DERIVED MIRROR BOX (canonical identification having been "
            "refused). No metric is certified here."),
        "cell": {"S": list(S), "g": g, "eps": list(eps)},
        "n_bridges": len(clipped),
        "theta_required_preregistered": THETA_REQUIRED,
        "regime_census": {str(k): v for k, v in reg_census.items()},
        "theta_census": {f"{k[0]}|theta={k[1]}": v
                         for k, v in th_census.items()},
        "upper_theta_pattern_census": {str(k): v
                                       for k, v in pat_census.items()},
        "bridge_ledger_census": {str(k): v
                                 for k, v in sheet_record_census.items()},
        "real_structure_finding": (
            "ON THE CORNER, THE ANTIHOLOMORPHIC INVOLUTION ACTS WITH "
            "MIXED SIGNS. R is REAL there: positive on the principal "
            "row (REAL root, fixed by conjugation) and "
            "negative on the two canonical rows (root i.sqrt(-R) "
            "PURELY IMAGINARY, hence NEGATED). The conjugate point and the "
            "point continued analytically across the corner are therefore "
            "TWO DISTINCT POINTS of projective space above the same base point; the "
            "gauge being normalised to Z_g = 1, independent signs per "
            "row are NOT a projective renormalisation. "
            "CONSEQUENCE FOR THE FIRST LEVER OF THE SCOPING NOTE: the "
            "conjugate neighbour is a legitimate atlas, but it is NOT the one "
            "that analytic continuation reaches. The theta pattern is "
            "PREDICTED row by row from the regime, and the prediction "
            "is verified on 64 of 64: a result, not a "
            "residual surprise."),
        "min_theta_margin": marg_min,
        "max_recentred_difference_sup": diff_max,
        "min_overlap_width": wid_min,
        "min_target_gauge_absmin": gmin,
        "upstream_chain": up,
        "deck": {**deck_algebra(DECK_D),
                 "identity": "Z_upper_conj = D . Z_bridge on the upper "
                             "overlap, certified AT THE COEFFICIENT level",
                 "measured_on_bridges": len(rows),
                 "why_the_algebra_is_not_the_content": (
                     "any diagonal of signs preserves the quadrics; "
                     "what singles out THIS D is that it is "
                     "the one, and the only one, linking the conjugate section "
                     "to the bridge section; the R1e negative control "
                     "carries the discriminating power, not the algebra")},
        "domain_intersection_graph": {
            "note": ("A GRAPH OF DOMAINS, NOT an atlas nerve: an "
                     "edge here is only an intersection of boxes of "
                     "positive width. Published under that name since a "
                     "review showed that v1 counted as edges "
                     "gluings explicitly in failure."),
            "n_bridge_side_geometric": 2 * len(clipped),
            "n_bridge_bridge_geometric": len(bb_geo),
            "n_edges": dom_edges},
        "nerve": {
            "note": ("NERVE: only the edges whose TRANSITION "
                     "is certified. The lower-to-lower edges "
                     "are IMPORTED from the atlas step (verified green), "
                     "not recomputed."),
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
            "note": ("STRATIFIED, corrected after review: calling the 64 charts "
                     "become ambient was an overreach."),
            "bilateral_in_both_im_directions": len(clipped),
            "fully_open_in_all_four_coordinates": fully_ambient,
            "still_relative_to_one_re_face": re_flush[1],
            "still_relative_to_two_re_faces": re_flush[2]},
        "negatives": {
            "F2e_i": "false bridge of half-width H: inclusion refused",
            "F2e_ii": ("component rule on the bridge: REFUSAL "
                       "`component_sigma_undetermined_on_bridge`, which is "
                       "the structural reason for the canonical regime"),
            "F2e_iii": ("an enlarged box where Re R straddles: REFUSAL "
                        "`re_R_straddles_zero_no_regime` (self-test)"),
            "F3d": ("NON SCALAR mutation D_bad = diag(1,-1,1,1,1,1) "
                    "(gauge and affine coordinates UNCHANGED, a single "
                    "canonical row negated) makes the gluing fall. v1 "
                    "negated ALL coordinates, yet Z and -Z are THE "
                    "SAME point of projective space: it tested only the "
                    "sensitivity of the normalisation Z_g = 1 to a "
                    "change of representative, not the detection of a "
                    "wrong POINT. A discriminance debt, raised by the "
                    "review and paid here."),
            "F3d_global_negation_control": (
                "the global negation stays the SAME projective point on "
                "64 of 64; published to say why v1 was worthless"),
            "R1e": ("a single sign changed in D makes the deck identity "
                    "BREAK: it is this negative control, and not the algebra of "
                    "the quadrics, that carries the discriminating power of the theorem")},
        "triple_cocycle_note": (
            "The 588 NEW triples are enumerated from the certified "
            "edges alone. Their cocycle theta_ij.theta_jk.theta_ki = +1 is "
            "IMPLIED: each theta equals +1 on a domain CONTAINING the "
            "triple box, and theta is discrete. It is therefore NOT an independent "
            "test and is not checked as one, with the same "
            "honesty as the atlas step on its own cocycle. What "
            "would remain to pay at the triple is METRIC congruence, "
            "which belongs to the metric step."),
        "per_bridge": rows,
        "not_paid_here": [
            "the METRIC of the bridge: Qmat, Weyl, lateral "
            "congruences: the full-run path requires the bilateral "
            "constructor, which now exists, but the run does not use it",
            "the canonical identification of the neighbour, which REFUSED; "
            "the upper atlas is DERIVED by conjugation, not enumerated",
            "gluing across the real faces, where the bridges stay "
            "RELATIVE charts (margin 0 against the face of the cell)",
            "the codimension-1 neighbours, the 895 other pairs, the later scaling"],
        "checks": checks, "checks_passed": npass, "checks_total": len(checks),
        "verdict": (
            "DELIVERED: 64 bilateral bridges, glued "
            "EXACTLY to the lower side (64) AND to the CONTINUED "
            "sheet D.Z_conj (64), 210 of 210 bridge-to-bridge transitions "
            "certified, a connected nerve of 380 nodes on edges "
            "certified ONLY. D = diag(+,-,+,+,+,-) is the "
            "deck transformation separating the conjugate from the continued sheet."
            if npass == len(checks) else
            f"RED: {len(checks) - npass} check(s) failed"),
        "provenance": {
            "git_head": head, "python": sys.version.split()[0],
            "platform": platform.platform(), "mp_prec": int(mp.prec),
            "tm_order": TM_ORDER, "unary_series_deg": UNARY_SERIES_DEG,
            "wall_s": round(time.time() - T0, 1), "n_workers": N_WORKERS,
            "preregistered": {"im_dirs": list(IM_DIRS),
                              "theta_required": THETA_REQUIRED},
            "inputs": {p.name: _sha(p) for p in
                       (COVER_JSON, ATLAS_JSON, C127E_JSON, PRELIMINARY_JSON,
                        F1_JSON)},
            "self_sha256": _sha(__file__)}}

    ART.write_text(json.dumps(out, indent=2, ensure_ascii=False),
                   encoding="utf-8")
    print("=" * 78)
    for k, v in checks.items():
        print(f"  {'OK  ' if v else 'FAIL'} {k}")
    print(f"\n{out['verdict']}")
    print(f"checks {npass}/{len(checks)} - artefact: {ART.name}")
    print("=" * 78)
    return npass == len(checks)


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(0 if selftest() else 1)
    sys.exit(0 if build() else 1)
