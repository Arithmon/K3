#!/usr/bin/env python3
"""
Closing the residual 1/64 of the volume by the THIRD determination.

The dyadic cover leaves 64 boxes (1/64 of the volume) without a chart: on
each, every candidate chart has at least one section row whose radicand
`R'` has `Re R' < 0` strictly and `Im R'` straddling 0 (a mechanism
measured by the chart criterion, on 29/29 probed rows). For those rows:

  . the PRINCIPAL branch is refused, `R'` meeting the cut on the
    non-positive reals;
  . the ROTATED branch of the rotated continuation has its guard
    SATISFIED, `R'` avoiding [0,+inf) since `Re R' < 0`, but its LABEL
    `sigma' = sign(Im R')` is UNDETERMINED (the enclosure of `Im R'`
    contains 0). The component rule forbids choosing it by trial, so the
    row was REFUSED.

THE THIRD DETERMINATION lifts this lock WITHOUT a new guard and WITHOUT
trial. Algebraic observation: in the section, the component and the sheet enter
only through their PRODUCT, `Z_s = sheet.component.i.sqrt_p(-R')`. One can therefore FIX
the component to +1 (canonical): `w = i.sqrt_p(-R')` is a legitimate continuous branch
on the plane cut along the non-negative reals (the guard applied to `-R'`, verbatim), and
the SHEET choice (`w` or `-w`) is carried by the sheet label, derived at the ANCHOR
by the STRICT SEPARATION of the transport step (margin > 0 serialised,
refusal if ambiguous) then verified with the SAME sign on the whole box.
The sheet certificate is the MARGIN: not a trial, not the component.

What is requalified and what is not:
  . the component-based sigma of the rotated continuation stays
    NECESSARY when one must
    coincide with the principal branch on an overlap; here the principal
    one DOES NOT EXIST on the row, so there is nothing to glue
    locally. Gluing BETWEEN tiles is the business of the cocycle,
    not of this script;
  . the SOURCE SECTION is NOT touched (the warning about a rotated cut
    and not relabelling the classes unilaterally
    concerns the source; here only the NATIVE section of the target chart
    of the residual tiles uses the third determination, and
    every such use is serialised).

Preregistered checks:
  R1 SERIALISED MECHANISM: every new tile publishes, per row, the
     determination, the provenance of the label, the enclosures of
     `Re R'` and `Im R'`, and the residual `w^2 - R'`.
  R2 UNIVERSAL IDENTITY: `w^2 - R'` contains 0 on EVERY row of EVERY new
     tile (widths published).
  R3 UNIVERSAL SHEET MARGIN: margin > 0 on every row.
  R4 NON-TAUTOLOGY: for each new tile, the chart criterion
     (imported, not reimplemented) REFUSES the same chart on the same box
     and the extended criterion ACCEPTS it; the third determination is
     used ONLY if both the principal and the certified-component rotated ones are refused
     (strict order, counted).
  R5 THE COVER CLOSES EXACTLY: dyadic addresses of the 252 tiles
     the transport step, the new tiles and the remaining residual:
     prefix-free, closed,
     Kraft equal to 1 (imported machinery). The number of boxes NOT
     closed is NOT pre-committed: it is published.
  R6 UNIVERSAL TRANSPORT ON THE NEW TILES: the same checks as
     the transport step: congruence, same projective point, relative
     residual at most 1e-5 (preregistered, unchanged), positive
     definiteness BY WEYL, zero failure filtered out.
  R7 SHEET NEGATIVE CONTROL: flipping the sheet on a row with a third
     determination BREAKS projective compatibility (probes).
  R8 CANONICAL COMPONENT NEGATIVE CONTROL: the rotated root at -1 is the EXACT NEGATION
     of `tm_sqrt_rotated(R, +1)`, bit for bit: (eps', sigma') enter only
     through their product, so a canonical sigma hides no degree of
     freedom.
  R9 NO SILENT CEILING: counts published and consistent.

Outputs: results/residual_closure{_pilot}.json
Usage  : residual_closure.py [--selftest]
Env    : K3_C127E_MODE    pilot (8 boxes, the default) | full (64)
         K3_C127E_WORKERS processes (default 4)
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
from multiprocessing import get_context
from pathlib import Path

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

_argv = sys.argv
sys.argv = [sys.argv[0]]
os.environ.setdefault("K3_TM_ORDER", "4")
os.environ.setdefault("K3_TM_SERIES", "4")
from mpmath import mp                                              # noqa: E402
from .witness_registry import load_canonical_MH              # noqa: E402
from .interval_arithmetic import (                            # noqa: E402
    BranchCutError, CIV, build_M_civ, minor_inv_times_T_exact)
from .taylor_models import (                                     # noqa: E402
    TMC, TM_ORDER, UNARY_SERIES_DEG, riv,
    rotated_sigma_from_coeffs, tm_sqrt_rotated)
from .owner_tiling import TRIPLES                    # noqa: E402
from .full_cell_charts import (                           # noqa: E402
    build_section, chart_certificate)
from .chart_selection_criterion import (                  # noqa: E402
    native_section_constructible, target_uv)
from .gram_congruence import GAMMA                     # noqa: E402
from .chart_transport import (                    # noqa: E402
    address_of, frontier_fractions, transport_hardened, tree_gates)
sys.argv = _argv

RES = Path(os.environ.get(
    "K3_RES_DIR", Path(__file__).resolve().parent / "data"))
COVER_JSON = RES / "dyadic_cover.json"
C127_JSON = RES / "chart_transport.json"
MODE = os.environ.get("K3_C127E_MODE", "pilot")
N_WORKERS = int(os.environ.get("K3_C127E_WORKERS", "4"))
ART = RES / ("residual_closure.json" if MODE == "full"
             else "residual_closure_pilot.json")

# --- PRÉ-ENREGISTRÉ ---------------------------------------------------------
DELTA_REL = 1e-5        # the same ceiling as the transport step, unchanged
N_PILOT_BOXES = 8       # boîtes résiduelles en mode pilot
N_PROBE_NEG = 4         # tiles probed by the negative controls R7/R8

T0 = time.time()
ALL_CHARTS = [(tuple(S2), g2) for S2 in TRIPLES
              for g2 in range(6) if g2 not in S2]


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def _sha(p):
    try:
        return hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
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
                              "n_pilot_boxes": N_PILOT_BOXES,
                              "n_probe_neg": N_PROBE_NEG},
            "inputs": {str(Path(x).name): _sha(x) for x in src},
            "self_sha256": _sha(__file__)}


def _rng(tm, im=False):
    x = tm.im_tm().to_iv() if im else tm.re_tm().to_iv()
    return (float(mp.mpf(x.a)), float(mp.mpf(x.b)))


def _exact_negation(wp, wm):
    """Is `wm` the EXACT NEGATION of `wp`? Checked coefficient by
    coefficient (the negated interval of [a,b] is [-b,-a]) and on the
    remainder, NOT on the enclosure of the sum, whose remainder adds
    `rem(wp) + rem(wm) > 0` and can never be exactly {0} on
    a true Taylor model (lesson from the pilot: check R8 v1 failed for that
    reason, on correct code)."""
    if len(wp.p) != len(wm.p):
        return False
    for c1, c2 in zip(wp.p, wm.p):
        for part in ("re", "im"):
            a1, b1 = getattr(c1, part).a, getattr(c1, part).b
            a2, b2 = getattr(c2, part).a, getattr(c2, part).b
            if not (mp.mpf(a2) == -mp.mpf(b1)
                    and mp.mpf(b2) == -mp.mpf(a1)):
                return False
    return (mp.mpf(wp.rem.a) == mp.mpf(wm.rem.a)
            and mp.mpf(wp.rem.b) == mp.mpf(wm.rem.b))


# ===========================================================================
#  The EXTENDED criterion, three steps in strict order
# ===========================================================================
def native_rows_ext(S2, g2, up, vp):
    """Mirror of `native_section_constructible` with the THIRD
    step. STRICT order per row:
      (1) plain;
      (2) rotated with a CERTIFIED sigma (enclosure of Im R', then the
          signs of the coefficients), never by trial;
      (3) CANONICAL rotated, sigma = +1, only if (1) AND (2) refused;
          the guard is that of the rotated root (applied to -R'),
          the sheet label is delegated to the sheet variable (margin).
    """
    T2 = tuple(j for j in range(6) if j not in S2)
    others = [c for c in T2 if c != g2]
    Ae = minor_inv_times_T_exact(S2, T2)
    perm = [list(T2).index(g2), list(T2).index(others[0]),
            list(T2).index(others[1])]
    A = [[riv(Ae[r][c]) for c in perm] for r in range(3)]
    ur, ui = _rng(up), _rng(up, True)
    vr, vi = _rng(vp), _rng(vp, True)
    u2, v2 = up * up, vp * vp
    rows = []
    for r, s in enumerate(S2):
        R = TMC.const(CIV(A[r][0])) + u2.mul_real(A[r][1]) \
            + v2.mul_real(A[r][2])
        rl, rh = _rng(R)
        il, ih = _rng(R, True)
        rec = {"row": r, "s_coord": int(s), "Re_R": [rl, rh],
               "Im_R": [il, ih], "determination": None, "sigma": None,
               "sigma_source": None}
        Zp = None
        try:
            Zp = R.sqrt_principal()
            rec["determination"] = "principal"
        except BranchCutError as exc:
            rec["principal_refused"] = exc.diag.get("guard")
            if il > 0:
                sg, how = 1, "Im_R_enclosure"
            elif ih < 0:
                sg, how = -1, "Im_R_enclosure"
            else:
                sg = rotated_sigma_from_coeffs(
                    Fraction(Ae[r][perm[1]]), Fraction(Ae[r][perm[2]]),
                    ur, ui, vr, vi)
                how = "coefficient_signs"
            if sg in (-1, 1):
                rec["sigma"], rec["sigma_source"] = sg, how
                try:
                    Zp = tm_sqrt_rotated(R, sg)
                    rec["determination"] = "rotated"
                except BranchCutError as exc2:
                    rec["rotated_refused"] = exc2.diag.get("guard")
            else:
                rec["rotated_refused"] = "rotated_component_undetermined"
                # (3) the THIRD determination: canonical sigma, sheet
                # déléguée à ε'. Garde inchangée (the cut guard sur −R').
                try:
                    Zp = tm_sqrt_rotated(R, 1)
                    rec["determination"] = "rotated_canonical"
                    rec["sigma"], rec["sigma_source"] = 1, "canonical"
                except BranchCutError as exc3:
                    rec["canonical_refused"] = exc3.diag.get("guard")
        if Zp is not None:
            d = Zp * Zp - R
            dr, di = _rng(d), _rng(d, True)
            rec["sq_residual_re"], rec["sq_residual_im"] = dr, di
            rec["sq_residual_contains_zero"] = bool(
                dr[0] <= 0 <= dr[1] and di[0] <= 0 <= di[1])
        rows.append(rec)
    return rows, all(x["determination"] is not None for x in rows)


def strong_chart_ok_ext(Z, dZ, S2, g2):
    """Chart criterion (domain, Jacobian, disjointness) with the
    constructibilité ÉTENDUE."""
    cert = chart_certificate(Z, dZ, S2, g2)
    if not (cert.get("admissible")
            and cert.get("disjoint_from_target_slice")):
        return False, cert, None
    uv = target_uv(Z, S2, g2)
    if uv is None:
        return False, cert, None
    rows, ok = native_rows_ext(S2, g2, uv[0], uv[1])
    return bool(ok), cert, rows


# ===========================================================================
#  Worker: chart search, then transport (fork, inherited state)
# ===========================================================================
_G = {}


def _init_worker(cell, M, c218, rw):
    _G["cell"], _G["M"], _G["c218"], _G["rw"] = cell, M, c218, rw


def _search_box(job):
    """Look for a chart for a residual box, extended criterion,
    seeds = charts of the cover (nearest first)."""
    bi, box, seeds = job
    S, g, eps = _G["cell"]
    c = [float.fromhex(x) for x in box["center_hex"]]
    h = float.fromhex(box["hw_hex"])
    Z, dZ, _rows = build_section(S, g, eps, c, h)
    if any(z is None for z in Z):
        return {"box_index": bi, "found": False,
                "reason": "source_section_incomplete"}
    order = seeds + [ch for ch in ALL_CHARTS if ch not in seeds]
    n_tried = 0
    for (S2, g2) in order:
        n_tried += 1
        ok, cert, rows = strong_chart_ok_ext(Z, dZ, S2, g2)
        if ok:
            # R4 non-tautology: the imported chart criterion must REFUSE
            uv = target_uv(Z, S2, g2)
            _old_rows, old_ok = native_section_constructible(
                S2, g2, uv[0], uv[1])
            return {
                "box_index": bi, "found": True,
                "center_hex": box["center_hex"],
                "hw_hex": box["hw_hex"], "depth": box["depth"],
                "chart": {"S": list(S2), "g": int(g2)},
                "n_charts_tried": n_tried,
                "rows": rows,
                "n_canonical_rows": sum(
                    1 for x in rows
                    if x["determination"] == "rotated_canonical"),
                "old_criterion_accepts": bool(old_ok),
                "gauge_absmin": cert.get("gauge_absmin"),
                "detJ_absmin": cert.get("detJ_absmin")}
    return {"box_index": bi, "found": False, "reason": "no_chart",
            "n_charts_tried": n_tried,
            "center_hex": box["center_hex"], "hw_hex": box["hw_hex"],
            "depth": box["depth"]}


def _transport_tile(job):
    t, kind, extra = job
    S, g, eps = _G["cell"]
    c = [float.fromhex(x) for x in t["center_hex"]]
    h = float.fromhex(t["hw_hex"])
    S2, g2 = tuple(t["chart"]["S"]), t["chart"]["g"]
    sigma2 = [x["sigma"] for x in t["rows"]]
    kw = {"fixed_sigma2": sigma2}
    if kind == "eps_flip":
        kw["fixed_eps2"] = extra["flipped_eps"]
    r = transport_hardened(S, g, eps, c, h, S2, g2,
                           _G["M"], _G["c218"], _G["rw"], **kw)
    r["box_index"], r["kind"] = t["box_index"], kind
    r["chart"] = t["chart"]
    r["center_hex"], r["hw_hex"] = t["center_hex"], t["hw_hex"]
    return r


# ===========================================================================
#  Build
# ===========================================================================
def build():
    print("=" * 78)
    print(f"residual closure BY THE THIRD DETERMINATION ({MODE}): "
          f"TM ({TM_ORDER},{UNARY_SERIES_DEG}), {N_WORKERS} workers, "
          f"relative delta = {DELTA_REL:.0e}")
    print("=" * 78)
    reg = load_canonical_MH()
    M = build_M_civ(reg["M_H_canonical"])
    c218 = reg["coeffs218"]
    rw = 1.0 - GAMMA
    cov = json.loads(COVER_JSON.read_text(encoding="utf-8"))
    tiles_old, residual = cov["tiles"], cov["residual"]
    cell = cov["cell"]
    S, g, eps = tuple(cell["S"]), cell["g"], tuple(cell["eps"])
    root_c = [float.fromhex(x) for x in cell["center_hex"]]
    root_h = float.fromhex(cell["hw_hex"])
    log(f"cover : {len(tiles_old)} tuiles, {len(residual)} boîtes "
        f"résiduelles · cellule S={list(S)} g={g}")

    boxes = residual if MODE == "full" else residual[:N_PILOT_BOXES]

    # seeds: the charts of the cover, sorted by proximity to each box
    def seeds_for(box):
        cb = [float.fromhex(x) for x in box["center_hex"]]
        d = {}
        for t in tiles_old:
            ct = [float.fromhex(x) for x in t["center_hex"]]
            key = (tuple(t["chart"]["S"]), t["chart"]["g"])
            dist = max(abs(cb[k] - ct[k]) for k in range(4))
            d[key] = min(d.get(key, math.inf), dist)
        return [k for k, _v in sorted(d.items(), key=lambda kv: kv[1])]

    mpctx = get_context("fork")
    jobs = [(i, b, seeds_for(b)) for i, b in enumerate(boxes)]
    with mpctx.Pool(N_WORKERS, initializer=_init_worker,
                    initargs=((S, g, eps), M, c218, rw)) as pool:
        found = pool.map(_search_box, jobs)
    new_tiles = [r for r in found if r["found"]]
    unclosed = [r for r in found if not r["found"]]
    n_canon = sum(t["n_canonical_rows"] for t in new_tiles)
    log(f"recherche : {len(new_tiles)}/{len(boxes)} boîtes fermées "
        f"({n_canon} lignes en 3ᵉ détermination), {len(unclosed)} non "
        f"fermées (PUBLIÉES)")

    # --- R6: transport on all the new tiles -----------------------------------
    tjobs = [(t, "transport", None) for t in new_tiles]
    with mpctx.Pool(N_WORKERS, initializer=_init_worker,
                    initargs=((S, g, eps), M, c218, rw)) as pool:
        transported = pool.map(_transport_tile, tjobs)
    tr_fail = [r for r in transported if r.get("failed")]
    tr_ok = [r for r in transported if not r.get("failed")]
    log(f"transport : {len(tr_ok)}/{len(new_tiles)} OK, "
        f"{len(tr_fail)} échecs (REFUSÉS, pas filtrés)")

    # --- R7: sheet negative control (sheet flipped on a canonical row) --------
    probe = [t for t in new_tiles if t["n_canonical_rows"] > 0]
    probe = probe[:N_PROBE_NEG]
    pjobs = []
    for t in probe:
        r0 = next(x for x in tr_ok if x["box_index"] == t["box_index"])
        canon_row = next(i for i, x in enumerate(t["rows"])
                         if x["determination"] == "rotated_canonical")
        fl = list(r0["eps_target"])
        fl[canon_row] = -fl[canon_row]
        pjobs.append((t, "eps_flip", {"flipped_eps": fl}))
    with mpctx.Pool(min(N_WORKERS, max(1, len(pjobs))),
                    initializer=_init_worker,
                    initargs=((S, g, eps), M, c218, rw)) as pool:
        negs = pool.map(_transport_tile, pjobs) if pjobs else []
    neg_breaks = [bool(n.get("failed"))
                  or not n["same_projective_point"] for n in negs]
    log(f"R7: a flipped eps' breaks compatibility {sum(neg_breaks)}"
        f"/{len(neg_breaks)}")

    # --- R8: a canonical sigma hides no degree of freedom -------------
    r8_ok = []
    for t in probe:
        c = [float.fromhex(x) for x in t["center_hex"]]
        h = float.fromhex(t["hw_hex"])
        Z, dZ, _r = build_section(S, g, eps, c, h)
        S2, g2 = tuple(t["chart"]["S"]), t["chart"]["g"]
        uv = target_uv(Z, S2, g2)
        T2 = tuple(j for j in range(6) if j not in S2)
        others = [x for x in T2 if x != g2]
        Ae = minor_inv_times_T_exact(S2, T2)
        perm = [list(T2).index(g2), list(T2).index(others[0]),
                list(T2).index(others[1])]
        rr = next(i for i, x in enumerate(t["rows"])
                  if x["determination"] == "rotated_canonical")
        A = [riv(Ae[rr][perm[j]]) for j in range(3)]
        u2, v2 = uv[0] * uv[0], uv[1] * uv[1]
        R = TMC.const(CIV(A[0])) + u2.mul_real(A[1]) + v2.mul_real(A[2])
        wp = tm_sqrt_rotated(R, 1)
        wm = tm_sqrt_rotated(R, -1)
        r8_ok.append(_exact_negation(wp, wm))
    log(f"R8: w(sigma=-1) is the EXACT negation of w(sigma=+1), "
        f"coefficient by coefficient: {sum(r8_ok)}/{len(r8_ok)}")

    # --- R5: the cover, recounted from the addresses ------------------
    leaves = tiles_old + new_tiles + unclosed
    addresses, addr_fail = [], 0
    residual_addr = []
    for lf in leaves:
        a = address_of(root_c, root_h,
                       [float.fromhex(x) for x in lf["center_hex"]],
                       float.fromhex(lf["hw_hex"]))
        if a is None:
            addr_fail += 1
        else:
            addresses.append(a)
            if lf in unclosed:
                residual_addr.append(a)
    # in pilot mode, the residual boxes NOT examined count towards the
    # residual so that the partition stays exact
    if MODE != "full":
        for b in residual[N_PILOT_BOXES:]:
            a = address_of(root_c, root_h,
                           [float.fromhex(x) for x in b["center_hex"]],
                           float.fromhex(b["hw_hex"]))
            if a is None:
                addr_fail += 1
            else:
                addresses.append(a)
                residual_addr.append(a)
    tg = tree_gates(addresses)
    fr = frontier_fractions(addresses)
    vol_res = sum(Fraction(1, 16 ** len(a)) for a in residual_addr)
    vol_cov = 1 - vol_res
    log(f"R5 : Kraft {tg['kraft_sum'][0]}/{tg['kraft_sum'][1]}, couvert "
        f"{float(100 * vol_cov):.4f} %, résidu {float(100 * vol_res):.4f} %")

    # --- Checks ----------------------------------------------------------------
    max_rel = max((r["residual_relative"] for r in tr_ok), default=None)
    eps_m = [m for r in tr_ok for m in (r["eps_margins"] or [])]
    id_rows = [x for t in new_tiles for x in t["rows"]]
    checks = {
        "R1_mechanism_serialised": bool(new_tiles) and all(
            x.get("Re_R") and x.get("Im_R")
            and x.get("determination") for x in id_rows),
        "R2_identity_universal": bool(id_rows) and all(
            x["sq_residual_contains_zero"] for x in id_rows),
        "R3_eps_margin_positive_universal": bool(eps_m)
        and all(m["margin"] > 0 for m in eps_m),
        "R4_not_tautological": bool(new_tiles) and all(
            (not t["old_criterion_accepts"])
            and t["n_canonical_rows"] > 0 for t in new_tiles),
        "R5_cover_exact_partition": bool(
            addr_fail == 0 and tg["unique"] and tg["prefix_free"]
            and tg["tree_closed"] and tg["kraft_is_one"]),
        "R6_transport_universal_on_new_tiles": bool(
            len(tr_ok) == len(new_tiles) == len(transported)
            and not tr_fail) and all(
            r["gauge_invariance_ok"] and r["congruence_contains_zero"]
            and r["same_projective_point"]
            and r["residual_relative"] is not None
            and r["residual_relative"] <= DELTA_REL
            and r["spectral"]["weyl_transport_ok"]
            and r["pd_source"]["is_PD"] and r["pd_target"]["is_PD"]
            for r in tr_ok),
        "R7_eps_flip_breaks": bool(neg_breaks) and all(neg_breaks),
        "R8_sigma_canonical_no_dof": bool(r8_ok) and all(r8_ok),
        "R9_no_silent_cap": bool(
            len(found) == len(boxes)
            and len(transported) == len(new_tiles)
            and len(negs) == len(pjobs)
            and len(new_tiles) + len(unclosed) == len(boxes))}

    n_pass = sum(1 for v in checks.values() if v)
    verdict = (
        "Residual closure (%s): the THIRD DETERMINATION closes %d of %d residual "
        "boxes examined (%d rows on the canonical component, sheet "
        "pinned by the sheet margin, minimum %.3e). Non-tautology: the "
        "plain criterion refuses and the extended one accepts on EACH new "
        "tile. The identity w^2 = R' is certified on every row. "
        "Complete transport on the new tiles: %d of %d, "
        "maximum relative residual %.3e at most delta = %.0e, Weyl everywhere. The cover "
        "now covers %.4f %% (residual %.4f %%, %d boxes left "
        "unclosed and published). Negative controls: a flipped sheet breaks %d of %d, "
        "w(+1)+w(-1) vanishes exactly %d of %d. NOT PAID: the atlas step "
        "(halos/overlaps/cocycle), the join BETWEEN determinations of "
        "neighbouring tiles (cocycle), the exact contract, globalisation, "
        "the later scaling." % (
            MODE, len(new_tiles), len(boxes), n_canon,
            min((m["margin"] for m in eps_m), default=float("nan")),
            len(tr_ok), len(new_tiles),
            max_rel if max_rel is not None else float("nan"), DELTA_REL,
            float(100 * vol_cov), float(100 * vol_res), len(unclosed),
            sum(neg_breaks), len(neg_breaks), sum(r8_ok), len(r8_ok)))

    art = {"artifact": ART.stem, "mode": MODE,
           "claim": ("Closing the residual of the dyadic cover by the "
                     "third determination (canonical component, sheet "
                     "by margin), with the complete transport."),
           "cell": cell,
           "n_boxes_examined": len(boxes),
           "new_tiles": new_tiles, "unclosed": unclosed,
           "n_canonical_rows": n_canon,
           "transports": transported,
           "eps_flip_negatives": negs,
           "r8_sigma_check": r8_ok,
           "tree_gates": tg, "frontier_fractions": fr,
           "covered_volume": [vol_cov.numerator, vol_cov.denominator],
           "residual_volume": [vol_res.numerator, vol_res.denominator],
           "max_residual_relative": max_rel,
           "delta_rel_preregistered": DELTA_REL,
           "source_section_untouched": (
               "the SOURCE section never uses the third "
               "determination; only the NATIVE section of the target chart "
               "of the residual tiles is concerned, and the class labelling of "
               "the source arc is unchanged"),
           "not_paid_here": ["the atlas step halos/overlaps/cocycle",
                             "the join between determinations of "
                             "neighbouring tiles (nerve cocycle)",
                             "the exact contract", "globalisation", "the later scaling"],
           "verdict": verdict, "checks": checks,
           "checks_passed": n_pass, "checks_total": len(checks),
           "provenance": provenance([COVER_JSON, C127_JSON],
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
#  Self-test: the third determination on control radicands
# ===========================================================================
def _selftest():
    ok = True

    def chk(name, cond):
        nonlocal ok
        print(f"  {'ok  ' if cond else 'FAIL'} {name}")
        ok = ok and bool(cond)

    def const_tmc(re, im=0.0):
        return TMC.const(CIV(riv(re), riv(im)))

    # T1: R = -4 exactly (Re < 0, Im identically 0: the residual configuration)
    #      principale refusée, composante indéterminée, CANONIQUE passe
    #      et w² = R
    R = const_tmc(-4.0)
    p_refused = False
    try:
        R.sqrt_principal()
    except BranchCutError:
        p_refused = True
    w = tm_sqrt_rotated(R, 1)
    d = w * w - R
    dr, di = _rng(d), _rng(d, True)
    chk("T1 R=−4 : principale refusée, canonique passe, w² = R",
        p_refused and dr[0] <= 0 <= dr[1] and di[0] <= 0 <= di[1])

    # T2: w(component=-1) is the EXACT NEGATION of w(component=+1); the sheet and
    #      component enter only through their product. Checked COEFFICIENT BY COEFFICIENT:
    #      the enclosure of the sum cannot serve as a test, because
    #      adding Taylor models adds the remainders (rem+rem > 0 on a
    #      true model, which is the defect of check R8 v1, found in the pilot).
    wm = tm_sqrt_rotated(R, -1)
    chk("T2 w(σ=−1) = négation exacte de w(σ=+1) (coefficients + reste)",
        _exact_negation(w, wm))

    # T2b NEGATIVE CONTROL: w is NOT its own negation (the test cannot
    #      be satisfied by a function that returns True)
    chk("T2b negative control: w is not the negation of w",
        not _exact_negation(w, w))

    # T3: w = i.sqrt_p(-R) equals 2i for R = -4 (the right branch, not +-2)
    wr, wi = _rng(w), _rng(w, True)
    chk("T3 w(−4, σ=+1) = 2i (enclosure serrée autour de (0, 2))",
        abs(wr[0]) < 1e-12 and abs(wr[1]) < 1e-12
        and abs(wi[0] - 2.0) < 1e-12 and abs(wi[1] - 2.0) < 1e-12)

    # T4 NEGATIVE CONTROL: R = +4, the principal branch passes, so the strict order
    #      must NEVER invoke the canonical one (no substitution
    #      silent substitution for a valid determination)
    R4 = const_tmc(4.0)
    p_ok = True
    try:
        R4.sqrt_principal()
    except BranchCutError:
        p_ok = False
    chk("T4 R=+4 : la principale passe (la canonique resterait inerte)",
        p_ok)

    # T5 NEGATIVE CONTROL: R whose enclosure CONTAINS 0: all three
    #      determinations must refuse (no branch exists)
    R0 = TMC.const(CIV(riv([-1.0, 1.0]), riv([-1.0, 1.0])))
    refuse = 0
    try:
        R0.sqrt_principal()
    except BranchCutError:
        refuse += 1
    for sg in (1, -1):
        try:
            tm_sqrt_rotated(R0, sg)
        except BranchCutError:
            refuse += 1
    chk("T5 negative control: R containing 0 makes all three determinations refuse",
        refuse == 3)

    # T6: R = 4i (Im > 0 strictly): the component is CERTIFIABLE, so
    #      step (2) must suffice; the canonical one must not be
    #      reached when a certified component exists (this is the strict order of
    #      native_rows_ext, vérifié ici sur sa brique)
    Ri = const_tmc(0.0, 4.0)
    w2 = tm_sqrt_rotated(Ri, 1)
    d2 = w2 * w2 - Ri
    d2r, d2i = _rng(d2), _rng(d2, True)
    chk("T6 R=4i : σ=+1 certifiable par l'enclosure, w² = R tient",
        d2r[0] <= 0 <= d2r[1] and d2i[0] <= 0 <= d2i[1])

    print("SELF-TEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_selftest() if "--selftest" in sys.argv else build())
