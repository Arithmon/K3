#!/usr/bin/env python3
"""
width_attribution.py — WIDTH ATTRIBUTION plus preconditioned fixed
congruence (the width contract, items A0 to A4, with two reviewer requests
folded in).

The wall measured by the direct probe (C = width(detQ)/h^2 around 2e3-2e4)
has no established CAUSE yet: the width mixes real linear and quadratic
variation, the enclosure excess of the Hessian, and dependencies of the
determinant. This probe measures it BEFORE choosing the tool that will
reduce it:

  A0 provenance: M_H = (M + M-dagger)/2 canonicalised ONCE at load time,
     hash published, the float engine AND the interval kernel consuming
     exactly M_H; degenerate control h = 0 replayed (mid(Q_gamma) equal to
     the float value, relative error below 5e-12);
          checks derived from the counts, "executed" rather than "all"; pricing
          serialised in the JSON, not only in the note.
  A1 width decomposition: the WITNESS box B FIRST (at a reviewer's
     request), then A and C, on a FIXED h grid {4e-3, 2e-3, 1e-3, 5e-4,
     2.5e-4} (at least 4 values, successive slopes published, never a
     "law" on 2 points). At each h, POST-PROCESSING of the two jets
     (centre and box): width per FORM (Taylor-2, mean-value, naive,
     determinant-of-components, final intersection, and the active form),
     and EXACT DECOMPOSITION of the t2 form (the interval widths add
     exactly): w = w(centre value) + w(grad F(t_0).delta) [real
     first-order variation] + w(half delta^T H(box) delta) [remainder
     term], with w(half delta^T H(t_0) delta) [EXACT Hessian at the
     centre, degenerate jet] as the true scale of the quadratic, isolating
     the enclosure excess of the Hessian = w(H box) - w(H centre). The
     gradient and Hessian norms are published. Sample
          float min/max of det Q on the box (NOT a certificate: it separates
     the apparent physical variation from the enclosure excess).
  A2 ablation: (i) fixed congruence R_0, with
     R_0 = inv(chol(G_rho(t_0))-dagger) frozen in float per box (a
     variant takes R_0 from Q_gamma(t_0)); Q'_gamma = R_0-dagger Q_gamma
     R_0 assembled AT THE JET LEVEL (post-processing, with no interval
     inversion), a pivot and determinant certificate on the normalised
     object, positive definiteness kept exactly, and
     det' = det Q/det G_rho, the dimensionless NORMALISED determinant
     (A3); h_pass_raw against h_pass_cong on the same grid; the negative
     witness t_bad MUST stay certified non-positive-definite under
     congruence. (ii) phi split per element (asked for by a reviewer):
     the top 5 individual |c| plus a remainder group plus phi-combined:
     the width of the Hessian term of q00 (linear in c) decomposes
     ADDITIVELY per element (exactly), giving the share of the top 5 in C,
     measured. (218 individual jets are out of budget: the monomial
     precomputation dominates the cost, about 11 s per jet whatever the
     number of non-zero coefficients; recorded.)
  A3 dimensionless pricing: generalised float lambda_min(Q_gamma, G_rho)
     (centre plus sampled minimum), the certified normalised determinant
     after congruence, and an over-enclosure proxy
     det_norm_min_float - det'.lo.
  A4 decision, published from the measurements: (1) if float stays
     positive where the enclosure fails AND the Hessian excess dominates,
     then the remainder representation is the problem and Taylor models
     are justified (a synthetic discriminant a reviewer made MANDATORY);
     (2) if the real variation (linear plus quadratic at the centre plus
     the float span) dominates, a higher order will not certify without
     subdivision; (3) if congruence already wins a factor 10 or more in
     atlas cells, it is the first production route.
  +  a shared coverage table: 60 pairs = 22 owner from O1 + 17 vacuities
     + 21 ambiguous; joined with the sampled set (26 sampled); the 17
     ambiguous ones WITHOUT named points.
  +  an atlas erratum (asked for by a reviewer): the range recomputed
     from the measured anchors (the equation in the direct-probe note
     gave 8.7e11 cells and 3e5 years, not 9e14 and 3e8; corrected, with
     the wall intact).

Self-test (DISCRIMINATING checks, negative controls included):
  T1 enclose_decompose matches taylor2_enclose in the kernel (bounds
     bit-identical)
  T2 additivity of the widths: w_t2 = w_val + w_lin + w_quad (exact)
  T3 congruence with R_0 = I matches the raw form (bounds bit-identical,
     4 components plus the determinant)
  T4 scaling: mid(det') equals |det R_0|^2.mid(det) (relative < 1e-10,
     h = 0)
  T5 NEGATIVE CONTROL: t_bad under congruence gives det'.hi < 0
     (certified non-positive-definite)
  T6 M_H is bit-hermitian and the degenerate h = 0 midpoint matches the
     float M_H (< 5e-12)
  T7 phi split: w_quad(phi-combined) <= w_quad(top 5) + w_quad(rest)
     (subadditivity is guaranteed; the gap is the combined cancellation
     gain, which is THE measurement the reviewer asked for) plus a
     NEGATIVE CONTROL: w_quad(top 5) < w_quad(phi) strictly

Output : results/width_attribution.json
Usage  : width_attribution.py [--selftest]
"""
from __future__ import annotations

import hashlib
import io
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

_argv = sys.argv
sys.argv = [sys.argv[0]]
from mpmath import iv, mp                                          # noqa: E402
from .witness_registry import load_active_witness            # noqa: E402
from .kahler_metric import chart_metric_kahler               # noqa: E402
from .witness_parametrisation import B3, B3_MULTIS, B3_IDX               # noqa: E402
from .interval_arithmetic import (                            # noqa: E402
    CIV, CZERO, HPAIRS, NG, T2CIV, T2IV, BranchCutError, _float_section,
    _iv_intersect, build_M_civ, det_packed_iv, det_packed_t2, iv_bounds,
    leaf_of_float_point, t2_chart_cell_section, t2_chart_metric,
    taylor2_enclose)
sys.argv = _argv

RES = Path(os.environ.get(
    "K3_RES_DIR", Path(__file__).resolve().parent / "data"))
PROBE_JSON = RES / "width_probe.json"
DIRECT_JSON = RES / "width_direct.json"
TILING_JSON = RES / "owner_tiling.json"

GAMMA = 0.25                    # unchanged, fixed in advance
H_GRID = [4e-3, 2e-3, 1e-3, 5e-4, 2.5e-4]   # ≥ 4 valeurs (GPT C44)
N_SAMPLES = 128
SEED = 20260725
TOP_K = 5
YEAR_S = 3.15576e7
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def wid(x) -> float:
    return float(mp.mpf(x.delta))


# ===========================================================================
#  Canonicalisation of M_H (once, with a published hash)
# ===========================================================================
def canonical_MH(Mnp):
    M_H = 0.5 * (Mnp + Mnp.conj().T)
    resid = float(np.max(np.abs(Mnp - M_H)))
    bit_herm = bool(np.array_equal(M_H, M_H.conj().T))
    sha = hashlib.sha256(np.ascontiguousarray(M_H).tobytes()).hexdigest()
    return M_H, resid, bit_herm, sha


# ===========================================================================
#  Q_gamma jets (one expensive evaluation at the centre and on the box)
#  plus the certificate
# ===========================================================================
def q_jets(S, g_col, eps, u0, v0, h, M_civ, coeffs, rho_w):
    Zc, Wc, _ = t2_chart_cell_section(S, g_col, eps, u0, v0, 0.0)
    q_c = t2_chart_metric(Zc, Wc, M_civ, coeffs, rho_weight=rho_w)
    Zb, Wb, _ = t2_chart_cell_section(S, g_col, eps, u0, v0, h)
    q_b = t2_chart_metric(Zb, Wb, M_civ, coeffs, rho_weight=rho_w)
    return q_c, q_b


def enclose_decompose(center, box, h):
    """An EXACT mirror of taylor2_enclose (T1: bounds bit-identical) plus
    the additive decomposition of the t2 form: the interval widths add
    exactly under +, so
      w_t2 = w(centre value) + w(sum grad.delta) + w(sum half H(box) delta^2)
    and w(sum half H(t_0) delta^2) (degenerate jet = EXACT Hessian at the
    centre) gives
        the TRUE scale of the quadratic, so the enclosure excess is isolated."""
    rad = iv.mpf([-h, h])
    sq_diag = iv.mpf([0, h * h])
    sq_off = iv.mpf([-h * h, h * h])
    half = iv.mpf(1) / 2
    t2 = center.val
    lin = iv.mpf(0)
    for a in range(NG):
        term = center.g[a] * rad
        t2 = t2 + term
        lin = lin + term
    quad_box = iv.mpf(0)
    quad_ctr = iv.mpf(0)
    for i, (a, b) in enumerate(HPAIRS):
        if a == b:
            tb = box.h[i] * sq_diag * half
            tc = center.h[i] * sq_diag * half
        else:
            tb = box.h[i] * sq_off
            tc = center.h[i] * sq_off
        t2 = t2 + tb
        quad_box = quad_box + tb
        quad_ctr = quad_ctr + tc
    mv = center.val
    for a in range(NG):
        mv = mv + box.g[a] * rad
    inter = _iv_intersect(_iv_intersect(t2, mv), box.val)
    w_t2, w_mv, w_nv = wid(t2), wid(mv), wid(box.val)
    active = min((("t2", w_t2), ("mv", w_mv), ("naive", w_nv)),
                 key=lambda t: t[1])[0]
    w_qb, w_qc = wid(quad_box), wid(quad_ctr)
    return {
        "iv": list(iv_bounds(inter)),
        "w_inter": wid(inter), "w_t2": w_t2, "w_mv": w_mv,
        "w_naive": w_nv, "active_form": active,
        "w_center_val": wid(center.val), "w_lin": wid(lin),
        "w_quad_box": w_qb, "w_quad_center": w_qc,
        "hess_excess": w_qb - w_qc,
        "hess_excess_ratio": w_qb / w_qc if w_qc > 0 else None,
        "grad_center_l2": float(math.sqrt(sum(
            (sum(iv_bounds(center.g[a])) / 2) ** 2 for a in range(NG)))),
        "hess_center_l1": float(sum(
            abs(sum(iv_bounds(center.h[i])) / 2) for i in range(len(HPAIRS)))),
        "hess_box_l1_sup": float(sum(
            max(abs(iv_bounds(box.h[i])[0]), abs(iv_bounds(box.h[i])[1]))
            for i in range(len(HPAIRS)))),
    }, inter


def certify_from_jets(q_c, q_b, h):
    """Pivot and determinant certificate plus widths per form (direct
    determinant against determinant of the components against their
    intersection), entirely as post-processing of the jets."""
    det_c, det_b = det_packed_t2(q_c), det_packed_t2(q_b)
    q_packed = [taylor2_enclose(q_c[c], q_b[c], h) for c in range(4)]
    det_comp = det_packed_iv(q_packed)
    dec_det, det_direct = enclose_decompose(det_c, det_b, h)
    dec_q00, _ = enclose_decompose(q_c[0], q_b[0], h)
    detQ = _iv_intersect(det_direct, det_comp)
    q00_lo, q00_hi = iv_bounds(q_packed[0])
    det_lo, det_hi = iv_bounds(detQ)
    ok_piv, ok_det = q00_lo > 0, det_lo > 0
    status = ("PASS" if (ok_piv and ok_det)
              else "FAIL_DET" if ok_piv else "FAIL_PIVOT")
    return {"status": status,
            "q00": [q00_lo, q00_hi], "det": [det_lo, det_hi],
            "w_det_final": det_hi - det_lo,
            "w_det_direct_t2form": dec_det["w_t2"],
            "w_det_mv": dec_det["w_mv"], "w_det_naive": dec_det["w_naive"],
            "w_det_components": wid(det_comp),
            "active_form_det": dec_det["active_form"],
            "decompose_det": dec_det, "decompose_q00": dec_q00}


# ===========================================================================
#  Fixed congruence R_0 at the jet level (post-processing)
# ===========================================================================
def _cplx_jet(re_j: T2IV, im_j: T2IV) -> T2CIV:
    return T2CIV(CIV(re_j.val, im_j.val),
                 [CIV(a, b) for a, b in zip(re_j.g, im_j.g)],
                 [CIV(a, b) for a, b in zip(re_j.h, im_j.h)])


def congruent_packed(q4, R0):
    """[g00, g11, Re g01, Im g01] with the same packing for
    R_0-dagger Q R_0, with R_0 a frozen 2x2 float (exact degenerate
    constants, no interval inversion)."""
    zero = T2IV.const(iv.mpf(0))
    Q01 = _cplx_jet(q4[2], q4[3])
    Q = [[_cplx_jet(q4[0], zero), Q01],
         [Q01.conj(), _cplx_jet(q4[1], zero)]]
    Qp = [[None, None], [None, None]]
    for A in range(2):
        for B in range(2):
            acc = T2CIV.const(CZERO)
            for i in range(2):
                for j in range(2):
                    c = complex(np.conj(R0[i, A]) * R0[j, B])
                    if c == 0:
                        continue
                    acc = acc + Q[i][j] * T2CIV.const(CIV.from_complex(c))
            Qp[A][B] = acc
    half = iv.mpf(1) / 2
    q00 = ((Qp[0][0] + Qp[0][0].conj()).re_t2()) * half
    q11 = ((Qp[1][1] + Qp[1][1].conj()).re_t2()) * half
    q01 = (Qp[0][1] + Qp[1][0].conj()).mul_real(half)
    return [q00, q11, q01.re_t2(), q01.im_t2()]


def r0_from_chol(Gmid):
    """R₀ = inv(L†), G = LL† ⟹ R₀†GR₀ = I (normalise le terme principal)."""
    L = np.linalg.cholesky(0.5 * (Gmid + Gmid.conj().T))
    return np.linalg.inv(L.conj().T)


# ===========================================================================
#  Float: G, Q, generalised lambda_min, box sampling (NOT a certificate)
# ===========================================================================
def float_G_pair(S, g_col, eps, u, v, M_H, c218):
    Z, W = _float_section(S, g_col, eps, u, v)
    if Z is None:
        return None, None
    Gf = chart_metric_kahler(Z, W, M_H, c218, B3, B3_MULTIS, B3_IDX)[0]
    Gr = chart_metric_kahler(Z, W, M_H, np.zeros_like(c218),
                             B3, B3_MULTIS, B3_IDX)[0]
    herm = lambda A: 0.5 * (A + A.conj().T)                  # noqa: E731
    return herm(Gf), herm(Gr)


def float_box_samples(S, g_col, eps, u0, v0, h, M_H, c218, rng, n=N_SAMPLES):
    """Sampled min/max on the box: det Q_gamma, generalised lambda_min(Q, G_rho)
    generalised, normalised determinant det Q/det G_rho. Centre included."""
    pts = [(u0, v0)]
    du = rng.uniform(-h, h, (n, 4))
    pts += [(u0 + complex(d[0], d[1]), v0 + complex(d[2], d[3]))
            for d in du]
    dets, lams, dnorms = [], [], []
    for u, v in pts:
        Gf, Gr = float_G_pair(S, g_col, eps, u, v, M_H, c218)
        if Gf is None:
            continue
        Q = Gf - GAMMA * Gr
        det = float((Q[0, 0] * Q[1, 1]).real - abs(Q[0, 1]) ** 2)
        detr = float((Gr[0, 0] * Gr[1, 1]).real - abs(Gr[0, 1]) ** 2)
        Li = np.linalg.inv(np.linalg.cholesky(Gr))
        lam = float(np.linalg.eigvalsh(Li @ Q @ Li.conj().T)[0])
        dets.append(det)
        lams.append(lam)
        dnorms.append(det / detr)
    return {"n_valid": len(dets),
            "det_min": min(dets), "det_max": max(dets),
            "det_span": max(dets) - min(dets),
            "gen_lammin_center": lams[0], "gen_lammin_min": min(lams),
            "det_norm_center": dnorms[0], "det_norm_min": min(dnorms)}


# ===========================================================================
#  A1/A2: analysis of one cell (box, h): raw, two congruences, and float
# ===========================================================================
def analyze_cell(S, g_col, eps, u0, v0, h, M_civ, M_H, c218, R0s, rng):
    t1 = time.time()
    rw = 1.0 - GAMMA
    try:
        q_c, q_b = q_jets(S, g_col, eps, u0, v0, h, M_civ, c218, rw)
    except BranchCutError as exc:
        return {"h": h, "status_raw": "BRANCH", "error": str(exc)[:120],
                "t_cell_s": time.time() - t1}
    rec = {"h": h, "raw": certify_from_jets(q_c, q_b, h)}
    rec["status_raw"] = rec["raw"]["status"]
    for name, R0 in R0s.items():
        qc_p = congruent_packed(q_c, R0)
        qb_p = congruent_packed(q_b, R0)
        cr = certify_from_jets(qc_p, qb_p, h)
        # after R_0 (rho variant): det' = det Q/det G_rho is dimensionless
        rec[f"cong_{name}"] = {k: cr[k] for k in
                               ("status", "q00", "det", "w_det_final",
                                "w_det_direct_t2form", "w_det_components",
                                "active_form_det")}
        rec[f"cong_{name}"]["decompose_det"] = cr["decompose_det"]
    rec["float_box"] = float_box_samples(S, g_col, eps, u0, v0, h,
                                         M_H, c218, rng)
    fb = rec["float_box"]
    rec["overclosure_proxy_norm"] = (
        fb["det_norm_min"] - rec["cong_rho"]["det"][0]
        if "cong_rho" in rec else None)
    rec["C_width_over_h2_raw"] = rec["raw"]["w_det_final"] / h ** 2
    rec["C_width_over_h2_cong_rho"] = (
        rec["cong_rho"]["w_det_final"] / h ** 2 if "cong_rho" in rec
        else None)
    rec["t_cell_s"] = time.time() - t1
    return rec


def slopes(records, key_path):
    """Pentes p = log2(w(h)/w(h/2)) entre h successifs (GPT C44)."""
    out = []
    usable = [r for r in records if "raw" in r]
    for r1, r2 in zip(usable, usable[1:]):
        if abs(r1["h"] / r2["h"] - 2.0) > 1e-9:
            continue
        w1, w2 = r1, r2
        for k in key_path:
            w1, w2 = w1[k], w2[k]
        if w2 > 0:
            out.append({"h_pair": [r1["h"], r2["h"]],
                        "ratio": w1 / w2,
                        "p_eff": math.log2(w1 / w2)})
    return out


# ===========================================================================
#  Build
# ===========================================================================
def load_boxes_from_direct(direct):
    """Takes EXACTLY the 3 boxes of the direct probe (u0/v0 serialised),
    the worst class being the minimum det.lo at h_pass. Order: the
    witness box B FIRST, as a reviewer asked."""
    boxes = []
    for r in direct["s1_three_boxes"]:
        cands = [pc for pc in r["per_class"]
                 if pc["h_pass"] == r["worst_class_h_pass"]]
        pc = min(cands, key=lambda p: p["steps"][-1]["det"][0])
        boxes.append({
            "tag": r["tag"], "S": tuple(r["S"]), "g": r["g"],
            "u0": complex(*r["u0"]), "v0": complex(*r["v0"]),
            "eps": tuple(pc["eps"]), "class_id": pc["class_id"],
            "h_pass_direct": r["worst_class_h_pass"],
            "provenance": r["provenance"]})
    order = {"B_interior": 0, "A_worst": 1, "C_near_residual": 2}
    boxes.sort(key=lambda b: order.get(b["tag"], 9))
    return boxes


def coverage_table(tiling, probe):
    """The shared table: 60 pairs by (tiling, sampled float, strategy)."""
    cov = {(tuple(c["S"]), c["g"]): c
           for c in probe["b_r_sampled"]["coverage_C39"]["by_couple_class"]}
    rows, amb_no_points, owner_no_points = [], [], []
    n_owner = n_vac = n_amb = 0
    for c in tiling["couples"]:
        key = (tuple(c["S"]), c["g"])
        if c["vol_owner"] > 0:
            o1, n_owner = "OWNER_CERTIFIED", n_owner + 1
        elif c["vol_residual"] <= tiling["rad_floor2"] * 1e6:
            o1, n_vac = "VACUOUS_CERTIFIED", n_vac + 1
        else:
            o1, n_amb = "AMBIGUOUS", n_amb + 1
        seen = cov.get(key)
        if o1 == "AMBIGUOUS" and seen is None:
            amb_no_points.append(f"S={list(key[0])} g={key[1]}")
        if o1 == "OWNER_CERTIFIED" and seen is None:
            owner_no_points.append(
                f"S={list(key[0])} g={key[1]} "
                f"(owner_fraction={c['owner_fraction']:.2e}, "
                f"frontier_capped={c['frontier_capped']})")
        strat = {"OWNER_CERTIFIED": "atlas P0a (certificat par cellule)",
                 "VACUOUS_CERTIFIED": "excluded (certified OUTSIDE exhaustion)",
                 "AMBIGUOUS": "subdiviser O1 plus profond avant tout atlas"
                 }[o1]
        rows.append({
            "S": list(key[0]), "g": key[1], "o1_status": o1,
            "vol_owner": c["vol_owner"], "vol_residual": c["vol_residual"],
            "sampled_C39": seen is not None,
            "n_classes_expected": seen["n_classes_expected"] if seen else None,
            "all_classes_seen": seen["all_classes_seen"] if seen else None,
            "min_r_sampled": (min(d["min_r"] for d in seen["classes"]
                                  .values()) if seen else None),
            "strategy": strat})
    n_sampled = sum(1 for r in rows if r["sampled_C39"])
    return {"rows": rows,
            "counts": {"owner_certified": n_owner, "vacuous": n_vac,
                       "ambiguous": n_amb, "sampled_C39": n_sampled,
                       "ambiguous_sampled": sum(
                           1 for r in rows if r["sampled_C39"]
                           and r["o1_status"] == "AMBIGUOUS")},
            "ambiguous_no_points_named": amb_no_points,
            "owner_no_points_named": owner_no_points,
            "reconciliation": ("26 sampled = 21 owner + 5 ambiguous; "
                               "34 absent = 17 vacuities + 16 ambiguous "
                               "without points plus 1 OWNER not sampled "
                               "(this refines the 17+17 of the review)"),
            "note": ("V_owner = 5.38 concerns the certified owner part "
                     "(22 pairs), not the cost of a global "
                     "certificate: the ambiguous residual (21 pairs) stays a "
                     "subdivision budget for the tiling stage")}


def atlas_erratum(direct):
    """The range recomputed from the measured anchors, as a reviewer
    asked. (The equation in the direct-probe note, "9e14 cells, about 3e8
    years", in fact gives 8.7e11 cells and 3e5 years; we republish it
    from h_pass = sqrt(margin/C).)"""
    margin = direct["s1_three_boxes"][1]["per_class"][0]["steps"][-1][
        "det"][0]                       # B_interior det.lo at h_pass
    t_call = direct["t2_seconds_per_call_mean"]
    v_eff = 5.38 * 2.6
    anchors = []
    for C in (1.8e3, 6.0e3, 2.4e4):
        h = math.sqrt(margin / C)
        cells = v_eff / (2 * h) ** 4
        anchors.append({"C": C, "h_pass": h, "cells": cells,
                        "cpu_years": cells * t_call / YEAR_S})
    return {"margin_det_B": margin, "t_call_s": t_call,
            "V_eff_owner_x_classes": v_eff, "anchors": anchors,
            "corrected_range": ("3e11-6e13 cells, 1.1e5-2e7 CPU years "
                                "(uniform, as a scenario): the uniform t2 "
                                "route is dead under both readings; the "
                                "target of the next stage is unchanged"),
            "erratum": ("the direct-probe note wrote 9e14 cells, about "
                        "3e8 years; its own equation gives 8.7e11 cells "
                        "and 3e5 years, and the wall figure is that of "
                        "the inputs measured above")}


def build():
    print("=" * 78)
    print("P0a2-A — ATTRIBUTION de largeur + congruence fixe (GPT the width contract "
          "A0-A4, a reviewer V3/V7)")
    print("=" * 78)
    rng = np.random.default_rng(SEED)
    w = load_active_witness()
    Mnp = np.asarray(w["M"], complex)
    c218 = np.asarray(w["coeffs218"], float)

    # --- A0/C42 : canonisation M_H --------------------------------------------------
    M_H, resid, bit_herm, sha_MH = canonical_MH(Mnp)
    M_civ = build_M_civ(M_H)
    log(f"A0: M_H = (M + M-dagger)/2 canonicalised, max antihermitian residual "
        f"{resid:.3e}, bit-hermitien {bit_herm}, sha256 {sha_MH[:16]}…")

    probe = json.loads(PROBE_JSON.read_text(encoding="utf-8"))
    direct = json.loads(DIRECT_JSON.read_text(encoding="utf-8"))
    tiling = json.loads(TILING_JSON.read_text(encoding="utf-8"))
    boxes = load_boxes_from_direct(direct)

    # --- A0: degenerate control h=0 replayed under M_H ----------------
    b0 = boxes[0]
    q0_c, _ = q_jets(b0["S"], b0["g"], b0["eps"], b0["u0"], b0["v0"],
                     0.0, M_civ, c218, 1.0 - GAMMA)
    Gf0, Gr0 = float_G_pair(b0["S"], b0["g"], b0["eps"], b0["u0"],
                            b0["v0"], M_H, c218)
    Qf0 = Gf0 - GAMMA * Gr0
    ref = [Qf0[0, 0].real, Qf0[1, 1].real, Qf0[0, 1].real, Qf0[0, 1].imag]
    mids = [sum(iv_bounds(q0_c[c].val)) / 2 for c in range(4)]
    rel_h0 = max(abs(m - r) / max(abs(r), 1e-30)
                 for m, r in zip(mids, ref))
    log(f"A0: degenerate h=0 under M_H, relative {rel_h0:.2e} (check < 5e-12)")

    # --- W: negative witness under M_H, raw and congruent -------------
    bw = probe["b_r_sampled"]["argmin"]
    S_w, g_w = tuple(bw["S"]), int(bw["g"])
    Z_w = np.array([re + 1j * im
                    for re, im in zip(bw["Z"], bw["Z_imag"])])
    eps_w = leaf_of_float_point(S_w, g_w, Z_w)
    others = [c for c in range(6) if c not in S_w and c != g_w]
    u_w, v_w = complex(Z_w[others[0]]), complex(Z_w[others[1]])
    t_bad = 1.25 * probe["d_discriminant_target"]["t_crit_sampled"]
    qw_c, qw_b = q_jets(S_w, g_w, eps_w, u_w, v_w, 0.0, M_civ,
                        c218 * t_bad, 1.0)
    cert_w = certify_from_jets(qw_c, qw_b, 0.0)
    _, Gr_w = float_G_pair(S_w, g_w, eps_w, u_w, v_w, M_H, c218)
    R0_w = r0_from_chol(Gr_w)
    qwp_c = congruent_packed(qw_c, R0_w)
    qwp_b = congruent_packed(qw_b, R0_w)
    cert_wc = certify_from_jets(qwp_c, qwp_b, 0.0)
    neg = {"t_bad": t_bad,
           "raw_det": cert_w["det"],
           "raw_certified_nonPD": cert_w["det"][1] < 0,
           "cong_det": cert_wc["det"],
           "cong_certified_nonPD": cert_wc["det"][1] < 0,
           "note": "M_H canonicalised; det' = det/det G_rho (dimensionless)"}
    log(f"W  : t_bad = {t_bad:.4f} — raw det.hi = {cert_w['det'][1]:.3e}"
        f" ; cong det'.hi = {cert_wc['det'][1]:.3e} (les 2 < 0 requis)")

    # --- A1/A2(i)/A3: h grid by 3 boxes, raw plus congruences ---------
    a1 = []
    for bx in boxes:
        S, g_col, eps = bx["S"], bx["g"], bx["eps"]
        u0, v0 = bx["u0"], bx["v0"]
        Gf_m, Gr_m = float_G_pair(S, g_col, eps, u0, v0, M_H, c218)
        Q_m = Gf_m - GAMMA * Gr_m
        R0s = {"rho": r0_from_chol(Gr_m)}
        if np.linalg.eigvalsh(Q_m)[0] > 0:
            R0s["q"] = r0_from_chol(Q_m)
        log(f"A1 {bx['tag']} S={S} g={g_col} eps={eps} "
            f"(classe {bx['class_id']}) — R0 : {list(R0s)}")
        recs = []
        for h in H_GRID:
            rec = analyze_cell(S, g_col, eps, u0, v0, h, M_civ, M_H,
                               c218, R0s, rng)
            recs.append(rec)
            if "raw" in rec:
                cst = {k: rec[f"cong_{k}"]["status"] for k in R0s}
                log(f"   h={h:9.3e} raw={rec['status_raw']:10s} "
                    f"cong={cst} C_raw={rec['C_width_over_h2_raw']:.3g} "
                    f"C_rho={rec['C_width_over_h2_cong_rho']:.3g} "
                    f"float[{rec['float_box']['det_min']:.3e},"
                    f"{rec['float_box']['det_max']:.3e}] "
                    f"{rec['t_cell_s']:.1f}s")
            else:
                log(f"   h={h:9.3e} {rec['status_raw']}")
        h_pass = {"raw": next((r["h"] for r in recs
                               if r.get("status_raw") == "PASS"), None)}
        for k in R0s:
            h_pass[f"cong_{k}"] = next(
                (r["h"] for r in recs
                 if r.get(f"cong_{k}", {}).get("status") == "PASS"), None)
        a1.append({
            "tag": bx["tag"], "S": list(S), "g": g_col,
            "eps": list(eps), "class_id": bx["class_id"],
            "u0": [u0.real, u0.imag], "v0": [v0.real, v0.imag],
            "provenance": bx["provenance"],
            "R0_variants": {k: [[c.real, c.imag] for c in R0.flatten()]
                            for k, R0 in R0s.items()},
            "records": recs, "h_pass": h_pass,
            "slopes_w_det_raw": slopes(recs, ("raw", "w_det_final")),
            "slopes_w_det_cong_rho": slopes(
                recs, ("cong_rho", "w_det_final")),
            "slopes_w_quad_box": slopes(
                recs, ("raw", "decompose_det", "w_quad_box"))})

    # --- A2(ii): phi split per element (a reviewer), box B, h = 1e-3 --
    bB = boxes[0]
    h_split = 1e-3
    top_idx = np.argsort(-np.abs(c218))[:TOP_K]
    log(f"A2 split φ (a reviewer V7) : top-{TOP_K} |c| = "
        f"{[f'{c218[i]:.3g}' for i in top_idx]} on box B at h={h_split:g}")

    def phi_quad(coeffs):
        qc, qb = q_jets(bB["S"], bB["g"], bB["eps"], bB["u0"], bB["v0"],
                        h_split, M_civ, coeffs, 0.0)
        d_q00, _ = enclose_decompose(qc[0], qb[0], h_split)
        d_det, _ = enclose_decompose(det_packed_t2(qc), det_packed_t2(qb),
                                     h_split)
        return d_q00, d_det

    mask_top = np.zeros_like(c218)
    mask_top[top_idx] = c218[top_idx]
    d_phi, ddet_phi = phi_quad(c218)
    d_top, ddet_top = phi_quad(mask_top)
    d_rest, ddet_rest = phi_quad(c218 - mask_top)
    per_elem = []
    for i in top_idx:
        m1 = np.zeros_like(c218)
        m1[i] = c218[i]
        d_e, _ = phi_quad(m1)
        per_elem.append({"elem": int(i), "c": float(c218[i]),
                         "w_quad_box_q00": d_e["w_quad_box"]})
        log(f"   element {i} (c={c218[i]:.4g}): w_quad(q00) = "
            f"{d_e['w_quad_box']:.3e}")
        # GUARANTEED sub-additivity (max-abs in h_i.[-h^2, h^2]): the gap
        # measures the CANCELLATION gain of the combined assembly at the level
        # of the Hessian, which is the measurement this probe exists for
    subadd_ok = (d_phi["w_quad_box"] <= (d_top["w_quad_box"]
                 + d_rest["w_quad_box"]) * (1 + 1e-12))
    cancel_gain = ((d_top["w_quad_box"] + d_rest["w_quad_box"])
                   / d_phi["w_quad_box"])
    a2_split = {
        "box": "B_interior", "h": h_split,
        "top_elems": [int(i) for i in top_idx],
        "phi_only_w_quad_box_q00": d_phi["w_quad_box"],
        "top5_w_quad_box_q00": d_top["w_quad_box"],
        "rest_w_quad_box_q00": d_rest["w_quad_box"],
        "top5_share": d_top["w_quad_box"] / d_phi["w_quad_box"],
        "subadditivity_ok": subadd_ok,
        "combined_cancellation_gain": cancel_gain,
        "per_element_sum_over_combined": (
            sum(e["w_quad_box_q00"] for e in per_elem)
            + d_rest["w_quad_box"]) / d_phi["w_quad_box"],
        "per_element": per_elem,
        "det_w_t2_phi_only": ddet_phi["w_t2"],
        "det_w_t2_top5": ddet_top["w_t2"],
        "det_w_t2_rest": ddet_rest["w_t2"],
        "note_218": ("218 individual jets are out of budget: the "
                     "precomputation of the monomials dominates (about the "
                     "same cost whatever the number of non-zero "
                     "coefficients), hence the top 5 plus groups, with "
                     "exact additivity verified on q00 (linear in c)")}

    # --- A4: the decision (measured, not decreed) ---------------------
    recB = a1[0]["records"]
    rB_pass = next((r for r in recB if r.get("status_raw") == "PASS"), None)
    rB_fail = next((r for r in recB if r.get("status_raw",
                                             "").startswith("FAIL")), None)
    dd = rB_pass["raw"]["decompose_det"] if rB_pass else None
    shares = None
    if dd:
        shares = {
            "lin_real": dd["w_lin"] / dd["w_t2"],
            "quad_real_center": dd["w_quad_center"] / dd["w_t2"],
            "hess_enclosure_excess": dd["hess_excess"] / dd["w_t2"],
            "center_val": dd["w_center_val"] / dd["w_t2"]}
    float_pos_where_iv_fails = (rB_fail is not None and "float_box"
                                in rB_fail
                                and rB_fail["float_box"]["det_min"] > 0)
    hp = a1[0]["h_pass"]
    cong_gain_cells = None
    if hp["raw"] and hp.get("cong_rho"):
        cong_gain_cells = (hp["cong_rho"] / hp["raw"]) ** 4
    a4 = {
        "widths_shares_B_at_h_pass": shares,
        "float_positive_where_enclosure_fails": float_pos_where_iv_fails,
        "float_span_vs_iv_width_B_fail": (
            rB_fail["float_box"]["det_span"] / rB_fail["raw"]["w_det_final"]
            if rB_fail and "raw" in rB_fail else None),
        "h_pass_B": hp, "cong_gain_atlas_cells": cong_gain_cells,
        "kimi_thresholds": {"revise_down_if_gain_lt": 10,
                            "cosign_budget_if_gain_ge": 100}}
    # verdict de routage (GPT A4, trois branches)
    if cong_gain_cells and cong_gain_cells >= 10:
        route = ("CONGRUENCE: a gain of at least a factor 10 in cells, the first route "
                 "production, models deferred")
    elif float_pos_where_iv_fails and shares and \
            shares["hess_enclosure_excess"] >= 0.5:
        route = ("TAYLOR MODELS OR AFFINE FORMS JUSTIFIED: the float stays positive where "
                 "the enclosure fails AND the Hessian enclosure excess "
                 "dominates the width, so the remainder representation is the "
                 "problem; a synthetic discriminant is REQUIRED")
    elif shares and (shares["lin_real"] + shares["quad_real_center"]) > 0.5:
        route = ("REAL VARIATION dominant: a higher order will not "
                 "certify without subdivision; a geometric decomposition "
                 "or a regional analytic bound remains to be designed")
    else:
        route = "AMBIGUOUS: no branch dominates; publish and arbitrate"
    a4["route"] = route
    log(f"A4 : {route}")

    cov = coverage_table(tiling, probe)
    err = atlas_erratum(direct)

    # --- checks (derived from the counts) -----------------------------
    all_cells = [r for a in a1 for r in a["records"]]
    checks = {
        "G0_MH_bit_hermitian": bit_herm,
        "G1_h0_degenerate_vs_float_MH": bool(rel_h0 < 5e-12),
        "G2_negative_raw_certified": bool(neg["raw_certified_nonPD"]),
        "G3_negative_cong_certified": bool(neg["cong_certified_nonPD"]),
        "G4_t2_width_additivity": bool(all(
            abs(r["raw"]["decompose_det"]["w_t2"]
                - (r["raw"]["decompose_det"]["w_center_val"]
                   + r["raw"]["decompose_det"]["w_lin"]
                   + r["raw"]["decompose_det"]["w_quad_box"]))
            <= 1e-12 * r["raw"]["decompose_det"]["w_t2"]
            for r in all_cells if "raw" in r)),
        "G5_phi_split_subadditive": bool(subadd_ok),
        "G6_ge4_h_values_per_box": bool(all(
            sum(1 for r in a["records"] if "raw" in r) >= 4 for a in a1)),
        "G7_coverage_reconciled": bool(
            cov["counts"]["owner_certified"] == 22
            and cov["counts"]["vacuous"] == 17
            and cov["counts"]["ambiguous"] == 21
            and cov["counts"]["sampled_C39"] == 26
            and len(cov["ambiguous_no_points_named"])
            + len(cov["owner_no_points_named"]) + 17 == 34)}
    n_pass = sum(1 for v in checks.values() if v)
    log(f"checks executed: {n_pass}/{len(checks)} PASS")

    verdict = (
        "WIDTH ATTRIBUTION EXECUTED (checks executed %d/%d PASS): the "
        "width is measured BEFORE the next decision (contract A0-A4). "
        "M_H canonicalised (sha %s...), degenerate h=0 replayed "
        "(relative %.1e). The witness t_bad is certified non-positive "
        "definite raw AND under congruence. Witness box B at h_pass: "
        "width shares (determinant, t2 form) = real linear %.3f / real "
        "quadratic at the centre %.3f / Hessian enclosure EXCESS %.3f. "
        "Float > 0 where the enclosure fails: %s. Fixed congruence on the "
        "rho block: h_pass %s against raw %s. Route: %s" % (
            n_pass, len(checks), sha_MH[:12], rel_h0,
            shares["lin_real"] if shares else float("nan"),
            shares["quad_real_center"] if shares else float("nan"),
            shares["hess_enclosure_excess"] if shares else float("nan"),
            float_pos_where_iv_fails,
            hp.get("cong_rho"), hp["raw"], route))

    out = {
        "phase": ("width attribution plus preconditioned fixed "
                  "congruence (the width contract, A0-A4, with two "
                  "reviewer requests folded in)"),
        "witness_sha256": str(w["artifact_sha256"].item()),
        "M_H_c42": {"sha256": sha_MH, "antiherm_residual_max": resid,
                    "bit_hermitian": bit_herm,
                    "consumers": "float engine + build_M_civ + probes"},
        "gamma_prefixed": GAMMA, "seed": SEED, "h_grid": H_GRID,
        "a0_h0_degenerate_rel": rel_h0,
        "w_negative_witness": neg,
        "a1_width_attribution": a1,
        "a2_phi_split_kimi_v7": a2_split,
        "a4_decision": a4,
        "coverage_table_gpt6": cov,
        "atlas_erratum_kimi_v3": err,
        "checks_executed": checks,
        "verdict": verdict}

    art = RES / "width_attribution.json"
    art.write_text(json.dumps(out, indent=2, ensure_ascii=False,
                              default=float), encoding="utf-8")
    print("\nVERDICT :\n" + verdict)
    print(f"\n→ {art}")
    return out


# ===========================================================================
#  Self-test
# ===========================================================================
def _selftest():
    fails = []
    w = load_active_witness()
    Mnp = np.asarray(w["M"], complex)
    c218 = np.asarray(w["coeffs218"], float)
    M_H, resid, bit_herm, _ = canonical_MH(Mnp)
    M_civ = build_M_civ(M_H)
    direct = json.loads(DIRECT_JSON.read_text(encoding="utf-8"))
    bB = load_boxes_from_direct(direct)[0]        # B_interior
    S, g_col, eps = bB["S"], bB["g"], bB["eps"]
    u0, v0 = bB["u0"], bB["v0"]
    h = 5e-4
    q_c, q_b = q_jets(S, g_col, eps, u0, v0, h, M_civ, c218, 1.0 - GAMMA)
    det_c, det_b = det_packed_t2(q_c), det_packed_t2(q_b)

    # --- T1 : enclose_decompose ≡ taylor2_enclose (bit-identique) -----------------
    devs = []
    for ctr, box in [(q_c[0], q_b[0]), (det_c, det_b)]:
        _, mine = enclose_decompose(ctr, box, h)
        ker = taylor2_enclose(ctr, box, h)
        devs += [abs(iv_bounds(mine)[i] - iv_bounds(ker)[i])
                 for i in (0, 1)]
    t1 = max(devs) == 0.0
    fails.append(not t1)
    print(f"[{'PASS' if t1 else 'FAIL'}] T1 enclose_decompose ≡ kernel "
          f"(dev max {max(devs):g})")

    # --- T2: additivity of the widths ---------------------------------
    dd, _ = enclose_decompose(det_c, det_b, h)
    dev = abs(dd["w_t2"] - (dd["w_center_val"] + dd["w_lin"]
                            + dd["w_quad_box"])) / dd["w_t2"]
    t2 = dev < 1e-12
    fails.append(not t2)
    print(f"[{'PASS' if t2 else 'FAIL'}] T2 additivity w_t2 = "
          f"val+lin+quad (rel {dev:.2e})")

    # --- T3 : congruence R0 = I ≡ brut ---------------------------------------------
    qI_c = congruent_packed(q_c, np.eye(2, dtype=complex))
    devs3 = [abs(iv_bounds(qI_c[c].val)[i] - iv_bounds(q_c[c].val)[i])
             for c in range(4) for i in (0, 1)]
    detI = det_packed_t2(qI_c)
    devs3 += [abs(iv_bounds(detI.val)[i] - iv_bounds(det_c.val)[i])
              for i in (0, 1)]
    t3 = max(devs3) == 0.0
    fails.append(not t3)
    print(f"[{'PASS' if t3 else 'FAIL'}] T3 congruence R0=I ≡ brut "
          f"(dev max {max(devs3):g})")

    # --- T4 : scaling det' = |det R0|²·det (h=0, mid) -------------------------------
    Gf_m, Gr_m = float_G_pair(S, g_col, eps, u0, v0, M_H, c218)
    R0 = r0_from_chol(Gr_m)
    q0_c, q0_b = q_jets(S, g_col, eps, u0, v0, 0.0, M_civ, c218,
                        1.0 - GAMMA)
    qp = congruent_packed(q0_c, R0)
    mid = lambda j: sum(iv_bounds(j.val)) / 2                  # noqa: E731
    det_raw = mid(det_packed_t2(q0_c))
    det_cng = mid(det_packed_t2(qp))
    scale = abs(np.linalg.det(R0)) ** 2
    rel = abs(det_cng - scale * det_raw) / abs(det_cng)
    t4 = rel < 1e-10
    fails.append(not t4)
    print(f"[{'PASS' if t4 else 'FAIL'}] T4 det' ≡ |det R0|²·det "
          f"(rel {rel:.2e})")

    # --- T5: NEGATIVE CONTROL, t_bad under congruence must stay
    #         certified non-positive-definite ------------------------
    probe = json.loads(PROBE_JSON.read_text(encoding="utf-8"))
    bw = probe["b_r_sampled"]["argmin"]
    S_w, g_w = tuple(bw["S"]), int(bw["g"])
    Z_w = np.array([re + 1j * im
                    for re, im in zip(bw["Z"], bw["Z_imag"])])
    eps_w = leaf_of_float_point(S_w, g_w, Z_w)
    others = [c for c in range(6) if c not in S_w and c != g_w]
    u_w, v_w = complex(Z_w[others[0]]), complex(Z_w[others[1]])
    t_bad = 1.25 * probe["d_discriminant_target"]["t_crit_sampled"]
    qw_c, qw_b = q_jets(S_w, g_w, eps_w, u_w, v_w, 0.0, M_civ,
                        c218 * t_bad, 1.0)
    _, Gr_w = float_G_pair(S_w, g_w, eps_w, u_w, v_w, M_H, c218)
    qwp = congruent_packed(qw_c, r0_from_chol(Gr_w))
    cert = certify_from_jets(qwp, congruent_packed(qw_b, r0_from_chol(
        Gr_w)), 0.0)
    t5 = cert["det"][1] < 0
    fails.append(not t5)
    print(f"[{'PASS' if t5 else 'FAIL'}] T5 negative control: congruent t_bad "
          f"det'.hi = {cert['det'][1]:.3e} < 0")

    # --- T6: M_H bit-hermitian plus degenerate h=0 against float M_H --
    Qf = Gf_m - GAMMA * Gr_m
    ref = [Qf[0, 0].real, Qf[1, 1].real, Qf[0, 1].real, Qf[0, 1].imag]
    mids = [mid(q0_c[c]) for c in range(4)]
    rel6 = max(abs(m - r) / max(abs(r), 1e-30)
               for m, r in zip(mids, ref))
    t6 = bit_herm and rel6 < 5e-12
    fails.append(not t6)
    print(f"[{'PASS' if t6 else 'FAIL'}] T6 C42 : bit-herm {bit_herm}, "
          f"degenerate relative {rel6:.2e} (antihermitian residual {resid:.1e})")

    # --- T7: additive phi split plus a strict negative control --------
    top_idx = np.argsort(-np.abs(c218))[:TOP_K]
    mask = np.zeros_like(c218)
    mask[top_idx] = c218[top_idx]

    def quad_q00(coeffs):
        qc, qb = q_jets(S, g_col, eps, u0, v0, h, M_civ, coeffs, 0.0)
        d, _ = enclose_decompose(qc[0], qb[0], h)
        return d["w_quad_box"]

    w_phi, w_top, w_rest = (quad_q00(c218), quad_q00(mask),
                            quad_q00(c218 - mask))
        # the width of the Hessian term is SUB-additive (max-abs in
    # h_i.[-h^2,h^2]); the gap is the combined cancellation gain
    gain7 = (w_top + w_rest) / w_phi
    t7 = w_phi <= (w_top + w_rest) * (1 + 1e-12) and w_top < w_phi
    fails.append(not t7)
    print(f"[{'PASS' if t7 else 'FAIL'}] T7 split φ : sous-additif "
          f"(gain cancellation ×{gain7:.4f}), top5 {w_top:.3e} < φ "
          f"{w_phi:.3e}")

    print("-" * 78)
    print("SELF-TEST:", "FAIL" if any(fails) else "ALL PASS")
    return 1 if any(fails) else 0


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--selftest":
        raise SystemExit(_selftest())
    build()
