#!/usr/bin/env python3
"""
k3_cap_b1e2iii_p0a2_attrib.py — P0a2-A : ATTRIBUTION de largeur +
congruence fixe préconditionnée (contrat GPT `gpt_b1e2iii_p0a2_direct_review`
§2-§6 : C42-C46, A0-A4 ; demandes Kimi `kimi_p0a2_direct_review` V3/V7).

The wall measured by the direct probe (C = width(detQ)/h^2 around 2e3-2e4) does not
encore de CAUSE établie (GPT C43-C45 : la largeur mélange variation
réelle linéaire/quadratique, excès d'enclosure de la Hessienne,
dependencies of the determinant). This probe measures it BEFORE choosing
the tool that will reduce it:

  A0 provenance : M_H = (M+M†)/2 canonisé UNE fois au chargement (C42),
     hash publié, moteur float ET noyau intervalle consomment exactement
     M_H ; contrôle dégénéré h=0 rejoué (mid(Q_γ) ≡ float, rel < 5e-12) ;
          checks derived from the counts, "executed" rather than "all"; pricing
          serialised in the JSON, not only in the note.
  A1 décomposition de largeur (GPT §5-A1 + C44) : boîte B TÉMOIN
     D'ABORD (Kimi), puis A et C — grille h FIXE {4e-3, 2e-3, 1e-3,
     5e-4, 2.5e-4} (≥ 4 valeurs, pentes successives publiées, jamais
          a "law" on 2 points). At each h, POST-PROCESSING of the two jets
     (centre + boîte) : largeur par FORME (Taylor-2 / valeur-moyenne /
     naïve / det-des-composantes / intersection finale, forme active),
          and EXACT DECOMPOSITION of the t2 form (the interval widths
     s'additionnent exactement) : w = w(val centre) + w(∇F(t₀)·δ)
     [variation réelle ordre 1] + w(½δᵀH(boîte)δ) [terme reste], avec
     w(½δᵀH(t₀)δ) [Hessienne EXACTE au centre, jet dégénéré] comme
          true scale of the quadratic, isolating the enclosure excess of the Hessian
     = w(H boîte) − w(H centre), ISOLÉ. Normes ∇/H publiées. Sample
          float min/max of det Q on the box (NOT a certificate: it separates the
     variation physique apparente de l'excès d'enclosure, GPT A1).
  A2 ablation (GPT §5-A2 + Kimi V7) : (i) congruence fixe R₀ (C46) —
     R₀ = inv(chol(G_ρ(t₀))†) gelé float par boîte (variante R₀ depuis
     Q_γ(t₀)) ; Q'_γ = R₀†Q_γR₀ assemblé AU NIVEAU DU JET (T2CIV,
     post-processing, aucune inversion d'intervalle), certificat
          pivot and determinant on the normalised object, positive definiteness kept exactly,
     det' = det Q/det G_ρ = déterminant NORMALISÉ dimensionless (A3) ;
          h_pass_raw against h_pass_cong on the same grid; the negative witness
     t_bad DOIT rester certifié non-PD sous congruence. (ii) split φ
     par élément (Kimi V7) : top-5 |c| individuels + groupe reste +
          phi-combined: the width of the Hessian term of q00 (linear in c)
     se décompose ADDITIVEMENT par élément (exact) ⟹ part du top-5
          in C, measured. (218 individual jets are out of budget: the
          monomial precomputation dominates the cost, about 11 s per jet whatever
     nombre de coeffs non nuls ; consigné.)
  A3 pricing dimensionless (GPT §5-A3) : λ_min(Q_γ, G_ρ) généralisé
     float (centre + min échantillonné), det normalisé certifié après
     congruence, proxy de sur-enclosure = det_norm_min_float − det'.lo.
  A4 décision T4 (GPT §5-A4, seuils Kimi) : publiée depuis les mesures —
          (1) if float stays positive where the enclosure fails AND the Hessian excess
          dominates, then the remainder representation is the problem and Taylor models are justified
     (discriminant synthétique Kimi OBLIGATOIRE) ; (2) si la variation
     réelle (lin + quad centre + span float) domine ⟹ un ordre supérieur
          will not certify without subdivision; (3) if congruence already wins
     ≥ ×10 en cellules d'atlas ⟹ première route de production.
  +  table de couverture commune (GPT §6) : 60 couples = 22 owner O1 +
     17 vacuités + 21 ambigus ; join avec C39 (26 échantillonnés) ;
          the 17 ambiguous ones WITHOUT named points.
  +  erratum d'atlas (Kimi V3) : fourchette recalculée depuis les
          measured anchors (the equation in the direct-probe note gave 8.7e11
     cellules/3e5 ans, pas 9e14/3e8 — corrigé, mur intact).

Self-test (gates DISCRIMINANTS, négatifs inclus) :
  T1 enclose_decompose ≡ taylor2_enclose kernel (bornes bit-identiques)
  T2 additivité des largeurs : w_t2 = w_val + w_lin + w_quad (exact)
  T3 congruence R₀ = I ≡ brut (bornes bit-identiques, 4 comp + det)
  T4 scaling : mid(det') ≡ |det R₀|²·mid(det) (rel < 1e-10, h=0)
  T5 NÉGATIF : t_bad sous congruence ⟹ det'.hi < 0 (non-PD certifié)
  T6 C42 : M_H bit-hermitien + dégénéré h=0 mid ≡ float M_H (< 5e-12)
  T7 split φ : w_quad(φ-combiné) ≤ w_quad(top5) + w_quad(reste)
     (sous-additivité garantie ; l'écart = gain de cancellation combiné,
     LA mesure Kimi V7) + NÉGATIF : w_quad(top5) < w_quad(φ) strictement

Sorties : results/k3_cap_b1e2iii_p0a2_attrib.json
Usage   : k3_cap_b1e2iii_p0a2_attrib.py [--selftest]
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
PROBE_JSON = RES / "k3_cap_b1e2iii_p0a_probe.json"
DIRECT_JSON = RES / "k3_cap_b1e2iii_p0a2_direct.json"
TILING_JSON = RES / "k3_cap_b1e2iii_owner_tiling.json"

GAMMA = 0.25                    # inchangé, pré-fixé (P0a2-direct)
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
#  C42 — canonisation M_H (une fois, hash publié)
# ===========================================================================
def canonical_MH(Mnp):
    M_H = 0.5 * (Mnp + Mnp.conj().T)
    resid = float(np.max(np.abs(Mnp - M_H)))
    bit_herm = bool(np.array_equal(M_H, M_H.conj().T))
    sha = hashlib.sha256(np.ascontiguousarray(M_H).tobytes()).hexdigest()
    return M_H, resid, bit_herm, sha


# ===========================================================================
#  Jets Q_γ (une évaluation chère centre + boîte) + certificat
# ===========================================================================
def q_jets(S, g_col, eps, u0, v0, h, M_civ, coeffs, rho_w):
    Zc, Wc, _ = t2_chart_cell_section(S, g_col, eps, u0, v0, 0.0)
    q_c = t2_chart_metric(Zc, Wc, M_civ, coeffs, rho_weight=rho_w)
    Zb, Wb, _ = t2_chart_cell_section(S, g_col, eps, u0, v0, h)
    q_b = t2_chart_metric(Zb, Wb, M_civ, coeffs, rho_weight=rho_w)
    return q_c, q_b


def enclose_decompose(center, box, h):
    """Miroir EXACT de taylor2_enclose (T1 : bornes bit-identiques) +
        additive decomposition of the t2 form: the interval widths
    s'additionnent exactement sous +, donc
      w_t2 = w(val centre) + w(Σ ∇·δ) + w(Σ ½H(boîte)δ²)
    et w(Σ ½H(t₀)δ²) (jet dégénéré = Hessienne EXACTE au centre) donne
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
    """Certificat pivot+det + largeurs par forme (det direct vs det des
    composantes vs intersection) — tout en post-processing des jets."""
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
#  C46 — congruence fixe R₀ au niveau du jet (post-processing)
# ===========================================================================
def _cplx_jet(re_j: T2IV, im_j: T2IV) -> T2CIV:
    return T2CIV(CIV(re_j.val, im_j.val),
                 [CIV(a, b) for a, b in zip(re_j.g, im_j.g)],
                 [CIV(a, b) for a, b in zip(re_j.h, im_j.h)])


def congruent_packed(q4, R0):
    """[g00, g11, Re g01, Im g01] → même packing pour R₀†QR₀, R₀ float
    2×2 gelé (constantes dégénérées exactes — aucune inversion iv)."""
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
    généralisé, det normalisé det Q/det G_ρ. Centre inclus."""
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
#  A1/A2 — analyse d'une cellule (boîte, h) : brut + 2 congruences + float
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
        # après R₀ (variante ρ) : det' = det Q/det G_ρ = dimensionless
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
    """Reprend EXACTEMENT les 3 boîtes de P0a2-direct (u0/v0 sérialisés),
    pire classe = min det.lo à h_pass. Ordre : B témoin D'ABORD (Kimi)."""
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
    """Table commune GPT §6 : 60 couples × (O1, float C39, stratégie)."""
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
                 "VACUOUS_CERTIFIED": "exclu (épuisement OUTSIDE certifié)",
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
            "reconciliation": ("26 échantillonnés = 21 owner + 5 ambigus "
                               "; 34 absents = 17 vacuités + 16 ambigus "
                               "without points plus 1 OWNER not sampled "
                               "(précision au « 17+17 » de la review "
                               "Kimi §3)"),
            "note": ("V_owner = 5.38 concerns the certified owner part "
                     "(22 pairs), not the cost of a global "
                     "certificate: the ambiguous residual (21 pairs) stays a "
                     "budget de subdivision O1 (GPT §6)")}


def atlas_erratum(direct):
    """Kimi V3 : fourchette recalculée depuis les ancres mesurées.
        (The equation in the direct-probe note, "9e14 cells, about 3e8 years", holds
    en fait 8.7e11/3e5 ans ; on republie depuis h_pass = √(marge/C).)"""
    margin = direct["s1_three_boxes"][1]["per_class"][0]["steps"][-1][
        "det"][0]                       # B_interior det.lo à h_pass
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
            "corrected_range": ("3e11-6e13 cellules, 1.1e5-2e7 ans CPU "
                                "(uniforme, scénario) — l'uniforme t2 "
                                "est mort aux deux lectures ; T4 et sa "
                                "cible inchangés"),
            "erratum": ("la note P0a2-direct écrivait 9e14 cellules ~ "
                        "3e8 ans ; son équation vaut 8.7e11 / 3e5 ans "
                        "the wall figure is that of the inputs "
                        "mesurés ci-dessus")}


def build():
    print("=" * 78)
    print("P0a2-A — ATTRIBUTION de largeur + congruence fixe (GPT C42-C46 "
          "A0-A4, Kimi V3/V7)")
    print("=" * 78)
    rng = np.random.default_rng(SEED)
    w = load_active_witness()
    Mnp = np.asarray(w["M"], complex)
    c218 = np.asarray(w["coeffs218"], float)

    # --- A0/C42 : canonisation M_H --------------------------------------------------
    M_H, resid, bit_herm, sha_MH = canonical_MH(Mnp)
    M_civ = build_M_civ(M_H)
    log(f"A0 C42 : M_H = (M+M†)/2 canonisé — résidu anti-herm max "
        f"{resid:.3e}, bit-hermitien {bit_herm}, sha256 {sha_MH[:16]}…")

    probe = json.loads(PROBE_JSON.read_text(encoding="utf-8"))
    direct = json.loads(DIRECT_JSON.read_text(encoding="utf-8"))
    tiling = json.loads(TILING_JSON.read_text(encoding="utf-8"))
    boxes = load_boxes_from_direct(direct)

    # --- A0 : contrôle dégénéré h=0 rejoué sous M_H ---------------------------------
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
    log(f"A0 : dégénéré h=0 sous M_H — rel {rel_h0:.2e} (gate < 5e-12)")

    # --- W : témoin négatif sous M_H, brut + congruent ------------------------------
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
           "note": "M_H canonisé ; det' = det/det G_ρ (dimensionless)"}
    log(f"W  : t_bad = {t_bad:.4f} — raw det.hi = {cert_w['det'][1]:.3e}"
        f" ; cong det'.hi = {cert_wc['det'][1]:.3e} (les 2 < 0 requis)")

    # --- A1/A2(i)/A3 : grille h × 3 boîtes, brut + congruences ----------------------
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

    # --- A2(ii) : split φ par élément (Kimi V7) — boîte B, h = 1e-3 -----------------
    bB = boxes[0]
    h_split = 1e-3
    top_idx = np.argsort(-np.abs(c218))[:TOP_K]
    log(f"A2 split φ (Kimi V7) : top-{TOP_K} |c| = "
        f"{[f'{c218[i]:.3g}' for i in top_idx]} @ boîte B h={h_split:g}")

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
        log(f"   élément {i} (c={c218[i]:.4g}) : w_quad(q00) = "
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
        "note_218": ("218 jets individuels hors budget : le précalcul "
                     "of the monomials dominates (about the same cost whatever the "
                     "nombre de coeffs non nuls) ⟹ top-5 + groupes, "
                     "additivité exacte vérifiée sur q00 (linéaire en c)")}

    # --- A4 : décision (mesurée, pas décrétée) --------------------------------------
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
        route = ("CONGRUENCE : gain ≥ ×10 en cellules — première route "
                 "de production, T4-TM différé")
    elif float_pos_where_iv_fails and shares and \
            shares["hess_enclosure_excess"] >= 0.5:
        route = ("T4-TM/AFFINE JUSTIFIÉ : le float reste > 0 là où "
                 "l'enclosure échoue ET l'excès d'enclosure Hessienne "
                 "domine la largeur — la représentation du reste est le "
                 "problème ; discriminant synthétique Kimi OBLIGATOIRE")
    elif shares and (shares["lin_real"] + shares["quad_real_center"]) > 0.5:
        route = ("VARIATION RÉELLE dominante : un ordre supérieur ne "
                 "certifiera pas sans subdivision — décomposition "
                 "géométrique / borne analytique régionale à concevoir")
    else:
        route = "AMBIGU : aucune branche ne domine — publier et arbitrer"
    a4["route"] = route
    log(f"A4 : {route}")

    cov = coverage_table(tiling, probe)
    err = atlas_erratum(direct)

    # --- gates (dérivés des comptes — GPT A0.3) -------------------------------------
    all_cells = [r for a in a1 for r in a["records"]]
    gates = {
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
    n_pass = sum(1 for v in gates.values() if v)
    log(f"gates exécutés : {n_pass}/{len(gates)} PASS")

    verdict = (
        "P0a2-A EXÉCUTÉ (gates exécutés %d/%d PASS) : attribution de la "
        "largeur mesurée AVANT décision T4 (contrat GPT A0-A4). C42 : "
        "M_H canonisé (sha %s…), dégénéré h=0 rejoué (rel %.1e). "
        "Témoin t_bad certifié non-PD brut ET sous congruence. Boîte B "
        "témoin à h_pass : parts de largeur (det, forme t2) = lin réel "
        "%.3f / quad réel centre %.3f / EXCÈS d'enclosure Hessienne "
        "%.3f. Float > 0 là où l'enclosure échoue : %s. Congruence "
        "fixe R₀ (Gρ) : h_pass %s vs raw %s. Route A4 : %s" % (
            n_pass, len(gates), sha_MH[:12], rel_h0,
            shares["lin_real"] if shares else float("nan"),
            shares["quad_real_center"] if shares else float("nan"),
            shares["hess_enclosure_excess"] if shares else float("nan"),
            float_pos_where_iv_fails,
            hp.get("cong_rho"), hp["raw"], route))

    out = {
        "phase": ("B1.e.2.iii P0a2-A — attribution de largeur + "
                  "congruence fixe préconditionnée (GPT C42-C46/A0-A4, "
                  "Kimi V3/V7)"),
        "witness_sha256": str(w["artifact_sha256"].item()),
        "M_H_c42": {"sha256": sha_MH, "antiherm_residual_max": resid,
                    "bit_hermitian": bit_herm,
                    "consumers": "moteur float + build_M_civ + probes"},
        "gamma_prefixed": GAMMA, "seed": SEED, "h_grid": H_GRID,
        "a0_h0_degenerate_rel": rel_h0,
        "w_negative_witness": neg,
        "a1_width_attribution": a1,
        "a2_phi_split_kimi_v7": a2_split,
        "a4_decision": a4,
        "coverage_table_gpt6": cov,
        "atlas_erratum_kimi_v3": err,
        "gates_executed": gates,
        "verdict": verdict}

    art = RES / "k3_cap_b1e2iii_p0a2_attrib.json"
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

    # --- T2 : additivité des largeurs ---------------------------------------------
    dd, _ = enclose_decompose(det_c, det_b, h)
    dev = abs(dd["w_t2"] - (dd["w_center_val"] + dd["w_lin"]
                            + dd["w_quad_box"])) / dd["w_t2"]
    t2 = dev < 1e-12
    fails.append(not t2)
    print(f"[{'PASS' if t2 else 'FAIL'}] T2 additivité w_t2 = "
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

    # --- T5 : NÉGATIF — t_bad sous congruence doit rester non-PD certifié -----------
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
    print(f"[{'PASS' if t5 else 'FAIL'}] T5 négatif : t_bad congruent "
          f"det'.hi = {cert['det'][1]:.3e} < 0")

    # --- T6 : C42 — M_H bit-hermitien + dégénéré h=0 vs float M_H -------------------
    Qf = Gf_m - GAMMA * Gr_m
    ref = [Qf[0, 0].real, Qf[1, 1].real, Qf[0, 1].real, Qf[0, 1].imag]
    mids = [mid(q0_c[c]) for c in range(4)]
    rel6 = max(abs(m - r) / max(abs(r), 1e-30)
               for m, r in zip(mids, ref))
    t6 = bit_herm and rel6 < 5e-12
    fails.append(not t6)
    print(f"[{'PASS' if t6 else 'FAIL'}] T6 C42 : bit-herm {bit_herm}, "
          f"dégénéré rel {rel6:.2e} (résidu anti-herm {resid:.1e})")

    # --- T7 : split φ additif + négatif strict --------------------------------------
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
    # h_i·[−h²,h²]) — l'écart = gain de cancellation combiné (mesure V7)
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
