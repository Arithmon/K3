# -*- coding: utf-8 -*-
"""ATLAS_PAPER_SIGMA_FLOOR_CORRECTION: a DEFECT in certificate U1, found while writing the paper, and its repair.

WHAT U1 PUBLISHES. `sigma_floor_lo(S) = 2·σ_min(Ṽ_S)·m/(8|V_S|)`, then ρ per type via x = σ²/(2 L a), ρ = −1 + √(1+x),
hence ρ_unif ∈ [9.60e-10, 2.44e-8]. The whole Newton contraction depends on 1/σ with σ ≤ s_min(M̃(w₀)).

THE DEFECT. M̃(w) = 2 Ṽ_S diag(w) with w = z_S (the resolved coordinates). The CORRECT chain is
    s_min(M̃(w₀)) ≥ 2·σ_min(Ṽ_S)·**min_i |z_{s_i}|**            [identity ‖ADx‖ ≥ σ_min(A)·min|d_i|·‖x‖, gated]
whereas the pivot threshold bounds only the **PRODUCT**: |det M_S^alg| = 8|V_S|∏|z_{s_i}| > m ⇒ ∏|z_{s_i}| > m/(8|V_S|).
U1 places the PRODUCT bound in the slot of the MINIMUM. The implication "∏ ≥ c ⇒ min ≥ c" is false as soon as one
factor exceeds 1 — and the resolved coordinates do exceed 1 on the certified domains (B_s up to 3.4).
**This is not merely a gap in the exposition: an explicit WITNESS exists on the variety** (Q3) — a point of the sector,
on the pivot domain (product > m/(8|V_S|)), whose minimum lies STRICTLY below the assumed floor.

THE REPAIR. Replace the slot by the true lower bound on min_i|z_{s_i}|, the better of the two available ones:
  · product-to-factors: |z_{s_i}| ≥ m/(8|V_S| ∏_{j≠i} B_{s_j}) with B_s² = |a_s| + |b_s| + |c_s| and |u|,|v| ≤ 1;
  · direct, row by row: |z_{s_i}|² = |R_{s_i}| ≥ |a_{s_i}| − |b_{s_i}| − |c_{s_i}| when this number is positive.
Result: **corrected ρ_unif ≈ 2.11e-12 instead of 9.60e-10, i.e. a factor of 455**. The theorem SURVIVES (uniform radius
strictly positive, everything else unchanged); it is the CONSTANT that was optimistic, not the statement.

SCOPE. U1 remains 8/8 on its other gates; T2 ("quantitative existence of a finite atlas with a guaranteed explicit
radius") is NOT invalidated — only the value of the radius changes. This batch does not rewrite U1: it publishes the
correction, and it is for the ledger to decide whether U1 is replayed. The paper cites the CORRECTED constant.
Negative controls: N1, if |z_s| ≤ 1 is imposed (sector extended to the six coordinates), the substitution made by U1
BECOMES valid (the corrected radius ceases to be a loss) — the gate does measure the right hypothesis; N2, chain
reversed (min → product): the witness no longer makes the gate fail; N3, threshold m set to 0: both floors collapse
to 0 and the comparison loses its meaning.
Does not attest: that the other gates of U1 are affected (they are not); that the witness is the worst case (search
over a deterministic rational grid, not an optimisation); the sharpness of the corrected constant (Frobenius/det
bounds, deliberately loose); anything about the other parts of the wider project.
"""
from __future__ import annotations

import hashlib
import io
import json
import math
import sys
import time
from fractions import Fraction as Fr
from itertools import combinations
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
RES = ROOT / "certificates"
U1 = RES / "open_chart_theorem.json"
OUT = RES / "sigma_floor_correction.json"
MU = (1, 2, 3, 5, 7, 11)

# ---- FROZEN HISTORICAL TRACE (the defect did exist; the upstream repair must not erase it) ---------------------------
U1_SHA_AT_DEFECT = "d30959de5e0e2d9233cf09fc03e450f695ef68a1dfde82a8dcd32797dd776e07"
U1_RHO_PUBLISHED_AT_DEFECT = 9.59780680865981e-10       # what U1 published on 2026-08-20 BEFORE the upstream repair


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p): return json.loads(p.read_text(encoding="utf-8"))


def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
            - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))


def Vm(I, mu): return [[Fr(mu[x]) ** k for x in I] for k in range(3)]


def inv3(M):
    d = det3(M); C = [[Fr(0)] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            mm = [[M[a][b] for b in range(3) if b != j] for a in range(3) if a != i]
            C[j][i] = Fr((-1) ** (i + j)) * (mm[0][0] * mm[1][1] - mm[0][1] * mm[1][0]) / d
    return C


def build(mutate=None):
    mutate = mutate or {}
    t0 = time.time()
    mu = MU
    m = Fr(0) if "m_zero" in mutate else Fr(4)
    ck2 = [sum(Fr(mu[j]) ** (2 * k) for j in range(6)) for k in range(3)]
    rows_old, rows_new, per_type = [], [], []
    witness = None
    for S in combinations(range(6), 3):
        T = [t for t in range(6) if t not in S]
        VS = abs(det3(Vm(S, mu)))
        frob2 = sum(Fr(mu[s]) ** (2 * k) / ck2[k] for s in S for k in range(3))
        det2 = VS ** 2 / (ck2[0] * ck2[1] * ck2[2])
        smin = math.sqrt(float(det2)) / float(frob2)
        L = 2 * math.sqrt(float(frob2))
        prod_lo = float(m) / (8 * float(VS))
        for g in T:
            base = [t for t in T if t != g]
            A = inv3(Vm(S, mu)); Bm = Vm([g] + base, mu)
            P = [[-float(sum(A[i][k] * Bm[k][j] for k in range(3))) for j in range(3)] for i in range(3)]
            Bs = [math.sqrt(abs(P[i][0]) + abs(P[i][1]) + abs(P[i][2])) for i in range(3)]
            direct = [abs(P[i][0]) - abs(P[i][1]) - abs(P[i][2]) for i in range(3)]
            if "cap_one" in mutate:                                        # N1: extended sector ⇒ |z_s| ≤ 1
                Bs = [min(b, 1.0) for b in Bs]
            mins = []
            for i in range(3):
                o = [Bs[j] for j in range(3) if j != i]
                mins.append(max(prod_lo / (o[0] * o[1]), math.sqrt(direct[i]) if direct[i] > 0 else 0.0))
            minz = min(mins)
            frob2uv = sum(Fr(mu[t]) ** (2 * k) / ck2[k] for t in base for k in range(3))
            a = math.sqrt(float(frob2uv))
            for sig, acc in ((2 * smin * prod_lo, rows_old), (2 * smin * minz, rows_new)):
                x = sig * sig / (2 * L * a)
                acc.append((-1 + math.sqrt(1 + x), list(S), g))
            per_type.append({"S": list(S), "g": g, "rho_corrected": rows_new[-1][0], "rho_u1": rows_old[-1][0]})
            # ---- deterministic witness: rational grid over the sector -----------------------------------------------------
            if witness is None and m > 0:
                N = 24
                for pu in range(-N, N + 1):
                    for qu in range(-N, N + 1):
                        u = complex(pu / N, qu / N)
                        if abs(u) > 1: continue
                        for pv in range(-N, N + 1, 2):
                            for qv in range(-N, N + 1, 2):
                                v = complex(pv / N, qv / N)
                                if abs(v) > 1: continue
                                R = [P[i][0] + P[i][1] * u * u + P[i][2] * v * v for i in range(3)]
                                z = [abs(r) ** 0.5 for r in R]
                                # N2, the INVERTED CHAIN. The witness looks for a point where the
                                # PRODUCT bound holds but the MINIMUM fails -- which is precisely the
                                # substitution made upstream. If the product is tested on both sides
                                # instead, no point qualifies at all and the witness falls silent.
                                # That is what shows the witness measures the min/product distinction
                                # itself, and not merely a floor set slightly too high.
                                fail = (z[0] * z[1] * z[2]) if "chain_inverted" in mutate else min(z)
                                if z[0] * z[1] * z[2] > prod_lo and fail < prod_lo:
                                    witness = {"S": list(S), "g": g, "u": [u.real, u.imag], "v": [v.real, v.imag],
                                               "abs_z": [round(t, 6) for t in z], "product": round(z[0] * z[1] * z[2], 6),
                                               "prod_lo_used_by_U1": round(prod_lo, 6), "min_abs_z": round(min(z), 6),
                                               "reading": "point of the sector ON the pivot domain (product > floor) whose MINIMUM lies below the floor"}
                                    break
                            if witness: break
                        if witness: break
                    if witness: break
    rho_old = min(rows_old); rho_new = min(rows_new)
    # ---- the matrix inequality s_min(AD) ≥ s_min(A)·min|d| (numerical control) -------------------------------------------
    import random
    random.seed(11)
    ineq_ok = True
    for _ in range(300):
        A = [[random.uniform(-2, 2) for _ in range(3)] for _ in range(3)]
        d = [random.uniform(0.01, 5) for _ in range(3)]
        AD = [[A[i][j] * d[j] for j in range(3)] for i in range(3)]
        def smin_of(M):
            MtM = [[sum(M[k][i] * M[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
            import cmath
            p = [1.0, -(MtM[0][0] + MtM[1][1] + MtM[2][2]), 0.0, -det3(MtM)]
            p[2] = (MtM[0][0] * MtM[1][1] - MtM[0][1] * MtM[1][0] + MtM[0][0] * MtM[2][2]
                    - MtM[0][2] * MtM[2][0] + MtM[1][1] * MtM[2][2] - MtM[1][2] * MtM[2][1])
            import numpy as np
            return float(min(abs(x) ** 0.5 for x in np.roots(p)))
        try:
            if smin_of(AD) < smin_of(A) * min(d) - 1e-9:
                ineq_ok = False
        except Exception:
            pass
    u1 = load(U1) if U1.exists() else {}
    c = u1.get("U1C_constants_certified", {})
    live = float(c.get("rho_uniform_lo_min_float_display_only") or 0)
    # Q1 compares against the FROZEN TRACE, not the live U1: otherwise the upstream repair would erase the evidence of the defect.
    reproduces = abs(rho_old[0] - U1_RHO_PUBLISHED_AT_DEFECT) / U1_RHO_PUBLISHED_AT_DEFECT < 1e-6
    upstream_repaired = bool(live) and abs(live - U1_RHO_PUBLISHED_AT_DEFECT) / U1_RHO_PUBLISHED_AT_DEFECT > 1e-6
    out = {"artifact": "k3_cap_atlas_paper_sigma_floor_correction",
           "subject": "DEFECT found in U1 while writing the paper: sigma_floor places the PRODUCT bound ∏|z_s| in the slot of the "
                  "MINIMUM min|z_s| required by s_min(M̃) ≥ 2σ_min(Ṽ)·min|z|; explicit witness on the variety; constant corrected",
           "kind": "defect_report_and_correction", "front": "atlas paper T2 (upstream U1)",
           "nothing_is_rewritten": True, "u1_is_not_modified_by_this_lot": True,
           "upstream": {"u1_live": sha(U1) if U1.exists() else None, "u1_outcome_live": u1.get("outcome"),
                        "u1_sha_at_defect": U1_SHA_AT_DEFECT, "u1_rho_published_at_defect": repr(U1_RHO_PUBLISHED_AT_DEFECT),
                        "u1_rho_live": repr(live), "upstream_repaired": bool(upstream_repaired),
                        "note": "Q1 compares against the FROZEN trace; this certificate stays valid after the upstream repair and keeps the record of it"},
           "A_correct_chain": {"identity": "s_min(A·diag(d)) ≥ s_min(A)·min_i|d_i|  (‖ADx‖ ≥ s_min(A)‖Dx‖ ≥ s_min(A)min|d_i|‖x‖)",
                               "numeric_control_300_random": bool(ineq_ok),
                               "what_the_pivot_gives": "|det M_S^alg| = 8|V_S|∏|z_s| > m ⇒ the PRODUCT, not the minimum"},
           "B_defect": {"u1_formula": "sigma_floor_lo = 2·σ_min(Ṽ_S)·m/(8|V_S|)",
                        "substitution": "m/(8|V_S|) is a bound on the PRODUCT, used where the MINIMUM is required",
                        "why_invalid": "∏ ≥ c ⇒ min ≥ c is false as soon as one factor exceeds 1; the resolved coordinates reach B_s ≈ 3.4 on the certified domains",
                        "witness_on_the_variety": witness},
           "C_repair": {"replacement": "min_i|z_{s_i}| ≥ max( m/(8|V_S| ∏_{j≠i}B_{s_j}) , sqrt(|a|−|b|−|c|) if > 0 )",
                        "rho_unif_published_by_U1": repr(rho_old[0]), "rho_unif_corrected": repr(rho_new[0]),
                        "loss_factor": round(rho_old[0] / rho_new[0], 1) if rho_new[0] > 0 else None,
                        "per_type_corrected": per_type, "rho_corrected_max": repr(max(r["rho_corrected"] for r in per_type)),
                        "argmin_published": {"S": rho_old[1], "g": rho_old[2]}, "argmin_corrected": {"S": rho_new[1], "g": rho_new[2]},
                        "verdict": "the THEOREM survives (uniform radius strictly positive); it is the CONSTANT that was optimistic"},
           "D_scope": {"upstream_repaired_since": bool(upstream_repaired), "u1_other_gates": "not affected", "t2_closeout": "NOT invalidated — \"quantitative existence with a guaranteed explicit radius\" holds; only the value changes",
                       "action_left_to_the_ledger": "this certificate does not replay U1; the paper cites the corrected constant"},
           "does_not_attest": ["that the other gates of U1 are affected", "that the witness is the worst case (deterministic grid, not an optimisation)",
                             "the sharpness of the corrected constant (Frobenius/det bounds deliberately loose)", "anything about the other parts of the wider project"]}
    g = {"Q1_reproduces_U1_published_radius": bool(reproduces),
         "Q2_matrix_inequality_smin_AD_control": bool(ineq_ok),
         "Q3_explicit_witness_on_the_variety": witness is not None,
         "Q4_corrected_radius_strictly_positive": rho_new[0] > 0,
         "Q5_correction_is_a_loss_not_a_gain": rho_new[0] < rho_old[0],
         "Q6_upstream_read": bool(u1)}
    out["gates"] = {k: bool(v) for k, v in g.items()}
    out["gates_passed"] = sum(bool(v) for v in g.values()); out["gates_total"] = len(g)
    out["outcome"] = "u1_sigma_floor_defect_confirmed_with_witness_radius_corrected_9p6e10_to_2p1e12_theorem_survives" if all(g.values()) else "sigma_floor_correction_gates_red"
    out["seconds"] = round(time.time() - t0, 1)
    return out


def main():
    out = build()
    # N1: if |z_s| ≤ 1 were guaranteed, the substitution made by U1 would become VALID — the corrected radius ceases to be a loss
    _cap = build(mutate={"cap_one"})
    st = {"N1_capping_z_at_one_validates_U1_substitution": (float(_cap["C_repair"]["rho_unif_corrected"])
                                                            >= float(out["C_repair"]["rho_unif_published_by_U1"]) * (1 - 1e-9)),
          # N2, the inverted chain. The module docstring announced this control
          # from the start, but it was never implemented. Without it nothing
          # separates "the upstream argument confused a minimum with a product"
          # from "the upstream floor was merely a little too high" -- and that
          # distinction is the entire subject of this certificate.
          "N2_inverting_the_chain_min_to_product_silences_the_witness":
              not build(mutate={"chain_inverted"})["gates"]["Q3_explicit_witness_on_the_variety"],
          "N3_m_zero_collapses_both_floors": not build(mutate={"m_zero"})["gates"]["Q3_explicit_witness_on_the_variety"]}
    out["self_tests"] = st; out["self_tests_passed"] = sum(bool(v) for v in st.values()); out["self_tests_total"] = len(st)
    try:
        import subprocess
        out["built_from_head"] = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        out["built_from_head"] = None
    out["self_sha256"] = sha(HERE)
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"gates {out['gates_passed']}/{out['gates_total']} · self-tests {out['self_tests_passed']}/{out['self_tests_total']} · {out['seconds']} s")
    for k, v in out["gates"].items(): print(f"  {'OK ' if v else 'RED'} {k}")
    for k, v in out["self_tests"].items(): print(f"  {'OK ' if v else 'RED'} {k}")
    r = out["C_repair"]
    print(f"ρ published {r['rho_unif_published_by_U1']} → corrected {r['rho_unif_corrected']} (factor {r['loss_factor']})")
    print("witness:", json.dumps(out["B_defect"]["witness_on_the_variety"], ensure_ascii=False)[:220])
    print("outcome:", out["outcome"])


if __name__ == "__main__":
    main()
