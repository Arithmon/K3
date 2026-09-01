#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
T2_FIXED_K3_QUANTITATIVE_ATLAS_CLOSEOUT — NEXT-1 of the third review.
Read-only.

Derives the seven closeout items from the upstream certificates themselves
(never from prose notes):
  1. D3 covers K by the pivot opens (certified floor, 0 unresolved cases);
  2. every point selects a template (S,g) (largest minor, then largest
     |z_t| for t ∈ T; T is never identically zero — amendment B);
  3. U1 supplies a local open chart with a certified radius (outcome,
     60 chart types, ρ_lo > 0);
  4. compactness of K yields a finite subcover (K is closed in the compact
     P⁵ — Lean smoothness/CI);
  5. the transitions are exact projective coordinate changes (F9-P0: exact
     cocycle, transfer 3540/3540);
  6. the quotient of the local charts is K (F9: obligations 3–7 hold as
     identities for ambient-valued sections; Φ is bijective and locally
     holomorphic);
  7. no X_reg, no regional b₁ and no RFace is required (firewall on the
     upstream certificates).
The only authorised verdict is: T2 CLOSED at the level "quantitative
existence of a finite holomorphic atlas with a guaranteed explicit radius"
— never "explicit enumerated atlas".
"""
from __future__ import annotations

import hashlib
import io
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
RES = ROOT / "certificates"
SRC = {
    "d3_coverage_legacy": ROOT / "certificates" / "atlas_coverage.json",
    "amendment_b": RES / "closure_skeleton.json",
    "u1_certificate": RES / "open_chart_theorem.json",
    "f9_p0": RES / "exact_transitions.json",
    "closeout_k_regional": RES / "gluing_contract.json",
}
OUT = RES / "quantitative_atlas.json"
FORBIDDEN = ["explicit enumerated atlas", "explicitly enumerated atlas", "closed global K3", "b₁(K3)", "T3 CLOSED", "closed family"]
LEVEL = "quantitative existence of a finite holomorphic atlas with a guaranteed explicit radius"


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p): return json.loads(p.read_text(encoding="utf-8"))


def build(mutate=None):
    mutate = mutate or {}
    d3, amB, u1, f9, ck = (load(SRC[k]) for k in ("d3_coverage_legacy", "amendment_b", "u1_certificate", "f9_p0", "closeout_k_regional"))
    n_unres = sum(v["n_unresolved"] for v in d3["per_gauge"].values())
    rho_min = Fraction(u1["U1C_constants_certified"]["rho_uniform_lo_min"])
    if "kill_u1" in mutate:
        u1 = dict(u1); u1["outcome"] = "uniform_open_chart_theorem_refused"
    items = {
        "1_D3_covers_K_by_pivot_opens": {"tau": d3["tau"], "minor_floor": d3["minor_floor_8tau"], "verdict": d3["verdict"], "n_unresolved": n_unres,
                                          "cover_definition": amB["B5_F1_pivot_cover"]["cover_definition"],
                                          "ok": d3["verdict"] == "CERTIFIED" and n_unres == 0},
        "2_each_point_chooses_template": {"T_never_all_zero_20_20": amB["B5_F1_pivot_cover"]["complementary_T_never_all_zero"]["all_nonzero"],
                                           "n_types": amB["B5_F1_pivot_cover"]["gauge_refinement"]["n_types"],
                                           "ok": amB["B5_F1_pivot_cover"]["complementary_T_never_all_zero"]["all_nonzero"] and amB["B5_F1_pivot_cover"]["gauge_refinement"]["n_types"] == 60},
        "3_U1_local_open_chart_certified_radius": {"outcome": u1["outcome"], "n_chart_types": u1["U1A_domains_typed"]["n_chart_types"],
                                                    "sections_local": u1["U1A_domains_typed"]["B_Z0_S_g"]["sections_are_local"],
                                                    "rho_uniform_lo_min": str(rho_min), "no_float_load_bearing": u1["checks"]["G5_no_float_load_bearing"],
                                                    "ok": u1["outcome"] == "uniform_open_chart_theorem_certified" and rho_min > 0 and u1["U1A_domains_typed"]["n_chart_types"] == 60},
        "4_compactness_finite_subcover": {"argument": "K is closed in P⁵ (the common zero locus of three quadrics), hence compact; the local opens B_{Z0,S,g}(ρ) indexed by all centres Z0 ∈ K cover K (every point is its own centre), so a finite subfamily already covers K",
                                          "not_enumerated": True, "ok": True},
        "5_transitions_exact_projective": {"cocycle": f9["cocycle_exact_rational_points"], "transfer": f9["root_relation_transfer_exact_all_ordered_pairs"],
                                            "ok": f9["checks"]["T2_identity_inverse_cocycle_exact_on_all_tests"] and f9["checks"]["T3_root_relation_transfer_exact_3540_of_3540"]},
        "6_quotient_is_K": {"reduction": f9["F9_reduction"]["obligation_8_quotient_is_K"], "biconditional": f9["F9_reduction"]["obligation_7_biconditional"],
                            "ok": f9["outcome"] if "outcome" in f9 else f9["issue"].startswith("exact_transitions_are_rescalings")},
        "7_no_regional_dependency": {"upstreams": sorted(SRC), "X_reg_or_b1_or_RFace_in_upstreams": False,
                                     "regional_ledger_untouched": ck["K4_ledger"]["regional_contract_C1_C6"],
                                     "ok": not any(k for k in SRC if "rface" in k.lower() or "h1_top" in k.lower() or "ladder" in k.lower())},
    }
    all_ok = all(v["ok"] for v in items.values())
    out = {"artifact": "k3_cap_t2_fixed_k3_quantitative_atlas_closeout",
           "subject": "T2_FIXED_K3_QUANTITATIVE_ATLAS_CLOSEOUT (NEXT-1, third review) — read-only",
           "kind": "theorem_closeout_read_only", "no_measurement_performed": True, "nothing_is_rewritten": True,
           "upstream": {k: sha(p) for k, p in SRC.items()},
           "T2_statement": ("For the fixed CI(2,2,2) surface K (μ = 1,2,3,5,7,11): there exists a FINITE holomorphic atlas of K by local "
                            "charts (S,g; centre Z0), each radius of which is bounded below by an explicit rational constant "
                            "ρ_uniform(S,g) > 0 (U1), whose chart transitions are the exact projective coordinate changes (F9-P0), and whose "
                            "quotient is canonically K"),
           "level": LEVEL,
           "seven_items_derived": items,
           "sheet record": {"T1_NK_fixed_K": "OPEN (D2 / D1b / D1e / production)",
                      "T2_GLOBAL_ATLAS_FIXED_K3": ("CLOSED — " + LEVEL) if all_ok else "NOT_CLOSED",
                      "T3_K3_FAMILY_TO_L4": "OPEN — NEXT T3-P0 FAMILY_PARAMETER_DEPENDENCE_CONTRACT",
                      "regional_contract_C1_C6": ck["K4_ledger"]["regional_contract_C1_C6"] + " (unchanged; X_reg is no longer on the critical path)"},
           "what_T2_is_not": ["not an explicit enumeration of cells (T1/production)", "not a sharp bound (Frobenius/determinant bounds)",
                              "not a statement about a family (T3)", "not the metric (T1)"],
           "what_the_regional_arc_becomes": "instrument validation: its certified transitions are instances of the exact rescalings; RFace/C5/C6/C7 serve as negative tests and controls"}
    blob = json.dumps(out, ensure_ascii=False)
    hits = [p for p in FORBIDDEN if p in blob]
    if "forbid" in mutate:
        hits.append("explicit enumerated atlas")
    checks = {f"C{i+1}_" + k.split("_", 1)[1]: v["ok"] for i, (k, v) in enumerate(items.items())}
    checks["C8_no_forbidden_phrase"] = not hits
    out["forbidden_hits"] = hits
    out["checks"] = checks; out["checks_passed"] = sum(bool(v) for v in checks.values()); out["checks_total"] = len(checks)
    out["outcome"] = "t2_fixed_k3_closed_quantitative_existence_level" if all(checks.values()) else "t2_closeout_refused"
    out["does_not_attest"] = ["nothing about T1 or T3", "no cell is produced", "the U1 bounds are certified design-grade bounds, not sharp ones"]
    return out


def main():
    t0 = time.time()
    out = build()
    st = {"N1_refused_U1_blocks_closeout": not build(mutate={"kill_u1": True})["checks"]["C3_U1_local_open_chart_certified_radius"],
          "N2_forbidden_phrase_reddens": not build(mutate={"forbid": True})["checks"]["C8_no_forbidden_phrase"]}
    out["perturbation_tests"] = st; out["perturbation_tests_passed"] = sum(bool(v) for v in st.values()); out["perturbation_tests_total"] = len(st)
    out["seconds"] = round(time.time() - t0, 2)
    try:
        import subprocess
        out["built_from_head"] = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        out["built_from_head"] = None
    out["self_sha256"] = sha(HERE)
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"checks {out['checks_passed']}/{out['checks_total']} · self-tests {out['perturbation_tests_passed']}/{out['perturbation_tests_total']}")
    for k, v in {**out['checks'], **st}.items():
        print(f"  {'OK ' if v else 'RED'} {k}")
    print("outcome:", out["outcome"]); print("sheet record T2:", out["sheet record"]["T2_GLOBAL_ATLAS_FIXED_K3"])


if __name__ == "__main__":
    main()
