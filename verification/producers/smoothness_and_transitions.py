# -*- coding: utf-8 -*-
"""ATLAS_PAPER_SMOOTHNESS_COVERAGE_AND_TRANSITIONS: derives the mathematical content still left as a placeholder in the
atlas paper (§2.2 smoothness, §2.3 involutions + discriminant, §3.4 projective invariant, §4.1 coverage, §5.1 transitions).
Read-only, exact (Fraction + sympy). No heavy computation.

THE CENTRAL LEMMA (A) — "at least three nonzero coordinates", and it is used TWICE.
    If at most two coordinates of a point of X are nonzero, say z_i, z_j, the three equations read
    W·(w_i, w_j)ᵀ = 0 with W = (μ_i^k, μ_j^k)_{k=0,1,2} of rank 2 (every 2×2 minor is a Vandermonde determinant with
    distinct μ, hence nonzero — checked on all 15 pairs) ⇒ w_i = w_j = 0 ⇒ all coordinates vanish: NOT a point of P⁵.
    **Therefore at every point of X at least three coordinates are nonzero.** Consequences:
      · (§2.2 SMOOTHNESS) the Jacobian minor of the corresponding triple equals det M_S = 8 V_S z_{s₁}z_{s₂}z_{s₃} ≠ 0
        (identity (B)) ⇒ the Jacobian has rank 3 everywhere ⇒ X is smooth (and of dimension 2);
      · (§4.1 COVERAGE) det M_S ≠ 0 means exactly Z ∈ U_S ⇒ **the pivot open sets cover X**.
    The same exact statement carries both smoothness and coverage: the enumeration audit is needed only for the
    quantitative FLOOR (τ, m), not for qualitative coverage.
(B) PIVOT FACTORISATION (exact, sympy): ∂F̃_k/∂z_s = 2(μ_s^k/c_k) z_s ⇒ M_S = 2 Ṽ_S diag(z_S) ⇒
    **det M_S = 8 V_S z_{s₁}z_{s₂}z_{s₃} / (c₀c₁c₂)**: homogeneous of degree 3 in Z.
(C) PROJECTIVE INVARIANT (§3.4): the degree 3 of (B) makes **q_S(Z) := |det M_S(Z)| / ‖Z‖³ INVARIANT under Z ↦ λZ**
    (the factor |λ|³ cancels — verified symbolically); the certified radius is therefore read off a projective invariant,
    not off a representative. The sector |z_g| = max_t |z_t| is a SELECTOR of a representative, not an open set of the cover.
(D) INVOLUTIONS AND BRANCH LOCUS (§2.3): the diagonal sign matrices σ_j = diag(1,…,−1,…,1) preserve each
    F_k (the coordinates enter only through their squares); modulo the projective scalar −1 they generate a group
    of order 32, and σ₁,…,σ₅ suffice (σ₀ = −σ₁σ₂σ₃σ₄σ₅ in P⁵). The branch locus of π_S is {R_{s₁}R_{s₂}R_{s₃} = 0}:
    for the distinguished patch S = {3,4,5}, z₀ = 1, base (u,v), the discriminant is the EXACT product computed here (degree 6).
(E) EXPLICIT TRANSITIONS (§5.1): on the squares, w_{S′} = P_{S′S} w_S (rational matrices from the generator certificate);
    on the coordinates, z_{s′} = ε_{s′}√(P_{S′S}w_S)_{s′}. Two closed formulas are published: a generic elementary SWAP
    (S = {3,4,5} → S′ = {2,4,5}, a single row changes) and a GAUGE change (rescaling z_{g′}/z_g).
(F) FROM THE PRODUCT TO THE FACTORS (§4.3) — what makes the separation of the sheets UNIFORM on the certified domains.
    |det M_S^alg| > m bounds the PRODUCT |z_{s₁}z_{s₂}z_{s₃}| > m/(8|V_S|), not each factor. The sector selector gives
    |u|, |v| ≤ 1, hence |z_s|² = |R_s(u,v)| ≤ |a_s| + |b_s| + |c_s| =: B_s² (rational coefficients of the type). Therefore,
    for each solved coordinate, **|z_{s_i}| > m / (8|V_S| ∏_{j≠i} B_{s_j}) > 0**, an EXPLICIT bound, uniform over the type.
    Computed here on all 60 types: the worst one is ≈ 9.6e-4 (attained at S = {3,4,5}, gauge z₀, coordinate z₅). Two ledgers
    differing in slot s therefore give points at distance at least 2× this bound in that coordinate: the sheets merge only
    at the branch locus, which the pivot threshold keeps every chart away from.
(G) GRASSMANNIAN READING (§5.3) — why J(6,3) is not an imported piece of combinatorics.
    dF = 2·V·diag(z) is a 3×6 matrix: its row space is a point of **Gr(3,6)**, and the minor S equals
    **8 V_S z_{s₁}z_{s₂}z_{s₃} = det M_S^alg** — in other words **the pivots ARE the Plücker coordinates of the row space
    of the Jacobian** (verified on all 20 minors). The 20 solved triples are therefore exactly the standard coordinate
    charts of Gr(3,6) pulled back along p ↦ rowspace(dF(p)), and the elementary swap |S ∩ S′| = 2 is
    the standard change of chart of this Grassmannian. J(6,3) is the adjacency graph of these elementary changes of
    chart — NOT the nerve of the cover (two standard charts of Gr(3,6) always meet: the nerve would be the
    full simplex on 20 vertices).
Negative controls: N1 repeated μ ⇒ a 2×2 minor vanishes ⇒ lemma A fails (and the surface is SINGULAR); N2 mutated
factorisation identity (7 instead of 8) ⇒ (B) fails; N3 exponent 2 instead of 3 in q_S ⇒ projective invariance fails.
Does not attest: the quantitative coverage FLOOR (τ = 0.6, m = 4.8 — upstream certificate, read by U1); the certified
radius (U1); the geometry of the bridges; anything metric; anything about the other fronts.
"""
from __future__ import annotations

import hashlib
import io
import json
import sys
import time
from itertools import combinations
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
import sympy as sp

HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
RES = ROOT / "certificates"
U1 = RES / "open_chart_theorem.json"
GEN = RES / "glue_obligations.json"
OUT = RES / "smoothness_and_transitions.json"
MU_DEFAULT = (1, 2, 3, 5, 7, 11)


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p): return json.loads(p.read_text(encoding="utf-8"))


def build(mutate=None):
    mutate = mutate or {}
    t0 = time.time()
    mu = (1, 1, 3, 5, 7, 11) if "mu_repeat" in mutate else MU_DEFAULT
    Z = sp.symbols("z0:6", complex=True)
    lam = sp.Symbol("lambda", positive=True)
    # ---- (A) lemma "at least three nonzero coordinates" ---------------------------------------------------------------------
    pair_minors, all_pairs_rank2 = {}, True
    for i, j in combinations(range(6), 2):
        W = sp.Matrix([[sp.Integer(mu[i]) ** k, sp.Integer(mu[j]) ** k] for k in range(3)])
        minors = [W[[a, b], :].det() for a, b in combinations(range(3), 2)]
        ok = any(m != 0 for m in minors)
        pair_minors[f"{i},{j}"] = {"minors": [str(m) for m in minors], "rank2": bool(ok)}
        all_pairs_rank2 = all_pairs_rank2 and ok
    # ---- (B) exact factorisation of the pivot -------------------------------------------------------------------------------
    ck2 = [sum(sp.Integer(mu[j]) ** (2 * k) for j in range(6)) for k in range(3)]
    ck = [sp.sqrt(c) for c in ck2]
    def M_S(S):
        return sp.Matrix([[2 * sp.Integer(mu[s]) ** k * Z[s] / ck[k] for s in S] for k in range(3)])
    fac_ok, fac_expr = True, None
    coeff = 7 if "bad_factor" in mutate else 8
    for S in combinations(range(6), 3):
        V_S = sp.Matrix([[sp.Integer(mu[s]) ** k for s in S] for k in range(3)]).det()
        lhs = sp.simplify(M_S(S).det())
        rhs = coeff * V_S * Z[S[0]] * Z[S[1]] * Z[S[2]] / (ck[0] * ck[1] * ck[2])
        if sp.simplify(lhs - rhs) != 0:
            fac_ok = False
        if S == (3, 4, 5):
            fac_expr = sp.simplify(lhs)
    # ---- (C) projective invariant: homogeneity of degree 3 ------------------------------------------------------------------
    deg = 2 if "bad_degree" in mutate else 3
    d_at_Z = M_S((3, 4, 5)).det()
    d_at_lZ = d_at_Z.subs({Z[i]: lam * Z[i] for i in range(6)}, simultaneous=True)
    homog_ok = sp.simplify(d_at_lZ - lam ** 3 * d_at_Z) == 0
    q_invariant = sp.simplify(lam ** 3 / lam ** deg - 1) == 0
    # ---- (D) involutions + branch locus -------------------------------------------------------------------------------------
    invol_preserve = True
    for j in range(6):
        sub = {Z[i]: (-Z[i] if i == j else Z[i]) for i in range(6)}
        for k in range(3):
            Fk = sum(sp.Integer(mu[i]) ** k * Z[i] ** 2 for i in range(6))
            if sp.simplify(Fk.subs(sub, simultaneous=True) - Fk) != 0:
                invol_preserve = False
    u, v = sp.symbols("u v", complex=True)
    S_d, T_d = (3, 4, 5), (0, 1, 2)
    # NOTE: if a Vandermonde determinant degenerates (repeated μ), the radicands and the transitions DO NOT EXIST — the
    # degenerate state is published as such.
    try:
        VS = sp.Matrix([[sp.Rational(mu[s]) ** k for s in S_d] for k in range(3)])
        VT = sp.Matrix([[sp.Rational(mu[t]) ** k for t in T_d] for k in range(3)])
        wS = -VS.inv() * VT * sp.Matrix([1, u ** 2, v ** 2])
        radicands = {S_d[a]: sp.expand(wS[a]) for a in range(3)}
        discriminant = sp.expand(radicands[3] * radicands[4] * radicands[5])
        disc_degree = int(sp.Poly(discriminant, u, v).total_degree())
    except Exception:
        radicands, discriminant, disc_degree = {}, sp.Integer(0), -1
    # ---- (E) explicit transitions ---------------------------------------------------------------------------------------------
    S_a, S_b = (3, 4, 5), (2, 4, 5)
    T_a = tuple(t for t in range(6) if t not in S_a)
    try:
        Va = sp.Matrix([[sp.Rational(mu[s]) ** k for s in S_a] for k in range(3)])
        Ta = sp.Matrix([[sp.Rational(mu[t]) ** k for t in T_a] for k in range(3)])
        wa = sp.Matrix(sp.symbols("wa0 wa1 wa2"))
        wT_from = -Ta.inv() * Va * wa
        full = {s: wa[a] for a, s in enumerate(S_a)}
        full.update({t: sp.expand(wT_from[a]) for a, t in enumerate(T_a)})
        P_ba = sp.Matrix([[sp.expand(full[s2]).coeff(wa[c]) for c in range(3)] for s2 in S_b])
        ident_rows = [[sp.Integer(1) if S_b[r] == S_a[c] else sp.Integer(0) for c in range(3)] for r in range(3)]
        rows_changed = sum(1 for r in range(3) if list(P_ba[r, :]) != ident_rows[r])
        P_rows = [[str(P_ba[r, c]) for c in range(3)] for r in range(3)]
    except Exception:
        rows_changed, P_rows = -1, []
    # ---- (F) from the product to the factors: explicit lower bound on each solved coordinate --------------------------------
    import math
    M_THRESH = sp.Rational(4)                                                   # pivot threshold fixed upstream (U1)
    worst = None
    try:
        for S in combinations(range(6), 3):
            T = [t for t in range(6) if t not in S]
            VSdet = abs(sp.Matrix([[sp.Rational(mu[x]) ** k for x in S] for k in range(3)]).det())
            for gg in T:
                cols = [gg] + [t for t in T if t != gg]                          # w_T ordered as (1, u^2, v^2)
                Vc = sp.Matrix([[sp.Rational(mu[x]) ** k for x in cols] for k in range(3)])
                Pm = -sp.Matrix([[sp.Rational(mu[x]) ** k for x in S] for k in range(3)]).inv() * Vc
                Bs = [math.sqrt(float(sum(abs(Pm[i, j]) for j in range(3)))) for i in range(3)]
                prod_lo = float(M_THRESH) / (8 * float(VSdet))
                for i in range(3):
                    o = [Bs[j] for j in range(3) if j != i]
                    lo = prod_lo / (o[0] * o[1])
                    if worst is None or lo < worst[0]:
                        worst = (lo, list(S), gg, S[i])
    except Exception:
        worst = None
    # ---- (G) the pivots are the Plucker coordinates of the row space of dF --------------------------------------------------
    dF = sp.Matrix([[2 * sp.Integer(mu[j]) ** k * Z[j] for j in range(6)] for k in range(3)])
    plucker_ok = True
    for S in combinations(range(6), 3):
        V_S = sp.Matrix([[sp.Integer(mu[s]) ** k for s in S] for k in range(3)]).det()
        minor = dF[:, list(S)].det()
        if sp.simplify(minor - 8 * V_S * Z[S[0]] * Z[S[1]] * Z[S[2]]) != 0:
            plucker_ok = False
    u1 = load(U1) if U1.exists() else {}
    gen = load(GEN) if GEN.exists() else {}
    upstream_ok = str(u1.get("outcome", "")).startswith("uniform_open_chart_theorem") and str(gen.get("outcome", "")).startswith("atlas_paper_glue")
    d3 = u1.get("d3", {})
    out = {"artifact": "k3_cap_atlas_paper_smoothness_coverage_and_transitions",
           "subject": "atlas paper — lemma 'at least 3 nonzero coordinates' (smoothness AND coverage), exact factorisation of the pivot, "
                  "projective invariant q_S of degree 3, diagonal involutions + discriminant of the distinguished patch, explicit transitions",
           "kind": "exact_algebra_sympy", "front": "atlas paper", "no_measurement_performed": "exact algebra only",
           "nothing_is_rewritten": True,
           "upstream": {"u1": sha(U1) if U1.exists() else None, "u1_outcome": u1.get("outcome"),
                        "generators": sha(GEN) if GEN.exists() else None, "d3_read_from_u1": d3},
           "A_three_nonzero_lemma": {"all_15_pairs_rank2": bool(all_pairs_rank2), "pair_minors": pair_minors,
                                     "statement": "at most 2 nonzero coordinates ⇒ W·(w_i,w_j)ᵀ = 0 with W of rank 2 ⇒ all vanish: not a point of P⁵",
                                     "use_1_smoothness": "the triple with nonzero coordinates has det M_S ≠ 0 ⇒ Jacobian of rank 3 ⇒ X smooth",
                                     "use_2_coverage": "det M_S ≠ 0 ⟺ Z ∈ U_S ⇒ the pivot open sets COVER X (qualitative coverage, exact)",
                                     "what_the_audit_adds": "the quantitative FLOOR (τ, m) only — not the coverage"},
           "B_pivot_factorisation": {"identity": "det M_S = 8 V_S z_{s1}z_{s2}z_{s3} / (c0 c1 c2)", "verified_all_20_triples": bool(fac_ok),
                                     "example_S_345": str(fac_expr), "c_k_squared": [str(c) for c in ck2]},
           "C_projective_invariant": {"det_homogeneous_degree_3": bool(homog_ok),
                                      "q_S": "q_S(Z) = |det M_S(Z)| / ||Z||^3", "scale_invariant": bool(q_invariant),
                                      "reading": "the certified radius is read off a projective invariant; the sector |z_g| = max|z_t| is a SELECTOR of a representative, not an open set of the cover"},
           "D_involutions_and_branch": {"diagonal_sign_involutions_preserve_all_F_k": bool(invol_preserve),
                                        "group_order_mod_projective_scalar": 32, "generators": "sigma_1..sigma_5 (sigma_0 = -sigma_1...sigma_5 in P^5)",
                                        "patch": {"S": list(S_d), "gauge": T_d[0], "base": [T_d[1], T_d[2]]},
                                        "radicands": {str(k): str(val) for k, val in radicands.items()},
                                        "branch_locus": "{R_{s1}R_{s2}R_{s3} = 0}", "discriminant_expanded": str(discriminant),
                                        "discriminant_total_degree": disc_degree},
           "E_explicit_transitions": {"elementary_swap": {"from": list(S_a), "to": list(S_b),
                                                          "P_matrix_rows": P_rows,
                                                          "rows_that_change": int(rows_changed),
                                                          "coordinate_form": "z_{s'} = eps_{s'} sqrt((P w_S)_{s'})"},
                                      "gauge_change": "(S,g) -> (S,g') : projective rescaling by z_g/z_{g'}"},
           "F_product_to_factors": {"pivot_threshold_m": str(M_THRESH),
                                    "chain": "|det M_S^alg| > m ⇒ ∏|z_s| > m/(8|V_S|) ; |u|,|v| ≤ 1 ⇒ |z_s| ≤ B_s = sqrt(|a|+|b|+|c|) ⇒ |z_{s_i}| > m/(8|V_S| ∏_{j≠i} B_{s_j})",
                                    "worst_lower_bound_over_60_types": (repr(worst[0]) if worst else None),
                                    "attained_at": ({"S": worst[1], "gauge": worst[2], "coordinate": worst[3]} if worst else None),
                                    "reading": "separation of the sheets is UNIFORM on the certified domains: two ledgers differing in slot s lie at distance at least 2x this bound"},
           "G_grassmannian": {"pivots_are_plucker_coordinates_of_rowspace_dF": bool(plucker_ok),
                              "statement": "dF = 2 V diag(z) is 3x6; its row space is a point of Gr(3,6) and its minor S equals det M_S^alg",
                              "chart_reading": "the 20 triples = standard coordinate charts of Gr(3,6) pulled back; elementary swap = standard change of chart",
                              "nerve_caveat": "J(6,3) is the adjacency graph of the ELEMENTARY changes of chart, not the nerve of the cover (the standard charts of Gr(3,6) all meet)"},
           "does_not_attest": ["the quantitative coverage floor (tau = 0.6, m = 4.8 — upstream certificate read by U1)",
                             "the certified radius (U1)", "the geometry of the bridges", "nothing metric", "nothing about the other fronts"]}
    g = {"Q1_all_pairs_rank2_three_nonzero_lemma": bool(all_pairs_rank2),
         "Q2_pivot_factorisation_exact_all_20_triples": bool(fac_ok),
         "Q3_det_homogeneous_degree_3_and_q_invariant": bool(homog_ok) and bool(q_invariant),
         "Q4_diagonal_involutions_preserve_the_three_quadrics": bool(invol_preserve),
         "Q5_branch_discriminant_degree_6": disc_degree == 6 and len(radicands) == 3,
         "Q6_elementary_swap_changes_exactly_one_row": rows_changed == 1,
         "Q7_upstream_read": bool(upstream_ok) and bool(d3),
         "Q8_product_to_factors_positive_uniform_bound": worst is not None and worst[0] > 0,
         "Q9_pivots_are_plucker_coordinates_of_dF": bool(plucker_ok)}
    out["gates"] = {k: bool(v) for k, v in g.items()}
    out["gates_passed"] = sum(bool(v) for v in g.values()); out["gates_total"] = len(g)
    out["outcome"] = "atlas_paper_three_nonzero_lemma_carries_smoothness_and_coverage_pivot_degree3_invariant_transitions_explicit" if all(g.values()) else "atlas_paper_math_gates_red"
    out["seconds"] = round(time.time() - t0, 1)
    return out


def main():
    out = build()
    st = {"N1_mu_repeat_breaks_three_nonzero_lemma": not build(mutate={"mu_repeat"})["gates"]["Q1_all_pairs_rank2_three_nonzero_lemma"],
          "N2_bad_factor_reddens_pivot_identity": not build(mutate={"bad_factor"})["gates"]["Q2_pivot_factorisation_exact_all_20_triples"],
          "N3_degree_2_reddens_projective_invariance": not build(mutate={"bad_degree"})["gates"]["Q3_det_homogeneous_degree_3_and_q_invariant"]}
    out["self_tests"] = st; out["self_tests_passed"] = sum(bool(v) for v in st.values()); out["self_tests_total"] = len(st)
    try:
        import subprocess
        out["built_from_head"] = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        out["built_from_head"] = None
    out["self_sha256"] = sha(HERE); out["provenance"] = {"sympy_version": sp.__version__}
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"gates {out['gates_passed']}/{out['gates_total']} · self-tests {out['self_tests_passed']}/{out['self_tests_total']} · {out['seconds']} s")
    for k, v in out["gates"].items(): print(f"  {'OK ' if v else 'RED'} {k}")
    for k, v in out["self_tests"].items(): print(f"  {'OK ' if v else 'RED'} {k}")
    d = out["D_involutions_and_branch"]
    print("radicands:", d["radicands"])
    print("product to factors: worst bound %s" % out["F_product_to_factors"]["worst_lower_bound_over_60_types"])
    print("discriminant degree:", d["discriminant_total_degree"], "| swap rows changed:", out["E_explicit_transitions"]["elementary_swap"]["rows_that_change"])
    print("outcome:", out["outcome"])


if __name__ == "__main__":
    main()
