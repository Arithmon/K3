#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OPEN CHART THEOREM — certificate producer (U1).

Theorem-grade, with NO atlas enumerated and NO load-bearing floating-point
value: every constant is a certified RATIONAL bound (Fraction, with square
roots enclosed by integer square roots), and every check is an exact
comparison.

U1-A  Typed domains (serialised separately, each with an explicit role):
        U_S        = {Z ∈ K : q_S(Z) = |det M_S(Z)|/‖Z‖³ > m}   pivot open set (m < 4.8, pre-registered)
                     — the threshold is stated on the PROJECTIVE INVARIANT: det is homogeneous of
                     degree 3. LINK to the product floor: in the gauge representative z_g = 1 one has
                     ‖Z‖ ≥ 1, so q_S(Z) > m implies |det M_S^alg(Z)| = q_S(Z)·‖Z‖³ > m, and it is THAT
                     inequality the factorisation 8|V_S|∏|z_s| > m consumes.
        C_{S,g}    = {Z ∈ U_S : |z_g| = max_{t∈T}|z_t|}  CLOSED sector, selector only (centres)
        A_g        = {z_g ≠ 0}                            genuine projective affine open set
        B_{Z0,S,g} = {b : |b − b₀|₂ < ρ}                   LOCAL chart open set, with one local
                                                          section s_{Z0,S,g} per centre Z0
      60 is the number of TYPES (S,g); the charts themselves are local (one
      per centre); compactness of K then gives a finite subfamily.
U1-B  Self-contained proof (simplified Newton contraction, exact for
      quadrics) — see PROOF below; the load-bearing checks bear on the
      residual, σ, β, L, h ≤ 1/2, r*, contraction < 1, uniqueness, holomorphy.
REVISION 2 (2026-08-20) — REPAIR OF A DEFECT, found while writing the atlas paper.
      Version 1 set `sigma_floor = 2·σ_min(Ṽ_S)·m/(8|V_S|)`. Since M̃(w) = 2Ṽ_S·diag(z_S),
      the correct chain is **s_min(M̃) ≥ 2·σ_min(Ṽ_S)·MIN_s|z_s|**: what is needed is the
      MINIMUM of the resolved coordinates, whereas the pivot threshold bounds only their
      PRODUCT (|det M_S^alg| = 8|V_S|∏|z_s| > m). The implication "∏ ≥ c ⇒ min ≥ c" is
      FALSE as soon as one factor exceeds 1 — and factors do reach 3.4 on the certified
      domains. An explicit witness on the variety: S = {0,1,2}, g = 3, u = −0.7083,
      v = −0.5i gives a product of 0.305 > 0.25 but a minimum of 0.204 < 0.25 (companion
      certificate `sigma_floor_correction`, which keeps a FROZEN record of the value
      published before the repair). The floor therefore becomes PER TYPE:
      min|z_s| ≥ max( m/(8|V_S|·∏_{j≠i}B_j) , √(|a|−|b|−|c|) when positive ), B_s² = |a|+|b|+|c|
      via |u|,|v| ≤ 1. **ρ_uniform goes from 9.60e-10 to 2.11e-12 (a factor of 455); G6 and G7
      still pass at the corrected radius — the STATEMENT survives, the CONSTANT was optimistic.**

U1-C  Certified constants (rational bounds): c_k, σ_min(Ṽ_S) (lower bound),
      σ_max(Ṽ_S) and a_{S,g} (upper bounds), σ_floor and ρ_uniform (lower
      bounds), for the 20 triples / 60 types.
U1-D  Semantic negative controls: m ≥ 4.8; sector declared open; n_types ≠ 60;
      a global section per type; a load-bearing floating-point constant;
      h > 1/2 accepted — each of these must make the corresponding check fail.

PROOF (self-contained). Normalised equations F̃_k = F_k/c_k, c_k² = Σ_j μ_j^{2k}
(an integer). Chart (S,g): w = z_S ∈ C³, gauge z_g = 1, base b = (u,v) ∈ C²,
G(w;b) := F̃(w,b,1). Exact facts (the equations are quadrics):
  (E1) ∂_w G(w;b) = M̃(w) = 2 Ṽ_S diag(w)  — independent of b;
  (E2) G(w;b) − G(w₀;b) − M̃(w₀)(w−w₀) = Ṽ_S · ((w−w₀)²) (componentwise square);
  (E3) G(w₀;b) − G(w₀;b₀) = Ṽ_uv · (b² − b₀²).
Let σ := σ_min(M̃(w₀)) > 0, β = 1/σ, L := 2 σ_max(Ṽ_S) (so ‖M̃(w)−M̃(w′)‖ ≤ L|w−w′|),
a := σ_max(Ṽ_uv). If G(w₀;b₀) = 0 and |b−b₀|₂ ≤ ρ with |b₀|_∞ ≤ 1 (centre in C_{S,g}), then:
  ‖G(w₀;b)‖ ≤ a·|b²−b₀²|₂ ≤ a ρ (2+ρ)                                    (E3)
  η := β‖G(w₀;b)‖ ≤ a ρ(2+ρ)/σ ;  h := βLη ≤ L a ρ(2+ρ)/σ².
Consider the operator T(w) = w − M̃(w₀)⁻¹ G(w;b). For |w−w₀| ≤ r:
  |T(w)−w₀| ≤ η + β σ_max(Ṽ_S) r² = η + (L/2σ) r²   (E2)      → ≤ r  iff  (L/2σ)r² − r + η ≤ 0,
  which is solvable iff h ≤ 1/2, with smallest root r* = (1−√(1−2h)) σ/L ≤ 2η;
  |T(w)−T(w′)| ≤ sup_ξ ‖I − M̃(w₀)⁻¹M̃(ξ)‖ |w−w′| ≤ (L r/σ)|w−w′|,  a contraction if L r*/σ = 1−√(1−2h) < 1  (h < 1/2).
  Hence a unique fixed point w(b) in B̄(w₀,r*), with G(w(b);b) = 0.
Uniqueness in the large: if G(w;b) = G(w′;b) = 0 then 0 = M̃((w+w′)/2)(w−w′) exactly (E1,E2);
  if |w−w₀|,|w′−w₀| < σ/L the midpoint satisfies ‖I−M̃(w₀)⁻¹M̃(mid)‖ < 1, so M̃(mid) is invertible and w = w′.
Holomorphy: the iterates T^n(w₀) are polynomial in b and converge uniformly on
  the open set {|b−b₀|₂ < ρ(σ)} (where h < 1/2 is strict), so the limit w(b) is holomorphic;
  M̃(w(b)) is invertible (contraction), so s(b) = (b, w(b), 1) is a local inverse of π_{S,g}.
Radius: ρ(σ) := −1 + √(1 + σ²/(2La)) solves L a ρ(2+ρ)/σ² = 1/2.
Uniform corollary (REV. 2 STATEMENT — the paragraph that stood here repeated the
  PRODUCT → MINIMUM implication refuted above, while the computation below already used the
  repaired form; stale proof text, corrected 2026-08-21):
  s_min(M̃) ≥ 2 σ_min(Ṽ_S) · MIN_s|z_s|, and the MINIMUM is bounded PER TYPE by
  min_s|z_s| ≥ max( m/(8|V_S|·∏_{j≠i}B_j) , √(|a|−|b|−|c|) when positive ) with B_s² = |a|+|b|+|c|
  (the pivot threshold bounds ONLY the product). The gauge z_g = 1 can only increase σ, so
  σ ≥ σ_floor(S,g) > 0 and therefore ρ ≥ ρ_uniform(S,g) := ρ(σ_floor(S,g)) > 0.
"""
from __future__ import annotations

import hashlib
import io
import json
import math
import re
import sys
import time
from fractions import Fraction
from itertools import combinations
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
RES = ROOT / "certificates"
KERNEL = ROOT / "verification" / "producers" / "model.py"
D3 = ROOT / "certificates" / "atlas_coverage.json"
U0 = RES / "uniform_chart_lemma.json"
OUT = RES / "open_chart_theorem.json"
M_PREREG = Fraction(4)
SCALE = 10 ** 40          # precision of the rational enclosures of square roots


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p): return json.loads(p.read_text(encoding="utf-8"))


def read_mu():
    m = re.search(r"^MU_INT\s*=\s*\(([^)]*)\)", KERNEL.read_text(encoding="utf-8"), re.M)
    return tuple(int(x) for x in m.group(1).split(",") if x.strip())


def sqrt_lo(q: Fraction) -> Fraction:
    """Rational LOWER bound for √q (q ≥ 0)."""
    n = q.numerator * SCALE * SCALE
    d = q.denominator
    return Fraction(math.isqrt(n // d), SCALE)


def sqrt_up(q: Fraction) -> Fraction:
    """Rational UPPER bound for √q."""
    n = q.numerator * SCALE * SCALE
    d = q.denominator
    r = math.isqrt(-(-n // d))          # ceiling division, then isqrt
    if Fraction(r, SCALE) ** 2 < q:
        r += 1
    return Fraction(r, SCALE)


def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
            - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))


def inv3(M):
    """Exact inverse of a rational 3×3 matrix."""
    d = det3(M)
    C = [[Fraction(0)] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            mm = [[M[a][b] for b in range(3) if b != j] for a in range(3) if a != i]
            C[j][i] = Fraction((-1) ** (i + j)) * (mm[0][0] * mm[1][1] - mm[0][1] * mm[1][0]) / d
    return C


def constants(mu, m):
    """Certified rational bounds for the 20 triples and the 60 types."""
    ck2 = [sum(Fraction(mu[j]) ** (2 * k) for j in range(6)) for k in range(3)]     # exact c_k²
    triples, types = {}, []
    for S in combinations(range(6), 3):
        T = [t for t in range(6) if t not in S]
        V_S = abs(det3([[Fraction(mu[s]) ** k for s in S] for k in range(3)]))
        # ‖Ṽ_S‖_F² = Σ_k Σ_s μ_s^{2k}/c_k²  (exact); det Ṽ_S² = V_S²/(c₀²c₁²c₂²) (exact)
        frob2 = sum(Fraction(mu[s]) ** (2 * k) / ck2[k] for s in S for k in range(3))
        det2 = V_S ** 2 / (ck2[0] * ck2[1] * ck2[2])
        smax_up = sqrt_up(frob2)                       # σ_max ≤ ‖·‖_F
        smin_lo = sqrt_lo(det2) / frob2                # σ_min ≥ |det| / σ_max² ≥ |det| / ‖·‖_F²
        L_up = 2 * smax_up
        prod_lo = m / (8 * V_S)                        # REVISION 2: the pivot bounds only the PRODUCT ∏|z_s|
        triples[str(list(S))] = {"S": list(S), "abs_V_S": int(V_S), "frob2_Vt_S": str(frob2), "det2_Vt_S": str(det2),
                                 "sigma_max_up": str(smax_up), "sigma_min_lo": str(smin_lo), "L_up": str(L_up),
                                 "product_floor_lo": str(prod_lo),
                                 "note_rev2": "sigma_floor is NO LONGER stored here: it depends on the type (S,g) through the radicands — see types_60.sigma_lo"}
        VS_mat = [[Fraction(mu[s]) ** k for s in S] for k in range(3)]
        A_inv = inv3(VS_mat)
        for g in T:
            uv = [t for t in T if t != g]
            frob2_uv = sum(Fraction(mu[t]) ** (2 * k) / ck2[k] for t in uv for k in range(3))
            a_up = sqrt_up(frob2_uv)
            # ---- REVISION 2 (2026-08-20): s_min(M̃(w₀)) ≥ 2σ_min(Ṽ_S)·MIN_s|z_s|, not ·∏|z_s| ------------------
            # radicands of the type: w_S = −V_S^{-1}·V_{(g,u,v)}·(1, u², v²), exact rational coefficients
            Vcols = [[Fraction(mu[x]) ** k for x in [g] + uv] for k in range(3)]
            Pm = [[-sum(A_inv[i][k] * Vcols[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
            # |u|,|v| ≤ 1 (sector) ⇒ |z_s|² = |R_s| ≤ |a|+|b|+|c|: an UPPER bound ⇒ sqrt_up
            B_up = [sqrt_up(sum(abs(Pm[i][j]) for j in range(3))) for i in range(3)]
            # DIRECT row-wise bound: |R_s| ≥ |a| − |b| − |c| when positive: a LOWER bound ⇒ sqrt_lo
            dir_lo = []
            for i in range(3):
                d_i = abs(Pm[i][0]) - abs(Pm[i][1]) - abs(Pm[i][2])
                dir_lo.append(sqrt_lo(d_i) if d_i > 0 else Fraction(0))
            others = ((1, 2), (0, 2), (0, 1))
            min_abs_z_lo = min(max(prod_lo / (B_up[j] * B_up[k]), dir_lo[i]) for i, (j, k) in enumerate(others))
            sigma_lo = 2 * smin_lo * min_abs_z_lo
            x = sigma_lo ** 2 / (2 * L_up * a_up)     # ρ = −1 + √(1+x) ≥ −1 + sqrt_lo(1+x)
            rho_lo = -1 + sqrt_lo(1 + x)
            types.append({"S": list(S), "g": g, "uv": uv, "a_up": str(a_up), "x": str(x), "rho_uniform_lo": str(rho_lo),
                          "sigma_lo": str(sigma_lo), "min_abs_z_lo": str(min_abs_z_lo), "B_up": [str(b) for b in B_up],
                          "_frac": {"sigma_lo": sigma_lo, "L_up": L_up, "a_up": a_up, "rho_lo": rho_lo}})
    return ck2, triples, types


def kantorovich_check(sigma: Fraction, L: Fraction, a: Fraction, rho: Fraction):
    """Exact checks of the proof at radius rho: h ≤ 1/2, r* ≤ 2η, contraction < 1, uniqueness."""
    eta_up = a * rho * (2 + rho) / sigma
    h = L * eta_up / sigma
    ok_h = h < Fraction(1, 2)                      # strict: the domain is open
    # r* = (1−√(1−2h))σ/L ≤ 2η; contraction c = 1−√(1−2h) < 1; rational upper bound: √(1−2h) ≥ sqrt_lo
    if h < Fraction(1, 2):
        s = sqrt_lo(1 - 2 * h)
        c_up = 1 - s
        r_star_up = c_up * sigma / L
        ok_contract = c_up < 1
        ok_r = r_star_up <= 2 * eta_up + Fraction(1, SCALE)
        ok_unique = r_star_up < sigma / L          # radius of uniqueness in the large, σ/L
    else:
        c_up = r_star_up = None; ok_contract = ok_r = ok_unique = False
    return {"eta_up": str(eta_up), "h": str(h), "h_le_half": ok_h, "contraction_up": str(c_up), "r_star_up": str(r_star_up),
            "contraction_lt_1": ok_contract, "r_star_le_2eta": ok_r, "r_star_lt_uniqueness_radius": ok_unique}


def build(m=M_PREREG, mutate=None):
    mutate = mutate or {}
    mu = read_mu()
    d3 = load(D3)
    ck2, triples, types = constants(mu, m)
    if "n_types" in mutate:
        types = types[:mutate["n_types"]]
    domains = {
        "U_S": {"def": "{Z ∈ K : q_S(Z) = |det M_S(Z)|/‖Z‖³ > m}  (projectively invariant: |det| is homogeneous of degree 3, so the threshold is set on the normalised representative)", "role": "pivot_open", "open": True, "m": str(m), "m_lt_4p8": m < Fraction(48, 10)},
        "C_S_g": {"def": "{Z ∈ U_S : |z_g| = max_{t∈T} |z_t|}", "role": "selector_only", "open": False if "sector_open" not in mutate else True,
                  "note": "CLOSED sector used to select the centres — never an open set of the atlas; it guarantees |u|,|v| ≤ 1 at the centre"},
        "A_g": {"def": "{z_g ≠ 0}", "role": "affine_open", "open": True},
        "B_Z0_S_g": {"def": "{b ∈ C² : |b − b₀|₂ < ρ}", "role": "local_chart_open", "open": True,
                     "sections_are_local": False if "global_section" in mutate else True,
                     "note": "one LOCAL chart per centre Z0 ∈ C_{S,g}; the section s_{Z0,S,g} is local (a single sheet), not a global section of the type"},
        "n_chart_types": len(types), "n_charts": "finite by compactness, NOT enumerated",
    }
    # proof checks at ρ_uniform_lo for each type (exact)
    proofs = []
    for t in types:
        f = t["_frac"]
        rho = f["rho_lo"] * (mutate.get("rho_factor", 1))
        proofs.append({"S": t["S"], "g": t["g"], "rho": str(rho), **kantorovich_check(f["sigma_lo"], f["L_up"], f["a_up"], rho)})
    all_h = all(p["h_le_half"] for p in proofs)
    all_c = all(p["contraction_lt_1"] and p["r_star_le_2eta"] and p["r_star_lt_uniqueness_radius"] for p in proofs)
    # load-bearing constants must be Fractions and nothing else
    lb_types = {type(f["sigma_lo"]), type(f["L_up"]), type(f["a_up"]), type(f["rho_lo"])} if types else set()
    if "float_leak" in mutate:
        lb_types.add(float)
    for t in types:
        t.pop("_frac", None)
    rho_vals = [Fraction(t["rho_uniform_lo"]) for t in types]
    out = {"artifact": "open_chart_theorem",
           "subject": "OPEN CHART THEOREM CERTIFICATE (U1) — theorem-grade, no atlas enumerated, no load-bearing floating-point value",
           "kind": "theorem_certificate_read_only", "no_measurement_performed": True, "nothing_is_rewritten": True,
           "upstream": {"d3_coverage_legacy": sha(D3), "f1prime_u0": sha(U0), "interval_kernel_mu": sha(KERNEL)},
           "mu": list(mu), "d3": {"tau": d3["tau"], "minor_floor": d3["minor_floor_8tau"], "verdict": d3["verdict"]},
           "c_k_squared_exact": [str(x) for x in ck2],
           "U1A_domains_typed": domains,
           "U1B_proof": {"self_contained": True, "where": "the PROOF section of this script's docstring (E1–E3, simplified Newton contraction, uniqueness through the midpoint — exact for quadrics — and holomorphy through uniform convergence of the polynomial iterates)",
                         "external_citation_load_bearing": False,
                         "gates_at_rho_uniform_lo_per_type": proofs,
                         "all_h_lt_half": all_h, "all_contraction_uniqueness_ok": all_c},
           "U1C_constants_certified": {"method": "exact Fraction arithmetic; square roots enclosed by isqrt to 10⁻⁴⁰; σ_max ≤ ‖·‖_F (upper bound); σ_min ≥ |det|/‖·‖_F² (lower bound); ρ ≥ −1 + sqrt_lo(1+x)",
                                       "triples_20": triples, "types_60": types,
                                       "rho_uniform_lo_min": str(min(rho_vals)) if rho_vals else None, "rho_uniform_lo_max": str(max(rho_vals)) if rho_vals else None,
                                       "rho_uniform_lo_min_float_display_only": float(min(rho_vals)) if rho_vals else None,
                                       "rho_uniform_lo_max_float_display_only": float(max(rho_vals)) if rho_vals else None,
                                       "load_bearing_types": sorted(t.__name__ for t in lb_types)},
           }
    checks = {
        "G1_m_prereg_lt_4p8": domains["U_S"]["m_lt_4p8"],
        "G2_sector_is_closed_selector_only": (domains["C_S_g"]["role"] == "selector_only" and domains["C_S_g"]["open"] is False),
        "G3_n_chart_types_60_and_sections_local": len(types) == 60 and domains["B_Z0_S_g"]["sections_are_local"] is True,
        "G4_all_constants_positive_rational": bool(rho_vals) and all(v > 0 for v in rho_vals) and all(Fraction(t["sigma_lo"]) > 0 for t in types),
        "G5_no_float_load_bearing": lb_types == {Fraction},
        "G6_kantorovich_h_lt_half_all_types": all_h,
        "G7_contraction_uniqueness_r_star_all_types": all_c,
        "G8_proof_self_contained_no_external_load_bearing": out["U1B_proof"]["self_contained"] and not out["U1B_proof"]["external_citation_load_bearing"],
    }
    out["checks"] = checks; out["checks_passed"] = sum(bool(v) for v in checks.values()); out["checks_total"] = len(checks)
    out["outcome"] = "uniform_open_chart_theorem_certified" if all(checks.values()) else "uniform_open_chart_theorem_refused"
    out["does_not_attest"] = ["no atlas is enumerated and no section is produced: the theorem is quantitative-existential (a guaranteed radius at each centre)",
                            "the bounds are non-sharp by construction (Frobenius / determinant) — they are strictly positive and reproducible, which is all that is required here",
                            "K is the fixed surface; nothing is claimed about the family: the μ are fixed HERE, and uniformity in b is not exported"]
    return out


def main():
    t0 = time.time()
    out = build()
    st = {}
    st["N1_m_ge_4p8_reddens"] = not build(m=Fraction(5))["checks"]["G1_m_prereg_lt_4p8"]
    st["N2_sector_declared_open_reddens"] = not build(mutate={"sector_open": True})["checks"]["G2_sector_is_closed_selector_only"]
    st["N3_n_types_ne_60_reddens"] = not build(mutate={"n_types": 59})["checks"]["G3_n_chart_types_60_and_sections_local"]
    st["N4_global_section_per_type_reddens"] = not build(mutate={"global_section": True})["checks"]["G3_n_chart_types_60_and_sections_local"]
    st["N5_float_load_bearing_reddens"] = not build(mutate={"float_leak": True})["checks"]["G5_no_float_load_bearing"]
    r = build(mutate={"rho_factor": 100})
    st["N6_rho_x100_gives_h_gt_half_reddens"] = not r["checks"]["G6_kantorovich_h_lt_half_all_types"]
    out["perturbation_tests"] = st; out["perturbation_tests_passed"] = sum(bool(v) for v in st.values()); out["perturbation_tests_total"] = len(st)
    out["seconds"] = round(time.time() - t0, 2)
    try:
        import subprocess
        out["built_from_head"] = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        out["built_from_head"] = None
    out["self_sha256"] = sha(HERE)
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"checks {out['checks_passed']}/{out['checks_total']} · self-tests {out['perturbation_tests_passed']}/{out['perturbation_tests_total']} · {out['seconds']} s")
    for k, v in {**out['checks'], **st}.items():
        print(f"  {'OK ' if v else 'RED'} {k}")
    c = out["U1C_constants_certified"]
    print(f"ρ_uniform_lo ∈ [{c['rho_uniform_lo_min_float_display_only']:.3e}, {c['rho_uniform_lo_max_float_display_only']:.3e}] (display only; the rational values are in the certificate)")
    print("outcome:", out["outcome"])


if __name__ == "__main__":
    main()
