#!/usr/bin/env python3
"""
k3_cap_refit_param.py — R0 du chantier refit witness v2 : paramétrisation
identifiable gelée. MODULE importable (R1/R2 consomment design_on_sample,
m_from_params9, params9_from_witness) + script de génération des artefacts.

Cadrage : k3_cap_refit_cadrage_2026_07_15.md §3 (+ §6b résultats).
L'ancien fit travaillait dans la base raw `basis_upto(3)[10:]` = 657 params
pour un espace de fonctions identifiable de 208 (440 directions plates
mesurées). Ce module construit et GÈLE la paramétrisation v2 :

  potentiel v2 :  K̃(Z) = log(Z†MZ) + Σ_j c_j · ψ_j(Z)

  - M = LL†, L Cholesky structurée par caractères Z₂³, JAUGE det M = 1
    (Σ u_i = 0) ⟹ 9 paramètres réels libres (pack/unpack ici) ;
  - ψ_j = Σ_e C[e,j]·q̃_e, q̃_e ∈ B₃ (base quotient degré 3, 218 éléments,
    d ≡ 3 homogène — consommable telle quelle par chart_metric_kahler) ;
  - C (218 × 208) : base ORTHONORMÉE L²(dV_FS, sample gelé seed=11) d'un
    complément de span(const, V₁) dans V₃, HIÉRARCHIQUE par blocs :
    48 colonnes = V₂⊖V₁ puis 160 colonnes = V₃⊖V₂.
    (V₁ est porté par M ; la constante est tuée par le centrage — NB la
    constante vit DANS V₁ : Σ_i |z_i|²/s = 1, d'où rang(V₁ centré) = 9.)

Construction (route équilibrée — jamais le Gram) : matrice de design
A = √w·(Q₃ − ⟨Q₃⟩_w), w = det G_FS·(16/N) (mesure de fit dV_FS, queue
légère — tail study), CGS2 + passe de raffinement, représentations en
coefficients exactes conservées.

Artefacts (GELÉS — données du witness v2, pas des recalculables) :
  canonical/results/k3_cap_refit_param_C.npz   (C, C1, E13, mean, spec)
  canonical/results/k3_cap_refit_param.json    (selftests + diagnostics)

Selftests bloquants (résultats du run gelé : note cadrage §6b, 8/8) :
  S1 rangs hiérarchiques exacts (9 / 48 / 160) ;
  S2 orthonormalité sur le sample gelé ;
  S3 complément : (A·C) ⊥ (A·C1), fuite constante < 1e-8 (∂∂̄const = 0) ;
  S4 G est AFFINE en c ;
  S5 invariance de rejauge du résidu Ricci r = log det G + 2 log|det M_S| ;
  S6 généralisation : cond(Gram) sur un sample FRAIS (mesuré : 5.9) ;
  S7 contraste : design raw 657 cond 1.88e19, rang(1e-8) = 217 = 208 + 9.

Usage : k3_cap_refit_param.py [N_DRAW=1000] [SEED_FROZEN=11]
"""
from __future__ import annotations

import io
import json
import sys
import time
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

from .spectral_basis import (                        # noqa: E402
    basis_upto, multis_of, basis_values, det2_herm,
    sample_all_charts, load_witness, MU,
)
from .kahler_metric import (                          # noqa: E402
    chart_metric_kahler, fs_pullback,
)
from .invariant_quotient_ring import basis_at_deg_quotient    # noqa: E402
from .nested_inclusion import build_E                 # noqa: E402

RES = Path(__file__).resolve().parent / "data"
CHUNK = 4000
RANK_TOL = 1e-8
DEFAULT_SEED_FROZEN = 11
DEFAULT_N_DRAW = 1000

# ===========================================================================
#  Base quotient V3 (module-level : consommée par R1/R2)
# ===========================================================================
B1 = basis_at_deg_quotient(1)
B2 = basis_at_deg_quotient(2)
B3 = basis_at_deg_quotient(3)
NB1, NB2, NB3 = len(B1), len(B2), len(B3)
B3_MULTIS, B3_IDX = multis_of(B3)


# ===========================================================================
#  M-bloc : L structurée, jauge det M = 1 (9 params libres)
# ===========================================================================
# Structure héritée du witness v1 (load_witness) : diag exp(u_0..u_5) +
# 2 sous-diagonales complexes (blocs de caractères A et B). Jauge :
# Σ u_i = 0 ⟺ det L réel = 1 ⟺ det M = 1 ; on paramètre les 6 u par
# 5 différences + les 4 réels off-diag ⟹ 9.
def m_from_params9(p9):
    """p9 = (du_1..du_5, bAr, bAi, bBr, bBi) → M = LL†, det M = 1.
    u_0 = −(du_1+..+du_5), u_i = du_i pour i ≥ 1 ⟹ Σ u = 0."""
    du = np.asarray(p9[:5], float)
    u = np.concatenate([[-du.sum()], du])
    bAr, bAi, bBr, bBi = p9[5:]
    L = np.zeros((6, 6), dtype=complex)
    for i in range(6):
        L[i, i] = np.exp(u[i])
    L[3, 2] = bAr + 1j * bAi
    L[5, 4] = bBr + 1j * bBi
    return L @ L.conj().T


def params9_from_witness():
    """Projette le M du witness v1 sur la jauge det M = 1 (M → M/det^{1/6}
    = rescale global de L par exp(−ū) ; le potentiel ne change que d'une
    constante, invisible pour ∂∂̄ et pour var(r))."""
    d = load_witness()["npz"]
    rho10 = d["params_full"][:10]
    u1, ut, uA1, bAr, bAi, uA2, uB1, bBr, bBi, uB2 = rho10
    u = np.array([u1, ut, uA1, uA2, uB1, uB2])
    ub = u.mean()
    du = (u - ub)[1:]
    s = np.exp(-ub)
    return np.concatenate([du, [bAr * s, bAi * s, bBr * s, bBi * s]])


# ===========================================================================
#  Design A = √w (Q3 − mean) sous dV_FS (réutilisé par la projection v1)
# ===========================================================================
def q3_values(Zs):
    """Valeurs des 218 éléments B₃ (K, 218)."""
    s_sl = (np.abs(Zs) ** 2).sum(axis=1)
    m = np.ones((Zs.shape[0], len(B3_MULTIS)), dtype=complex)
    for i, I in enumerate(B3_MULTIS):
        vv = np.ones(Zs.shape[0], dtype=complex)
        for o in I:
            vv = vv * Zs[:, o]
        m[:, i] = vv
    return basis_values(B3, m, s_sl, B3_IDX)


def design_on_sample(seed, n_draw):
    """Sample propriétaire + design pondéré-centré. Retourne
    (A, mean, w, Zr, Wr)."""
    rng = np.random.default_rng(seed)
    _, Zr, Wr, _ = sample_all_charts(rng, n_draw)
    K = Zr.shape[0]
    w = np.empty(K)
    Q3 = np.empty((K, NB3))
    for i0 in range(0, K, CHUNK):
        sl = slice(i0, min(i0 + CHUNK, K))
        w[sl] = det2_herm(fs_pullback(Zr[sl], Wr[sl]))
        Q3[sl] = q3_values(Zr[sl])
    w *= 16.0 / n_draw
    mean = (w @ Q3) / w.sum()
    A = np.sqrt(w)[:, None] * (Q3 - mean[None, :])
    return A, mean, w, Zr, Wr


def orth_block(A, cols, prev_C_list, tol=RANK_TOL):
    """Orthonormalise A·cols contre les blocs précédents (A·C_prev déjà
    orthonormés) puis en interne (SVD). CGS2 + raffinement. Toutes les
    opérations restent des combinaisons de colonnes ⟹ représentation
    coefficient EXACTE conservée. Retourne (C_blk, sv)."""
    cols = cols.copy()
    B = A @ cols
    for _ in range(2):                        # CGS2 : « twice is enough »
        for Cp in prev_C_list:
            Qp = A @ Cp
            coef = Qp.T @ B
            cols = cols - Cp @ coef
            B = B - Qp @ coef
    U, S, Vt = np.linalg.svd(B, full_matrices=False)
    k = int((S > tol * S[0]).sum())
    C_blk = cols @ (Vt[:k].T / S[:k][None, :])
    Bk = A @ C_blk                            # passe de raffinement
    for Cp in prev_C_list:
        Qp = A @ Cp
        coef = Qp.T @ Bk
        C_blk = C_blk - Cp @ coef
        Bk = Bk - Qp @ coef
    U2, S2, Vt2 = np.linalg.svd(Bk, full_matrices=False)
    C_blk = C_blk @ (Vt2.T / S2[None, :])
    return C_blk, S


def load_param_artifact():
    """Charge l'artefact gelé (C, C1, mean, p9_v1, spec). Pour R1/R2."""
    d = np.load(RES / "k3_cap_refit_param_C.npz")
    return {k: d[k] for k in d.files}


# ===========================================================================
#  Script de génération + selftests
# ===========================================================================
def main():
    T0 = time.time()
    n_draw = int(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_N_DRAW
    seed_frozen = int(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_SEED_FROZEN
    seed_fresh = seed_frozen + 1
    results = {"phase": "K3 CAP refit R0 — parametrisation identifiable",
               "N_DRAW": n_draw, "seed_frozen": seed_frozen,
               "seed_fresh": seed_fresh}
    checks = {}

    def log(msg):
        print(f"[{time.time() - T0:6.1f}s] {msg}")

    def check(name, ok, detail):
        checks[name] = bool(ok)
        log(f"  [{'PASS' if ok else 'FAIL'}] {name} — {detail}")

    p9w = params9_from_witness()
    M9 = m_from_params9(p9w)
    check("gauge_detM", abs(np.linalg.det(M9).real - 1.0) < 1e-12,
          f"det M(p9_witness) = {np.linalg.det(M9).real:.15f}")

    log(f"bases quotient : V1 = {NB1}, V2 = {NB2}, V3 = {NB3}")
    E12, _, _ = build_E(1)
    E23, _, _ = build_E(2)
    E13 = E23 @ E12
    V1cols = E13[1:, 1:]
    V2cols = E23[1:, 1:]
    results["dims"] = {"V1": NB1, "V2": NB2, "V3": NB3,
                       "rank_E13": int(np.linalg.matrix_rank(E13))}

    log("sample gelé + design...")
    A, mean_frozen, w_frozen, Zr, Wr = design_on_sample(seed_frozen, n_draw)
    K_PTS = A.shape[0]
    log(f"design gelé : {K_PTS} pts × {NB3}")

    log("orthonormalisation hiérarchique...")
    C1, S1v = orth_block(A, V1cols, [])
    C2, S2v = orth_block(A, V2cols, [C1])
    C3, S3v = orth_block(A, np.eye(NB3), [C1, C2])
    k1, k2, k3 = C1.shape[1], C2.shape[1], C3.shape[1]
    C = np.concatenate([C2, C3], axis=1)
    results["ranks"] = {"V1_centered": k1, "V2_minus_V1": k2,
                        "V3_minus_V2": k3, "C_total": int(C.shape[1])}
    results["block_sv_gaps"] = {
        "V1": [float(S1v[k1 - 1] / S1v[0]),
               float(S1v[k1] / S1v[0]) if k1 < len(S1v) else 0.0],
        "V2": [float(S2v[k2 - 1] / S2v[0]),
               float(S2v[k2] / S2v[0]) if k2 < len(S2v) else 0.0],
        "V3": [float(S3v[k3 - 1] / S3v[0]),
               float(S3v[k3] / S3v[0]) if k3 < len(S3v) else 0.0]}
    check("S1_ranks", (k1, k2, k3) == (9, 48, 160),
          f"rangs (V1c, V2⊖V1, V3⊖V2) = ({k1}, {k2}, {k3}) attendu "
          f"(9, 48, 160) ; gaps sv = {results['block_sv_gaps']}")

    AC = A @ C
    gram = AC.T @ AC
    dev = float(np.linalg.norm(gram - np.eye(C.shape[1]))
                / np.sqrt(C.shape[1]))
    check("S2_orthonormal_frozen", dev < 1e-10, f"‖Gram − I‖/√n = {dev:.2e}")
    AC1 = A @ C1
    cross = float(np.abs(AC1.T @ AC).max())
    const_leak = float(np.abs((np.sqrt(w_frozen) @ AC)).max()
                       / np.sqrt(w_frozen.sum()))
    # seuil fuite constante : 1e-8 — bruit float de centrage amplifié par
    # ‖C‖ (petits sv ~2e-4 du bloc V3), et direction NULLE de la métrique
    # (∂∂̄ const = 0 : une composante constante de c ne change pas G)
    check("S3_complement", cross < 1e-10 and const_leak < 1e-8,
          f"max|⟨V1, ψ⟩| = {cross:.2e}, fuite constante = {const_leak:.2e} "
          f"(inoffensive : ∂∂̄const = 0)")

    # S4 — G affine en c
    rng_t = np.random.default_rng(99)
    ca = 1e-2 * rng_t.standard_normal(C.shape[1])
    cb = 1e-2 * rng_t.standard_normal(C.shape[1])
    Zt, Wt = Zr[:200], Wr[:200]

    def G_of(c208):
        return chart_metric_kahler(Zt, Wt, M9, C @ c208, B3,
                                   B3_MULTIS, B3_IDX)

    G0 = G_of(np.zeros(C.shape[1]))
    Ga = G_of(ca)
    Gb = G_of(cb)
    Gab = G_of(ca + cb)
    aff = float(np.abs(Gab - Ga - Gb + G0).max() / np.abs(Gab).max())
    check("S4_affine_in_c", aff < 1e-12,
          f"‖G(a+b) − G(a) − G(b) + G(0)‖/‖G‖ = {aff:.2e}")

    # S5 — invariance de rejauge du résidu Ricci (c aléatoire)
    def residual_r(Z, W, M, coeffs, det_MS):
        G = chart_metric_kahler(Z, W, M, coeffs, B3, B3_MULTIS, B3_IDX)
        return (np.log(np.abs(det2_herm(G)))
                + 2.0 * np.log(np.abs(det_MS)))

    rngr = np.random.default_rng(5)
    blocks_r, _, _, _ = sample_all_charts(rngr, 200)
    b = blocks_r[0]
    T = tuple(j for j in range(6) if j not in b["S"])
    o1 = [c for c in T if c != b["g"]][0]
    Zb, Wb = b["Z"][:80], b["W"][:80]
    Zg = Zb[:, o1]
    keepm = np.abs(Zg) > 0.3
    Zb, Wb, Zg = Zb[keepm], Wb[keepm], Zg[keepm]
    a1 = b["g"]
    a2 = [c for c in T if c not in (b["g"], o1)][0]
    D = (Wb * Zg[:, None, None]
         - Zb[:, :, None] * Wb[:, o1, None, :]) / (Zg ** 2)[:, None, None]
    J = np.stack([D[:, a1, :], D[:, a2, :]], axis=1)
    Zp = Zb / Zg[:, None]
    Wp = np.einsum("kaA,kAB->kaB", D, np.linalg.inv(J))
    i_, j_, k_ = b["S"]
    VS = (MU[j_] - MU[i_]) * (MU[k_] - MU[i_]) * (MU[k_] - MU[j_])
    dMS = 8.0 * Zb[:, i_] * Zb[:, j_] * Zb[:, k_] * VS
    dMSp = 8.0 * Zp[:, i_] * Zp[:, j_] * Zp[:, k_] * VS
    cr = C @ (1e-2 * rng_t.standard_normal(C.shape[1]))
    rA = residual_r(Zb, Wb, M9, cr, dMS)
    rB = residual_r(Zp, Wp, M9, cr, dMSp)
    rej = float(np.abs(rA - rB).max())
    check("S5_r_gauge_invariance", rej < 1e-10,
          f"max|Δr| rejauge = {rej:.2e} (n = {len(Zb)}, c aléatoire)")

    # S6 — généralisation sample frais
    log("sample frais (généralisation)...")
    Af, _, _, _, _ = design_on_sample(seed_fresh, n_draw)
    ACf = Af @ C
    gram_f = ACf.T @ ACf
    ev = np.linalg.eigvalsh(0.5 * (gram_f + gram_f.T))
    cond_fresh = float(ev[-1] / ev[0])
    results["fresh_gram"] = {"cond": cond_fresh,
                             "ev_min": float(ev[0]), "ev_max": float(ev[-1])}
    check("S6_fresh_conditioning", cond_fresh < 1e3,
          f"cond(Gram frais) = {cond_fresh:.1f} "
          f"(λ ∈ [{ev[0]:.3f}, {ev[-1]:.3f}])")

    # S7 — contraste raw 657
    log("contraste base raw 657...")
    RAW = basis_upto(3)[10:]
    raw_multis, raw_idx = multis_of(RAW)
    Qraw = np.empty((K_PTS, len(RAW)))
    for i0 in range(0, K_PTS, CHUNK):
        sl = slice(i0, min(i0 + CHUNK, K_PTS))
        Zs = Zr[sl]
        s_sl = (np.abs(Zs) ** 2).sum(axis=1)
        m = np.ones((Zs.shape[0], len(raw_multis)), dtype=complex)
        for i, I in enumerate(raw_multis):
            vv = np.ones(Zs.shape[0], dtype=complex)
            for o in I:
                vv = vv * Zs[:, o]
            m[:, i] = vv
        Qraw[sl] = basis_values(RAW, m, s_sl, raw_idx)
    mean_raw = (w_frozen @ Qraw) / w_frozen.sum()
    Araw = np.sqrt(w_frozen)[:, None] * (Qraw - mean_raw[None, :])
    nrm = np.linalg.norm(Araw, axis=0)
    Sraw = np.linalg.svd(Araw / np.where(nrm == 0, 1.0, nrm)[None, :],
                         compute_uv=False)
    rank_raw = int((Sraw > RANK_TOL * Sraw[0]).sum())
    results["raw_contrast"] = {
        "n_params_raw": len(RAW),
        "cond_raw_design": float(Sraw[0] / Sraw[-1]),
        "rank_raw_at_1e-8": rank_raw,
        "n_flat_raw": len(RAW) - rank_raw,
        "cond_C_design": 1.0}
    check("S7_raw_contrast", rank_raw <= C.shape[1] + k1,
          f"raw 657 : cond = {Sraw[0] / Sraw[-1]:.2e}, rang(1e-8) = "
          f"{rank_raw} ⟹ {len(RAW) - rank_raw} directions plates ; "
          f"design C : cond = 1 (orthonormée par construction)")

    all_pass = all(checks.values())
    np.savez_compressed(
        RES / "k3_cap_refit_param_C.npz",
        C=C, C1=C1, E13=E13, V1cols=V1cols, V2cols=V2cols,
        mean_frozen=mean_frozen,
        p9_witness_v1=p9w,
        spec=np.array([seed_frozen, n_draw, K_PTS]),
        block_dims=np.array([k1, k2, k3]))
    results["checks"] = checks
    results["all_pass"] = bool(all_pass)
    results["artifact"] = "k3_cap_refit_param_C.npz"
    results["elapsed_seconds"] = time.time() - T0
    (RES / "k3_cap_refit_param.json").write_text(
        json.dumps(results, indent=2, default=float))
    print("\n" + "=" * 74)
    print(f"BILAN R0 : {sum(checks.values())}/{len(checks)} PASS — "
          f"paramétrisation v2 = 9 (M) + {C.shape[1]} (c) = "
          f"{9 + C.shape[1]}")
    print("=" * 74)
    for kx, v in checks.items():
        print(f"  {'PASS' if v else 'FAIL'}  {kx}")
    print(f"\nnpz + JSON → results/   ({time.time() - T0:.0f}s)")


if __name__ == "__main__":
    main()
