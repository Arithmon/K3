#!/usr/bin/env python3
"""
witness_parametrisation.py — stage R0 of the witness v2 refit: the frozen
identifiable parametrisation. An importable MODULE (R1 and R2 consume
design_on_sample, m_from_params9, params9_from_witness) and a script that
generates the artefacts.

Framing: the fit scoping note. The earlier fit worked in the raw basis
`basis_upto(3)[10:]` = 657 parameters for an identifiable function space
of dimension 208 (440 flat directions measured). This module builds and
FREEZES the v2 parametrisation:

  potential v2:  K(Z) = log(Z-dagger M Z) + sum_j c_j . psi_j(Z)

  - M = L L-dagger, with L a Cholesky factor structured by the Z_2^3
    characters, under the GAUGE det M = 1 (sum u_i = 0), giving 9 free
    real parameters (packed and unpacked here);
  - psi_j = sum_e C[e,j].q_e, with q_e in the degree-3 quotient basis (218
    elements,
    homogeneous of degree 3, consumable as such by
    chart_metric_kahler);
  - C (218 x 208): an L^2-ORTHONORMAL basis (dV_FS, frozen sample with
    seed 11) of a complement of span(const, V_1) inside V_3, HIERARCHICAL
    by blocks: 48 columns for V_2 minus V_1, then 160 columns for V_3
    minus V_2. (V_1 is carried by M; the constant is killed by centring.
    Note that the constant lives IN V_1, since the sum of |z_i|^2/s is 1,
    so the centred V_1 has rank 9.)

Construction (a balanced route that never forms the Gram matrix): the
design matrix A = sqrt(w).(Q_3 - <Q_3>_w) with w = det G_FS.(16/N) (the
fit measure dV_FS, a light tail; see the tail study), CGS2 plus a
refinement pass, with exact coefficient representations preserved.

Artefacts (FROZEN: witness data, not recomputables):
  canonical/results/witness_parameters_C.npz   (C, C1, E13, mean, spec)
  canonical/results/witness_parametrisation.json    (selftests + diagnostics)

Blocking self-tests (results of the frozen run, 8/8):
  S1 exact hierarchical ranks (9 / 48 / 160);
  S2 orthonormality on the frozen sample;
  S3 complement: (A.C) orthogonal to (A.C1), constant leak below 1e-8
     (the complex Hessian of a constant vanishes);
  S4 G is AFFINE in c;
  S5 regauging invariance of the Ricci residual
     r = log det G + 2 log|det M_S|;
  S6 generalisation: cond(Gram) on a FRESH sample (measured: 5.9);
  S7 contrast: the raw 657 design has cond 1.88e19 and rank(1e-8) = 217 =
     208 + 9.

Usage : witness_parametrisation.py [N_DRAW=1000] [SEED_FROZEN=11]
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
#  Degree-3 quotient basis (module level: consumed by stages R1 and R2)
# ===========================================================================
B1 = basis_at_deg_quotient(1)
B2 = basis_at_deg_quotient(2)
B3 = basis_at_deg_quotient(3)
NB1, NB2, NB3 = len(B1), len(B2), len(B3)
B3_MULTIS, B3_IDX = multis_of(B3)


# ===========================================================================
#  M block: structured L, gauge det M = 1 (9 free parameters)
# ===========================================================================
# Structure inherited from the witness v1 (load_witness): diagonal
# exp(u_0..u_5) plus 2 complex subdiagonals (character blocks A and B).
# Gauge: sum u_i = 0, equivalently det L real equal to 1, equivalently
# det M = 1; the 6 values u are parametrised by 5 differences plus the 4
# real off-diagonal entries, giving 9.
def m_from_params9(p9):
    """p9 = (du_1..du_5, bAr, bAi, bBr, bBi) gives M = L L-dagger with
    det M = 1, where u_0 = -(du_1 + ... + du_5) and u_i = du_i for
    i >= 1, so that sum u = 0."""
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
    """Project the M of the v1 witness onto the gauge det M = 1 (M -> M/det^{1/6}),
    a global rescaling of L by exp(-u); the potential changes only by an
    additive constant, invisible to the complex Hessian and to var(r))."""
    d = load_witness()["npz"]
    rho10 = d["params_full"][:10]
    u1, ut, uA1, bAr, bAi, uA2, uB1, bBr, bBi, uB2 = rho10
    u = np.array([u1, ut, uA1, uA2, uB1, uB2])
    ub = u.mean()
    du = (u - ub)[1:]
    s = np.exp(-ub)
    return np.concatenate([du, [bAr * s, bAi * s, bBr * s, bBi * s]])


# ===========================================================================
#  Design A = sqrt(w).(Q3 - mean) under dV_FS (reused by the v1 projection)
# ===========================================================================
def q3_values(Zs):
    """Values of the 218 elements of the degree-3 basis (K, 218)."""
    s_sl = (np.abs(Zs) ** 2).sum(axis=1)
    m = np.ones((Zs.shape[0], len(B3_MULTIS)), dtype=complex)
    for i, I in enumerate(B3_MULTIS):
        vv = np.ones(Zs.shape[0], dtype=complex)
        for o in I:
            vv = vv * Zs[:, o]
        m[:, i] = vv
    return basis_values(B3, m, s_sl, B3_IDX)


def design_on_sample(seed, n_draw):
    """Ownership sample plus weighted, centred design. Returns
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
    """Orthonormalise A.cols against the previous blocks (A.C_prev
    already orthonormal), then internally (SVD). CGS2 plus refinement.
    Every operation stays a combination of columns, so the EXACT
    coefficient representation is preserved. Returns (C_blk, sv)."""
    cols = cols.copy()
    B = A @ cols
    for _ in range(2):                        # CGS2: "twice is enough"
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
    """Load the frozen artefact (C, C1, mean, p9_v1, spec)."""
    d = np.load(RES / "witness_parameters_C.npz")
    return {k: d[k] for k in d.files}


# ===========================================================================
#  Generation script plus self-tests
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

    log("frozen sample plus design...")
    A, mean_frozen, w_frozen, Zr, Wr = design_on_sample(seed_frozen, n_draw)
    K_PTS = A.shape[0]
    log(f"frozen design: {K_PTS} pts by {NB3}")

    log("hierarchical orthonormalisation...")
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
    # constant-leak threshold 1e-8: the float noise of centring amplified by
        # norm of C (small singular values around 2e-4 of the V3 block), and the NULL
        # direction of the metric (the Hessian of a constant vanishes: a constant
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
          f"||G(a+b) - G(a) - G(b) + G(0)||/||G|| = {aff:.2e}")

    # S5: regauging invariance of the Ricci residual (random c)
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
          f"max regauging deviation of r = {rej:.2e} (n = {len(Zb)}, random c)")

    # S6: generalisation on a fresh sample
    log("fresh sample (generalisation)...")
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
          f"raw 657: cond = {Sraw[0] / Sraw[-1]:.2e}, rank(1e-8) = "
          f"{rank_raw}, so {len(RAW) - rank_raw} flat directions; "
          f"design C: cond = 1 (orthonormal by construction)")

    all_pass = all(checks.values())
    np.savez_compressed(
        RES / "witness_parameters_C.npz",
        C=C, C1=C1, E13=E13, V1cols=V1cols, V2cols=V2cols,
        mean_frozen=mean_frozen,
        p9_witness_v1=p9w,
        spec=np.array([seed_frozen, n_draw, K_PTS]),
        block_dims=np.array([k1, k2, k3]))
    results["checks"] = checks
    results["all_pass"] = bool(all_pass)
    results["artifact"] = "witness_parameters_C.npz"
    results["elapsed_seconds"] = time.time() - T0
    (RES / "witness_parametrisation.json").write_text(
        json.dumps(results, indent=2, default=float))
    print("\n" + "=" * 74)
    print(f"BILAN R0 : {sum(checks.values())}/{len(checks)} PASS — "
          f"v2 parametrisation = 9 (M) + {C.shape[1]} (c) = "
          f"{9 + C.shape[1]}")
    print("=" * 74)
    for kx, v in checks.items():
        print(f"  {'PASS' if v else 'FAIL'}  {kx}")
    print(f"\nnpz + JSON → results/   ({time.time() - T0:.0f}s)")


if __name__ == "__main__":
    main()
