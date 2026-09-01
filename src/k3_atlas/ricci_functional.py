#!/usr/bin/env python3
"""
ricci_functional.py — R1 du chantier refit witness v2 : fonctionnelle
Ricci + gradients analytiques + baselines. MODULE importable (R2 = optimiseur).

Cadrage : the fit scoping note.md §4-5 ; paramétrisation R0
(`witness_parameters_C.npz`, 9 + 208 params).

  Fonctionnelle :  F(p9, c) = Var_w[ r ],   r = log det G + 2 log|det M_S|
  Mesure de fit  :  w = det G_FS · 16/n_c  (dV_FS, queue légère — tail study)

Structure exploitée : G(p9, c) = G_ρ(M(p9)) + Σ_e (C·c)_e T_e(x) est AFFINE
in c, so the 218 tensor fields T_e are PRECOMPUTED once (packed
hermitiens 4-réels [g00, g11, Re g01, Im g01]) et

  ∂r/∂c_j = tr(G⁻¹ T̃_j),  T̃ = T·C   (analytique exact, un GEMM)
  ∂F/∂c = 2·⟨ŵ (r − r̄), ∂r/∂c⟩ ;    ∂F/∂p9 : FD central (9 params, ρ-bloc seul)

Frozen STRATIFIED fit sample: a uniform budget per chart, boosted by 4 on the
2 charts porteurs de queue (tail study : S=(0,2,5) g=1 ~40 %, S=(1,3,5) g=2
(about 25 percent), exact because the charts are independent strata (weight 16/n_c
per chart). The min|R| stratum is DEFERRED to a later revision if needed: the
fit measure has a light tail (Hill 7-9), so the band boost is only
required for the diagnostics (the monitors will say so).

Baselines (script __main__) :
  B0  FS pur (p9 = 0, c = 0) ;
  B1  v1 projeté : M(p9_v1) + projection L²(dV_FS, sample R0) du potentiel
      raw 657 of v1 onto the 208 psi (starting point for the continuation) and decomposition
      (part V₁ absorbée par M, part const, résidu).

Sortie : canonical/results/ricci_functional.json
Usage : ricci_functional.py [N_BASE=500] [SEED_FIT=21]
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
    basis_upto, multis_of, basis_values, det2_herm, load_witness,
    sample_chart, detMS_on_block, TRIPLES,
)
from .kahler_metric import (                          # noqa: E402
    chart_metric_kahler, fs_pullback, holomorphic_grads, dlog_s,
)
from .witness_parametrisation import (                            # noqa: E402
    B3, B3_MULTIS, B3_IDX, NB3, m_from_params9, params9_from_witness,
    design_on_sample, load_param_artifact, q3_values,
)

RES = Path(__file__).resolve().parent / "data"
CHUNK = 4000
HEAVY_CHARTS = {((0, 2, 5), 1), ((1, 3, 5), 2)}    # tail study 07-15
BOOST = 4
D3 = 3.0                                            # degré homogène de B3


def pack_herm(G):
    """(K,2,2) hermitien → (K,4) réels [g00, g11, Re g01, Im g01]."""
    return np.stack([G[:, 0, 0].real, G[:, 1, 1].real,
                     G[:, 0, 1].real, G[:, 0, 1].imag], axis=1)


def det_packed(g):
    return g[:, 0] * g[:, 1] - g[:, 2] ** 2 - g[:, 3] ** 2


def tr_inv_packed(g, t, det):
    """tr(G⁻¹T) pour G (K,4), T (K,4) ou (K,E,4) hermitiens packés."""
    sh = (slice(None),) + (None,) * (t.ndim - 2)
    g0, g1, g2, g3 = (g[:, i][sh] for i in range(4))
    return (g1 * t[..., 0] + g0 * t[..., 1]
            - 2.0 * (g2 * t[..., 2] + g3 * t[..., 3])) / det[sh]


# ===========================================================================
#  Sample de fit stratifié (charts = strates indépendantes, poids 16/n_c)
# ===========================================================================
def stratified_fit_sample(seed, n_base, boost=BOOST):
    rng = np.random.default_rng(seed)
    Zs, Ws, wfacs, dMSs, labels = [], [], [], [], []
    for S in TRIPLES:
        T = tuple(j for j in range(6) if j not in S)
        for g_col in T:
            n_c = n_base * (boost if (S, g_col) in HEAVY_CHARTS else 1)
            r = sample_chart(rng, S, g_col, n_c)
            if r is None:
                continue
            Z, W, _ = r
            blk = {"S": S, "g": g_col, "Z": Z}
            Zs.append(Z)
            Ws.append(W)
            dMSs.append(detMS_on_block(blk))
            wfacs.append(np.full(Z.shape[0], 16.0 / n_c))
            labels.append((S, g_col, Z.shape[0]))
    Z = np.concatenate(Zs)
    W = np.concatenate(Ws)
    wfac = np.concatenate(wfacs)
    det_MS = np.concatenate(dMSs)
    w = np.empty(Z.shape[0])
    for i0 in range(0, Z.shape[0], CHUNK):
        sl = slice(i0, min(i0 + CHUNK, Z.shape[0]))
        w[sl] = det2_herm(fs_pullback(Z[sl], W[sl]))
    w *= wfac
    return Z, W, w, det_MS, labels


# ===========================================================================
#  Tenseurs par élément T_e (packés) — miroir du bloc φ du moteur
# ===========================================================================
def element_tensors_packed(Z, W):
    """(K, 218, 4) : T_e = ∂∂̄(q̃_e) pullback, par élément B₃ (d = 3).
    EXACT mirror of the phi block of chart_metric_kahler, verified by the
    selftest T1 (contraction vs moteur à ~1e-14)."""
    K = Z.shape[0]
    s = (np.abs(Z) ** 2).sum(axis=1)
    zW = dlog_s(Z, W)
    WtWb = np.einsum("kaA,kaB->kAB", W, W.conj())
    m, p = holomorphic_grads(Z, W, B3_MULTIS)
    zWc = zW.conj()
    base_id = (12.0 * s ** (-5.0))[:, None, None] \
        * np.einsum("kA,kB->kAB", zW, zWc) \
        - (3.0 * s ** (-4.0))[:, None, None] * WtWb
    sd = s ** (-D3)
    sd1 = s ** (-D3 - 1.0)
    Tp = np.empty((K, NB3, 4))
    for e, be in enumerate(B3):
        I, Kk, typ = be["ij"], be["kl"], be["type"]
        mI = m[:, B3_IDX[I]]
        mK = m[:, B3_IDX[Kk]]
        pI = p[:, B3_IDX[I], :]
        pK = p[:, B3_IDX[Kk], :]
        mv = mI * np.conj(mK)
        Tm = np.einsum("kA,kB->kAB", pI, pK.conj())
        if typ == "self":
            lead = Tm
            g1 = pI * np.conj(mK)[:, None]
            phi = mv.real
        else:
            Ts = np.einsum("kA,kB->kAB", pK, pI.conj())
            if typ == "real_pair":
                lead = Tm + Ts
                g1 = pI * np.conj(mK)[:, None] + pK * np.conj(mI)[:, None]
                phi = 2.0 * mv.real
            else:
                lead = 1j * (Tm - Ts)
                g1 = 1j * (pI * np.conj(mK)[:, None]
                           - pK * np.conj(mI)[:, None])
                phi = -2.0 * mv.imag
        Te = sd[:, None, None] * lead \
            - (D3 * sd1)[:, None, None] * (
                np.einsum("kA,kB->kAB", g1, zWc)
                + np.einsum("kA,kB->kAB", zW, g1.conj())) \
            + phi[:, None, None] * base_id
        Te = 0.5 * (Te + np.conj(np.transpose(Te, (0, 2, 1))))
        Tp[:, e, :] = pack_herm(Te)
    return Tp


# ===========================================================================
#  Problème de fit
# ===========================================================================
class FitProblem:
    """F(p9, c208) = Var_w(r) sur sample stratifié gelé. Gradient c
    analytique exact ; gradient p9 par FD central (ρ-bloc seul)."""

    def __init__(self, seed=21, n_base=500, log=print):
        t0 = time.time()
        art = load_param_artifact()
        self.C = art["C"]
        self.p9_v1 = art["p9_witness_v1"]
        Z, W, w, det_MS, labels = stratified_fit_sample(seed, n_base)
        self.Z, self.W = Z, W
        self.w = w
        self.what = w / w.sum()
        self.log_dMS2 = 2.0 * np.log(np.abs(det_MS))
        self.labels = labels
        self.K = Z.shape[0]
        log(f"  fit sample stratifié : {self.K} pts "
            f"({len(labels)} charts, boost ×{BOOST} sur "
            f"{len(HEAVY_CHARTS)})")
        Tps = []
        for i0 in range(0, self.K, CHUNK):
            sl = slice(i0, min(i0 + CHUNK, self.K))
            Tps.append(element_tensors_packed(Z[sl], W[sl]))
        self.Tp = np.concatenate(Tps)                  # (K, 218, 4)
        log(f"  T_e packés : {self.Tp.shape} "
            f"({self.Tp.nbytes / 1e6:.0f} MB) en {time.time() - t0:.1f}s")

    # --- blocs -----------------------------------------------------------
    def rho_packed(self, p9):
        M = m_from_params9(p9)
        rows = []
        for i0 in range(0, self.K, CHUNK):
            sl = slice(i0, min(i0 + CHUNK, self.K))
            rows.append(pack_herm(chart_metric_kahler(
                self.Z[sl], self.W[sl], M, np.zeros(0), [], [], {})))
        return np.concatenate(rows)

    def G_packed(self, rho_p, c208):
        return rho_p + np.einsum("keq,e->kq", self.Tp, self.C @ c208)

    # --- résidu et fonctionnelle ------------------------------------------
    def r_of(self, gp):
        det = det_packed(gp)
        pd_mask = (det > 0) & (gp[:, 0] > 0)
        r = np.where(pd_mask, np.log(np.abs(det) + 1e-300), 0.0) \
            + self.log_dMS2
        return r, det, pd_mask

    def F_and_grad_c(self, rho_p, c208):
        """F, grad_c (208,), + diagnostics. Var pondérée sur points PD."""
        gp = self.G_packed(rho_p, c208)
        r, det, pd = self.r_of(gp)
        wm = self.what * pd
        wm = wm / wm.sum()
        rbar = wm @ r
        dr = (r - rbar) * pd
        F = float(wm @ dr ** 2)
        # ∂r/∂c_e = tr(G⁻¹ T_e) ; grad = 2 Σ ŵ (r−r̄) ∂r  (points PD)
        Dr218 = tr_inv_packed(gp, self.Tp, np.where(pd, det, 1.0))
        gvec = 2.0 * ((wm * dr) @ Dr218) @ self.C
        return F, gvec, {"pd_frac": float(pd.mean()),
                         "r_mean": float(rbar), "det": det, "r": r,
                         "pd": pd}

    def F_only(self, p9, c208):
        F, _, diag = self.F_and_grad_c(self.rho_packed(p9), c208)
        return F, diag

    def grad_p9_fd(self, p9, c208, h=1e-5):
        g = np.zeros(9)
        for i in range(9):
            pp = p9.copy()
            pp[i] += h
            Fp, _ = self.F_only(pp, c208)
            pp[i] -= 2 * h
            Fm, _ = self.F_only(pp, c208)
            g[i] = (Fp - Fm) / (2 * h)
        return g


# ===========================================================================
#  Projection of the v1 potential onto the 208 psi (start of the continuation)
# ===========================================================================
def project_v1(log=print):
    """Project f_v1 = sum coeffs_raw onto (psi_j) through the design
    (seed 11: the SAME sample as the orthonormalisation of C). Returns
    c208_v1 + décomposition en normes²."""
    art = load_param_artifact()
    C, C1 = art["C"], art["C1"]
    seed, n_draw = int(art["spec"][0]), int(art["spec"][1])
    A, _, w, Zr, _ = design_on_sample(seed, n_draw)
    wit = load_witness()
    RAW = basis_upto(3)[10:]
    raw_multis, raw_idx = multis_of(RAW)
    K = Zr.shape[0]
    f = np.empty(K)
    for i0 in range(0, K, CHUNK):
        sl = slice(i0, min(i0 + CHUNK, K))
        Zs = Zr[sl]
        s_sl = (np.abs(Zs) ** 2).sum(axis=1)
        m = np.ones((Zs.shape[0], len(raw_multis)), dtype=complex)
        for i, I in enumerate(raw_multis):
            vv = np.ones(Zs.shape[0], dtype=complex)
            for o in I:
                vv = vv * Zs[:, o]
            m[:, i] = vv
        f[sl] = basis_values(RAW, m, s_sl, raw_idx) @ wit["coeffs"]
    fbar = (w @ f) / w.sum()
    ftil = np.sqrt(w) * (f - fbar)
    c208 = (A @ C).T @ ftil
    cV1 = (A @ C1).T @ ftil
    n2 = float(ftil @ ftil)
    res2 = n2 - float(c208 @ c208) - float(cV1 @ cV1)
    dec = {"norm2_total": n2, "norm2_c208": float(c208 @ c208),
           "norm2_V1_centered": float(cV1 @ cV1),
           "norm2_residual": res2, "const_part_mean": float(fbar)}
    log(f"  projection v1 : ‖f̃‖² = {n2:.4f} = ψ {dec['norm2_c208']:.4f} "
        f"+ V1 {dec['norm2_V1_centered']:.4f} "
        f"+ résidu {res2:.2e} (attendu ~0 : f ∈ V₃)")
    return c208, dec


# ===========================================================================
#  Script : selftests + baselines
# ===========================================================================
def main():
    T0 = time.time()
    n_base = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    seed_fit = int(sys.argv[2]) if len(sys.argv) > 2 else 21
    results = {"phase": "K3 CAP refit R1 — fonctionnelle + baselines",
               "n_base": n_base, "seed_fit": seed_fit,
               "heavy_charts": [[list(S), g] for S, g in HEAVY_CHARTS],
               "boost": BOOST}
    checks = {}

    def log(msg):
        print(f"[{time.time() - T0:6.1f}s] {msg}")

    def check(name, ok, detail):
        checks[name] = bool(ok)
        log(f"  [{'PASS' if ok else 'FAIL'}] {name} — {detail}")

    log("construction du problème de fit...")
    prob = FitProblem(seed=seed_fit, n_base=n_base, log=log)
    results["K_pts"] = prob.K
    p9_v1 = prob.p9_v1
    rng = np.random.default_rng(123)

    # T1: contraction T.c against the engine (is the mirror exact?)
    c_test = 1e-2 * rng.standard_normal(208)
    rho_p = prob.rho_packed(p9_v1)
    gp = prob.G_packed(rho_p, c_test)
    Zt = prob.Z[:1500]
    Wt = prob.W[:1500]
    G_eng = chart_metric_kahler(Zt, Wt, m_from_params9(p9_v1),
                                prob.C @ c_test, B3, B3_MULTIS, B3_IDX)
    dev = float(np.abs(gp[:1500] - pack_herm(G_eng)).max()
                / np.abs(pack_herm(G_eng)).max())
    check("T1_tensors_vs_engine", dev < 1e-12,
          f"G packé (ρ + T·Cc) vs moteur : rel max = {dev:.2e}")

    # T2 — gradient c analytique vs FD (au point FS + petit c, PD partout)
    c0 = 1e-3 * rng.standard_normal(208)
    rho_fs = prob.rho_packed(np.zeros(9))
    F0, g_ana, diag0 = prob.F_and_grad_c(rho_fs, c0)
    h = 1e-6
    errs = []
    for j in rng.choice(208, 6, replace=False):
        cp = c0.copy()
        cp[j] += h
        Fp, _, _ = prob.F_and_grad_c(rho_fs, cp)
        cp[j] -= 2 * h
        Fm, _, _ = prob.F_and_grad_c(rho_fs, cp)
        gfd = (Fp - Fm) / (2 * h)
        errs.append(abs(gfd - g_ana[j]) / max(abs(gfd), 1e-300))
    check("T2_grad_c_vs_fd", max(errs) < 1e-5,
          f"6 coords : max rel = {max(errs):.2e} (PD = "
          f"{diag0['pd_frac']:.4f})")

    # B0 — baseline FS pur
    F_fs, g_fs, d_fs = prob.F_and_grad_c(rho_fs, np.zeros(208))
    gp9_fs = prob.grad_p9_fd(np.zeros(9), np.zeros(208))
    results["baseline_FS"] = {
        "var_r": F_fs, "pd_frac": d_fs["pd_frac"],
        "grad_c_norm": float(np.linalg.norm(g_fs)),
        "grad_p9_norm": float(np.linalg.norm(gp9_fs))}
    log(f"  B0 FS pur : var(r) = {F_fs:.4f}, PD = {d_fs['pd_frac']:.4f}, "
        f"‖∇c‖ = {np.linalg.norm(g_fs):.3f}, "
        f"‖∇p9‖ = {np.linalg.norm(gp9_fs):.3f}")

    # B1 — baseline v1 projeté
    log("projection v1 (design R0)...")
    c208_v1, dec = project_v1(log=log)
    results["v1_projection"] = dec
    F_v1, g_v1, d_v1 = prob.F_and_grad_c(rho_p, c208_v1)
    gp9_v1 = prob.grad_p9_fd(p9_v1, c208_v1)
    results["baseline_v1proj"] = {
        "var_r": F_v1, "pd_frac": d_v1["pd_frac"],
        "grad_c_norm": float(np.linalg.norm(g_v1)),
        "grad_p9_norm": float(np.linalg.norm(gp9_v1))}
    log(f"  B1 v1 projeté : var(r) = {F_v1:.4f}, PD = "
        f"{d_v1['pd_frac']:.4f}, ‖∇c‖ = {np.linalg.norm(g_v1):.3f}, "
        f"‖∇p9‖ = {np.linalg.norm(gp9_v1):.3f}")

    # contrôle de cohérence : var(r) v1 plein (coeffs raw, moteur) ≈ 0.84
    # mesurée au cadrage sur seed 7 — ici v1 PROJETÉ (V₁ tronquée) diffère
    check("B_baselines_finite",
          np.isfinite(F_fs) and np.isfinite(F_v1),
          f"var_FS = {F_fs:.4f}, var_v1proj = {F_v1:.4f}")

    results["checks"] = checks
    results["all_pass"] = bool(all(checks.values()))
    results["elapsed_seconds"] = time.time() - T0
    (RES / "ricci_functional.json").write_text(
        json.dumps(results, indent=2, default=float))
    print("\n" + "=" * 74)
    print(f"BILAN R1 : {sum(checks.values())}/{len(checks)} PASS — "
          f"baselines : FS var = {F_fs:.4f} (PD {d_fs['pd_frac']:.2%}) ; "
          f"v1proj var = {F_v1:.4f} (PD {d_v1['pd_frac']:.2%})")
    print("=" * 74)
    for kx, v in checks.items():
        print(f"  {'PASS' if v else 'FAIL'}  {kx}")
    print(f"\nJSON → ricci_functional.json   "
          f"({time.time() - T0:.0f}s)")


if __name__ == "__main__":
    main()
