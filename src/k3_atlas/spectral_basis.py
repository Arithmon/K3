#!/usr/bin/env python3
"""
spectral_basis.py — the float engine of the spectral dry run.

Design: see the spectral dry-run design note.
Convention post-mortem: see the convention fix note.

Radical charts (gauge Z_g = 1, holomorphic coordinates w = (u,v),
holomorphic sheet). Bridge between the chart and the ambient witness:

    z = Z/sqrt(s)                (sphere representative, |z|^2 = 1)
    U = (I - z z-dagger) W/sqrt(s)   (horizontal chart frame,
                                      z_bar^T U = 0)
    U = V.A     (V the orthonormal horizontal SVD frame of the witness;
                 A = V-dagger U invertible)
    g_chart = A-dagger . g_V . A, so positivity is equivalent

RETRACTED 2026-07-13: the metric-dependent path in this module uses the
historical convention below. It is retained for audit only. ``load_witness``
refuses the archived v1 artifact unless the caller explicitly opts in.

HISTORICAL metric convention (inherited from the polish route):
    g_{alpha, beta-bar} = sum_{i,j} conj(W[i,alpha]) . W[j,beta] . gI[i]
                          . conj(gK[j])
                        = A[alpha] . B[beta],  A = W-dagger . gI,
                          B = sum conj(gK) . W

  This is NOT the standard holomorphic chain-rule pullback (which would give
  W^T.gI rather than W-dagger.gI), but the convention against which the witness
  was fitted AND certified. Both forms are valid Hermitian ones;
  only the polish convention is consistent with epsilon = 2.824e-4 and the
  mini-cover 2048/2048.

Projective scalar residual (invariant under the frame U against V):
    r_0 = log det g - log det(U-dagger U) + log jjh,
    with jjh = det(J_Q . J_Q-dagger)

Volume: dV_0 = det g_chart . d^4(u,v) uniform, calibrated so that the
integral is 4 pi^2.

Ownership sampling: (u,v) drawn uniformly on [-1,1]^4 per chart (S,
gauge), 8 reconstructed sheets; a point is kept exactly when THIS (S, g)
owns it canonically (argmax pivot |Z_i Z_j Z_k|^2 . V_S, argmax gauge
|Z_g|).
"""
from __future__ import annotations

import io
import sys
from collections import Counter
from itertools import combinations, combinations_with_replacement
from pathlib import Path

import numpy as np

from .witness_registry import (  # noqa: E402
    WitnessArtifactError,
    load_witness_artifact,
)

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

ROOT = Path(__file__).resolve().parent
NPZ = ROOT / "data" / "retracted" / "k3_closedform_witness_kahler_v1.npz"

MU = np.array([1., 2., 3., 5., 7., 11.])
LAMBDA = np.vstack([np.ones(6), MU, MU ** 2])
COORD_CHARS = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0],
                        [0, 1, 0], [0, 0, 1], [0, 0, 1]], dtype=int)
TRIPLES = list(combinations(range(6), 3))
V2 = np.array([((MU[j] - MU[i]) * (MU[k] - MU[i]) * (MU[k] - MU[j])) ** 2
               for i, j, k in TRIPLES])
VOL_TARGET = 4.0 * np.pi ** 2


# ===========================================================================
#  Witness
# ===========================================================================
def load_witness(*, allow_retracted=False):
    """Load the historical witness only for an explicit audit replay."""
    try:
        d = load_witness_artifact(NPZ, allow_retracted=allow_retracted)
    except WitnessArtifactError as exc:
        raise WitnessArtifactError(
            "the metric-dependent spectral engine is retracted; "
            "no active v2 adapter exists") from exc
    params = d["params_full"]
    rho10 = params[:10]
    coeffs = params[10:]
    u1, ut, uA1, bAr, bAi, uA2, uB1, bBr, bBi, uB2 = rho10
    L = np.zeros((6, 6), dtype=complex)
    L[0, 0] = np.exp(u1)
    L[1, 1] = np.exp(ut)
    L[2, 2] = np.exp(uA1)
    L[3, 2] = bAr + 1j * bAi
    L[3, 3] = np.exp(uA2)
    L[4, 4] = np.exp(uB1)
    L[5, 4] = bBr + 1j * bBi
    L[5, 5] = np.exp(uB2)
    M = L @ L.conj().T
    return {"M": M, "coeffs": np.asarray(coeffs, float), "npz": d}


# ===========================================================================
#  Bases (Z2^3-invariant sector, homogenised real elements phi/s^d)
# ===========================================================================
def char_sum(idx):
    c = np.zeros(3, dtype=int)
    for i in idx:
        c += COORD_CHARS[i]
    return tuple((c % 2).tolist())


def enumerate_sector(degree):
    hol = list(combinations_with_replacement(range(6), degree))
    out = []
    for I in hol:
        for K in hol:
            if char_sum([*I, *K]) != (0, 0, 0):
                continue
            if I < K:
                out.append({"type": "real_pair", "ij": I, "kl": K})
                out.append({"type": "imag_pair", "ij": I, "kl": K})
            elif I == K:
                out.append({"type": "self", "ij": I, "kl": K})
    return out


def basis_upto(deg_max):
    b = []
    for dg in range(1, deg_max + 1):
        b += enumerate_sector(dg)
    return b


def multis_of(basis):
    ms = sorted({b["ij"] for b in basis} | {b["kl"] for b in basis})
    return ms, {m: i for i, m in enumerate(ms)}


# ===========================================================================
#  Ownership sampling by radical charts
# ===========================================================================
def minor_inv_times_T_float(S, T):
    VS = LAMBDA[:, list(S)]
    VT = LAMBDA[:, list(T)]
    return -np.linalg.solve(VS, VT)          # (3,3) real


def owner_scores(Z):
    """(K,20) pivot scores |Z_i Z_j Z_k|^2 V2; the owner is the argmax."""
    m2 = np.abs(Z) ** 2
    sc = np.empty(Z.shape[:-1] + (len(TRIPLES),))
    for t, (i, j, k) in enumerate(TRIPLES):
        sc[..., t] = m2[..., i] * m2[..., j] * m2[..., k] * V2[t]
    return sc


def sample_chart(rng, S, g_col, n_draw, uv_offset=(0., 0., 0., 0.)):
    """Draw n_draw points (u,v) uniformly on [-1,1]^4 in the chart (S, g_col),
    reconstruct the 8 sheets, and keep the points owned by (S, g_col).
    Returns Z (K,6), W (K,6,2), uv (K,2).

    uv_offset is an optional shift of the centre of the (u, v) box. The
    default (0, 0, 0, 0) is backward compatible and gives a box centred on
    [-1,1]^4, so the Z_2 symmetries u -> -u and v -> -v hold. Choosing a
    NON-symmetric offset (non-zero, unequal components) breaks these Z_2
    in the sample, exposing to the analysis basis the whole predicted
    quotient dimension (at zero shift, 103 of 611 dimensions are
    invisible at V_4)."""
    T = tuple(j for j in range(6) if j not in S)
    others = [c for c in T if c != g_col]
    o1, o2 = others
    A = minor_inv_times_T_float(S, T)         # colonnes ordre T
    perm = [list(T).index(g_col), list(T).index(o1), list(T).index(o2)]
    A = A[:, perm]                            # (a, b, c) pour (1, u^2, v^2)
    du_r, du_i, dv_r, dv_i = uv_offset
    u = (du_r + rng.uniform(-1, 1, n_draw)) + 1j * (
        du_i + rng.uniform(-1, 1, n_draw))
    v = (dv_r + rng.uniform(-1, 1, n_draw)) + 1j * (
        dv_i + rng.uniform(-1, 1, n_draw))
    R = A[:, 0][None, :] + np.outer(u * u, A[:, 1]) + np.outer(v * v, A[:, 2])
    w0 = np.sqrt(R + 0j)                      # (n,3) branche principale
    Zs, Ws, UVs = [], [], []
    for eps_id in range(8):
        eps = np.array([1 if (eps_id >> b) & 1 == 0 else -1
                        for b in range(3)], dtype=float)
        Z = np.zeros((n_draw, 6), dtype=complex)
        Z[:, g_col] = 1.0
        Z[:, o1] = u
        Z[:, o2] = v
        Z[:, list(S)] = eps[None, :] * w0
        # ownership: argmax pivot = S AND argmax gauge = g_col
        sc = owner_scores(Z)
        own_S = np.array(TRIPLES)[np.argmax(sc, axis=1)]
        ok_S = (own_S == np.array(S)[None, :]).all(axis=1)
        absT = np.abs(Z[:, list(T)])
        ok_g = np.array(list(T))[np.argmax(absT, axis=1)] == g_col
        keep = ok_S & ok_g & np.isfinite(Z).all(axis=1)
        # exclut les radicands quasi nuls (mesure nulle, robustesse float)
        keep &= np.abs(R).min(axis=1) > 1e-12
        if not keep.any():
            continue
        Zk = Z[keep]
        # tangentes radicales dZ_s/du = b_s u / Z_s, dZ_s/dv = c_s v / Z_s
        W = np.zeros((Zk.shape[0], 6, 2), dtype=complex)
        W[:, o1, 0] = 1.0
        W[:, o2, 1] = 1.0
        Zsk = Zk[:, list(S)]
        W[:, list(S), 0] = (A[:, 1][None, :] * u[keep][:, None]) / Zsk
        W[:, list(S), 1] = (A[:, 2][None, :] * v[keep][:, None]) / Zsk
        Zs.append(Zk)
        Ws.append(W)
        UVs.append(np.stack([u[keep], v[keep]], axis=1))
    if not Zs:
        return None
    return np.concatenate(Zs), np.concatenate(Ws), np.concatenate(UVs)


def sample_all_charts(rng, n_draw_per_chart, uv_offset=(0., 0., 0., 0.)):
    """Boucle (S, g) : 20 triples x 3 jauges. Retourne dict de blocs par
    chart plus concatenated arrays. uv_offset is passed to sample_chart."""
    blocks = []
    for S in TRIPLES:
        T = tuple(j for j in range(6) if j not in S)
        for g_col in T:
            r = sample_chart(rng, S, g_col, n_draw_per_chart, uv_offset)
            if r is None:
                continue
            Z, W, UV = r
            blocks.append({"S": S, "g": g_col, "Z": Z, "W": W, "UV": UV})
    Z = np.concatenate([b["Z"] for b in blocks])
    W = np.concatenate([b["W"] for b in blocks])
    det_MS = np.concatenate([detMS_on_block(b) for b in blocks])
    return blocks, Z, W, det_MS


def detMS_on_block(b):
    """det M_S = 8 Z_i Z_j Z_k V_S on the block (for the density of Omega)."""
    i, j, k = b["S"]
    VS = ((MU[j] - MU[i]) * (MU[k] - MU[i]) * (MU[k] - MU[j]))
    return 8.0 * b["Z"][:, i] * b["Z"][:, j] * b["Z"][:, k] * VS


# ===========================================================================
#  Per-point primitives: monomials, projected gradients, chart metric
# ===========================================================================
def multi_values_and_projected_grads(Z, W, multis):
    """m (K,nm) holds the values z^I; p (K,nm,2) holds the ambient
    gradients contracted with W-dagger (the polish convention):
        p[k, I, α] = Σ_a conj(W[k, a, α]) · (∂z^I/∂Z_a)(Z[k])
                   = (W†·∂z^I)[k, α]
    This convention is NOT the standard pullback (W^T.dz^I through the pure
    holomorphic chain rule), but the one against which the witness was
    fitted (in the polish route: A = einsum('nab,nb->na',
    Vh, gI) with Vh = conj(V.T)), and the one that certifies positivity
    (mini-cover / interval_metric_kahler)."""
    K = Z.shape[0]
    nm = len(multis)
    m = np.ones((K, nm), dtype=complex)
    p = np.zeros((K, nm, 2), dtype=complex)
    Wc = W.conj()
    for mi, I in enumerate(multis):
        cnt = Counter(I)
        val = np.ones(K, dtype=complex)
        for o, mo in cnt.items():
            val = val * Z[:, o] ** mo
        m[:, mi] = val
        for a, ma in cnt.items():
            gv = np.full(K, complex(ma))
            for o, mo in cnt.items():
                e = mo - 1 if o == a else mo
                if e:
                    gv = gv * Z[:, o] ** e
            p[:, mi, 0] += gv * Wc[:, a, 0]
            p[:, mi, 1] += gv * Wc[:, a, 1]
    return m, p


def sphere_horizontal_frame(Z, W):
    """Bridge between the radical chart and the sphere/SVD witness.

    The witness was fitted on the unit sphere (Z_sph with |Z| = 1) with V
    orthonormal and horizontal (SVD of [J_Q ; Z_bar]) and a correction
    -f.I_2. This format is THE convention in which c_a is meaningful.

    In the radical chart one has (Z, W) with |Z|^2 = s and gauge Z_g = 1;
    W = dZ/dw is neither horizontal nor orthonormal. Frame covariance
    gives the equivalent representative directly:

        z = Z/sqrt(s)          (sphere representative of the same
                                projective point)
        P_z = I - z.z-dagger   (horizontal projector at z)
        U = P_z . W / sqrt(s)  (non-orthonormal horizontal frame of (u,v))

    Properties verified:
      z_bar^T U = 0  (horizontal)
      J_Q(z).U = 0   (tangent to the K3, inherited from J_Q(Z).W = 0 and
                      Q(Z) = 0)
        U spans the SAME horizontal plane as V, so there is an invertible A = V-dagger U with
      que U = V·A  ⟹  g_chart = A† · g_V · A (covariance).

    The homogenised correction becomes exactly -f.U-dagger U (instead of
    -f.I_2
    which assumes an orthonormal frame). Returns (z, U)."""
    s = (np.abs(Z) ** 2).sum(axis=1)
    sqrt_s = np.sqrt(s)
    z = Z / sqrt_s[:, None]
    zbW = np.einsum("ka,kaA->kA", z.conj(), W)         # (K, 2) = z̄ᵀW
    U = (W - np.einsum("ka,kA->kaA", z, zbW)) / sqrt_s[:, None, None]
    return z, U


def chart_metric(Z, W, M, coeffs, basis, multis, midx, want_element_data=False):
    """g_chart (K,2,2) is the pullback of the complex Hessian of K in the
    (u,v) coordinates. EXPECTED as input: (Z, W) = (z, U) already passed
    through sphere_horizontal_frame; the caller converts once. Ambient
    formula: the rho block plus sum c.T_raw minus f.U-dagger U.
    Positivity is guaranteed (the witness certifies 4000/4000 on this
    frame by covariance U = V.A). The volume convention
    dV_0 = det(g_chart).d(Re u, Im u, Re v, Im v) is inherited from the
    chart."""
    K = Z.shape[0]
    s = (np.abs(Z) ** 2).sum(axis=1)                  # 1 by construction
    m, p = multi_values_and_projected_grads(Z, W, multis)
    # rho block: W-dagger (M/rho) W - (W-dagger M Z_bar)(...)-dagger/rho^2
    MZc = Z.conj() @ M.T
    rho = np.einsum("ki,ki->k", Z, MZc).real
    WHM = np.einsum("kaA,ab->kAb", W.conj(), M)     # (K,2,6)
    G = np.einsum("kAb,kbB->kAB", WHM, W) / rho[:, None, None]
    wv = np.einsum("kaA,ka->kA", W.conj(), MZc)     # (K,2)
    G = G - np.einsum("kA,kB->kAB", wv, wv.conj()) / rho[:, None, None] ** 2
    # homogenised phi blocks: s^{-d}.H_bare - c.d.phi.s^{-d-1}.W_h-dagger W_h
    WHW = np.einsum("kaA,kaB->kAB", W.conj(), W)    # (K,2,2)
    ident_coef = np.zeros(K)
    for e, be in enumerate(basis):
        c = coeffs[e]
        if c == 0.0:
            continue
        I, Kk, typ = be["ij"], be["kl"], be["type"]
        d = float(len(I))
        sc = c * s ** (-d)
        pI = p[:, midx[I], :]
        pK = p[:, midx[Kk], :]
        mv = m[:, midx[I]] * np.conj(m[:, midx[Kk]])
        Tm = np.einsum("kA,kB->kAB", pI, pK.conj())
        if typ == "self":
            G += sc[:, None, None] * Tm
            phi = mv.real
        else:
            Ts = np.einsum("kA,kB->kAB", pK, pI.conj())
            if typ == "real_pair":
                G += sc[:, None, None] * (Tm + Ts)
                phi = 2.0 * mv.real
            else:
                G += sc[:, None, None] * (1j * (Tm - Ts))
                phi = -2.0 * mv.imag
        ident_coef += c * d * phi * s ** (-d - 1.0)
    G = G - ident_coef[:, None, None] * WHW
    G = 0.5 * (G + np.conj(np.transpose(G, (0, 2, 1))))
    if want_element_data:
        return G, (m, p, s)
    return G


def det2_herm(G):
    return (G[:, 0, 0].real * G[:, 1, 1].real
            - (G[:, 0, 1].real ** 2 + G[:, 0, 1].imag ** 2))


def F0_pointwise(G, jjh, GFS):
    """r_0 = log det G_U − log det(U†U) + log jjh   [scalaire projectif].

    G_U and G_FS_U = U-dagger U pulled back in a frame U = V.A that is not orthonormal
    (chart radical horizontal). Sous changement V → U : det G_U = |det A|²·
    det G_V and det(U-dagger U) = |det A|^2. The factor
    |Omega|^2(U_1,U_2) = |det A|^2 .
    |Omega|^2(V_1,V_2) is absorbed by the volume form dV_0 = det G_U.d^4(u,v).
    so r_0 is frame-invariant: r_U = r_V to machine precision.
    In the polish convention (V orthonormal), r_V = log det G_V + log jjh
    (det(V-dagger V) = 1 disappears). Here -log det(U-dagger U) must be
    kept explicit."""
    return np.log(det2_herm(G)) - np.log(det2_herm(GFS)) + np.log(jjh)


def jjh_on_sphere(z):
    """jjh(z) = det(J_Q(z) · J_Q(z)†), J_Q_k(z) = 2 μ^{k-1} z, k=1..3.
    Requires z on the sphere (|z|^2 = 1) to match the polish convention."""
    Jq = 2.0 * LAMBDA[None, :, :] * z[:, None, :]      # (K, 3, 6)
    G3 = np.einsum("kai,kbi->kab", Jq, Jq.conj())      # (K, 3, 3) Hermitien
    return np.linalg.det(G3).real


# ===========================================================================
#  Fonctions de base : valeurs, grads de chart, Hessiens de chart
# ===========================================================================
def basis_values(basis, m, s, midx):
    """q_e = phi_e/s^d, shape (K, nb), real."""
    K = m.shape[0]
    nb = len(basis)
    Q = np.empty((K, nb))
    for e, be in enumerate(basis):
        I, Kk, typ = be["ij"], be["kl"], be["type"]
        d = float(len(I))
        mv = m[:, midx[I]] * np.conj(m[:, midx[Kk]])
        if typ == "self":
            v = mv.real
        elif typ == "real_pair":
            v = 2.0 * mv.real
        else:
            v = -2.0 * mv.imag
        Q[:, e] = v * s ** (-d)
    return Q


def basis_chart_derivs(basis, m, p, s, midx, Z, W):
    """First derivatives (K,nb,2) and chart Hessian (K,nb,2,2) of the
    basis functions.

    Warning: DIAGNOSTIC ONLY. This routine builds the Laplacian in STRONG
    FORM (the complex Hessian by chain rule plus corrections) and is used
    only for the self-adjointness test. The core spectral pivot now goes
    through the weak form (matrices M, K): the ownership partition creates
    artificial fluxes that the weak form eliminates.

    Convention: p is expected to be W-dagger times the ambient gradient
    (consistent with multi_values_and_projected_grads after the fix). The
    homogenisation corrections s^{-d} times raw minus
    d.phi.s^{-(d+1)}.(W-dagger W) are written for the ambient polish
    convention."""
    K = m.shape[0]
    nb = len(basis)
    dQ = np.empty((K, nb, 2), dtype=complex)
    HQ = np.empty((K, nb, 2, 2), dtype=complex)
    zbW = np.einsum("ka,kaA->kA", Z.conj(), W)      # (K,2) = Z̄ᵀW
    WHW = np.einsum("kaA,kaB->kAB", W.conj(), W)    # (K,2,2)
    for e, be in enumerate(basis):
        I, Kk, typ = be["ij"], be["kl"], be["type"]
        d = float(len(I))
        sd = s ** (-d)
        mI = m[:, midx[I]]
        mK = m[:, midx[Kk]]
        pI = p[:, midx[I], :]
        pK = p[:, midx[Kk], :]
        mv = mI * np.conj(mK)
        Tm = np.einsum("kA,kB->kAB", pI, pK.conj())
        if typ == "self":
            g1 = pI * np.conj(mK)[:, None]
            H = Tm
            phi = mv.real
        else:
            Ts = np.einsum("kA,kB->kAB", pK, pI.conj())
            if typ == "real_pair":
                g1 = pI * np.conj(mK)[:, None] + pK * np.conj(mI)[:, None]
                H = Tm + Ts
                phi = 2.0 * mv.real
            else:
                g1 = 1j * (pI * np.conj(mK)[:, None]
                           - pK * np.conj(mI)[:, None])
                H = 1j * (Tm - Ts)
                phi = -2.0 * mv.imag
        corr1 = (d * phi / s)[:, None] * zbW
        dQ[:, e, :] = (g1 - corr1) * sd[:, None]
        HQ[:, e, :, :] = (H * sd[:, None, None]
                          - (d * phi * s ** (-d - 1.0))[:, None, None] * WHW)
    return dQ, HQ


def laplacian_of_basis(G, HQ):
    """The Laplacian of q_e: 2 g^{ab} Hess[a,b], shape (K, nb), real."""
    Ginv = np.linalg.inv(G)                       # (K,2,2)
    lap = 2.0 * np.einsum("kBA,keAB->ke", Ginv, HQ)
    return lap.real


def dirichlet_pairing(G, dQa, dQb=None):
    """E(a,b) pointwise = 2 Re[g^{αβ̄} ∂_αa conj(∂_βb)] ; (K, na, nb)."""
    Ginv = np.linalg.inv(G)
    if dQb is None:
        dQb = dQa
    E = np.einsum("kBA,keA,kfB->kef", Ginv, dQa, dQb.conj())
    return 2.0 * E.real
