#!/usr/bin/env python3
"""
kahler_metric.py — COHERENT Kaehler primitives (holomorphic convention).

Rebuilt after the Ritz wall, and recorded in the design note of the wall
review: the historical engine `spectral_basis` assembles its blocks with
W-dagger (the polish convention), so G_code is NOT the complex Hessian of
K (E1 median 0.39-0.57), the form is not closed (E2 median 0.187), the
congruence of K is violated by a factor 31.7 (B), and the volume is 71.76
against 4 pi^2 (E3). This module rewrites EVERY primitive in the single
holomorphic chain-rule convention; it does NOT touch the older engine,
which stays the artefact against which the witness v1 was certified.

=== CONVENTION LOCK ============================================================
One block of conventions: the design notes CITE this block instead of
of re-deriving them.

 1. Input: the RAW HOLOMORPHIC section of the radical chart (Z, W),
    gauge Z_g = 1, W[a,alpha] = dZ_a/dw^alpha (w = (u,v)). NEVER the
    sphere frame
    (z, U) from sphere_horizontal_frame: z = Z/sqrt(s) is not holomorphic
    in w, so the plain chain rule does not apply there.
 2. Potential: K(Z) = log(Z-dagger M Z) + sum_e c_e phi_e(Z)/s^{d_e},
    with s = |Z|^2 and M = L L-dagger Hermitian. Projective invariance:
    Z becoming lambda(w).Z (lambda holomorphic without zeros) adds a
    pluriharmonic log|lambda|^2, so G is unchanged (tested: check G5,
    chart transition).
 3. Holomorphic gradient: p_I = W^T . grad z^I, WITHOUT conjugation
    (the former W-dagger convention is the root cause of the Ritz wall).
 4. d_alpha s = zW_alpha = sum_a conj(Z_a).W[a,alpha]  (s is not holomorphic: a real term).
 5. Metric: G_{alpha, beta-bar} = the second derivative of K composed
    with Z, equal to the sum over a, b of
    W[a,alpha].conj(W[b,beta]).(second derivative of K in Z_a, Z_b-bar):
    a pure chain rule, with NO first-order term (Z(w) being holomorphic).
 6. Measure: dV = det G . d^4(Re u, Im u, Re v, Im v).
    the volume integral equals 4 pi^2 for EVERY vector c (the phi/s^d are functions
    global, hence exact; log(rho/s) is globally smooth, so the class is
    that of the Fubini-Study form and is
    rigid). This is a verification check (G7), not an adjustable datum: the free
    scale of the older note was an artefact.
================================================================================

Derivation of the blocks (Wirtinger calculus):

 rho block: r_a = drho/dZ_a = (Z-dagger M)_a;
   H_ρ[a,b] = M_ba/ρ − r_a·conj(r_b)/ρ²
   G_ρ = (WᵀMᵀW̄)/ρ − (Wᵀr)(Wᵀr)†/ρ²

 phi block (Hermitian pair c.z^J.conj(z^L).s^{-d}):
   contracted H = c.[ s^{-d}.pJ (x) conj(pL)
                    − d·s^{-d-1}·( conj(m_L)·pJ ⊗ z̄W + m_J·zW ⊗ p̄L )
                    + m_J·conj(m_L)·( d(d+1)·s^{-d-2}·zW⊗z̄W
                                      − d·s^{-d-1}·WᵀW̄ ) ]
   The real elements (self / real_pair / imag_pair) group their
   pairs; the cross term factors as g1 (x) conj(z)W + zW (x) conj(g1) with
   g1 the holomorphic gradient of the real element.

 first derivatives (weak form): q = phi/s^d, so
   ∂_α q̃ = s^{-d}·( g1_α − d·φ·zW_α/s )

Fubini-Study control (M = I, c = 0): G_FS = (W^T conj(W))/s - zW (x) conj(zW)/s^2.
Note: conj(G_FS) equals U-dagger U of the older engine in the sphere frame, consistent with
the E0 diagnostic check (conj matches to 5.6e-08).
"""
from __future__ import annotations

import io
import sys
from collections import Counter

import numpy as np

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

# shared structures (bases, sampling, values), convention-neutral
from .spectral_basis import (  # noqa: F401  (deliberate re-exports)
    MU, LAMBDA, TRIPLES, V2, VOL_TARGET,
    load_witness, basis_upto, multis_of, basis_values,
    sample_all_charts, sample_chart, minor_inv_times_T_float,
    det2_herm,
)


# ===========================================================================
#  Root primitive: projected holomorphic gradients (W^T, NO conjugation)
# ===========================================================================
def holomorphic_grads(Z, W, multis):
    """m (K,nm) holds the values z^I; p (K,nm,2) the HOLOMORPHIC chart
    gradients:

        p[k, I, α] = Σ_a W[k, a, α] · (∂z^I/∂Z_a)(Z[k])  =  (Wᵀ·∇z^I)[k, α]

    This is the scalar chain rule for a holomorphic chart Z(w):
    d_alpha(z^I of Z) = p[., I, alpha] EXACTLY (not one convention among
    among others). It replaces multi_values_and_projected_grads
    (W-dagger), whose "gradients" were the derivatives of no function at
    all (congruence violated by a factor 31.7)."""
    K = Z.shape[0]
    nm = len(multis)
    m = np.ones((K, nm), dtype=complex)
    p = np.zeros((K, nm, 2), dtype=complex)
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
            p[:, mi, 0] += gv * W[:, a, 0]
            p[:, mi, 1] += gv * W[:, a, 1]
    return m, p


def dlog_s(Z, W):
    """zW_α = ∂_α s = Σ_a Z̄_a W[a,α]  (K,2)."""
    return np.einsum("ka,kaA->kA", Z.conj(), W)


def fs_pullback(Z, W):
    """The complex Hessian of log s in a chart (a control on the machinery)."""
    s = (np.abs(Z) ** 2).sum(axis=1)
    zW = dlog_s(Z, W)
    WtWb = np.einsum("kaA,kaB->kAB", W, W.conj())
    return (WtWb / s[:, None, None]
            - np.einsum("kA,kB->kAB", zW, zW.conj()) / (s ** 2)[:, None, None])


# ===========================================================================
#  Potential (THE definition: the finite differences of the check derive exactly this)
# ===========================================================================
def potential_value(Z, M, coeffs, basis, multis, midx):
    """K(Z) = log(Z-dagger M Z) + sum_e c_e phi_e(Z)/s^{d_e}, real of shape (K,)."""
    rho = np.einsum("ki,ij,kj->k", Z.conj(), M, Z).real
    s = (np.abs(Z) ** 2).sum(axis=1)
    m = np.ones((Z.shape[0], len(multis)), dtype=complex)
    for i, I in enumerate(multis):
        vv = np.ones(Z.shape[0], dtype=complex)
        for o in I:
            vv = vv * Z[:, o]
        m[:, i] = vv
    phi = basis_values(basis, m, s, midx)
    return np.log(rho) + phi @ np.asarray(coeffs, float)


# ===========================================================================
#  Chart metric: G is the pullback of the complex Hessian of K (all blocks)
# ===========================================================================
def chart_metric_kahler(Z, W, M, coeffs, basis, multis, midx,
                        want_element_data=False):
    """g_chart (K,2,2) = W^T . (second derivative of K in Z, conj(Z)) . conj(W) in the (u,v) coordinates.

    INPUT: (Z, W), the RAW holomorphic section of the chart (gauge
    Z_g = 1), as produced by sample_chart or reconstruct. No frame is
    assumed, and no normalisation to the sphere. The formula is GENERAL
    (zW terms included): valid for any holomorphic section, not only a
    horizontal one; gauge covariance is tested by check G5."""
    K = Z.shape[0]
    s = (np.abs(Z) ** 2).sum(axis=1)
    zW = dlog_s(Z, W)
    WtWb = np.einsum("kaA,kaB->kAB", W, W.conj())
    m, p = holomorphic_grads(Z, W, multis)
    # --- rho block ----------------------------------------------------
    rho = np.einsum("ki,ij,kj->k", Z.conj(), M, Z).real
    r = Z.conj() @ M                                   # r_a = (Z†M)_a
    T1 = np.einsum("kaA,ba,kbB->kAB", W, M, W.conj())  # WᵀMᵀW̄
    v = np.einsum("kaA,ka->kA", W, r)                  # Wᵀr
    G = (T1 / rho[:, None, None]
         - np.einsum("kA,kB->kAB", v, v.conj()) / (rho ** 2)[:, None, None])
    # --- blocs phi ------------------------------------------------------
    cross = np.zeros((K, 2), dtype=complex)   # Σ c·d·s^{-d-1}·g1_e
    id1 = np.zeros(K)                          # Σ c·d·φ_e·s^{-d-1}
    id2 = np.zeros(K)                          # Σ c·d(d+1)·φ_e·s^{-d-2}
    coeffs = np.asarray(coeffs, float)
    for e, be in enumerate(basis):
        c = coeffs[e]
        if c == 0.0:
            continue
        I, Kk, typ = be["ij"], be["kl"], be["type"]
        d = float(len(I))
        sd1 = s ** (-d - 1.0)
        mI = m[:, midx[I]]
        mK = m[:, midx[Kk]]
        pI = p[:, midx[I], :]
        pK = p[:, midx[Kk], :]
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
        G += (c * s ** (-d))[:, None, None] * lead
        cross += (c * d * sd1)[:, None] * g1
        id1 += c * d * phi * sd1
        id2 += c * d * (d + 1.0) * phi * s ** (-d - 2.0)
    G -= (np.einsum("kA,kB->kAB", cross, zW.conj())
          + np.einsum("kA,kB->kAB", zW, cross.conj()))
    G += (id2[:, None, None] * np.einsum("kA,kB->kAB", zW, zW.conj())
          - id1[:, None, None] * WtWb)
    G = 0.5 * (G + np.conj(np.transpose(G, (0, 2, 1))))
    if want_element_data:
        return G, (m, p, s, zW)
    return G


# ===========================================================================
#  Fonctions d'analyse : gradients de chart (forme faible)
# ===========================================================================
def basis_chart_grads(basis, m, p, s, midx, zW):
    """∂_α q̃_e (K,nb,2), q̃_e = φ_e/s^d, convention holomorphe :

        ∂_α q̃ = s^{-d}·( ∂_α φ − d·φ·zW_α/s )

    with d_alpha phi assembled from the holomorphic p (same groupings
    self/real/imag que chart_metric_kahler). p DOIT venir de
    holomorphic_grads on the same section (Z, W) as zW."""
    K = m.shape[0]
    nb = len(basis)
    dQ = np.empty((K, nb, 2), dtype=complex)
    for e, be in enumerate(basis):
        I, Kk, typ = be["ij"], be["kl"], be["type"]
        d = float(len(I))
        sd = s ** (-d)
        mI = m[:, midx[I]]
        mK = m[:, midx[Kk]]
        pI = p[:, midx[I], :]
        pK = p[:, midx[Kk], :]
        mv = mI * np.conj(mK)
        if typ == "self":
            g1 = pI * np.conj(mK)[:, None]
            phi = mv.real
        elif typ == "real_pair":
            g1 = pI * np.conj(mK)[:, None] + pK * np.conj(mI)[:, None]
            phi = 2.0 * mv.real
        else:
            g1 = 1j * (pI * np.conj(mK)[:, None] - pK * np.conj(mI)[:, None])
            phi = -2.0 * mv.imag
        dQ[:, e, :] = (g1 - (d * phi / s)[:, None] * zW) * sd[:, None]
    return dQ


def dirichlet_pairing(G, dQa, dQb=None):
    """E(a,b) pointwise = 2 Re[g^{αβ̄} ∂_α a · conj(∂_β b)] ; (K,na,nb).
    Same shape as the older engine, but G and dQ must come
    from the primitives of THIS module (mixing conventions is the bug
    that this module removes)."""
    Ginv = np.linalg.inv(G)
    if dQb is None:
        dQb = dQa
    E = np.einsum("kBA,keA,kfB->kef", Ginv, dQa, dQb.conj())
    return 2.0 * E.real
