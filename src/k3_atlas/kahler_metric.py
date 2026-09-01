#!/usr/bin/env python3
"""
k3_cap_kahler_engine.py — primitives Kähler COHÉRENTES (convention holomorphe).

Rebuilt after the Ritz wall, and recorded in the design note
k3_cap_ritz_wall_reviews_confrontation_2026_07_15.md §7 : le moteur
historique `k3_cap_spectral_engine` assemble ses blocs avec W†
(polish convention), so G_code is NOT the complex Hessian of K (E1 median 0.39-0.57),
dω̃ ≠ 0 (E2 med 0.187), congruence K violée ×31.7 (B), volume 71.76 vs
4 pi^2 (E3). This module rewrites EVERY primitive in the single
chain-rule holomorphe (GPT §Route : « primitive séparée, orientation
convention; it does NOT touch the older engine, which stays the artefact
contre lequel le witness v1 a été certifié.

=== CONVENTION LOCK ============================================================
One block of conventions: the design notes CITE this block instead of
de redériver (reco R2, friction_report_2026-07-15).

 1. Entrée : section HOLOMORPHE brute du chart radical (Z, W), jauge
    Z_g = 1, W[a,α] = ∂Z_a/∂w^α (w = (u,v)). JAMAIS la frame sphère
    (z, U) from sphere_horizontal_frame: z = Z/sqrt(s) is not holomorphic
    in w, so the plain chain rule does not apply there.
 2. Potentiel : K̃(Z) = log(Z†MZ) + Σ_e c_e φ_e(Z)/s^{d_e},
    s = |Z|², M = LL† hermitienne. Invariance projective : Z → λ(w)·Z
    (λ holomorphe sans zéro) ajoute log|λ|² pluriharmonique ⟹ G inchangé
    (testé : check G5, transition de charts).
 3. Gradient holomorphe : p_I = Wᵀ·∇z^I — SANS conjugaison
    (the former W-dagger convention is the root cause of the Ritz wall).
 4. d_alpha s = zW_alpha = sum_a conj(Z_a).W[a,alpha]  (s is not holomorphic: a real term).
 5. Métrique : G_αβ̄ = ∂²(K̃∘Z)/∂w^α∂w̄^β
             = Σ_{a,b} W[a,α]·conj(W[b,β])·(∂²K̃/∂Z_a∂Z̄_b)
    — chain rule pure, AUCUN terme du premier ordre (Z(w) holomorphe).
 6. Mesure : dV = det G · d⁴(Re u, Im u, Re v, Im v).
    the volume integral equals 4 pi^2 for EVERY vector c (the phi/s^d are functions
    globales ⟹ ∂∂̄-exactes ; log(ρ/s) global lisse ⟹ [ω̃] = [ω_FS]
    rigid). This is a verification check (G7), not an adjustable datum: the free
    scale of the older note was an artefact.
================================================================================

Dérivation des blocs (Wirtinger, ∂̄_b s = Z_b, ∂̄_b Z̄_a = δ_ab) :

 bloc ρ : r_a = ∂ρ/∂Z_a = (Z†M)_a ;
   H_ρ[a,b] = M_ba/ρ − r_a·conj(r_b)/ρ²
   G_ρ = (WᵀMᵀW̄)/ρ − (Wᵀr)(Wᵀr)†/ρ²

 bloc φ (paire hermitienne c·z^J·conj(z^L)·s^{-d}) :
   H contracté = c·[ s^{-d}·pJ⊗p̄L
                    − d·s^{-d-1}·( conj(m_L)·pJ ⊗ z̄W + m_J·zW ⊗ p̄L )
                    + m_J·conj(m_L)·( d(d+1)·s^{-d-2}·zW⊗z̄W
                                      − d·s^{-d-1}·WᵀW̄ ) ]
   Les éléments réels (self / real_pair / imag_pair) regroupent leurs
   pairs; the cross term factors as g1 (x) conj(z)W + zW (x) conj(g1) with
   g1 = ∂φ (gradient holomorphe de l'élément réel).

 dérivée premières (forme faible) : q̃ = φ/s^d ⟹
   ∂_α q̃ = s^{-d}·( g1_α − d·φ·zW_α/s )

Contrôle FS (M = I, c = 0) : G_FS = (WᵀW̄)/s − zW⊗z̄W/s² = ∂∂̄ log s.
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

# structures partagées (bases, sampling, valeurs) — convention-neutres
from .spectral_basis import (  # noqa: F401  (ré-exports volontaires)
    MU, LAMBDA, TRIPLES, V2, VOL_TARGET,
    load_witness, basis_upto, multis_of, basis_values,
    sample_all_charts, sample_chart, minor_inv_times_T_float,
    det2_herm,
)


# ===========================================================================
#  Primitive racine : gradients holomorphes projetés (Wᵀ, SANS conjugaison)
# ===========================================================================
def holomorphic_grads(Z, W, multis):
    """m (K,nm) valeurs z^I ; p (K,nm,2) gradients de chart HOLOMORPHES :

        p[k, I, α] = Σ_a W[k, a, α] · (∂z^I/∂Z_a)(Z[k])  =  (Wᵀ·∇z^I)[k, α]

    This is the scalar chain rule for a holomorphic chart Z(w):
    d_alpha(z^I of Z) = p[., I, alpha] EXACTLY (not one convention among
    d'autres). Remplace multi_values_and_projected_grads (W†) dont les
    « gradients » n'étaient les dérivées d'aucune fonction (congruence K
    violée ×31.7, diagnostic 07-15 §B)."""
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
    """∂∂̄ log s en chart : (WᵀW̄)/s − zW⊗z̄W/s²  (contrôle machinerie)."""
    s = (np.abs(Z) ** 2).sum(axis=1)
    zW = dlog_s(Z, W)
    WtWb = np.einsum("kaA,kaB->kAB", W, W.conj())
    return (WtWb / s[:, None, None]
            - np.einsum("kA,kB->kAB", zW, zW.conj()) / (s ** 2)[:, None, None])


# ===========================================================================
#  Potential (THE definition: the finite differences of the check derive exactly this)
# ===========================================================================
def potential_value(Z, M, coeffs, basis, multis, midx):
    """K̃(Z) = log(Z†MZ) + Σ_e c_e φ_e(Z)/s^{d_e}  (K,) réel."""
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
#  Métrique de chart : G = ∂∂̄K̃ pullback (tous blocs, convention unique)
# ===========================================================================
def chart_metric_kahler(Z, W, M, coeffs, basis, multis, midx,
                        want_element_data=False):
    """g_chart (K,2,2) = W^T . (second derivative of K in Z, conj(Z)) . conj(W) in the (u,v) coordinates.

    ENTRÉE : (Z, W) section holomorphe BRUTE du chart (jauge Z_g = 1),
    as produced by sample_chart or reconstruct. No frame is assumed: any
    sphère. Formule GÉNÉRALE (termes zW inclus) : valide pour toute
    holomorphic section, not only a horizontal one; gauge covariance is
    tested by check G5."""
    K = Z.shape[0]
    s = (np.abs(Z) ** 2).sum(axis=1)
    zW = dlog_s(Z, W)
    WtWb = np.einsum("kaA,kaB->kAB", W, W.conj())
    m, p = holomorphic_grads(Z, W, multis)
    # --- bloc rho -----------------------------------------------------
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
