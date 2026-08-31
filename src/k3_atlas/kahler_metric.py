#!/usr/bin/env python3
"""
k3_cap_kahler_engine.py — primitives Kähler COHÉRENTES (convention holomorphe).

Chantier refit post-mur-Ritz, acté dans la note
k3_cap_ritz_wall_reviews_confrontation_2026_07_15.md §7 : le moteur
historique `k3_cap_spectral_engine` assemble ses blocs avec W†
(convention refit_polish) ⟹ G_code n'est PAS ∂∂̄K̃ (E1 med 0.39-0.57),
dω̃ ≠ 0 (E2 med 0.187), congruence K violée ×31.7 (B), volume 71.76 vs
4π² (E3). Ce module réécrit TOUTES les primitives dans la convention
chain-rule holomorphe (GPT §Route : « primitive séparée, orientation
unique ») ; il ne modifie PAS l'ancien moteur, qui reste l'artefact
contre lequel le witness v1 a été certifié.

=== CONVENTION LOCK ============================================================
Bloc unique de conventions — les notes de design CITENT ce bloc au lieu
de redériver (reco R2, friction_report_2026-07-15).

 1. Entrée : section HOLOMORPHE brute du chart radical (Z, W), jauge
    Z_g = 1, W[a,α] = ∂Z_a/∂w^α (w = (u,v)). JAMAIS la frame sphère
    (z, U) de sphere_horizontal_frame : z = Z/√s n'est pas holomorphe
    en w, la chain rule pure ne s'y applique pas.
 2. Potentiel : K̃(Z) = log(Z†MZ) + Σ_e c_e φ_e(Z)/s^{d_e},
    s = |Z|², M = LL† hermitienne. Invariance projective : Z → λ(w)·Z
    (λ holomorphe sans zéro) ajoute log|λ|² pluriharmonique ⟹ G inchangé
    (testé : gate G5, transition de charts).
 3. Gradient holomorphe : p_I = Wᵀ·∇z^I — SANS conjugaison
    (l'ex-convention W† est le bug racine du mur Ritz).
 4. ∂_α s = zW_α = Σ_a Z̄_a·W[a,α]   (s n'est pas holomorphe : terme réel).
 5. Métrique : G_αβ̄ = ∂²(K̃∘Z)/∂w^α∂w̄^β
             = Σ_{a,b} W[a,α]·conj(W[b,β])·(∂²K̃/∂Z_a∂Z̄_b)
    — chain rule pure, AUCUN terme du premier ordre (Z(w) holomorphe).
 6. Mesure : dV = det G · d⁴(Re u, Im u, Re v, Im v).
    ∫_{K3} dV = 4π² pour TOUT vecteur c (les φ/s^d sont des fonctions
    globales ⟹ ∂∂̄-exactes ; log(ρ/s) global lisse ⟹ [ω̃] = [ω_FS]
    rigide). C'est un GATE (G7), pas un datum ajustable — la « scale
    libre » de l'ancienne note T2 était un artefact.
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
   paires ; le terme croisé se factorise en g1⊗z̄W + zW⊗conj(g1) avec
   g1 = ∂φ (gradient holomorphe de l'élément réel).

 dérivée premières (forme faible) : q̃ = φ/s^d ⟹
   ∂_α q̃ = s^{-d}·( g1_α − d·φ·zW_α/s )

Contrôle FS (M = I, c = 0) : G_FS = (WᵀW̄)/s − zW⊗z̄W/s² = ∂∂̄ log s.
NB : conj(G_FS) = U†U de l'ancien moteur en frame sphère — cohérent avec
le contrôle E0 du diagnostic (conj matche à 5.6e-08).
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

    C'est la règle de chaîne scalaire pour un chart holomorphe Z(w) —
    ∂_α(z^I∘Z) = p[·, I, α] EXACTEMENT (pas une convention parmi
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
#  Potentiel (LA définition — les FD du gate dérivent exactement ceci)
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
    """g_chart (K,2,2) = Wᵀ·(∂²K̃/∂Z∂Z̄)·W̄ dans les coords (u,v).

    ENTRÉE : (Z, W) section holomorphe BRUTE du chart (jauge Z_g = 1),
    telle que produite par sample_chart / reconstruct. Pas de frame
    sphère. Formule GÉNÉRALE (termes zW inclus) : valide pour toute
    section holomorphe, pas seulement horizontale — la covariance de
    jauge est testée par la gate G5."""
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

    avec ∂_α φ assemblé des p holomorphes (mêmes regroupements
    self/real/imag que chart_metric_kahler). p DOIT venir de
    holomorphic_grads sur la même section (Z, W) que zW."""
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
    Identique en forme à l'ancien moteur — mais G et dQ doivent venir
    des primitives de CE module (le mélange des conventions est le bug
    que ce module supprime)."""
    Ginv = np.linalg.inv(G)
    if dQb is None:
        dQb = dQa
    E = np.einsum("kBA,keA,kfB->kef", Ginv, dQa, dQb.conj())
    return 2.0 * E.real
