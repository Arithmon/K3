#!/usr/bin/env python3
"""
k3_cap_b1e2iii_owner_tiling.py — B1.e.2.iii-owner : le noyau O0
(classification OWNER d'une BOÎTE par arithmétique d'intervalle) et le
pilote O1 (tuilage par branch-and-bound depuis les 60 couples candidats).
Exécute les gates O0/O1 du contrat GPT
`gpt_b1e3c3_global_owner_pilot_review_2026_07_24.md` §3.1-3.2, en réponse
aux corrections C17 (owner certifié au centre seulement), C18 (disjonction
intra-chart seulement) et C19 (28 charts = échantillonnage, pas preuve).

-------------------------------------------------------------------------
O0 — classification d'une boîte (u₀±h)×(v₀±h) pour un couple (S, g)
-------------------------------------------------------------------------
Critère propriétaire du moteur (miroir EXACT de sample_chart) :
  (1) argmax_t sc_t = S,  sc_t = Π_{i∈t}|Z_i|²·V2[t]  (20 triples) ;
  (2) argmax_{c∈T} |Z_c| = g  (3 jauges du complément) ;
  (3) min_s |R_s| > 1e-12  (radicands, robustesse branche).
Sur la boîte, |Z_g|² = 1, |Z_o1|² = |u|², |Z_o2|² = |v|², |Z_s|² = |R_s|
avec R_s = a_s + b_s·u² + c_s·v² — tout est calculable en intervalle
RÉEL (modules de complexes par re²+im²), SANS évaluer la métrique, et la
classification est INDÉPENDANTE de la feuille ε (les 8 feuilles d'un
(u,v) propriétaire le sont ensemble).

Verdicts (stricts, déterministes) :
  OWNER     sc_S.lo > max_{t≠S} sc_t.hi  ET  |Z_g|².lo > |Z_o|².hi (×2)
            ET  |R_s|².lo > (1e-12)²  ∀s ;
  OUTSIDE   ∃t≠S : sc_t.lo > sc_S.hi  OU  ∃o : |Z_o|².lo > |Z_g|².hi ;
  BRANCH    ∃s : |R_s|² rencontre [0, (1e-12)²]  (et pas OUTSIDE) ;
  AMBIGUOUS sinon — subdivisée par le tuilage, JAMAIS intégrée.

-------------------------------------------------------------------------
O1 — pilote de tuilage (branch-and-bound, budget publié)
-------------------------------------------------------------------------
Pour CHAQUE couple des 60 (S, g) candidats (C19) : grille initiale
N0⁴ sur [-1,1]⁴, subdivision 2⁴ des AMBIGUOUS/BRANCH jusqu'à D_MAX.
Sorties par couple : volumes paramétriques OWNER / OUTSIDE / résiduel
ambigu (Σ = 16), comptes par profondeur. Le volume est PARAMÉTRIQUE
(coordonnées (u,v) du couple) — la masse (∫detg, fermeture 4π², gate O2)
est l'étage suivant, qui ne consommera que des boîtes OWNER certifiées.

Self-test (gates DISCRIMINANTS, dont les tests négatifs obligatoires O0) :
  S1 point possédé (MC) ⟹ boîte dégénérée OWNER ; même (u,v) sous une
     JAUGE fausse ⟹ OUTSIDE ; sous un TRIPLE faux ⟹ OUTSIDE
  S2 boîte volontairement à cheval (domaine entier) ⟹ ni OWNER ni
     OUTSIDE (AMBIGUOUS/BRANCH)
  S3 radicand nul construit (v=0, u² = −a/b) ⟹ BRANCH sur un voisinage
  S4 partition ponctuelle : un point possédé est OWNER pour EXACTEMENT
     un couple sur les 60 (unicité = disjonction, réponse C18)
  S5 validation MC : la fraction possédée mesurée par sample_chart tombe
     dans [V_owner, V_owner + V_ambigu]/16 (couple témoin, profondeur 2)

Sorties : results/k3_cap_b1e2iii_owner_tiling.json
Usage   : k3_cap_b1e2iii_owner_tiling.py [--selftest]
"""
from __future__ import annotations

import io
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

_argv = sys.argv
sys.argv = [sys.argv[0]]
from .spectral_basis import (                               # noqa: E402
    TRIPLES, V2, minor_inv_times_T_float, owner_scores, sample_chart)
from .interval_arithmetic import iv                           # noqa: E402
sys.argv = _argv

RES = Path(os.environ.get(
    "K3_RES_DIR", Path(__file__).resolve().parent / "data"))

N0 = 4                            # grille initiale N0⁴ par couple
D_MAX = 4                         # profondeur de subdivision du pilote
RAD_FLOOR2 = 1e-24                # (1e-12)² — même seuil que le moteur
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:6.1f}s] {msg}", flush=True)


def _iv(lo, hi):
    return iv.mpf([lo, hi])


# ===========================================================================
#  O0 — noyau de classification
# ===========================================================================
def couple_setup(S, g_col):
    T = tuple(j for j in range(6) if j not in S)
    others = [c for c in T if c != g_col]
    o1, o2 = others
    A = minor_inv_times_T_float(S, T)
    perm = [list(T).index(g_col), list(T).index(o1), list(T).index(o2)]
    return T, o1, o2, A[:, perm]


def classify_box(S, g_col, box, setup=None, with_flags=False):
    """box = (ur_lo, ur_hi, ui_lo, ui_hi, vr_lo, vr_hi, vi_lo, vi_hi).
    Retourne 'OWNER' | 'OUTSIDE' | 'BRANCH' | 'AMBIGUOUS'.
    OUTSIDE prime sur BRANCH (sound : une boîte certifiée hors du
    propriétaire n'a pas besoin de traitement de branche) ; with_flags
    expose le drapeau branch interne (gate S3)."""
    T, o1, o2, A = setup if setup is not None else couple_setup(S, g_col)
    ur, ui = _iv(box[0], box[1]), _iv(box[2], box[3])
    vr, vi = _iv(box[4], box[5]), _iv(box[6], box[7])
    ur2, ui2, vr2, vi2 = ur ** 2, ui ** 2, vr ** 2, vi ** 2
    u2r, u2i = ur2 - ui2, 2 * ur * ui                # u² (complexe)
    v2r, v2i = vr2 - vi2, 2 * vr * vi
    m2u = ur2 + ui2                                  # |u|², |v|²
    m2v = vr2 + vi2
    # modules² des 6 coordonnées ; pour s ∈ S : |Z_s|⁴ = |R_s|², on
    # travaille en m2 via sqrt(|R_s|²) — iv.sqrt est enclosure-valide
    m2 = [None] * 6
    m2[g_col] = _iv(1, 1)
    m2[o1], m2[o2] = m2u, m2v
    branch = False
    radicands = []                                   # C101
    for si, s in enumerate(S):
        Rr = A[si, 0] + A[si, 1] * u2r + A[si, 2] * v2r
        Ri = A[si, 1] * u2i + A[si, 2] * v2i
        R2 = Rr ** 2 + Ri ** 2                       # |R_s|² ≥ 0
        if float(R2.a) <= RAD_FLOOR2:
            branch = True
        # C101 : l'enclosure du radicande de la racine de section est
        # calculée ICI, par O0, sur la boîte ENTIÈRE. C'est la seule
        # mesure « range-aware » déjà versionnée de l'arc, et c'est
        # l'objet que la garde √ du TM remplace par un disque
        # isotrope — donc ce qu'un artefact de branche doit publier.
        radicands.append({
            "s_coord": int(s), "row": int(si),
            "Rr": [float(Rr.a), float(Rr.b)],
            "Ri": [float(Ri.a), float(Ri.b)],
            "R2": [float(R2.a), float(R2.b)],
            "Rr_contains_zero": float(Rr.a) <= 0.0 <= float(Rr.b),
            "Ri_contains_zero": float(Ri.a) <= 0.0 <= float(Ri.b),
            "below_rad_floor2": float(R2.a) <= RAD_FLOOR2})
        m2[s] = iv.sqrt(R2)                          # |Z_s|² = |R_s|
    # (1) scores des 20 triples
    sc = [m2[i] * m2[j] * m2[k] * float(V2[t])
          for t, (i, j, k) in enumerate(TRIPLES)]
    t_S = TRIPLES.index(tuple(S)) if isinstance(S, tuple) \
        else TRIPLES.index(tuple(S))
    lo_S, hi_S = float(sc[t_S].a), float(sc[t_S].b)
    others_hi = max(float(sc[t].b) for t in range(len(TRIPLES))
                    if t != t_S)
    pivot_owner = lo_S > others_hi
    pivot_out = any(float(sc[t].a) > hi_S for t in range(len(TRIPLES))
                    if t != t_S)
    # (2) jauge : |Z_g|² strictement max sur T
    g_lo, g_hi = float(m2[g_col].a), float(m2[g_col].b)
    gauge_owner = all(g_lo > float(m2[o].b) for o in (o1, o2))
    gauge_out = any(float(m2[o].a) > g_hi for o in (o1, o2))
    if pivot_out or gauge_out:
        verdict = "OUTSIDE"
    elif branch:
        verdict = "BRANCH"
    elif pivot_owner and gauge_owner:
        verdict = "OWNER"
    else:
        verdict = "AMBIGUOUS"
    if with_flags:
        return verdict, {"branch": bool(branch),
                         "pivot_owner": bool(pivot_owner),
                         "pivot_out": bool(pivot_out),
                         "gauge_owner": bool(gauge_owner),
                         "gauge_out": bool(gauge_out),
                         "radicands": radicands,          # C101
                         "rad_floor2": RAD_FLOOR2}
    return verdict


# ===========================================================================
#  O0-fast — même classification, intervalle numpy VECTORISÉ, arrondi
#  dirigé émulé par nextafter (précédent B1.a ; statut : design-grade,
#  l'oracle mpmath ci-dessus reste la référence — gate S6 croise les deux)
# ===========================================================================
NEG, POS = -np.inf, np.inf


def _dn(x):
    return np.nextafter(x, NEG)


def _up(x):
    return np.nextafter(x, POS)


def i_add(al, ah, bl, bh):
    return _dn(al + bl), _up(ah + bh)


def i_sub(al, ah, bl, bh):
    return _dn(al - bh), _up(ah - bl)


def i_mul(al, ah, bl, bh):
    p = np.stack([al * bl, al * bh, ah * bl, ah * bh])
    return _dn(p.min(axis=0)), _up(p.max(axis=0))


def i_sq(al, ah):
    lo = np.where((al <= 0) & (ah >= 0), 0.0,
                  np.minimum(al * al, ah * ah))
    return _dn(np.where(lo < 0, 0.0, lo)), _up(np.maximum(al * al,
                                                          ah * ah))


def i_scal(c, al, ah):
    if c >= 0:
        return _dn(c * al), _up(c * ah)
    return _dn(c * ah), _up(c * al)


def i_sqrt(al, ah):
    return _dn(np.sqrt(np.maximum(al, 0.0))), _up(np.sqrt(ah))


def classify_boxes_np(S, g_col, boxes, setup=None):
    """boxes (N, 8) → codes (N,) : 0=OWNER 1=OUTSIDE 2=BRANCH 3=AMBIGUOUS.
    Miroir vectorisé de classify_box (mêmes inégalités strictes)."""
    T, o1, o2, A = setup if setup is not None else couple_setup(S, g_col)
    b = np.asarray(boxes, float)
    isq = i_sq
    url, urh, uil, uih = b[:, 0], b[:, 1], b[:, 2], b[:, 3]
    vrl, vrh, vil, vih = b[:, 4], b[:, 5], b[:, 6], b[:, 7]
    ur2 = isq(url, urh)
    ui2 = isq(uil, uih)
    vr2 = isq(vrl, vrh)
    vi2 = isq(vil, vih)
    u2r = i_sub(*ur2, *ui2)
    u2i = i_scal(2.0, *i_mul(url, urh, uil, uih))
    v2r = i_sub(*vr2, *vi2)
    v2i = i_scal(2.0, *i_mul(vrl, vrh, vil, vih))
    m2 = [None] * 6
    ones = np.ones(b.shape[0])
    m2[g_col] = (ones, ones)
    m2[o1] = i_add(*ur2, *ui2)
    m2[o2] = i_add(*vr2, *vi2)
    branch = np.zeros(b.shape[0], bool)
    for si, s in enumerate(S):
        Rr = i_add(np.full_like(ones, _dn(A[si, 0])),
                   np.full_like(ones, _up(A[si, 0])),
                   *i_add(*i_scal(A[si, 1], *u2r),
                          *i_scal(A[si, 2], *v2r)))
        Ri = i_add(*i_scal(A[si, 1], *u2i), *i_scal(A[si, 2], *v2i))
        R2 = i_add(*isq(*Rr), *isq(*Ri))
        branch |= R2[0] <= RAD_FLOOR2
        m2[s] = i_sqrt(*R2)
    sc_lo = np.empty((len(TRIPLES), b.shape[0]))
    sc_hi = np.empty((len(TRIPLES), b.shape[0]))
    for t, (i, j, k) in enumerate(TRIPLES):
        lo, hi = i_mul(*i_mul(*m2[i], *m2[j]), *m2[k])
        sc_lo[t], sc_hi[t] = i_scal(float(V2[t]), lo, hi)
    t_S = TRIPLES.index(tuple(S))
    oth = [t for t in range(len(TRIPLES)) if t != t_S]
    pivot_owner = sc_lo[t_S] > sc_hi[oth].max(axis=0)
    pivot_out = (sc_lo[oth] > sc_hi[t_S]).any(axis=0)
    g_lo, g_hi = m2[g_col]
    gauge_owner = (g_lo > m2[o1][1]) & (g_lo > m2[o2][1])
    gauge_out = (m2[o1][0] > g_hi) | (m2[o2][0] > g_hi)
    code = np.full(b.shape[0], 3, np.int8)
    code[branch] = 2
    code[pivot_owner & gauge_owner & ~branch] = 0
    code[pivot_out | gauge_out] = 1
    return code


# ===========================================================================
#  O1 — tuilage branch-and-bound d'un couple
# ===========================================================================
def split_boxes_np(boxes):
    """(N, 8) → (16N, 8) : subdivision dyadique des 4 dimensions."""
    N = boxes.shape[0]
    mids = (boxes[:, 0::2] + boxes[:, 1::2]) / 2
    out = np.empty((N, 16, 8))
    for m in range(16):
        for i in range(4):
            lo = (m >> i) & 1
            out[:, m, 2 * i] = boxes[:, 2 * i] if lo == 0 else mids[:, i]
            out[:, m, 2 * i + 1] = mids[:, i] if lo == 0 \
                else boxes[:, 2 * i + 1]
    return out.reshape(16 * N, 8)


MAX_FRONTIER = 400_000            # budget par profondeur (PUBLIÉ si atteint)


def tile_couple(S, g_col, n0=N0, d_max=D_MAX):
    """Branch-and-bound vectorisé. Retourne (stats volumes, counts par
    profondeur, capped). Aucune coupe silencieuse : si le budget de
    frontière est atteint, l'excédent part en résiduel et capped=True."""
    setup = couple_setup(S, g_col)
    edges = np.linspace(-1.0, 1.0, n0 + 1)
    g = np.arange(n0)
    a, b, c, d = np.meshgrid(g, g, g, g, indexing="ij")
    idx = np.stack([x.ravel() for x in (a, b, c, d)], axis=1)
    frontier = np.empty((n0 ** 4, 8))
    for i in range(4):
        frontier[:, 2 * i] = edges[idx[:, i]]
        frontier[:, 2 * i + 1] = edges[idx[:, i] + 1]
    stats = {"OWNER": 0.0, "OUTSIDE": 0.0, "residual": 0.0}
    counts, capped = [], False
    for depth in range(d_max + 1):
        vol = ((2.0 / n0) / 2 ** depth) ** 4
        code = classify_boxes_np(S, g_col, frontier, setup)
        cnt = {"OWNER": int((code == 0).sum()),
               "OUTSIDE": int((code == 1).sum()),
               "BRANCH": int((code == 2).sum()),
               "AMBIGUOUS": int((code == 3).sum())}
        counts.append({"depth": depth, **cnt})
        stats["OWNER"] += vol * cnt["OWNER"]
        stats["OUTSIDE"] += vol * cnt["OUTSIDE"]
        todo = frontier[code >= 2]
        if depth == d_max:
            stats["residual"] += vol * todo.shape[0]
            break
        if todo.shape[0] > MAX_FRONTIER // 16:
            keep = MAX_FRONTIER // 16
            stats["residual"] += vol * (todo.shape[0] - keep)
            todo = todo[:keep]
            capped = True
        if todo.shape[0] == 0:
            break
        frontier = split_boxes_np(todo)
    return stats, counts, capped


# ===========================================================================
#  Build
# ===========================================================================
def build():
    print("=" * 78)
    print("B1.e.2.iii-owner — O0 (owner par boîte, intervalle) + O1 "
          "(tuilage 60 couples)")
    print("=" * 78)
    log(f"grille N0={N0}⁴, D_MAX={D_MAX} (branch-and-bound, "
        f"AMBIGUOUS/BRANCH subdivisées, résiduel publié)")
    rows = []
    tot = {"OWNER": 0.0, "OUTSIDE": 0.0, "residual": 0.0}
    n_boxes, n_capped = 0, 0
    for S in TRIPLES:
        T = tuple(j for j in range(6) if j not in S)
        for g_col in T:
            stats, counts, capped = tile_couple(S, g_col)
            nb = sum(sum(c[k] for k in
                         ("OWNER", "OUTSIDE", "AMBIGUOUS", "BRANCH"))
                     for c in counts)
            n_boxes += nb
            n_capped += int(capped)
            for k in tot:
                tot[k] += stats[k]
            rows.append({"S": list(S), "g": g_col,
                         "vol_owner": stats["OWNER"],
                         "vol_outside": stats["OUTSIDE"],
                         "vol_residual": stats["residual"],
                         "owner_fraction": stats["OWNER"] / 16.0,
                         "frontier_capped": bool(capped),
                         "counts": counts})
            log(f"  S={tuple(S)} g={g_col} : owner {stats['OWNER']:8.4f} "
                f"({stats['OWNER'] / 16:6.2%}) · outside "
                f"{stats['OUTSIDE']:8.4f} · résiduel "
                f"{stats['residual']:8.4f} · {nb} boîtes"
                + (" · CAP frontière" if capped else ""))
    n_owner_couples = sum(1 for r in rows if r["vol_owner"] > 0)
    n_empty_cert = sum(1 for r in rows
                       if r["vol_owner"] == 0 and r["vol_residual"] == 0)
    n_empty_amb = sum(1 for r in rows
                      if r["vol_owner"] == 0 and r["vol_residual"] > 0)
    v_tot = 16.0 * len(rows)
    verdict = (
        "B1.e.2.iii-owner O0+O1 LIVRÉS : classification OWNER par BOÎTE "
        "en intervalle (miroir exact du critère du moteur — scores 20 "
        "triples + jauge + radicands ; répond C17), tuilage "
        "branch-and-bound depuis les 60 couples candidats (répond C19). "
        "À D_MAX=%d : %d couples ont un volume OWNER CERTIFIÉ ; les %d "
        "autres se répartissent en %d VACUITÉS CERTIFIÉES (épuisement "
        "OUTSIDE) + %d à résiduel ambigu non tranché à cette profondeur "
        "; volumes paramétriques globaux : owner %.2f%%, "
        "outside %.2f%%, résiduel ambigu %.2f%% (%d boîtes classées). Le "
        "résiduel ambigu est le budget de subdivision d'O2/O3 — aucune "
        "boîte ambiguë n'entre jamais dans un moment. Étage suivant : O2 "
        "fermeture de masse scalaire (Z sur boîtes OWNER + majorant du "
        "résiduel, cible 4π² ∈ [Z], tolérance fixée avant run)." % (
            D_MAX, n_owner_couples,
            60 - n_owner_couples, n_empty_cert, n_empty_amb,
            100 * tot["OWNER"] / v_tot, 100 * tot["OUTSIDE"] / v_tot,
            100 * tot["residual"] / v_tot, n_boxes))

    out = {
        "phase": ("B1.e.2.iii-owner — O0 classification owner par boîte "
                  "(intervalle) + O1 tuilage 60 couples (contrat GPT "
                  "e.3c.4 §3.1-3.2, corrections C17/C18/C19)"),
        "n0": N0, "d_max": D_MAX, "rad_floor2": RAD_FLOOR2,
        "max_frontier": MAX_FRONTIER, "n_couples_capped": n_capped,
        "fast_kernel": ("numpy + arrondi dirigé émulé nextafter "
                        "(design-grade ; oracle mpmath = référence, "
                        "gate S6)"),
        "n_boxes_classified": n_boxes,
        "n_couples_with_certified_owner": n_owner_couples,
        "n_couples_vacuous_certified": n_empty_cert,
        "n_couples_vacuous_ambiguous": n_empty_amb,
        "totals_fraction": {k: v / v_tot for k, v in tot.items()},
        "couples": rows,
        "verdict": verdict}

    print("\nVERDICT :\n" + verdict)
    art = RES / "k3_cap_b1e2iii_owner_tiling.json"
    art.write_text(json.dumps(out, indent=2, ensure_ascii=False,
                              default=float), encoding="utf-8")
    print(f"\n→ {art}")
    return out


# ===========================================================================
#  Self-test
# ===========================================================================
def _selftest():
    fails = []
    S0, g0 = (0, 1, 5), 2
    rng = np.random.default_rng(7)
    res = sample_chart(rng, S0, g0, 2000)
    Z, W, UV = res
    u0, v0 = complex(UV[0, 0]), complex(UV[0, 1])
    log(f"self-test : point possédé témoin (S={S0}, g={g0}), "
        f"u={u0:.3f}, v={v0:.3f}")

    def degen(u, v):
        return (u.real, u.real, u.imag, u.imag,
                v.real, v.real, v.imag, v.imag)

    # --- S1 : dégénéré OWNER / jauge fausse / triple faux ----------------------
    v_ok = classify_box(S0, g0, degen(u0, v0))
    T0_ = tuple(j for j in range(6) if j not in S0)
    g_bad = [c for c in T0_ if c != g0][0]
    v_gbad = classify_box(S0, g_bad, degen(u0, v0))
    S_bad = next(t for t in TRIPLES if t != S0)
    g_for_Sbad = [j for j in range(6) if j not in S_bad][0]
    v_sbad = classify_box(S_bad, g_for_Sbad, degen(u0, v0))
    s1 = v_ok == "OWNER" and v_gbad == "OUTSIDE" and v_sbad == "OUTSIDE"
    fails.append(not s1)
    print(f"[{'PASS' if s1 else 'FAIL'}] S1 négatifs O0 : possédé → "
          f"{v_ok} ; jauge fausse → {v_gbad} ; triple faux → {v_sbad}")

    # --- S2 : boîte à cheval ------------------------------------------------------
    v_big = classify_box(S0, g0, (-1, 1, -1, 1, -1, 1, -1, 1))
    s2 = v_big in ("AMBIGUOUS", "BRANCH")
    fails.append(not s2)
    print(f"[{'PASS' if s2 else 'FAIL'}] S2 boîte domaine entier "
          f"(frontières incluses) → {v_big} (ni OWNER ni OUTSIDE)")

    # --- S3 : radicand nul construit -----------------------------------------------
    s3, msg3 = False, "aucun zéro de radicand trouvé (60 couples, 2 axes)"
    h3 = 1e-6
    for S in TRIPLES:
        for g_col in (j for j in range(6) if j not in S):
            _, _, _, A = couple_setup(S, g_col)
            for si in range(3):
                for col, name in ((1, "u"), (2, "v")):
                    a_, b_ = A[si, 0], A[si, col]
                    if b_ == 0:
                        continue
                    # zéro sur l'axe réel (x² = −a/b) ou imaginaire
                    # pur (x = i·t, x² = −t² ⟹ t² = a/b)
                    for ratio, imag in ((-a_ / b_, False),
                                        (a_ / b_, True)):
                        if ratio <= 0 or ratio > 1:
                            continue
                        r = float(np.sqrt(ratio))
                        re_lo, re_hi = ((-h3, h3) if imag
                                        else (r - h3, r + h3))
                        im_lo, im_hi = ((r - h3, r + h3) if imag
                                        else (-h3, h3))
                        box = ((re_lo, re_hi, im_lo, im_hi,
                                -h3, h3, -h3, h3) if col == 1 else
                               (-h3, h3, -h3, h3,
                                re_lo, re_hi, im_lo, im_hi))
                        v_rad, fl = classify_box(S, g_col, box,
                                                 with_flags=True)
                        # au zéro de R_s : |Z_s| = 0 ⟹ score du triple
                        # s'effondre ⟹ OUTSIDE prime légitimement sur
                        # BRANCH ; le gate exige (a) drapeau branch levé
                        # (b) JAMAIS OWNER
                        if fl["branch"] and v_rad != "OWNER":
                            s3 = True
                            msg3 = (f"S={S} g={g_col} : R_{si} = 0 à "
                                    f"{name} = "
                                    f"{'i·' if imag else ''}{r:.4f} "
                                    f"→ {v_rad}, branch flag levé")
                            break
                    if s3:
                        break
                if s3:
                    break
            if s3:
                break
        if s3:
            break
    fails.append(not s3)
    print(f"[{'PASS' if s3 else 'FAIL'}] S3 radicand ~0 : {msg3}")

    # --- S4 : unicité du propriétaire au POINT PROJECTIF (60 couples) ---------------
    # même point de K3 transporté dans les coordonnées de chaque couple :
    # u' = Z_{o1'}/Z_{g'}, v' = Z_{o2'}/Z_{g'} (ownership projectivement
    # invariant — les scores scalent en |λ|⁶, l'argmax est inchangé)
    Z0 = Z[0]
    owners = []
    for S in TRIPLES:
        T_ = tuple(j for j in range(6) if j not in S)
        for g_col in T_:
            if abs(Z0[g_col]) < 1e-12:
                continue
            o1_, o2_ = [c for c in T_ if c != g_col]
            up, vp = Z0[o1_] / Z0[g_col], Z0[o2_] / Z0[g_col]
            if classify_box(S, g_col, degen(up, vp)) == "OWNER":
                owners.append((tuple(S), g_col))
    s4 = owners == [(S0, g0)]
    fails.append(not s4)
    print(f"[{'PASS' if s4 else 'FAIL'}] S4 unicité projective : "
          f"{len(owners)} couple(s) OWNER = {owners}")

    # --- S5 : validation MC vs volume intervalle -----------------------------------
    stats, counts, _ = tile_couple(S0, g0, n0=N0, d_max=4)
    n_mc = 4000
    res2 = sample_chart(np.random.default_rng(42), S0, g0, n_mc)
    freq = res2[0].shape[0] / 8 / n_mc if res2 is not None else 0.0
    lo, hi = stats["OWNER"] / 16, (stats["OWNER"] + stats["residual"]) / 16
    sig = 3 * np.sqrt(freq * (1 - freq) / n_mc)
    s5 = lo - sig <= freq <= hi + sig and stats["OWNER"] > 0
    fails.append(not s5)
    print(f"[{'PASS' if s5 else 'FAIL'}] S5 MC : freq possédée = "
          f"{freq:.4f} ∈ [{lo:.4f}, {hi:.4f}] ± {sig:.4f} "
          f"(owner/résiduel à D=4 ; volume OWNER certifié > 0 : "
          f"{stats['OWNER'] > 0})")

    # --- S6 : cohérence noyau rapide (numpy/nextafter) vs oracle mpmath -------------
    rng6 = np.random.default_rng(99)
    n6 = 120
    ctr = rng6.uniform(-0.9, 0.9, (n6, 4))
    hw = 10 ** rng6.uniform(-4, -0.5, (n6, 1))
    bx = np.empty((n6, 8))
    bx[:, 0::2] = ctr - hw
    bx[:, 1::2] = ctr + hw
    codes_np = classify_boxes_np(S0, g0, bx)
    names = {0: "OWNER", 1: "OUTSIDE", 2: "BRANCH", 3: "AMBIGUOUS"}
    n_mismatch = 0
    for i in range(n6):
        ref = classify_box(S0, g0, tuple(bx[i]))
        fast = names[int(codes_np[i])]
        # le rapide ne doit JAMAIS certifier (OWNER/OUTSIDE) contre
        # l'oracle ; un déclassement certifié→AMBIGUOUS serait toléré
        if fast != ref and (fast in ("OWNER", "OUTSIDE")
                            or ref in ("OWNER", "OUTSIDE")):
            n_mismatch += 1
    s6 = n_mismatch == 0
    fails.append(not s6)
    print(f"[{'PASS' if s6 else 'FAIL'}] S6 cross-check np vs mpmath : "
          f"{n6} boîtes aléatoires, {n_mismatch} désaccord(s) certifié(s)")

    print("-" * 78)
    print("SELF-TEST:", "FAIL" if any(fails) else "ALL PASS")
    return 1 if any(fails) else 0


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--selftest":
        raise SystemExit(_selftest())
    build()
