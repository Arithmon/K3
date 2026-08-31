#!/usr/bin/env python3
"""
k3_cap_b1e2iii_rface_p0_b1_leaf.py — RFace-P0 / B1 : LA CONTINUATION
TRAVERSE LA FACE Re, ET LA FEUILLE ATTEINTE EST DÉRIVÉE, PAS SUPPOSÉE.

CE QUE CE SCRIPT PAIE — la phase B1 du contrat RFace-P0 (13ᵉ revue §7.3-B,
amendée par la 14ᵉ revue §1/§4/§5) sur le témoin déterministe
`w_one_neighbor` (tuile 295, flush face Re u HAUTE, où le scout v2 a
adressé une cellule voisine canonique portant DEUX classes candidates).

L'ARCHITECTURE EST CELLE DU §5 DE LA 14ᵉ REVUE, DANS L'ORDRE :
  1. le PONT-Re est construit sur la géométrie 2H — la face Re u est
     traversée avec demi-largeur `w₀ = 2H` de chaque côté, les autres
     coordonnées héritées du pont F0 (déjà bilatéral en Im) ;
  2. le chart n'est PAS hérité : la jauge du chart cible est re-minorée
     sur TOUTE la boîte élargie, les radicandes sur toute la boîte ;
  3. la section du pont-Re est construite DIRECTEMENT (constructeur
     bilatéral : régime par signe CERTIFIÉ de `Re R`, refus si
     straddle), avec le ledger F0 CERTIFIÉ `(1,−1,−1)` du témoin ;
  4. sa RESTRICTION au pont F0 est comparée à la section F0 certifiée —
     recentrage anisotrope exact, soustraction coefficient par
     coefficient, `θ = +1` exigé sur les SIX coordonnées ;
  5. sur l'overlap voisin `W_cell = pont-Re ∩ cellule voisine` (qui est
     un CUBE de demi-largeur H — l'anisotropie disparaît, mesuré et
     gaté), le ledger côté voisin est DÉRIVÉ ligne par ligne contre la
     section continuée, puis RECONSTRUIT et re-vérifié (le motif
     F2d-bis) ;
  6. les DEUX candidates enregistrées (classe 0, ε=(1,1,1) ; classe 1,
     ε=(1,1,−1)) sont construites sur `W_cell` et comparées CHACUNE à la
     continuation, sous le gate à TROIS BRANCHES PRÉ-ENREGISTRÉES de la
     14ᵉ revue §1-B1 — ET LE RÉSULTAT MESURÉ EST LA TROISIÈME BRANCHE :
     **AUCUNE étiquette enregistrée ne ferme naïvement** (`none_closes`,
     sérialisé tel quel, PAS corrigé en douce) ;
  7. ce résultat est ÉLUCIDÉ PAR UN CONTRÔLE, pas par un argument : la
     MÊME expérience sur NOTRE PROPRE CÔTÉ (l'étiquette registre de
     NOTRE cellule, (1,1,1), construite sur le cube symétrique
     `W_ours = miroir-Re(W_cell)`) ne ferme PAS non plus, avec le MÊME
     motif de θ. La fermeture naïve d'étiquette n'identifie donc AUCUNE
     paire — pas même la nôtre : l'étiquette de registre vit dans la
     convention du registre, le ledger du pont dans la convention des
     régimes signe-certifiés (rotated→canonique absorbe σ), et la
     conversion `c` par ligne est MESURÉE par le contrôle ;
  8. l'IDENTIFICATION est donc RELATIVE AU CONTRÔLE (gate B1f) : la
     feuille atteinte est la classe dont le motif de θ côté voisin est
     IDENTIQUE au motif du contrôle côté source — même étiquette, même
     comportement à travers la face. La perdante doit différer du
     contrôle EXACTEMENT sur la coordonnée 5 prédite, marge O(1).

PRÉDICTIONS PRÉ-ENREGISTRÉES (falsifiables) : (i) les deux candidates
ne diffèrent que par le signe de la ligne 2 (`s = 5`) ⟹ si les régimes
sur `W_cell` restent `(principal, canonical, canonical)`, leurs motifs
de θ ne peuvent différer QUE sur la coordonnée 5 ; (ii) le motif du
contrôle doit relier étiquette et ledger dérivé par
`motif[s] = ε_dérivé[s]·étiquette[s]` ligne à ligne.

CONVENTIONS, DÉCLARÉES ET PAS ÉCRASÉES : le gate naïf et son
`none_closes` restent sérialisés (`naive_closure_outcome`) À CÔTÉ de
l'identification relative et du ledger dérivé. Si l'identification
relative et la dérivation désignaient des étiquettes incompatibles, ce
serait une DISCORDANCE PUBLIÉE (13ᵉ revue §B : « toute discordance
devient un résultat, pas une correction automatique »).

CE QUE CE SCRIPT NE PAIE PAS : le raccord du témoin 315 (B3 : deux
ordres, commutation — la diagonale n'est adressée qu'À L'ADRESSE) ; les
transitions d'atlas et le nerf étendu (C) ; la métrique (E) ; l'atlas de
la cellule voisine (le pilote certifie un raccord LOCAL à une cellule
adressée, pas son atlas) ; les faces basses (frontière d'énumération) ;
les 895 autres paires.

GATES
  B1a  amont vert et modes full là où le champ existe (scout v2 8/8,
       F0 scout 7/7, F2/F3 v2 17/17 full, atlas C127-D 14/14 full) —
       et témoin/face/boîtes IMPORTÉS des artefacts, jamais recodés ;
  B1b  géométrie : `w₀ = 2H` exact (la face F0 est affleurante), le
       pont-Re contient le pont F0 ET son miroir-Re, bornes ATTEINTES ;
       tous les coins dyadiques (float(x) exact, gaté) ;
       NÉGATIF : la boîte de demi-largeur `w₀/2` échoue DES DEUX CÔTÉS ;
  B1c  section sur la boîte ÉLARGIE : 3 régimes assignés, IDENTIQUES aux
       régimes F0 sérialisés, radicandes minorés > 0, jauge du chart
       cible minorée > 0 sur TOUTE la boîte (chart re-validé) ;
  B1d  restriction : `θ = +1` sur 6/6 coordonnées contre la section F0
       (elle-même re-construite et confrontée aux régimes sérialisés),
       marges de séparation publiées ;
  B1e  ledger voisin DÉRIVÉ puis RECONSTRUIT : `θ = +1` sur 6/6 après
       reconstruction — et W_cell est le CUBE attendu (4 largeurs
       égales, ancre centrale STRICTEMENT intérieure, Re u > face) ;
  B1f′ identification PAR CONVERSION EXPLICITE (v2, 15ᵉ revue §3) :
       B1f1 κ = ε_canonique_source ⊙ ε_registre_source, dérivé et
            confirmé par le motif MESURÉ du contrôle (6 coordonnées) ;
       B1f2 κ concorde avec les lois de détermination : conversion
            CONFINÉE aux lignes rotated (core C127-E) → canonique
            (pont), κ = +1 sur la ligne principale ;
       B1f3 les DEUX étiquettes converties `κ ⊙ ε_reg` sérialisées ;
       B1f4 les deux sections converties RECONSTRUITES : EXACTEMENT
            UNE ferme (et elle == le ledger dérivé de B1e) ;
       B1f5 la perdante convertie casse UNIQUEMENT sur la coordonnée 5
            prédite, marge O(1) ;
       B1f6 le gate naïf `none_closes` reste publié en diagnostic ;
  B1g  NÉGATIF : muter UNE ligne du ledger dérivé ⟹ le raccord tombe
       sur cette ligne (θ = −1), marge O(1) ;
  B1h  aucun filtrage silencieux : 6 coordonnées partout, tout refus
       publié.

Sortie : results/k3_cap_b1e2iii_rface_p0_b1_leaf.json
Usage  : k3_cap_b1e2iii_rface_p0_b1_leaf.py [--selftest]
"""
from __future__ import annotations

import hashlib
import io
import json
import os
import platform
import subprocess
import sys
import time
from fractions import Fraction
from pathlib import Path

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

_argv = sys.argv
sys.argv = [sys.argv[0]]
os.environ.setdefault("K3_TM_ORDER", "4")
os.environ.setdefault("K3_TM_SERIES", "4")
from mpmath import mp                                              # noqa: E402
from .taylor_models import CIV, civ_absmin                       # noqa: E402
from .atlas_assembly import _f_down, _f_up              # noqa: E402
from .bridge_continuation import (               # noqa: E402
    build_section_bilateral, center_hw, contains, inter,
    theta_lines)
sys.argv = _argv

RES = Path(os.environ.get(
    "K3_RES_DIR", Path(__file__).resolve().parent / "data"))
COVER_JSON = RES / "k3_cap_b1e2iii_d5_6_dyadic_cover.json"
ATLAS_JSON = RES / "k3_cap_b1e2iii_c127d_atlas.json"
C127E_JSON = RES / "k3_cap_b1e2iii_c127e_residual.json"
SCOUT0_JSON = RES / "k3_cap_b1e2iii_c129f_bridge_scout.json"
F2F3_JSON = RES / "k3_cap_b1e2iii_c129f_f2f3_bridge_atlas.json"
RSCOUT_JSON = RES / "k3_cap_b1e2iii_rface_p0_scout.json"
ART = RES / "k3_cap_b1e2iii_rface_p0_b1_leaf.json"

# --- PRÉ-ENREGISTRÉ -------------------------------------------------------
RE_DIR = 0                     # la direction traversée : Re u
PREDICTED_DIFF_COORD = 5       # les candidates ne diffèrent que sur s=5
THETA_REQUIRED = 1
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:6.1f}s] {msg}", flush=True)


def _sha(p):
    try:
        return hashlib.sha256(Path(p).read_bytes()).hexdigest()
    except OSError:
        return None


def _q(x):
    return [int(x.numerator), int(x.denominator)]


def _qbox(bx):
    return [[_q(a), _q(b)] for a, b in bx]


def _frombox(qb):
    return [(Fraction(*a), Fraction(*b)) for a, b in qb]


def mirror_re(bx, face, k=RE_DIR):
    """Le miroir à travers l'hyperplan `x_k = face` — Fractions exactes."""
    return [((2 * face - bx[j][1], 2 * face - bx[j][0]) if j == k
             else bx[j]) for j in range(4)]


def widen_across(bx, face, w, k=RE_DIR):
    """La boîte élargie à `[face − w, face + w]` dans la direction k."""
    return [((face - w, face + w) if j == k else bx[j])
            for j in range(4)]


def all_dyadic_floats(bx):
    """Chaque borne doit convertir en float SANS PERTE (dyadique) : les
    constructeurs TM reçoivent des floats, l'exactitude est un gate."""
    return all(Fraction(float(x)) == x for ab in bx for x in ab)


def classify_candidates(closed_flags):
    """Le gate à trois branches pré-enregistré de la 14ᵉ revue §1-B1."""
    n = sum(1 for v in closed_flags.values() if v)
    if n == 1:
        return "exactly_one_closes"
    return "both_close_not_discriminated" if n >= 2 else "none_closes"


def section_or_none(S, g, eps, box):
    """Construit la section bilatérale sur `box` ; None si incomplète."""
    c, h = center_hw(box)
    Z, dZ, rows = build_section_bilateral(
        S, g, eps, [float(x) for x in c], [float(x) for x in h])
    return (Z if all(z is not None for z in Z) else None), rows


def theta_summary(res, keys):
    """Résumé d'un theta_lines : tout +1 ?, marge min, motif par coord."""
    if "refused" in res:
        return {"closed": False, "refused": res["refused"]}
    pattern = {str(k): res[k].get("theta") for k in keys}
    closed = all(res[k].get("theta") == THETA_REQUIRED for k in keys)
    margs = [res[k]["margin"] for k in keys
             if res[k].get("theta") in (1, -1)
             and res[k].get("margin") is not None]
    return {"closed": bool(closed), "pattern": pattern,
            "min_margin": min(margs) if margs else None}


# ===========================================================================
#  Self-test
# ===========================================================================
def selftest():
    ok, tot = 0, 0

    def chk(name, cond):
        nonlocal ok, tot
        tot += 1
        ok += bool(cond)
        print(f"  {'OK  ' if cond else 'FAIL'} T{tot} {name}")

    f = Fraction(-183, 256)
    bx = [(f - Fraction(1, 512), f)] + [(Fraction(-1, 512),
                                         Fraction(1, 512))] * 3
    m = mirror_re(bx, f)
    chk("mirror_re réfléchit la seule direction Re u",
        m[0] == (f, f + Fraction(1, 512)) and m[1:] == bx[1:])
    chk("mirror_re est involutif", mirror_re(m, f) == bx)
    w = widen_across(bx, f, Fraction(1, 512))
    chk("widen_across pose [face−w, face+w]",
        w[0] == (f - Fraction(1, 512), f + Fraction(1, 512)))
    ok_in, _ = contains(w, bx, strict_dirs=())
    ok_mi, _ = contains(w, m, strict_dirs=())
    chk("le pont-Re contient le côté ET son miroir", ok_in and ok_mi)
    narrow = widen_across(bx, f, Fraction(1, 1024))
    chk("NÉGATIF : la boîte étroite échoue DES DEUX CÔTÉS",
        not contains(narrow, bx, strict_dirs=())[0]
        and not contains(narrow, m, strict_dirs=())[0])
    chk("classify : une seule ferme",
        classify_candidates({"0": True, "1": False})
        == "exactly_one_closes")
    chk("classify : les deux ⟹ refus non discriminé",
        classify_candidates({"0": True, "1": True})
        == "both_close_not_discriminated")
    chk("classify : aucune ⟹ résultat structurel",
        classify_candidates({"0": False, "1": False}) == "none_closes")
    chk("all_dyadic_floats accepte le dyadique, refuse 1/3",
        all_dyadic_floats(bx)
        and not all_dyadic_floats([(Fraction(1, 3), Fraction(1, 2))] * 4))
    print(f"\nself-test {ok}/{tot}")
    return ok == tot


# ===========================================================================
#  Build
# ===========================================================================
def build():
    print("=" * 78)
    print("RFace-P0 / B1 — la feuille atteinte à travers la face Re u "
          "haute, DÉRIVÉE")
    print("=" * 78)
    cov = json.loads(COVER_JSON.read_text(encoding="utf-8"))
    atl = json.loads(ATLAS_JSON.read_text(encoding="utf-8"))
    c127e = json.loads(C127E_JSON.read_text(encoding="utf-8"))
    scout0 = json.loads(SCOUT0_JSON.read_text(encoding="utf-8"))
    f23 = json.loads(F2F3_JSON.read_text(encoding="utf-8"))
    rs = json.loads(RSCOUT_JSON.read_text(encoding="utf-8"))

    # --- B1a : amont — C127-E inclus (15ᵉ revue §5 : il fournit le
    # chart cible S₂/g₂, dépendance load-bearing, donc GATÉE) ----------
    up = {}
    for name, blob in (("rface_scout_v2", rs), ("f0_scout", scout0),
                       ("f2f3_v2", f23), ("c127d", atl),
                       ("c127e_residual", c127e)):
        gp, gt = blob.get("gates_passed"), blob.get("gates_total")
        mode = blob.get("mode")
        up[name] = {"gates": f"{gp}/{gt}", "mode": mode,
                    "green": bool(gp == gt and gt
                                  and (mode == "full"
                                       if mode is not None else True))}
    up["c127e_residual"]["green"] = bool(
        up["c127e_residual"]["green"]
        and len(c127e.get("new_tiles", [])) == 64)
    b1a = all(v["green"] for v in up.values())
    log("B1a : amont — " + " ; ".join(
        f"{k} {v['gates']}" for k, v in up.items()) + f" ⟹ {b1a}")

    cell = cov["cell"]
    S, g = tuple(cell["S"]), cell["g"]
    tile = rs["witnesses"]["w_one_neighbor"]["tile"]
    face = Fraction(*rs["faces"]["u_hi"]["face_value"])
    nb_box = _frombox(rs["address_square"]["u_neighbor"]["box"])
    cands = {c: rs["address_square"]["u_neighbor"]["classes"][str(c)]
             ["eps_registered_c118"] for c in rs["candidate_classes"]}
    f0_box = [(Fraction(*b[0]), Fraction(*b[1]))
              for b in [t for t in scout0["per_tile"]
                        if t["tile"] == tile][0]
              ["bridge_corrected_bounds"]]
    pb = [x for x in f23["per_bridge"] if x["tile"] == tile][0]
    eps_f0 = tuple(pb["F2d_ledger_derived"])
    regimes_f0 = [r["regime"] for r in pb["bridge_rows"]]
    chart = c127e["new_tiles"][tile - 252]["chart"]
    S2, g2 = tuple(chart["S"]), chart["g"]
    log(f"témoin {tile} (importé du scout v2), face Re u = {face} "
        f"({float(face):+.9f}), ledger F0 {eps_f0}, chart cible "
        f"S₂={S2} g₂={g2}")

    # --- B1b : géométrie 2H ------------------------------------------
    w0 = face - f0_box[RE_DIR][0]
    flush = f0_box[RE_DIR][1] == face
    # 15ᵉ revue §5 : si le texte appelle la largeur « 2H », le symbole
    # doit être GATÉ contre le H sérialisé de l'atlas, pas affirmé.
    H_halo = Fraction(float.fromhex(
        {h["index"]: h["record"] for h in atl["halos"]
         if h["ok"]}[tile]["H_hex"]))
    w0_is_2H = (w0 == 2 * H_halo)
    bridge = widen_across(f0_box, face, w0)
    mirror = mirror_re(f0_box, face)
    inc_f0, marg_f0 = contains(bridge, f0_box, strict_dirs=())
    inc_mi, marg_mi = contains(bridge, mirror, strict_dirs=())
    attained = (bridge[RE_DIR][0] == f0_box[RE_DIR][0]
                and bridge[RE_DIR][1] == mirror[RE_DIR][1])
    narrow = widen_across(f0_box, face, w0 / 2)
    neg_narrow = (not contains(narrow, f0_box, strict_dirs=())[0]
                  and not contains(narrow, mirror, strict_dirs=())[0])
    dyadic = all_dyadic_floats(bridge) and all_dyadic_floats(f0_box)
    b1b = bool(flush and inc_f0 and inc_mi and attained and neg_narrow
               and dyadic and w0 > 0 and w0_is_2H)
    log(f"B1b : w₀ = {w0} == 2·H sérialisé ({w0_is_2H}, H = {H_halo}), "
        f"pont-Re ⊇ F0 ∪ miroir (bornes atteintes {attained}), négatif "
        f"étroit {neg_narrow}, dyadique {dyadic} ⟹ {b1b}")

    # --- B1c : section sur la boîte élargie --------------------------
    keys = list(range(6))
    Zw, rows_w = section_or_none(S, g, eps_f0, bridge)
    regimes_w = [r.get("regime") for r in rows_w]
    rad_min = min((r["radicand_absmin"] for r in rows_w
                   if r.get("radicand_absmin") is not None),
                  default=None)
    gmin = None
    if Zw is not None:
        gr, gi = Zw[g2].to_iv_pair()
        gmin = _f_down(mp.mpf(civ_absmin(CIV(gr, gi)).a))
    b1c = bool(Zw is not None and regimes_w == regimes_f0
               and rad_min is not None and rad_min > 0
               and gmin is not None and gmin > 0)
    log(f"B1c : régimes élargis {regimes_w} == F0 sérialisés "
        f"{regimes_f0} ; radicande min {rad_min} ; jauge cible min "
        f"{gmin} ⟹ {b1c}")

    # --- B1d : restriction au côté source ----------------------------
    Zf0, rows_f0 = section_or_none(S, g, eps_f0, f0_box)
    b1d = False
    th_src = None
    if Zw is not None and Zf0 is not None:
        regimes_f0_re = [r.get("regime") for r in rows_f0]
        res = theta_lines(Zw, bridge, Zf0, f0_box, f0_box, keys)
        th_src = theta_summary(res, keys)
        b1d = bool(th_src["closed"]
                   and regimes_f0_re == regimes_f0)
    log(f"B1d : restriction au pont F0 — {th_src} ⟹ {b1d}")

    # --- B1e : W_cell et ledger dérivé -------------------------------
    W = inter(bridge, nb_box)
    b1e = False
    eps_der = None
    th_der = None
    W_widths = None
    derive_th = None
    if W is not None:
        cW, hW = center_hw(W)
        W_widths = [b - a for a, b in W]
        is_cube = len(set(W_widths)) == 1
        anchor_interior = cW[RE_DIR] > face
        Zp, rows_p = section_or_none(S, g, (1, 1, 1), W)
        if Zp is not None and Zw is not None:
            res_p = theta_lines(Zp, W, Zw, bridge, W, keys)
            derive_th = theta_summary(res_p, keys)
            per_line = ({int(s): res_p[int(s)].get("theta") for s in S}
                        if "refused" not in res_p else {})
            if all(t in (1, -1) for t in per_line.values()):
                eps_der = tuple(per_line[int(s)] for s in S)
                Zd, rows_d = section_or_none(S, g, eps_der, W)
                if Zd is not None:
                    res_d = theta_lines(Zd, W, Zw, bridge, W, keys)
                    th_der = theta_summary(res_d, keys)
                    b1e = bool(th_der["closed"] and is_cube
                               and anchor_interior)
    log(f"B1e : W_cell largeurs {[str(x) for x in (W_widths or [])]} "
        f"(cube), ledger DÉRIVÉ {eps_der}, reconstruit ⟹ "
        f"{th_der} ⟹ {b1e}")

    # --- B1f : les deux candidates -----------------------------------
    cand_out = {}
    closed = {}
    if W is not None and Zw is not None:
        for cls, eps_c in sorted(cands.items()):
            Zc, rows_c = section_or_none(S, g, tuple(eps_c), W)
            if Zc is None:
                cand_out[str(cls)] = {
                    "eps_registered": eps_c,
                    "refused": "candidate_section_incomplete",
                    "rows": [{"s_coord": r.get("s_coord"),
                              "refused": r.get("refused")}
                             for r in rows_c]}
                closed[str(cls)] = False
                continue
            res_c = theta_lines(Zc, W, Zw, bridge, W, keys)
            s = theta_summary(res_c, keys)
            cand_out[str(cls)] = {
                "eps_registered": eps_c,
                "regimes": [r.get("regime") for r in rows_c], **s}
            closed[str(cls)] = s["closed"]
    outcome = classify_candidates(closed) if closed else "not_run"

    # --- B1f : le CONTRÔLE, puis l'identification RELATIVE -----------
    # Le contrôle : NOTRE étiquette de registre, sur NOTRE côté du pont
    # (le cube symétrique W_ours = miroir-Re de W_cell). Si elle non
    # plus ne ferme pas naïvement — et avec quel motif — c'est mesuré,
    # pas argumenté.
    control = None
    ctrl_pattern = None
    W_ours = mirror_re(W, face) if W is not None else None
    if W_ours is not None and Zw is not None:
        Zo, rows_o = section_or_none(S, g, tuple(cell["eps"]), W_ours)
        if Zo is not None:
            res_o = theta_lines(Zo, W_ours, Zw, bridge, W_ours, keys)
            control = theta_summary(res_o, keys)
            control["eps_registered_ours"] = list(cell["eps"])
            control["regimes"] = [r.get("regime") for r in rows_o]
            ctrl_pattern = control.get("pattern")
    # B1f1 — κ DÉRIVÉ du contrôle source (15ᵉ revue §3) : le ledger
    # canonique source est le ledger F0 certifié ; κ = ε_can ⊙ ε_reg,
    # et le motif MESURÉ du contrôle doit le confirmer ligne à ligne.
    kappa = tuple(int(eps_f0[r]) * int(cell["eps"][r])
                  for r in range(3))
    b1f1 = bool(
        ctrl_pattern is not None
        and all(ctrl_pattern.get(str(int(s))) == kappa[r]
                for r, s in enumerate(S))
        and all(ctrl_pattern.get(str(k)) == 1
                for k in keys if k not in [int(s) for s in S]))
    # B1f2 — κ CONCORDE avec les lois de détermination, et le gate
    # PORTE la loi qu'il nomme (durci, 16ᵉ revue B3-D1) :
    #     principal → principal : κ = +1
    #     rotated   → canonique : κ = −1
    # Les déterminations du CORE (C127-E) sont importées, pas récitées.
    core_dets = [t for t in c127e["transports"]
                 if t.get("box_index") == tile - 252]
    core_src = (core_dets[0].get("source_determinations")
                if core_dets else None)
    b1f2 = bool(
        core_src is not None
        and all((regimes_w[r] == "principal" and kappa[r] == 1
                 and core_src[r] == "principal")
                or (regimes_w[r] == "canonical" and kappa[r] == -1
                    and core_src[r] == "rotated")
                for r in range(3)))
    # B1f3 — les DEUX étiquettes CONVERTIES explicitement -------------
    eps_conv = {c: tuple(kappa[r] * int(e[r]) for r in range(3))
                for c, e in cands.items()}
    # B1f4/B1f5 — les candidates converties sont RECONSTRUITES --------
    conv_out = {}
    conv_closed = {}
    if W is not None and Zw is not None:
        for cls, e_c in sorted(eps_conv.items()):
            Zc, rows_c = section_or_none(S, g, e_c, W)
            if Zc is None:
                conv_out[str(cls)] = {"eps_converted": list(e_c),
                                      "refused": "section_incomplete"}
                conv_closed[str(cls)] = False
                continue
            s_ = theta_summary(theta_lines(Zc, W, Zw, bridge, W, keys),
                               keys)
            conv_out[str(cls)] = {"eps_converted": list(e_c), **s_}
            conv_closed[str(cls)] = s_["closed"]
    n_conv = sum(1 for v in conv_closed.values() if v)
    winner = ([c for c, v in conv_closed.items() if v][0]
              if n_conv == 1 else None)
    loser = ([c for c, v in conv_closed.items() if not v][0]
             if n_conv == 1 and len(conv_closed) == 2 else None)
    b1f4 = bool(n_conv == 1)
    loser_pred = False
    if loser is not None and "pattern" in conv_out[loser]:
        pat = conv_out[loser]["pattern"]
        bad = [k for k, v in pat.items() if v != 1]
        loser_pred = bad == [str(PREDICTED_DIFF_COORD)]
    b1f5 = bool(loser_pred)
    # cohérence : la gagnante convertie doit être le ledger DÉRIVÉ
    winner_is_derived = bool(
        winner is not None and eps_der is not None
        and eps_conv[int(winner)] == eps_der)
    b1f = bool(b1f1 and b1f2 and b1f4 and b1f5 and winner_is_derived)
    log(f"B1f′ : κ = {kappa} (dérivé du contrôle : {b1f1} ; confiné "
        f"aux lignes rotated→canonique du core : {b1f2}) ; converties "
        f"{ {c: list(v) for c, v in eps_conv.items()} } ; EXACTEMENT "
        f"UNE ferme : classe {winner} ({b1f4}, == ledger dérivé : "
        f"{winner_is_derived}) ; la perdante casse UNIQUEMENT sur "
        f"{PREDICTED_DIFF_COORD} : {b1f5} ; gate naïf = {outcome} "
        f"(diagnostic conservé) ⟹ {b1f}")

    # --- B1g : NÉGATIF — muter UNE ligne du ledger dérivé ------------
    b1g = False
    neg_mut = None
    if eps_der is not None and W is not None and Zw is not None:
        e_bad = list(eps_der)
        e_bad[1] = -e_bad[1]
        Zb_, _rows = section_or_none(S, g, tuple(e_bad), W)
        if Zb_ is not None:
            res_b = theta_lines(Zb_, W, Zw, bridge, W, keys)
            s_bad = theta_summary(res_b, keys)
            mut_coord = int(S[1])
            neg_mut = {"mutated_line_s": mut_coord,
                       "eps_mutated": e_bad, **s_bad}
            b1g = bool(not s_bad["closed"]
                       and s_bad.get("pattern", {}).get(str(mut_coord))
                       == -1)
    log(f"B1g : mutation d'une ligne du ledger ⟹ {neg_mut} ⟹ {b1g}")

    b1h = bool(Zw is not None and Zf0 is not None and W is not None
               and len(closed) == len(cands) == 2)
    log(f"B1h : 6 coordonnées partout, 2 candidates traitées ⟹ {b1h}")

    gates = {
        "B1a_upstream_green_full_modes": bool(b1a),
        "B1b_geometry_2H_exact_negative_narrow": bool(b1b),
        "B1c_section_on_widened_box_chart_revalidated": bool(b1c),
        "B1d_restriction_equals_certified_f0_section": bool(b1d),
        "B1e_neighbor_ledger_derived_and_rebuilt": bool(b1e),
        "B1f_explicit_label_to_ledger_conversion_identifies": bool(b1f),
        "B1g_single_line_mutation_breaks": bool(b1g),
        "B1h_no_silent_filtering": bool(b1h)}
    npass = sum(1 for v in gates.values() if v)

    try:
        head = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=Path(__file__).resolve().parent, capture_output=True,
            text=True, timeout=10).stdout.strip() or None
    except Exception:
        head = None

    out = {
        "artifact": "k3_cap_b1e2iii_rface_p0_b1_leaf",
        "mode": "full",
        "contract": ("phase B1 du contrat RFace-P0 — 13ᵉ revue §7.3-B, "
                     "amendée 14ᵉ revue §1/§4/§5 ; témoin unique "
                     "w_one_neighbor"),
        "claim": (
            "Sur le témoin 295 (flush face Re u haute), la section du "
            "pont-Re 2H — construite avec le ledger F0 certifié, régimes "
            "identiques, chart re-validé sur la boîte élargie — se "
            "restreint EXACTEMENT à la section F0 certifiée (θ=+1, 6/6) "
            "et TRAVERSE la face : sur le cube W_cell = pont ∩ cellule "
            "voisine, le ledger dérivé ligne à ligne se reconstruit et "
            "referme (6/6). AUCUNE étiquette enregistrée ne ferme "
            "naïvement (`none_closes`, y compris la NÔTRE sur notre "
            "côté : le contrôle le mesure) — l'identification est donc "
            "PAR CONVERSION EXPLICITE (B1f′, 15ᵉ revue §3) : "
            "κ = ε_canonique ⊙ ε_registre est dérivé du contrôle "
            "source, confiné aux lignes rotated→canonique du core, les "
            "DEUX étiquettes sont converties et RECONSTRUITES — "
            "exactement UNE ferme (classe 0, == ledger dérivé), la "
            "classe 1 convertie casse uniquement sur la coordonnée 5 "
            "prédite. Après conversion, AUCUNE transformation de deck "
            "SUPPLÉMENTAIRE n'est observée à cette face. Ceci est une "
            "PREMIÈRE CONTINUATION ANALYTIQUE LOCALE entre deux "
            "cellules canoniquement adressées — pas un atlas de la "
            "voisine : ni transition d'atlas (C), ni métrique (E), ni "
            "315/B3 ne sont payés ici ; le gate naïf `none_closes` "
            "reste publié comme diagnostic séparé."),
        "witness": {"tile": tile, "face_re_u": _q(face),
                    "chart_target": {"S": list(S2), "g": g2}},
        "cell": {"S": list(S), "g": g, "eps": list(cell["eps"])},
        "geometry": {
            "w0": _q(w0), "bridge_re": _qbox(bridge),
            "f0_bridge": _qbox(f0_box), "mirror_re_f0": _qbox(mirror),
            "W_cell": _qbox(W) if W else None,
            "W_cell_widths": [_q(x) for x in (W_widths or [])],
            "W_cell_is_cube": bool(W_widths
                                   and len(set(W_widths)) == 1)},
        "ledger_f0": list(eps_f0),
        "regimes_widened": regimes_w,
        "restriction_to_f0": th_src,
        "neighbor_ledger": {
            "derivation_probe_eps_plus": derive_th,
            "eps_derived": list(eps_der) if eps_der else None,
            "rebuilt": th_der},
        "candidates": cand_out,
        "naive_closure_outcome": outcome,
        "naive_closure_note": (
            "AUCUNE étiquette enregistrée ne ferme naïvement — Y "
            "COMPRIS la nôtre sur notre propre côté (voir `control`) : "
            "la fermeture naïve d'étiquette n'est pas un critère "
            "d'identification, c'est le contrôle qui le prouve. "
            "L'étiquette registre vit dans la convention du registre ; "
            "le ledger du pont dans celle des régimes signe-certifiés "
            "(rotated→canonique absorbe σ). Le motif du contrôle "
            "MESURE la conversion, et l'identification est RELATIVE."),
        "control": control,
        "conversion": {
            "kappa": list(kappa),
            "kappa_derived_from_control_pattern": bool(b1f1),
            "kappa_confined_to_regime_changed_lines": bool(b1f2),
            "core_source_determinations": core_src,
            "eps_converted_by_class": {str(c): list(v)
                                       for c, v in eps_conv.items()},
            "note": (
                "κ = ε_canonique_source ⊙ ε_registre_source, ligne à "
                "ligne (15ᵉ revue §3). La conversion est CONFINÉE aux "
                "lignes dont la détermination a changé (rotated au "
                "core → canonique sur le pont) ; κ = +1 sur la ligne "
                "restée principale.")},
        "candidates_converted": conv_out,
        "identified_class": int(winner) if winner is not None else None,
        "winner_converted_equals_derived_ledger": bool(winner_is_derived),
        "deck_statement": (
            "Après conversion explicite de l'étiquette de registre "
            "vers le ledger canonique, la continuation conserve la "
            "classe 0 et AUCUNE transformation de deck SUPPLÉMENTAIRE "
            "n'est observée à cette face. Le motif BRUT du contrôle, "
            "étendu aux six coordonnées, est (+1,−1,+1,+1,+1,−1) — le "
            "même diagonal que D de C129-F0 : ici il corrige "
            "l'injection naïve d'une étiquette d'une autre convention, "
            "il n'est PAS un changement de feuille après conversion. "
            "La phrase non qualifiée « pas de deck à cette face » "
            "était trop forte (15ᵉ revue §4)."),
        "predicted_diff_coord": PREDICTED_DIFF_COORD,
        "negative_mutation": neg_mut,
        "conventions_note": (
            "Les ε enregistrés (C118) sont des étiquettes de classe "
            "dans la convention du registre ; les candidates sont "
            "construites sous la règle de régime signe-certifié "
            "(convention C127-E des régions clippées, celle de l'arc "
            "F0). L'identification (B1f) et la dérivation (B1e) sont "
            "publiées côte à côte ; toute discordance serait un "
            "résultat, pas une correction."),
        "not_paid_here": [
            "B3 : témoin 315, deux ordres, commutation",
            "C : transitions d'atlas et nerf étendu",
            "E : métrique (transport_hardened, δ = 1e-5)",
            "l'atlas de la cellule voisine",
            "les faces basses (frontière d'énumération)",
            "les 895 autres paires"],
        "gates": gates, "gates_passed": npass, "gates_total": len(gates),
        "verdict": (
            f"B1 v2 LIVRÉ — première continuation analytique locale "
            f"entre deux cellules canoniquement adressées : la feuille "
            f"atteinte est la classe {winner}, identifiée PAR "
            f"CONVERSION EXPLICITE (κ dérivé, candidates converties "
            f"reconstruites, une seule ferme) ; aucune deck "
            f"SUPPLÉMENTAIRE après conversion."
            if npass == len(gates) else
            f"ROUGE — {len(gates) - npass} gate(s) en échec"),
        "provenance": {
            "git_head": head, "python": sys.version.split()[0],
            "platform": platform.platform(), "mp_prec": int(mp.prec),
            "wall_s": round(time.time() - T0, 1),
            "inputs": {p.name: _sha(p) for p in
                       (COVER_JSON, ATLAS_JSON, C127E_JSON, SCOUT0_JSON,
                        F2F3_JSON, RSCOUT_JSON)},
            "self_sha256": _sha(__file__)}}

    ART.write_text(json.dumps(out, indent=2, ensure_ascii=False),
                   encoding="utf-8")
    print("=" * 78)
    for k, v in gates.items():
        print(f"  {'OK  ' if v else 'FAIL'} {k}")
    print(f"\n{out['verdict']}")
    print(f"gates {npass}/{len(gates)} — artefact : {ART.name}")
    print("=" * 78)
    return npass == len(gates)


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        sys.exit(0 if selftest() else 1)
    sys.exit(0 if build() else 1)
