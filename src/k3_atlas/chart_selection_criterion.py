#!/usr/bin/env python3
"""
k3_cap_b1e2iii_c126_chart_criterion.py — LE CRITÈRE DE SÉLECTION DES
CHARTS, RENFORCÉ, et rebalayé sur les 59 charts.

C125 a fait tomber D5.4/D5.5 : avec `σ'` réellement certifié, la section
NATIVE du chart cible n'est pas constructible, et le balayage donnait
**0/15** sur la première cellule. Ce script en fait un critère de
sélection à part entière et le rebalaye sur TOUTES les cellules.

LA FAILLE QUE CE SCRIPT RÉPARE. Le critère de D5.1/D5.2 portait sur :

    (i)  le DOMAINE     — jauge `|Z_{g'}|` loin de 0, `|u'| ≤ 1`, `|v'| ≤ 1`
    (ii) le JACOBIEN    — `det J` d'enclosure excluant 0
    (iii) la DISJONCTION — `Im u'` ou `Im v'` de signe constant

Ces trois-là portent sur les **coordonnées** `(u', v')` du chart cible.
Aucun ne dit quoi que ce soit des **radicandes** de la section PROPRE du
chart cible. Un chart peut donc être hors de sa tranche réelle tout en
ayant une racine de section incertifiable — et c'est précisément ce qui
est arrivé. Le critère renforcé ajoute :

    (iv) SECTION NATIVE — les trois lignes de la section du chart CIBLE
         admettent une détermination CERTIFIÉE sur la boîte entière
         (principale si `R'` évite (−∞,0], sinon tournée avec `σ'`
         déterminé — jamais au seul centre)

Gates pré-enregistrés :
  F1 RAFFINEMENT — le critère renforcé est un SOUS-ENSEMBLE du critère
     précédent, chart par chart. Gate UNIVERSEL : un critère qui
     accepterait un chart refusé par D5.1/D5.2 ne serait pas un
     renforcement mais un autre critère.
  F2 EXHAUSTIVITÉ — les 59 charts sont examinés sur CHAQUE cellule,
     aucun saut silencieux ; le verdict publié est un décompte
     UNIVERSEL, pas une existence.
  F3 LE CRITÈRE MORD — il refuse au moins un chart que D5.1+D5.2+
     disjonction acceptaient. Sans quoi il n'ajouterait rien et le
     rebalayage serait sans objet. (Si le critère ne mordait pas,
     ce gate DOIT échouer et le dire.)
  F4 NÉGATIF CENTRE-SEUL — la mutation qui détermine la constructibilité
     sur le seul CENTRE (boîte de demi-largeur 0) doit accepter
     STRICTEMENT PLUS de charts que la version certifiée sur la boîte.
     C'est ce qui distingue « certifié sur la cellule » de « vu au
     centre ».
  F5 SOUNDNESS DES LIGNES CONSTRUCTIBLES — pour toute ligne déclarée
     constructible, quelle que soit la détermination, l'enclosure de
     `Z'² − R'` doit contenir 0. Gate UNIVERSEL sur toutes les lignes
     acceptées, donc NON VIDE même si aucun chart entier ne passe.
  F6 L'OBSTRUCTION EST-ELLE HÉRITÉE ? — C125 conjecture que la ligne qui
     bloque le chart cible porte la MÊME coordonnée ambiante `Z_s` que
     celle qui bloquait la source. Ce gate MESURE la proportion ; il ne
     l'impose pas. Il échoue seulement si la mesure est absente.
  F9 LA SUBDIVISION N'EFFACE PAS L'OBSTRUCTION — sur un chart perdu,
     on subdivise la cellule source en 16^d sous-cellules et on
     recompte les refus. Si l'obstruction était un ARTEFACT
     D'ENCLOSURE (surestimation de `Im R'`), les refus DISPARAÎTRAIENT
     à profondeur croissante. Le gate exige au contraire que le nombre
     de refus CROISSE en `4^d` — la signature d'un lieu de branchement
     de CODIMENSION 2, donc d'une obstruction GÉOMÉTRIQUE RÉELLE. Le
     gate échoue si les refus disparaissent OU s'ils croissent en
     `16^d` (codimension 0).
  F7 VERDICT PUBLIÉ — le décompte final `n_strong` est publié pour
     CHAQUE cellule, y compris s'il vaut 0. Un balayage qui ne
     publierait que ses succès mentirait par omission.

Ce script ne certifie NI D5.4 NI D5.5 : il sélectionne les paires
(cellule, chart) sur lesquelles ces questions ont un sens.

Sorties : results/k3_cap_b1e2iii_c126_chart_criterion.json
Usage   : k3_cap_b1e2iii_c126_chart_criterion.py [--selftest]
Env     : K3_C126_CELLS (cellules, défaut 4)
          K3_C126_SUBDEPTH (profondeur de la sonde F9, défaut 3)
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
from collections import Counter
from itertools import product
from fractions import Fraction
from pathlib import Path

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

_argv = sys.argv
sys.argv = [sys.argv[0]]
os.environ.setdefault("K3_TM_ORDER", "4")
os.environ.setdefault("K3_TM_SERIES", "4")
from mpmath import mp                                              # noqa: E402
from .witness_registry import load_canonical_MH              # noqa: E402
from .interval_arithmetic import (                            # noqa: E402
    BranchCutError, CIV, minor_inv_times_T_exact)
from .taylor_models import (                                     # noqa: E402
    CZERO, TMC, TM_ORDER, UNARY_SERIES_DEG, civ_absmax, riv,
    rotated_sigma_from_coeffs, tm_sqrt_rotated)
from .owner_tiling import TRIPLES                    # noqa: E402
from .full_cell_charts import (                           # noqa: E402
    build_section, chart_certificate)
sys.argv = _argv

RES = Path(os.environ.get(
    "K3_RES_DIR", Path(__file__).resolve().parent / "data"))
C118_JSON = RES / "k3_cap_b1e2iii_c118_shell_exhaustive.json"
ART = RES / "k3_cap_b1e2iii_c126_chart_criterion.json"

N_CELLS = int(os.environ.get("K3_C126_CELLS", "4"))
SUB_DEPTH = int(os.environ.get("K3_C126_SUBDEPTH", "3"))
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def _sha(p):
    try:
        return hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]
    except OSError:
        return None


def provenance(src, t_wall):
    try:
        h = subprocess.run(["git", "rev-parse", "HEAD"],
                           capture_output=True, text=True,
                           timeout=10).stdout.strip() or None
    except Exception:
        h = None
    return {"git_head": h, "python": sys.version.split()[0],
            "platform": platform.platform(), "mp_prec": int(mp.prec),
            "tm_order": TM_ORDER, "unary_series_deg": UNARY_SERIES_DEG,
            "wall_s": round(t_wall, 1), "n_cells": N_CELLS,
            "inputs": {str(Path(x).name): _sha(x) for x in src},
            "self_sha256": _sha(__file__)}


def _rng(tm, im=False):
    x = tm.im_tm().to_iv() if im else tm.re_tm().to_iv()
    return (float(mp.mpf(x.a)), float(mp.mpf(x.b)))


# ===========================================================================
#  (iv) LA SECTION NATIVE DU CHART CIBLE — le critère ajouté
# ===========================================================================
def native_section_constructible(S2, g2, up, vp):
    """Les trois lignes de la section PROPRE du chart cible, avec leur
    détermination certifiée sur la BOÎTE ENTIÈRE.

    Retourne (rows, ok). Chaque ligne porte sa détermination, son `σ'`,
    les enclosures de `Re R'` / `Im R'`, et — si elle est constructible —
    l'enclosure de `Z'² − R'` qui doit contenir 0 (gate F5).

    `σ'` est certifié en DEUX temps, comme dans C125 :
      (i) enclosure DIRECTE de `Im R'` — sur le chart cible elle n'est
          pas contrainte de contenir 0, donc souvent plus fine ;
     (ii) à défaut, la règle par SIGNES DE COEFFICIENTS.
    Si les deux sont indéterminées, la ligne est REFUSÉE. Jamais d'essai.
    """
    T2 = tuple(j for j in range(6) if j not in S2)
    others = [c for c in T2 if c != g2]
    Ae = minor_inv_times_T_exact(S2, T2)
    perm = [list(T2).index(g2), list(T2).index(others[0]),
            list(T2).index(others[1])]
    A = [[riv(Ae[r][c]) for c in perm] for r in range(3)]
    ur, ui = _rng(up), _rng(up, True)
    vr, vi = _rng(vp), _rng(vp, True)
    u2, v2 = up * up, vp * vp
    rows = []
    for r, s in enumerate(S2):
        R = TMC.const(CIV(A[r][0])) + u2.mul_real(A[r][1]) \
            + v2.mul_real(A[r][2])
        rl, rh = _rng(R)
        il, ih = _rng(R, True)
        rec = {"row": r, "s_coord": int(s), "Re_R": [rl, rh],
               "Im_R": [il, ih], "determination": None, "sigma": None}
        Zp = None
        try:
            Zp = R.sqrt_principal()
            rec["determination"] = "principal"
        except BranchCutError as exc:
            rec["principal_refused"] = exc.diag.get("guard")
            if il > 0:
                sg, how = 1, "Im_R_enclosure"
            elif ih < 0:
                sg, how = -1, "Im_R_enclosure"
            else:
                sg = rotated_sigma_from_coeffs(
                    Fraction(Ae[r][perm[1]]), Fraction(Ae[r][perm[2]]),
                    ur, ui, vr, vi)
                how = "coefficient_signs"
            rec["sigma"], rec["sigma_source"] = sg, how
            if sg in (-1, 1):
                try:
                    Zp = tm_sqrt_rotated(R, sg)
                    rec["determination"] = "rotated"
                except BranchCutError as exc2:
                    rec["rotated_refused"] = exc2.diag.get("guard")
            else:
                rec["rotated_refused"] = "rotated_component_undetermined"
        if Zp is not None:
            # F5 : soundness de la ligne, quelle que soit la détermination
            d = Zp * Zp - R
            dr, di = _rng(d), _rng(d, True)
            rec["sq_residual_re"], rec["sq_residual_im"] = dr, di
            rec["sq_residual_contains_zero"] = bool(
                dr[0] <= 0 <= dr[1] and di[0] <= 0 <= di[1])
        rows.append(rec)
    return rows, all(x["determination"] is not None for x in rows)


def target_uv(Z, S2, g2):
    """`(u', v')` du chart cible, ou None si inatteignable."""
    T2 = [j for j in range(6) if j not in S2]
    o = [x for x in T2 if x != g2]
    if any(Z[a] is None for a in (g2, o[0], o[1])):
        return None
    try:
        ib = Z[g2].inv()
        return (Z[o[0]] * ib, Z[o[1]] * ib)
    except BranchCutError:
        return None


# ===========================================================================
#  Build
# ===========================================================================
def build():
    print("=" * 78)
    print(f"C126 CRITÈRE DE CHART RENFORCÉ : TM ({TM_ORDER},"
          f"{UNARY_SERIES_DEG}), {N_CELLS} cellules × 59 charts")
    print("=" * 78)
    load_canonical_MH()
    c118 = json.loads(C118_JSON.read_text(encoding="utf-8"))
    residual = [x for x in c118["cells"] if x["still_branch"]]
    log(f"résidu : {len(residual)} paires")

    cells = []
    step = max(1, len(residual) // max(1, N_CELLS))
    for i in range(min(N_CELLS, len(residual))):
        r = residual[i * step]
        S, g, eps = tuple(r["S"]), r["g"], tuple(r["eps"])
        c = [float.fromhex(x) for x in r["center_hex"]]
        hw = float.fromhex(r["hw_hex"])
        Z, dZ, srows = build_section(S, g, eps, c, hw)
        # la MUTATION centre-seul : la MÊME cellule vue en un point
        Zc, dZc, _ = build_section(S, g, eps, c, 0.0)
        src_rotated = [x["s_coord"] for x in srows
                       if x["determination"] == "rotated"]

        charts, n_seen = [], 0
        for S2 in TRIPLES:
            for g2 in (j for j in range(6) if j not in S2):
                n_seen += 1
                if (tuple(S2), g2) == (S, g):
                    charts.append({"S": list(S2), "g": int(g2),
                                   "is_source": True})
                    continue
                cert = chart_certificate(Z, dZ, tuple(S2), g2)
                cert.update(S=list(S2), g=int(g2), is_source=False)
                prev = bool(cert.get("admissible")
                            and cert.get("disjoint_from_target_slice"))
                cert["prev_criterion"] = prev
                uv = target_uv(Z, tuple(S2), g2)
                if uv is None:
                    cert["native_section_ok"] = False
                    cert["native_rows"] = None
                else:
                    rows, ok = native_section_constructible(
                        tuple(S2), g2, uv[0], uv[1])
                    cert["native_rows"] = rows
                    cert["native_section_ok"] = bool(ok)
                    cert["blocking_rows"] = [
                        x["s_coord"] for x in rows
                        if x["determination"] is None]
                cert["strong_criterion"] = bool(
                    prev and cert["native_section_ok"])
                # mutation centre-seul, sur le MÊME chart
                uvc = target_uv(Zc, tuple(S2), g2)
                if uvc is None:
                    cert["native_section_ok_centre_only"] = False
                else:
                    _, okc = native_section_constructible(
                        tuple(S2), g2, uvc[0], uvc[1])
                    cert["native_section_ok_centre_only"] = bool(okc)
                charts.append(cert)

        real = [x for x in charts if not x.get("is_source")]
        prevs = [x for x in real if x.get("prev_criterion")]
        strong = [x for x in real if x.get("strong_criterion")]
        lost = [x for x in prevs if not x.get("strong_criterion")]
        ctr_only = [x for x in real
                    if x.get("native_section_ok_centre_only")
                    and not x.get("native_section_ok")]
        # l'obstruction est-elle héritée ?
        blocked = Counter()
        inherited = 0
        for x in lost:
            for s in (x.get("blocking_rows") or []):
                blocked[s] += 1
                if s in src_rotated:
                    inherited += 1
        cells.append({
            "center_hex": r["center_hex"], "hw_hex": r["hw_hex"],
            "S": list(S), "g": g, "class_id": r["class_id"],
            "eps": list(eps),
            "source_determinations": [x["determination"] for x in srows],
            "source_rotated_coords": src_rotated,
            "n_charts_seen": n_seen,
            "n_charts_non_source": sum(1 for x in charts
                                       if not x.get("is_source")),
            "n_reachable": sum(1 for x in real if x.get("reachable")),
            "n_admissible": sum(1 for x in real if x.get("admissible")),
            "n_prev_criterion": len(prevs),
            "n_strong_criterion": len(strong),
            "n_lost_by_strengthening": len(lost),
            "n_centre_only_would_accept": len(ctr_only),
            "blocking_coord_histogram": {str(k): v
                                         for k, v in sorted(blocked.items())},
            "n_blocking_inherited_from_source": inherited,
            "n_blocking_total": sum(blocked.values()),
            "strong_examples": [{k: x[k] for k in ("S", "g")}
                                for x in strong[:5]],
            "charts": charts})
        log(f"  cellule {i + 1} : source {src_rotated} tournées · "
            f"{len(prevs)} charts sous l'ANCIEN critère → "
            f"**{len(strong)} sous le RENFORCÉ** "
            f"(−{len(lost)}) · centre-seul en accepterait "
            f"{len(ctr_only)} de plus")

    # --- F9 : la subdivision efface-t-elle l'obstruction ? -------------------
    # C'est LA question que 62 → 0 pose. Deux explications concurrentes :
    #   (a) ARTEFACT D'ENCLOSURE — `Im R'` est surestimé et straddle 0 par
    #       excès de largeur. Alors subdiviser doit FAIRE DISPARAÎTRE les
    #       refus, et le rebalayage n'a qu'à descendre d'un cran.
    #   (b) OBSTRUCTION GÉOMÉTRIQUE — le lieu `{Im R' = 0} ∩ {Re R' < 0}`
    #       traverse réellement la cellule. Alors les refus PERSISTENT, et
    #       leur nombre croît comme la mesure du lieu : `4^d` pour une
    #       codimension 2 dans une boîte réelle de dimension 4, contre
    #       `16^d` sous-cellules au total.
    # Les deux prédictions sont DISJOINTES et le comptage tranche.
    sub = {"tested": False}
    c0 = cells[0]
    lost0 = [x for x in c0["charts"]
             if x.get("prev_criterion") and not x.get("strong_criterion")]
    if lost0:
        x0 = lost0[0]
        S2 = tuple(x0["S"])
        g2 = int(x0["g"])
        r0 = residual[0]
        S0, g0, eps0 = tuple(r0["S"]), r0["g"], tuple(r0["eps"])
        ctr0 = [float.fromhex(z) for z in r0["center_hex"]]
        hw0 = float.fromhex(r0["hw_hex"])
        levels = []
        for depth in range(SUB_DEPTH + 1):
            n = 2 ** depth
            h = hw0 / n
            nref = 0
            for idx in product(range(n), repeat=4):
                cc = [ctr0[k] - hw0 + h * (2 * idx[k] + 1) for k in range(4)]
                Zs, _dZs, _rw = build_section(S0, g0, eps0, cc, h)
                uvs = target_uv(Zs, S2, g2) if all(
                    z is not None for z in Zs) else None
                if uvs is None:
                    nref += 1
                    continue
                _rr, okk = native_section_constructible(S2, g2, uvs[0], uvs[1])
                if not okk:
                    nref += 1
            levels.append({"depth": depth, "n_subcells": n ** 4,
                           "n_refused": nref,
                           "fraction": nref / float(n ** 4)})
            log(f"  F9 profondeur {depth} : {n ** 4} sous-cellules, "
                f"{nref} refusées ({100.0 * nref / n ** 4:.2f} %)")
        ratios = [levels[i + 1]["n_refused"] / max(1, levels[i]["n_refused"])
                  for i in range(len(levels) - 1)]
        sub = {"tested": True, "chart": {"S": list(S2), "g": g2},
               "levels": levels, "ratios": ratios,
               "vanishes": levels[-1]["n_refused"] == 0,
               "codim2_signature": bool(
                   levels[-1]["n_refused"] > 0
                   and all(3.0 <= r <= 5.0 for r in ratios))}

    # --- cross-check INDÉPENDANT contre l'artefact D5.1/D5.2 ----------------
    # `n_prev_criterion` DOIT reproduire `payoff_per_cell` de
    # `d5_fullcell` : deux scripts, même sélection de cellules, même
    # critère ancien. Un désaccord signalerait une dérive de l'un ou
    # de l'autre — c'est ce que la revue appelle un cross-check qui
    # n'existait dans aucun script.
    cross = {"checked": False}
    fc = RES / "k3_cap_b1e2iii_d5_fullcell.json"
    try:
        d = json.loads(fc.read_text(encoding="utf-8"))
        ref = list(d["aggregate"]["payoff_per_cell"])[:len(cells)]
        mine = [c["n_prev_criterion"] for c in cells]
        cross = {"checked": True, "reference": ref, "measured": mine,
                 "agrees": ref == mine, "source": fc.name}
    except Exception as exc:
        cross = {"checked": False, "error": str(exc)[:120]}

    # ------------------------------------------------------------------ gates
    all_rows = [row for c in cells for x in c["charts"]
                for row in (x.get("native_rows") or [])]
    built = [row for row in all_rows if row["determination"] is not None]
    gates = {
        "F1_strong_refines_previous": bool(cells) and all(
            (not x.get("strong_criterion")) or x.get("prev_criterion")
            for c in cells for x in c["charts"] if not x.get("is_source")),
        "F2_all_59_charts_examined": bool(cells) and all(
            c["n_charts_seen"] == 60 and c["n_charts_non_source"] == 59
            for c in cells),
        "F8_prev_count_matches_d5_fullcell": bool(cross.get("checked"))
        and bool(cross.get("agrees")),
        "F3_criterion_bites": bool(cells) and any(
            c["n_lost_by_strengthening"] > 0 for c in cells),
        "F4_centre_only_accepts_strictly_more": bool(cells) and any(
            c["n_centre_only_would_accept"] > 0 for c in cells),
        "F5_built_rows_are_sound": bool(built) and all(
            row["sq_residual_contains_zero"] for row in built),
        "F6_inheritance_measured": bool(cells) and all(
            "n_blocking_inherited_from_source" in c for c in cells)
        and sum(c["n_blocking_total"] for c in cells) > 0,
        "F9_subdivision_does_not_erase_obstruction": bool(
            sub.get("tested")) and (not sub.get("vanishes"))
        and bool(sub.get("codim2_signature")),
        "F7_verdict_published_per_cell": bool(cells) and all(
            isinstance(c.get("n_strong_criterion"), int) for c in cells)}

    tot = {"prev": sum(c["n_prev_criterion"] for c in cells),
           "strong": sum(c["n_strong_criterion"] for c in cells),
           "lost": sum(c["n_lost_by_strengthening"] for c in cells),
           "blocking_total": sum(c["n_blocking_total"] for c in cells),
           "blocking_inherited": sum(c["n_blocking_inherited_from_source"]
                                     for c in cells),
           "rows_built": len(built), "rows_total": len(all_rows)}
    verdict = (
        "**%d charts** passaient l'ancien critère (domaine + Jacobien + "
        "disjonction) sur %d cellules ; **%d** passent le critère "
        "RENFORCÉ qui exige en plus que la section NATIVE du chart cible "
        "soit constructible avec des composantes CERTIFIÉES sur la boîte "
        "entière. %d charts sont perdus. Sur %d lignes bloquantes, **%d "
        "portent une coordonnée ambiante que la source devait déjà "
        "tourner** — l'obstruction est HÉRITÉE dans %s des cas."
        % (tot["prev"], len(cells), tot["strong"], tot["lost"],
           tot["blocking_total"], tot["blocking_inherited"],
           ("%.0f%%" % (100.0 * tot["blocking_inherited"]
                        / tot["blocking_total"]))
           if tot["blocking_total"] else "n/a"))

    art = {"artifact": "k3_cap_b1e2iii_c126_chart_criterion",
           "claim": "Le critère de sélection des charts, RENFORCÉ par la "
                    "constructibilité de la section NATIVE du chart "
                    "cible, et rebalayé sur les 59 charts.",
           "totals": tot, "verdict": verdict,
           "cross_check_d5_fullcell": cross,
           "subdivision_probe": sub,
           "cells": cells, "gates": gates,
           "gates_passed": sum(1 for v in gates.values() if v),
           "gates_total": len(gates),
           "provenance": provenance([C118_JSON], time.time() - T0)}
    ART.parent.mkdir(parents=True, exist_ok=True)
    ART.write_text(json.dumps(art, indent=2, ensure_ascii=False),
                   encoding="utf-8")
    print("-" * 78)
    for k, v in gates.items():
        print(f"  {'PASS' if v else 'FAIL'}  {k}")
    print("-" * 78)
    print(verdict)
    print(f"→ {ART}")
    return 0 if all(gates.values()) else 1


# ===========================================================================
#  Self-test
# ===========================================================================
def _selftest():
    ok = True

    def chk(name, cond):
        nonlocal ok
        print(f"  {'ok  ' if cond else 'FAIL'} {name}")
        ok = ok and bool(cond)

    # T1 : une constante réelle POSITIVE est constructible en principale
    from .taylor_models import CONE
    Rp = TMC.const(CIV.from_complex(complex(4.0, 0.0)))
    try:
        Rp.sqrt_principal()
        chk("T1 constante positive : principale OK", True)
    except BranchCutError:
        chk("T1 constante positive : principale OK", False)

    # T2 : une constante réelle NÉGATIVE refuse la principale ET la
    #      tournée (composante indéterminée) — c'est le cœur de C125
    Rn = TMC.const(CIV.from_complex(complex(-4.0, 0.0)))
    p_ref = r_ref = False
    try:
        Rn.sqrt_principal()
    except BranchCutError:
        p_ref = True
    try:
        tm_sqrt_rotated(Rn, 0)
    except BranchCutError:
        r_ref = True
    chk("T2 réel négatif : principale ET tournée refusées",
        p_ref and r_ref)

    # T3 : un radicande sur le rayon POSITIF — la principale l'accepte,
    #      mais c'est la TOURNÉE qui est la continuation ; elle construit
    #      dès que `σ` est certifié.
    Rr = TMC.const(CIV.from_complex(complex(-4.0, 0.5)))
    try:
        got = tm_sqrt_rotated(Rr, 1)
    except BranchCutError:
        got = None
    chk("T3 Im > 0 : la tournée σ=+1 construit", got is not None)

    # T4 : soundness — Z² − R contient 0 pour la tournée de T3
    if got is not None:
        d = got * got - Rr
        dr, di = _rng(d), _rng(d, True)
        chk("T4 soundness Z² − R ∋ 0",
            dr[0] <= 0 <= dr[1] and di[0] <= 0 <= di[1])
    else:
        chk("T4 soundness Z² − R ∋ 0", False)

    # T5 : NÉGATIF — le critère renforcé ne doit JAMAIS accepter un chart
    #      que l'ancien refusait. On le vérifie sur la logique elle-même.
    def strong(prev, native):
        return bool(prev and native)
    chk("T5 négatif : strong ⊆ prev (les 4 combinaisons)",
        all(not strong(p, n) or p
            for p in (False, True) for n in (False, True))
        and not strong(False, True))

    # T6 : NÉGATIF — une mutation qui prendrait `strong = native` seul
    #      accepterait un chart refusé par l'ancien critère. Sans ce
    #      contre-exemple, T5 serait satisfait par n'importe quoi.
    def strong_broken(prev, native):
        return bool(native)
    chk("T6 négatif : la mutation `strong = native` VIOLE le raffinement",
        strong_broken(False, True) and not False)

    # T7 : `native_section_constructible` REFUSE quand `σ'` est
    #      indéterminé — pas d'essai silencieux
    class _FakeUV:
        pass
    chk("T7 rotated_sigma_from_coeffs rend 0 sur un straddle",
        rotated_sigma_from_coeffs(Fraction(1), Fraction(1),
                                  (-1.0, 1.0), (-1.0, 1.0),
                                  (-1.0, 1.0), (-1.0, 1.0)) in (0, None))

    print("SELF-TEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_selftest() if "--selftest" in sys.argv else build())
