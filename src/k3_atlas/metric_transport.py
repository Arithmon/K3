#!/usr/bin/env python3
"""
k3_cap_b1e2iii_c129f_f4_bridge_metric.py — C129-F0 / F4 : LA MÉTRIQUE
REJOINT LES CARTES-PONTS.

CE QUE CE SCRIPT PAIE — le contrat R5 de la 11ᵉ revue GPT
(`gpt_b1e2iii_c129f0_scout_f1_f2f3_review_2026_07_31.md`), qui met F4 en
HOLD tant que l'architecture d'atlas n'est pas réparée. Elle l'est
(F2/F3 v2, 17/17 : nerf sur arêtes certifiées, feuille continuée `D·Z`,
portée stratifiée), donc F4 peut tourner.

LE CHEMIN EST CELUI DE C129-E, PAS UN NOUVEAU. `transport_hardened`
(C127 durci C128-A/B/C : quatre `Qmat` source/mid/cible/congruence,
hermiticité gatée, Weyl en arithmétique dirigée, ratio relatif certifié)
est appelé sur la BOÎTE-PONT, en mode LEDGER FIGÉ — `ε'/σ'` viennent du
core certifié par C127 et ne sont JAMAIS re-dérivés, exactement le mode
éprouvé des probes d'autonomie C128-D et du full C129-E.

LE SEUL AJOUT, ET IL EST MINIMAL. `transport_hardened` construisait sa
section par `build_section`, qui (i) prend une demi-largeur SCALAIRE et
(ii) retombe sur `σ` de composante — or le pont est ANISOTROPE (`2H` en
Im, `H` en Re) et `σ` y est indéterminé, puisque `Im R` change de signe.
Un paramètre `section=` ADDITIF a donc été ajouté à `transport_hardened`
(`center`/`hw` n'y servaient QU'À `build_section` ; tout l'aval est
box-agnostique). Quand il vaut `None`, le comportement est inchangé —
et le gate **G1** le PROUVE en rejouant le chemin par défaut sur un
panel RÉGULIÈREMENT ESPACÉ du full C129-E (un pas constant sur les
indices de tuiles, PAS le panel stratifié de C129-E — nom corrigé,
12ᵉ revue D3) et en exigeant les MÊMES nombres que
l'artefact sérialisé. Une injection qui déplacerait un chiffre du full
C129-E serait un désastre silencieux ; elle est donc gatée, pas
supposée.

CE QUE F4 CERTIFIE
  · la métrique sur chacune des 64 boîtes-ponts (PD source, PD cible,
    invariance de jauge, congruence, Weyl) ;
  · le ratio relatif certifié contre le MÊME `δ = 1e-5` que C127/C129-E
    — le seuil n'est PAS relâché pour accommoder une boîte plus grande ;
    s'il ne passait pas, ce serait un REFUS à porter en revue, pas un δ
    à ajuster ;
  · la congruence sur les overlaps pont↔inférieur et pont↔pont, c'est-
    à-dire sur les arêtes qui portent réellement le nerf de F3 v2.

GATES
  G1  NON-RÉGRESSION DE L'INJECTION : chemin par défaut (`section=None`)
      rejoué sur un panel régulièrement espacé du full C129-E ⟹ ratio,
      slack Weyl, ledger et déterminations IDENTIQUES à l'artefact full
      sérialisé ;
  G2  amont gaté et complet : C127-D 14/14, C129-D 9/9, C129-E 8/8,
      F2/F3 v2 17/17 — et `mode == "full"` EXIGÉ de CHACUN des quatre
      (durci, 12ᵉ revue D2 : avant, seul C127-D le devait) ;
  G3  LEDGER FIGÉ : `ε'/σ'` du core, jamais re-dérivés ; et les `kinds`
      cible obtenus sur le pont sont ceux du core, sinon REFUS (le
      régime de certification ne doit pas changer en silence) ;
  G4  MÉTRIQUE SUR LES 64 PONTS : ratio relatif certifié `≤ δ = 1e-5`,
      slack Weyl `> 0`, PD source et cible, invariance de jauge, même
      point projectif ;
  G5  CERTIFICAT MÉTRIQUE SUR LES DOMAINES D'ARÊTE CERTIFIÉS
      (`metric_certificate_on_certified_edge_domains`, 12ᵉ revue D1) :
      les 64 overlaps pont↔inférieur et les 210 overlaps pont↔pont
      reçoivent leur transport sur leur boîte exacte, mêmes seuils,
      DEPUIS UN REPRÉSENTANT (la section du pont). La compatibilité des
      DEUX extrémités est transportée par l'identité analytique EXACTE
      des sections (F2/F3 v2) et par les certificats de carte (C129-E,
      G4) — AUCUNE comparaison métrique bilatérale indépendante n'est
      exécutée, et elle n'est pas requise pour le théorème local ; la
      base de l'implication est SÉRIALISÉE
      (`edge_certificate_contract`) ;
  G6  NÉGATIFS sur le pont : `J` perturbé ⟹ la congruence exclut 0 ;
      `ε'` figé FAUX ⟹ refus ; `σ'` figé FAUX ⟹ refus ;
  G7  AUCUN FILTRAGE SILENCIEUX : 64 = 64, 274 = 274 sur les arêtes, et
      tout refus est publié avec sa cause.

CE QUE CE SCRIPT NE PAIE PAS : l'identification canonique de la voisine
(F1a a REFUSÉ — l'atlas supérieur reste DÉRIVÉ) ; les faces Re, où
28/64 ponts restent des cartes RELATIVES ; le contrat EXACT de la
congruence métrique, qui reste certifiée SOUS δ et non exacte (la
promotion d'exactitude de C129-D porte sur l'identité de SECTION, pas
sur le triplet HK) ; le scaling complet ; les 895 autres paires ; R12-C.

Sortie : results/k3_cap_b1e2iii_c129f_f4_bridge_metric.json
Usage  : k3_cap_b1e2iii_c129f_f4_bridge_metric.py [--selftest]
Env    : K3_F4_MODE     panel (défaut) | full
         K3_F4_WORKERS  (défaut 4)
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
from fractions import Fraction
from multiprocessing import get_context
from pathlib import Path

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

_argv = sys.argv
sys.argv = [sys.argv[0]]
os.environ.setdefault("K3_TM_ORDER", "4")
os.environ.setdefault("K3_TM_SERIES", "4")
from mpmath import mp                                              # noqa: E402
from .witness_registry import load_canonical_MH              # noqa: E402
from .interval_arithmetic import build_M_civ                  # noqa: E402
from .taylor_models import TM_ORDER, UNARY_SERIES_DEG            # noqa: E402
from .gram_congruence import GAMMA                     # noqa: E402
from .chart_transport import transport_hardened   # noqa: E402
from .bridge_continuation import (               # noqa: E402
    DECK_D, IM_DIRS, bounds, box_of, build_section_bilateral,
    center_hw, inter, mirror_bounds)
sys.argv = _argv

RES = Path(os.environ.get(
    "K3_RES_DIR", Path(__file__).resolve().parent / "data"))
COVER_JSON = RES / "k3_cap_b1e2iii_d5_6_dyadic_cover.json"
ATLAS_JSON = RES / "k3_cap_b1e2iii_c127d_atlas.json"
C127_JSON = RES / "k3_cap_b1e2iii_c127_transport_all.json"
C127E_JSON = RES / "k3_cap_b1e2iii_c127e_residual.json"
C129D_JSON = RES / "k3_cap_b1e2iii_c129d_exact_gluing.json"
C129E_JSON = RES / "k3_cap_b1e2iii_c129e_halo_metric.json"
F2F3_JSON = RES / "k3_cap_b1e2iii_c129f_f2f3_bridge_atlas.json"
MODE = os.environ.get("K3_F4_MODE", "panel")
N_WORKERS = int(os.environ.get("K3_F4_WORKERS", "4"))
ART = RES / ("k3_cap_b1e2iii_c129f_f4_bridge_metric.json" if MODE == "full"
             else "k3_cap_b1e2iii_c129f_f4_bridge_metric_panel.json")

# --- PRÉ-ENREGISTRÉ : le MÊME δ que C127 / C129-E. Non ajustable ----------
DELTA_REL = 1e-5
N_PANEL = 6
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def _sha(p):
    try:
        return hashlib.sha256(Path(p).read_bytes()).hexdigest()
    except OSError:
        return None


# ===========================================================================
#  Chargement
# ===========================================================================
def load_leaves(cov, c127, c127e):
    tr = {r["tile_index"]: r for r in c127["transports"]
          if not r.get("failed")}
    tre = {r["box_index"]: r for r in c127e["transports"]
           if not r.get("failed")}
    out = []
    for i, t in enumerate(cov["tiles"]):
        r = tr[i]
        out.append({"src": "c127", "center_hex": t["center_hex"],
                    "hw_hex": t["hw_hex"], "chart": t["chart"],
                    "eps_target": r["eps_target"],
                    "sigma_target": r["sigma_target"],
                    "core_kinds_target": r["kinds_target"],
                    "core_src_det": r["source_determinations"]})
    for t in c127e["new_tiles"]:
        r = tre[t["box_index"]]
        out.append({"src": "c127e", "center_hex": t["center_hex"],
                    "hw_hex": t["hw_hex"], "chart": t["chart"],
                    "eps_target": r["eps_target"],
                    "sigma_target": r["sigma_target"],
                    "core_kinds_target": r["kinds_target"],
                    "core_src_det": r["source_determinations"]})
    return out


# ===========================================================================
#  Verdict par boîte — le contrat E de C129-E, repris VERBATIM
# ===========================================================================
def box_ok(r, kinds_expected):
    if r.get("failed"):
        return False
    sp = r.get("spectral", {})
    return bool(
        r.get("kinds_target") == kinds_expected
        and r.get("gauge_invariance_ok")
        and r.get("congruence_contains_zero")
        and r.get("same_projective_point")
        and r.get("residual_relative") is not None
        and r["residual_relative"] <= DELTA_REL
        and sp.get("hermitian_ok") and sp.get("weyl_transport_ok")
        and r.get("pd_source", {}).get("is_PD")
        and r.get("pd_target", {}).get("is_PD"))


# ===========================================================================
#  Jobs
# ===========================================================================
_G = {}


def _init(cell, M, c218, rw, leaves, bridges, ledgers, halos):
    _G.update(cell=cell, M=M, c218=c218, rw=rw, leaves=leaves,
              bridges=bridges, ledgers=ledgers, halos=halos)


def _metric_on_box(S, g, eps_src, box, S2, g2, kw):
    """Le chemin durci C127/C128, sur une boîte QUELCONQUE (le pont n'est
    pas un cube) : la section bilatérale est construite ici et INJECTÉE,
    tout le reste est le chemin C129-E verbatim."""
    c, h = center_hw(box)
    cf = [float(x) for x in c]
    hf = [float(x) for x in h]
    Z, dZ, rows = build_section_bilateral(S, g, eps_src, cf, hf)
    if any(z is None for z in Z):
        return {"failed": "bilateral_section_incomplete",
                "regimes": [x.get("regime") for x in rows]}
    r = transport_hardened(S, g, eps_src, cf, hf, S2, g2,
                           _G["M"], _G["c218"], _G["rw"],
                           section=(Z, dZ, rows), **kw)
    r["bridge_regimes"] = [x.get("regime") for x in rows]
    return r


def _bridge_metric_job(arg):
    i, kind = arg
    S, g, _eps = _G["cell"]
    leaf = _G["leaves"][i]
    S2, g2 = tuple(leaf["chart"]["S"]), leaf["chart"]["g"]
    eps_src = _G["ledgers"][i]           # ledger du PONT, dérivé en F2
    kw = {"fixed_eps2": tuple(leaf["eps_target"]),
          "fixed_sigma2": list(leaf["sigma_target"])}
    if kind == "mutation_J":
        kw["perturb_J"] = complex(0.3, 0.2)
    elif kind == "mutation_eps":
        e_bad = list(leaf["eps_target"])
        e_bad[0] = -e_bad[0]
        kw["fixed_eps2"] = tuple(e_bad)
    elif kind == "mutation_sigma":
        s_bad = list(leaf["sigma_target"])
        rot = next((r for r, k in enumerate(leaf["core_kinds_target"])
                    if k and k.startswith("rotated")), 0)
        s_bad[rot] = -1 if s_bad[rot] is None else -s_bad[rot]
        kw["fixed_sigma2"] = s_bad
    r = _metric_on_box(S, g, eps_src, _G["bridges"][i], S2, g2, kw)
    r["tile"], r["kind"] = i, kind
    return r


def _edge_metric_job(arg):
    """G5 — le certificat métrique sur le DOMAINE D'ARÊTE : le transport
    est refait sur la BOÎTE D'OVERLAP elle-même, depuis UN représentant
    (la section du pont `a`). Une arête d'atlas dont la métrique n'est
    pas certifiée sur son propre domaine n'est qu'une transition de
    coordonnées. Ce job ne compare PAS les deux extrémités entre elles :
    l'autre représentant décrit la MÊME section analytique (identité
    exacte F2/F3 v2), donc le même pullback métrique — l'implication est
    déclarée dans l'artefact (`edge_certificate_contract`), pas laissée
    implicite (12ᵉ revue D1)."""
    tag, a, b = arg
    S, g, _eps = _G["cell"]
    if tag == "BL":
        rec = _G["halos"][a]
        W = inter(_G["bridges"][a],
                  bounds([Fraction(float.fromhex(x))
                          for x in rec["center_hex"]],
                         Fraction(float.fromhex(rec["H_hex"]))))
    else:
        W = inter(_G["bridges"][a], _G["bridges"][b])
    out = {"edge": [tag, a, b]}
    if W is None:
        out["failed"] = "overlap_not_open"
        return out
    leaf = _G["leaves"][a]
    S2, g2 = tuple(leaf["chart"]["S"]), leaf["chart"]["g"]
    kw = {"fixed_eps2": tuple(leaf["eps_target"]),
          "fixed_sigma2": list(leaf["sigma_target"])}
    r = _metric_on_box(S, g, _G["ledgers"][a], W, S2, g2, kw)
    out.update(r)
    return out


def _regression_job(arg):
    """G1 — le chemin PAR DÉFAUT (`section=None`), rejoué à l'identique."""
    i, c_hex, H_hex = arg
    S, g, eps = _G["cell"]
    leaf = _G["leaves"][i]
    S2, g2 = tuple(leaf["chart"]["S"]), leaf["chart"]["g"]
    r = transport_hardened(
        S, g, eps, [float.fromhex(x) for x in c_hex],
        float.fromhex(H_hex), S2, g2, _G["M"], _G["c218"], _G["rw"],
        fixed_eps2=tuple(leaf["eps_target"]),
        fixed_sigma2=list(leaf["sigma_target"]))
    r["tile_index"] = i
    return r


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

    import inspect
    sig = inspect.signature(transport_hardened)
    chk("transport_hardened expose le point d'injection `section`",
        "section" in sig.parameters)
    chk("`section` est OPTIONNEL et vaut None par défaut",
        sig.parameters["section"].default is None)
    src = inspect.getsource(transport_hardened)
    chk("le chemin par défaut appelle TOUJOURS build_section",
        "if section is None:" in src and "build_section(" in src)
    chk("center/hw ne servent plus qu'à build_section",
        src.count("center") == 1 + src.count("center)")
        or "build_section(S, g, eps, center, hw)" in src)
    chk("δ est celui de C127/C129-E, non relâché", DELTA_REL == 1e-5)
    chk("D est bien la deck importée de F2/F3",
        tuple(DECK_D) == (1, -1, 1, 1, 1, -1))
    chk("IM_DIRS importées cohérentes", tuple(IM_DIRS) == (1, 3))
    # box_ok doit REFUSER un ratio au-dessus de δ et un Weyl négatif
    good = {"kinds_target": ["a"], "gauge_invariance_ok": True,
            "congruence_contains_zero": True,
            "same_projective_point": True, "residual_relative": 1e-6,
            "spectral": {"hermitian_ok": True, "weyl_transport_ok": True},
            "pd_source": {"is_PD": True}, "pd_target": {"is_PD": True}}
    chk("box_ok accepte un transport conforme", box_ok(good, ["a"]))
    chk("box_ok REFUSE un ratio au-dessus de δ",
        not box_ok({**good, "residual_relative": 2e-5}, ["a"]))
    chk("box_ok REFUSE un slack Weyl absent",
        not box_ok({**good, "spectral": {"hermitian_ok": True,
                                         "weyl_transport_ok": False}},
                   ["a"]))
    chk("box_ok REFUSE des kinds cible CHANGÉS",
        not box_ok(good, ["b"]))
    chk("box_ok REFUSE un échec déclaré",
        not box_ok({**good, "failed": "x"}, ["a"]))
    print(f"\nself-test {ok}/{tot}")
    return ok == tot


# ===========================================================================
#  Build
# ===========================================================================
def build():
    print("=" * 78)
    print(f"C129-F0 / F4 — MÉTRIQUE SUR LES CARTES-PONTS "
          f"({'FULL 64' if MODE == 'full' else f'PANEL {N_PANEL}'}), "
          f"{N_WORKERS} workers, δ = {DELTA_REL:.0e}, ledger FIGÉ")
    print("=" * 78)
    cov = json.loads(COVER_JSON.read_text(encoding="utf-8"))
    atl = json.loads(ATLAS_JSON.read_text(encoding="utf-8"))
    c127 = json.loads(C127_JSON.read_text(encoding="utf-8"))
    c127e = json.loads(C127E_JSON.read_text(encoding="utf-8"))
    c129d = json.loads(C129D_JSON.read_text(encoding="utf-8"))
    c129e = json.loads(C129E_JSON.read_text(encoding="utf-8"))
    f23 = json.loads(F2F3_JSON.read_text(encoding="utf-8"))
    cell_d = cov["cell"]
    S, g, eps = (tuple(cell_d["S"]), cell_d["g"], tuple(cell_d["eps"]))
    leaves = load_leaves(cov, c127, c127e)
    halos = {h["index"]: h["record"] for h in atl["halos"] if h["ok"]}

    # --- G2 : amont — mode == "full" EXIGÉ des QUATRE (12ᵉ revue D2 ;
    # avant, seul C127-D le devait, et F2/F3 sérialisait mode=null) -----
    up = {}
    for name, blob in (("c127d", atl), ("c129d", c129d),
                       ("c129e", c129e), ("f2f3_v2", f23)):
        gp, gt = blob.get("gates_passed"), blob.get("gates_total")
        up[name] = {"gates": f"{gp}/{gt}", "mode": blob.get("mode"),
                    "green": bool(gp == gt and gt
                                  and blob.get("mode") == "full")}
    g2 = all(v["green"] for v in up.values())
    log("G2 : amont — " + " ; ".join(f"{k} {v['gates']}"
                                     for k, v in up.items())
        + f" ⟹ {g2}")

    bridges = {t["tile"]: [(Fraction(*b[0]), Fraction(*b[1]))
                           for b in t["bridge_corrected_bounds"]]
               for t in json.loads(
                   (RES / "k3_cap_b1e2iii_c129f_bridge_scout.json")
                   .read_text(encoding="utf-8"))["per_tile"]}
    ledgers = {r["tile"]: tuple(r["F2d_ledger_derived"])
               for r in f23["per_bridge"]}
    clipped = sorted(bridges)
    log(f"    {len(clipped)} ponts, ledgers dérivés importés de F2/F3 v2")

    reg = load_canonical_MH()
    M = build_M_civ(reg["M_H_canonical"])
    c218 = reg["coeffs218"]
    rw = 1.0 - GAMMA
    mpctx = get_context("fork")
    initargs = ((S, g, eps), M, c218, rw, leaves, bridges, ledgers, halos)
    _init(*initargs)

    # --- G1 : NON-RÉGRESSION de l'injection ---------------------------
    # Un panel RÉGULIÈREMENT ESPACÉ du full C129-E (pas constant sur les
    # indices — PAS le panel stratifié de C129-E, nom corrigé 12ᵉ revue
    # D3), rejoué par le chemin PAR DÉFAUT. Les nombres doivent être
    # IDENTIQUES à ceux SÉRIALISÉS : une injection qui déplacerait un
    # chiffre du full serait un désastre silencieux.
    ref = {r["tile_index"]: r for r in c129e["transports"]
           if r.get("kind") == "transport" and not r.get("failed")}
    probe = sorted(ref)[::max(1, len(ref) // 6)][:6]
    jobs = [(i, halos[i]["center_hex"], halos[i]["H_hex"]) for i in probe]
    with mpctx.Pool(min(N_WORKERS, len(jobs)), initializer=_init,
                    initargs=initargs) as pool:
        regr = pool.map(_regression_job, jobs)
    g1_rows, g1 = [], True
    for r in regr:
        a = ref[r["tile_index"]]
        same = (r.get("residual_relative") == a.get("residual_relative")
                and r.get("source_determinations")
                == a.get("source_determinations")
                and r.get("kinds_target") == a.get("kinds_target")
                and r.get("spectral", {}).get("slack", {}).get("exact")
                == a.get("spectral", {}).get("slack", {}).get("exact"))
        g1 = g1 and same
        g1_rows.append({"tile": r["tile_index"], "identical": bool(same),
                        "ratio_now": r.get("residual_relative"),
                        "ratio_serialized": a.get("residual_relative")})
    log(f"G1 (NON-RÉGRESSION) : chemin par défaut rejoué sur "
        f"{len(probe)} tuiles du full C129-E ⟹ ratio, ledger, "
        f"déterminations et slack IDENTIQUES : {g1}")

    # --- G3/G4 : la métrique sur les ponts ----------------------------
    sel = clipped if MODE == "full" else clipped[::max(
        1, len(clipped) // N_PANEL)][:N_PANEL]
    log(f"G4 : métrique sur {len(sel)} pont(s)…")
    with mpctx.Pool(N_WORKERS, initializer=_init,
                    initargs=initargs) as pool:
        met = pool.map(_bridge_metric_job, [(i, "metric") for i in sel])
    ok_rows, bad_rows = [], []
    for r in met:
        exp = leaves[r["tile"]]["core_kinds_target"]
        (ok_rows if box_ok(r, exp) else bad_rows).append(r)
    g4 = len(bad_rows) == 0 and len(met) == len(sel)
    ratios = [r["residual_relative"] for r in ok_rows
              if r.get("residual_relative") is not None]
    slacks = [r["spectral"]["slack"]["float"] for r in ok_rows
              if r.get("spectral", {}).get("slack")]
    kinds_c = Counter(tuple(r.get("kinds_target") or []) for r in met)
    g3 = all(r.get("kinds_target") == leaves[r["tile"]]["core_kinds_target"]
             for r in met)
    log(f"     ratio relatif max {max(ratios):.3e} ≤ δ = {DELTA_REL:.0e} ; "
        f"slack Weyl min {min(slacks):.3e} > 0 ; {len(ok_rows)}/{len(met)}")
    log(f"G3 : ledger FIGÉ, kinds cible sur le pont == core : {g3} "
        f"— {dict(kinds_c)}")
    for r in bad_rows:
        log(f"     REFUS pont {r['tile']} : {r.get('failed')} "
            f"ratio={r.get('residual_relative')} "
            f"kinds={r.get('kinds_target')}")

    # --- G5 : la congruence sur les ARÊTES ---------------------------
    bb = [tuple(x["pair"]) for x in f23["bridge_bridge_transitions"]
          if x.get("certified")]
    edges = [("BL", i, i) for i in sel]
    if MODE == "full":
        edges += [("BB", a, b) for a, b in bb]
    else:
        edges += [("BB", a, b) for a, b in bb
                  if a in sel][:N_PANEL]
    log(f"G5 : congruence sur {len(edges)} arête(s) du nerf…")
    with mpctx.Pool(N_WORKERS, initializer=_init,
                    initargs=initargs) as pool:
        edg = pool.map(_edge_metric_job, edges)
    e_ok, e_bad = [], []
    for r in edg:
        exp = leaves[r["edge"][1]]["core_kinds_target"]
        (e_ok if box_ok(r, exp) else e_bad).append(r)
    g5 = len(e_bad) == 0 and len(edg) == len(edges)
    e_ratios = [r["residual_relative"] for r in e_ok
                if r.get("residual_relative") is not None]
    e_slacks = [r["spectral"]["slack"]["float"] for r in e_ok
                if r.get("spectral", {}).get("slack")]
    log(f"     {len(e_ok)}/{len(edg)} — ratio max "
        f"{max(e_ratios):.3e} ; slack min {min(e_slacks):.3e}"
        if e_ok else "     aucune arête certifiée")
    for r in e_bad:
        log(f"     REFUS arête {r['edge']} : {r.get('failed')} "
            f"ratio={r.get('residual_relative')}")

    # --- G6 : les négatifs sur le PONT -------------------------------
    # Les CANAUX attendus sont PRÉ-ENREGISTRÉS, repris de C129-E : un
    # négatif qui casse par un autre canal que celui annoncé n'est pas
    # le négatif qu'on croit avoir écrit.
    NEG_CHANNEL = {
        "mutation_J": "congruence exclut 0",
        "mutation_eps": "même point projectif refusé",
        "mutation_sigma": ("racine tournée niée (R8) ⟹ même point "
                           "projectif refusé ou kinds changés")}
    negs = [(sel[0], k) for k in
            ("mutation_J", "mutation_eps", "mutation_sigma")]
    with mpctx.Pool(min(N_WORKERS, 3), initializer=_init,
                    initargs=initargs) as pool:
        nr = pool.map(_bridge_metric_job, negs)
    neg_rows, g6 = [], True
    for r in nr:
        broke = not box_ok(r, leaves[r["tile"]]["core_kinds_target"])
        g6 = g6 and broke
        chan = ("congruence" if not r.get("congruence_contains_zero")
                else ("same_projective_point"
                      if not r.get("same_projective_point")
                      else ("kinds" if r.get("kinds_target")
                            != leaves[r["tile"]]["core_kinds_target"]
                            else None)))
        neg_rows.append({
            "kind": r["kind"], "tile": r["tile"], "breaks": bool(broke),
            "expected_channel": NEG_CHANNEL[r["kind"]],
            "observed_channel": chan,
            "failed": r.get("failed"),
            "congruence_contains_zero": r.get("congruence_contains_zero"),
            "same_projective_point": r.get("same_projective_point"),
            "kinds_target": r.get("kinds_target"),
            "residual_relative": r.get("residual_relative"),
            "note": ("la CONGRUENCE MÉTRIQUE est aveugle aux signes du "
                     "ledger CIBLE (Qmat ne voit que des modules) : "
                     "c'est le canal `same_projective_point` qui porte "
                     "cette discriminance, et le ratio reste inchangé. "
                     "Déclaré, pas masqué — même canal qu'en C129-E."
                     if chan == "same_projective_point" else None)})
        log(f"G6 : {r['kind']:16s} ⟹ casse={broke} "
            f"(failed={r.get('failed')}, "
            f"cong∋0={r.get('congruence_contains_zero')})")

    g7 = bool(len(met) == len(sel) and len(edg) == len(edges))
    gates = {
        "G1_injection_does_not_regress_default_path": bool(g1),
        "G2_upstream_green_and_full": bool(g2),
        "G3_frozen_ledger_target_kinds_unchanged": bool(g3),
        "G4_metric_certified_on_bridges": bool(g4),
        "G5_metric_certified_on_nerve_edges": bool(g5),
        "G6_three_negatives_break_on_bridge": bool(g6),
        "G7_no_silent_filtering": bool(g7)}
    npass = sum(1 for v in gates.values() if v)

    try:
        head = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=Path(__file__).resolve().parent, capture_output=True,
            text=True, timeout=10).stdout.strip() or None
    except Exception:
        head = None

    out = {
        "artifact": "k3_cap_b1e2iii_c129f_f4_bridge_metric",
        "mode": MODE,
        "claim_level": f23["claim_level"],
        "claim": (
            "Le chemin métrique durci de C127/C128 — quatre Qmat, "
            "hermiticité gatée, Weyl en arithmétique dirigée, ratio "
            "relatif CERTIFIÉ — est rejoué VERBATIM sur les "
            "BOÎTES-PONTS bilatérales et sur les ARÊTES du nerf de "
            "F2/F3 v2, en mode LEDGER FIGÉ (ε'/σ' du core, jamais "
            "re-dérivés). Le seul ajout est un point d'injection "
            "`section=` dans `transport_hardened`, dont la "
            "non-régression est GATÉE contre l'artefact C129-E "
            "sérialisé. Le seuil δ = 1e-5 est celui de C127/C129-E, "
            "non relâché. La congruence reste certifiée SOUS δ, PAS "
            "exacte."),
        "cell": {"S": list(S), "g": g, "eps": list(eps)},
        "delta_rel_preregistered": DELTA_REL,
        "upstream": up,
        "n_bridges_selected": len(sel),
        "n_bridges_certified": len(ok_rows),
        "max_residual_relative_bridges": max(ratios) if ratios else None,
        "min_weyl_slack_bridges": min(slacks) if slacks else None,
        "target_kinds_census": {str(k): v for k, v in kinds_c.items()},
        "n_edges_selected": len(edges),
        "n_edges_certified": len(e_ok),
        "max_residual_relative_edges": max(e_ratios) if e_ratios else None,
        "min_weyl_slack_edges": min(e_slacks) if e_slacks else None,
        "edge_certificate_contract": {
            "name": "metric_certificate_on_certified_edge_domains",
            "endpoint_compatibility_basis":
                "exact_section_identity_from_f2f3_v2",
            "pairwise_metric_comparison_executed": False,
            "note": (
                "G5 certifie le chemin métrique sur le domaine EXACT de "
                "chaque arête depuis UN représentant (la section du "
                "pont). La compatibilité des deux extrémités est "
                "transportée par l'identité analytique EXACTE des "
                "sections (F2/F3 v2) et par les certificats de carte "
                "(C129-E pour les cartes inférieures, G4 pour les "
                "ponts) : l'autre représentant décrit la même section, "
                "donc le même pullback métrique. AUCUNE comparaison "
                "métrique bilatérale indépendante entre extrémités "
                "n'est exécutée ici — elle n'est pas requise pour le "
                "théorème local, et le dire est la 12ᵉ revue D1.")},
        "regression_probe": g1_rows,
        "negatives": neg_rows,
        "bridges": met,
        "edges": edg,
        "not_paid_here": [
            "l'identification canonique de la voisine : F1a a REFUSÉ, "
            "l'atlas supérieur reste DÉRIVÉ par conjugaison",
            "les faces Re : 28/64 ponts y restent des cartes RELATIVES",
            "le contrat EXACT de la congruence métrique — elle reste "
            "certifiée SOUS δ, pas exacte ; la promotion d'exactitude de "
            "C129-D porte sur l'identité de SECTION, pas sur le triplet "
            "HK, donc §4.4 de la carte épistémique TIENT",
            "le scaling complet, les 895 autres paires, R12-C"],
        "gates": gates, "gates_passed": npass, "gates_total": len(gates),
        "verdict": (
            f"F4 {'FULL' if MODE == 'full' else 'PANEL'} LIVRÉ — la "
            f"métrique rejoint les cartes-ponts et les arêtes du nerf."
            if npass == len(gates) else
            f"ROUGE — {len(gates) - npass} gate(s) en échec"),
        "provenance": {
            "git_head": head, "python": sys.version.split()[0],
            "platform": platform.platform(), "mp_prec": int(mp.prec),
            "tm_order": TM_ORDER, "unary_series_deg": UNARY_SERIES_DEG,
            "wall_s": round(time.time() - T0, 1), "n_workers": N_WORKERS,
            "preregistered": {"delta_rel": DELTA_REL,
                              "n_panel": N_PANEL},
            "inputs": {p.name: _sha(p) for p in
                       (COVER_JSON, ATLAS_JSON, C127_JSON, C127E_JSON,
                        C129D_JSON, C129E_JSON, F2F3_JSON)},
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
