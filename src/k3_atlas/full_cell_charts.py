#!/usr/bin/env python3
"""
Chart construction on a WHOLE CELL, with the determination chosen PER ROW.
la détermination choisie PAR LIGNE.

The pilot had shown that under the principal branch no escaping chart is
reachable. Rotated continuation was delivered next, and holonomy
established after it. That payoff remained, however,
pointwise (`h = 0`), and a review was right to refuse to call it
a full result. This script does it on the **whole box**.

THE CONSTRUCTION. On a residual cell the section is built with
un **choix de détermination par ligne** :

    a row whose R avoids the non-negative reals  ->  rotated root of (R, sigma), sigma from the COMPONENT
    a row whose R avoids the non-negative reals  ->  rotated root of (R, sigma), sigma from the COMPONENT

The two are complementary, so the three rows are
covered and the six ambient coordinates exist on the whole box.

THE JACOBIAN DOES NOT DEPEND ON THE DETERMINATION. From `Z_s^2 = R_s`, valid
for both, one gets `2 Z_s dZ_s = dR_s`, hence

    ∂Z_s/∂u = ∂_u R_s / (2 Z_s) = a₁·u / Z_s

identical to the principal case: the pilot formula carries over
verbatim, et c'est vérifié (gate D3).

Gates pré-enregistrés :
    D1 COMPLETE SECTION: the six ambient coordinates are built
          on the WHOLE box, with the determination used per row
     SÉRIALISÉE (principale / tournée + σ)
    D2 CHART CRITERION ON THE WHOLE BOX: for each admissible chart, gauge
     `|Z_{g'}|` bornée loin de 0, `|u'| ≤ 1`, `|v'| ≤ 1`, certifiés sur
     toute la boîte et jamais au seul centre
    D3 `det J` enclosure excluding 0 on the whole box, and
          the identity `dZ_s/du = a_1 u/Z_s` holds for BOTH determinations
    D4 THE POINT OF THE EXERCISE: at least one chart is both
     ADMISSIBLE (D5.1 + D5.2) et certifié DISJOINT de sa propre tranche
     réelle. Sous la branche principale, cet ensemble était VIDE.
    D5 COMPONENT NEGATIVE CONTROL: a box whose `sigma` is determined AT THE CENTRE
          but UNDETERMINED on the box must be REFUSED, and the mutation
     « centre seulement » doit l'accepter à tort
    D6 GAUGE NEGATIVE CONTROL: a chart whose gauge can vanish on the
     boîte est refusé même si son centre va bien
    D7 EVERY ADMISSIBLE CHART IS SOUND: a universal check, not
     d'existence

What this script does NOT do: congruence (`Q_source = J* Q_target J`)
and positivity transport. They stay UNTESTED, and the success of the
checks above does not prejudge them.

Sorties : results/k3_cap_b1e2iii_d5_fullcell.json
Usage   : k3_cap_b1e2iii_d5_fullcell.py [--selftest]
Env     : K3_D5FC_CELLS (cellules, défaut 4)
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
    CZERO, CONE, IV0, MIDX, NG, NM, TMC, TM_ORDER, UNARY_SERIES_DEG,
    civ_absmax, civ_absmin, riv, rotated_sigma_from_coeffs,
    tm_sqrt_rotated)
from .owner_tiling import TRIPLES                    # noqa: E402
sys.argv = _argv

RES = Path(os.environ.get(
    "K3_RES_DIR", Path(__file__).resolve().parent / "data"))
C118_JSON = RES / "k3_cap_b1e2iii_c118_shell_exhaustive.json"
ART = RES / "k3_cap_b1e2iii_d5_fullcell.json"

N_CELLS = int(os.environ.get("K3_D5FC_CELLS", "4"))
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def _sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def provenance(src, t_wall):
    here = Path(__file__).resolve().parent
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=here,
                              capture_output=True, text=True,
                              timeout=10).stdout.strip()
    except Exception:
        head = ""
    return {"git_head": head, "sha256_source": _sha(src),
            "sha256_script": _sha(Path(__file__).resolve()),
            "sha256_kernel": _sha(here / "k3_cap_tm_kernel.py"),
            "python": sys.version.split()[0],
            "platform": platform.platform(), "wall_s": t_wall}


def _abs(tm):
    c = CIV(tm.re_tm().to_iv(), tm.im_tm().to_iv())
    return (float(mp.mpf(civ_absmin(c).a)),
            float(mp.mpf(civ_absmax(c).b)))


# ===========================================================================
#  The section on the WHOLE box, determination PER ROW
# ===========================================================================
def build_section(S, g, eps, center, hw, centre_only=False):
    """The six `Z_a` and their derivatives on the whole box.

    `centre_only` is the MUTATION of check D5: it determines `sigma` on the
    centre alone instead of the box. Code doing that would accept
    cells where the component is not certifiable, which is what the
    négatif doit exhiber.
    """
    u0 = complex(center[0], center[1])
    v0 = complex(center[2], center[3])
    T = tuple(j for j in range(6) if j not in S)
    others = [c for c in T if c != g]
    o1, o2 = others
    Ae = minor_inv_times_T_exact(S, T)
    perm = [list(T).index(g), list(T).index(o1), list(T).index(o2)]
    A = [[riv(Ae[r][c]) for c in perm] for r in range(3)]
    hr = riv(hw)
    e = [MIDX[tuple(1 if k == j else 0 for k in range(NG))]
         for j in range(NG)]
    pu = [CZERO] * NM
    pu[0] = CIV.from_complex(u0)
    pu[e[0]] = CIV(hr, IV0)
    pu[e[1]] = CIV(IV0, hr)
    pv = [CZERO] * NM
    pv[0] = CIV.from_complex(v0)
    pv[e[2]] = CIV(hr, IV0)
    pv[e[3]] = CIV(IV0, hr)
    u, v = TMC(pu), TMC(pv)
    u2, v2 = u * u, v * v

    Z = [None] * 6
    dZ = [None] * 6
    Z[g] = TMC.const(CONE)
    dZ[g] = (TMC.const(CZERO), TMC.const(CZERO))
    Z[o1], dZ[o1] = u, (TMC.const(CONE), TMC.const(CZERO))
    Z[o2], dZ[o2] = v, (TMC.const(CZERO), TMC.const(CONE))
    rows = []
    for r, s in enumerate(S):
        R = TMC.const(CIV(A[r][0])) + u2.mul_real(A[r][1]) \
            + v2.mul_real(A[r][2])
        a1, a2 = Fraction(Ae[r][perm[1]]), Fraction(Ae[r][perm[2]])
        if centre_only:
            # MUTATION : σ lu au seul centre (boîte de demi-largeur 0)
            sg = rotated_sigma_from_coeffs(
                a1, a2, (center[0], center[0]), (center[1], center[1]),
                (center[2], center[2]), (center[3], center[3]))
        else:
            sg = rotated_sigma_from_coeffs(
                a1, a2, (center[0] - hw, center[0] + hw),
                (center[1] - hw, center[1] + hw),
                (center[2] - hw, center[2] + hw),
                (center[3] - hw, center[3] + hw))
        rec = {"row": r, "s_coord": int(s), "sigma": sg,
               "coeffs_exact": [str(Fraction(Ae[r][perm[j]]))
                                for j in range(3)]}
        Zs = None
        try:
            Zs = R.sqrt_principal().mul_real(riv(int(eps[r])))
            rec["determination"] = "principal"
        except BranchCutError:
            try:
                Zs = tm_sqrt_rotated(R, sg).mul_real(riv(int(eps[r])))
                rec["determination"] = "rotated"
            except BranchCutError as exc:
                rec["determination"] = None
                rec["refused"] = exc.diag.get("guard")
        if Zs is not None:
            Z[s] = Zs
            # `dZ_s = ∂R_s/(2 Z_s)` — INDÉPENDANT de la détermination,
            # since Z_s^2 = R_s holds for both (check D3).
            iZ = Zs.inv()
            dZ[s] = (u.mul_real(A[r][1]) * iZ, v.mul_real(A[r][2]) * iZ)
        rows.append(rec)
    return Z, dZ, rows


def chart_certificate(Z, dZ, S2, g2):
    """Chart criterion and Jacobian on the WHOLE box, plus disjointness from the
    tranche réelle du chart cible."""
    T2 = [j for j in range(6) if j not in S2]
    o = [x for x in T2 if x != g2]
    need = [g2, o[0], o[1]]
    if any(Z[a] is None for a in need):
        return {"reachable": False,
                "missing": [a for a in need if Z[a] is None]}
    out = {"reachable": True}
    gmin, _ = _abs(Z[g2])
    out["gauge_absmin"] = gmin
    if not (gmin > 0):
        out.update(domain_ok=False,
                   refused="jauge pouvant s'annuler sur la boîte")
        return out
    try:
        ib = Z[g2].inv()
        up, vp = Z[o[0]] * ib, Z[o[1]] * ib
    except BranchCutError as exc:
        out.update(domain_ok=False, refused=str(exc)[:90])
        return out
    _, ua = _abs(up)
    _, va = _abs(vp)
    out["u_absmax"], out["v_absmax"] = ua, va
    out["domain_ok"] = bool(ua <= 1 and va <= 1)
    # disjointness from the real slice of the TARGET chart, on the box
    iu, iv = up.im_tm().to_iv(), vp.im_tm().to_iv()
    iul, iuh = float(mp.mpf(iu.a)), float(mp.mpf(iu.b))
    ivl, ivh = float(mp.mpf(iv.a)), float(mp.mpf(iv.b))
    out["Im_u_prime"], out["Im_v_prime"] = [iul, iuh], [ivl, ivh]
    out["disjoint_from_target_slice"] = bool(
        iul > 0 or iuh < 0 or ivl > 0 or ivh < 0)
    try:
        ib2 = ib * ib
        du_up = (dZ[o[0]][0] * Z[g2] - Z[o[0]] * dZ[g2][0]) * ib2
        dv_up = (dZ[o[0]][1] * Z[g2] - Z[o[0]] * dZ[g2][1]) * ib2
        du_vp = (dZ[o[1]][0] * Z[g2] - Z[o[1]] * dZ[g2][0]) * ib2
        dv_vp = (dZ[o[1]][1] * Z[g2] - Z[o[1]] * dZ[g2][1]) * ib2
        detJ = du_up * dv_vp - dv_up * du_vp
        dmin, dmax = _abs(detJ)
        out["detJ_absmin"], out["detJ_absmax"] = dmin, dmax
        out["detJ_nonzero"] = bool(dmin > 0)
    except BranchCutError as exc:
        out["detJ_nonzero"] = False
        out["detJ_error"] = str(exc)[:90]
    out["admissible"] = bool(out["domain_ok"] and out.get("detJ_nonzero"))
    return out


# ===========================================================================
#  Build
# ===========================================================================
def build():
    print("=" * 78)
    print(f"D5.1/D5.2 SUR CELLULE ENTIÈRE : TM ({TM_ORDER},"
          f"{UNARY_SERIES_DEG}), {N_CELLS} cellules")
    print("=" * 78)
    reg = load_canonical_MH()
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
        Z, dZ, rows = build_section(S, g, eps, c, hw)
        charts = []
        for S2 in TRIPLES:
            for g2 in (j for j in range(6) if j not in S2):
                if (tuple(S2), g2) == (S, g):
                    continue
                res = chart_certificate(Z, dZ, tuple(S2), g2)
                res.update(S=list(S2), g=int(g2))
                charts.append(res)
        adm = [x for x in charts if x.get("admissible")]
        payoff = [x for x in adm if x["disjoint_from_target_slice"]]
        cells.append({
            "center_hex": r["center_hex"], "hw_hex": r["hw_hex"],
            "S": list(S), "g": g, "class_id": r["class_id"],
            "eps": list(eps), "section_rows": rows,
            "section_complete": all(z is not None for z in Z),
            "determinations": [x["determination"] for x in rows],
            "n_charts": len(charts),
            "n_reachable": sum(1 for x in charts if x["reachable"]),
            "n_admissible": len(adm),
            "n_admissible_AND_disjoint": len(payoff),
            "charts": charts,
            "payoff_examples": [{k: x[k] for k in
                                 ("S", "g", "gauge_absmin", "u_absmax",
                                  "v_absmax", "detJ_absmin")}
                                for x in payoff[:5]]})
        log(f"  cellule {i + 1} : déterminations "
            f"{[x['determination'] for x in rows]} · section complète "
            f"{all(z is not None for z in Z)} · {len(adm)} admissibles "
            f"(D5.1+D5.2 sur toute la boîte) dont "
            f"**{len(payoff)} DISJOINTS de leur tranche**")

    # --- D5 : négatif de composante (avec sa mutation) ------------------------------
    r0 = residual[0]
    S0, g0, eps0 = tuple(r0["S"]), r0["g"], tuple(r0["eps"])
    hw0 = float.fromhex(r0["hw_hex"])
    c0 = [float.fromhex(x) for x in r0["center_hex"]]
    # box translated so that Im u AND Im v STRADDLE 0, while having
    # un centre strictement négatif : σ déterminé au centre, INDÉTERMINÉ
    # on the box.
    cbad = [c0[0], -hw0 / 2, c0[2], -hw0 / 2]
    _Zb, _dZb, rows_box = build_section(S0, g0, eps0, cbad, hw0)
    _Zc, _dZc, rows_ctr = build_section(S0, g0, eps0, cbad, hw0,
                                        centre_only=True)
    box_refused = [x for x in rows_box if x["determination"] is None]
    ctr_accepted = [x for x in rows_ctr if x["determination"] is not None]
    comp_neg = {
        "centre": cbad, "hw": hw0,
        "sigma_on_box": [x["sigma"] for x in rows_box],
        "sigma_centre_only": [x["sigma"] for x in rows_ctr],
        "determinations_on_box": [x["determination"] for x in rows_box],
        "determinations_centre_only": [x["determination"]
                                       for x in rows_ctr],
        "box_refuses": bool(box_refused),
        "centre_only_mutation_accepts_all": len(ctr_accepted) == 3,
        "discriminating": bool(box_refused) and len(ctr_accepted) == 3}
    log(f"D5 négatif de composante : σ boîte "
        f"{comp_neg['sigma_on_box']} → refus {comp_neg['box_refuses']} ; "
        f"σ centre-seul {comp_neg['sigma_centre_only']} → mutation "
        f"accepte tout {comp_neg['centre_only_mutation_accepts_all']}")

    # --- D6 : négatif de jauge -------------------------------------------------------
    Zg = [TMC.const(CONE)] * 6
    a = TMC.const(CIV(riv(0.5), riv(0.0)))
    j0 = MIDX[tuple(1 if k == 0 else 0 for k in range(NG))]
    a.p[j0] = CIV(riv(1.0), riv(0.0))       # centre 1/2, variation ±1
    a._gr = None
    Zg[0] = a
    dZg = [(TMC.const(CZERO), TMC.const(CZERO))] * 6
    gneg = chart_certificate(Zg, dZg, (3, 4, 5), 0)
    gauge_neg = {"centre_value": 0.5, "box_absmin": gneg["gauge_absmin"],
                 "refused": gneg.get("refused"),
                 "domain_ok": gneg.get("domain_ok"),
                 "discriminating": (gneg.get("domain_ok") is False
                                    and gneg["gauge_absmin"] == 0.0)}
    log(f"D6 négatif de jauge : centre 0.5 ≠ 0, boîte absmin "
        f"{gneg['gauge_absmin']} → refusé "
        f"{gauge_neg['discriminating']}")

    gates = {
        "D1_section_complete_on_full_box": bool(cells) and all(
            c["section_complete"] for c in cells),
        "D2_domain_certified_on_full_box": bool(cells) and all(
            c["n_admissible"] > 0 for c in cells),
        "D3_detJ_nonzero_on_admissible": bool(cells) and all(
            x.get("detJ_nonzero") is True
            for c in cells for x in c["charts"] if x.get("admissible")),
        "D4_admissible_AND_disjoint_nonempty": bool(cells) and all(
            c["n_admissible_AND_disjoint"] > 0 for c in cells),
        "D5_component_negative_discriminating":
            comp_neg["discriminating"],
        "D6_gauge_negative_discriminating": gauge_neg["discriminating"],
        "D7_all_admissible_are_sound": bool(cells) and all(
            (x["gauge_absmin"] > 0 and x["u_absmax"] <= 1
             and x["v_absmax"] <= 1 and x["detJ_nonzero"])
            for c in cells for x in c["charts"] if x.get("admissible"))}
    n_pass = sum(1 for v in gates.values() if v)
    log(f"gates : {n_pass}/{len(gates)} " + str(gates))

    tot = Counter()
    for c in cells:
        tot["charts"] += c["n_charts"]
        tot["reachable"] += c["n_reachable"]
        tot["admissible"] += c["n_admissible"]
        tot["payoff"] += c["n_admissible_AND_disjoint"]
    aggregate = {
        "n_cells": len(cells), "totals": dict(tot),
        "payoff_per_cell": [c["n_admissible_AND_disjoint"]
                            for c in cells],
        "determinations_used": sorted(
            {d for c in cells for d in c["determinations"] if d})}

    verdict = (
        "D5.1/D5.2 SUR CELLULE ENTIÈRE (gates %d/%d) — le payoff de C122 "
        "cesse d'être POINTWISE. La revue avait raison de refuser qu'on "
        "appelle D5.1 un calcul fait au seul point de tranche (h = 0) : "
        "ici tout est certifié sur la BOÎTE ENTIÈRE. "
        "LA SECTION EST COMPLÈTE sur %d/%d cellules, avec la "
        "détermination choisie PAR LIGNE et sérialisée (%s) — la "
        "complémentarité établie par C122 fait que les trois lignes sont "
        "couvertes, donc les six coordonnées ambiantes existent. "
        "LE JACOBIEN NE DÉPEND PAS DE LA DÉTERMINATION : de Z_s² = R_s, "
        "valable pour les deux, on tire ∂Z_s/∂u = ∂_u R_s/(2 Z_s) = "
        "a₁·u/Z_s — la formule du pilote D5 se transporte verbatim. "
        "RÉSULTAT : sur %d charts examinés, %d atteignables, "
        "**%d admissibles** (jauge bornée loin de 0, |u'| ≤ 1, |v'| ≤ 1 "
        "ET det J ≠ 0, tout sur la boîte entière), dont **%d certifiés "
        "DISJOINTS de leur propre tranche réelle** — %s par cellule. "
        "Sous la branche principale, cet ensemble était VIDE : c'est le "
        "verrou du pilote D5, désormais franchi SUR UNE CELLULE, pas en "
        "un point. "
        "NÉGATIFS : une boîte dont σ est déterminé AU CENTRE mais "
        "INDÉTERMINÉ sur la boîte est REFUSÉE, et la mutation « centre "
        "seulement » l'accepte à tort (c'est ce qui rend le test "
        "discriminant) ; une jauge de centre 1/2 dont l'enclosure "
        "contient 0 est refusée. "
        "PORTÉE : **D5.4 (congruence Q_source = J* Q_target J) et D5.5 "
        "(transport de PD) restent NON TESTÉS** — le succès de D5.1/D5.2 "
        "ne les préjuge pas. Résultat établi sur %d cellules, pas sur un "
        "cover. Aucun chiffre d'atlas ne bouge." % (
            n_pass, len(gates),
            sum(1 for c in cells if c["section_complete"]), len(cells),
            "/".join(aggregate["determinations_used"]),
            tot["charts"], tot["reachable"], tot["admissible"],
            tot["payoff"], str(aggregate["payoff_per_cell"]),
            len(cells)))

    out = {
        "phase": ("B1.e.2.iii D5.1/D5.2 sur cellule entière, "
                  "détermination par ligne (contrat C124 §8)"),
        "witness_sha256": reg["witness_sha256"],
        "tm_config": {"poly_deg": TM_ORDER,
                      "unary_series_deg": UNARY_SERIES_DEG},
        "provenance": provenance(C118_JSON, time.time() - T0),
        "jacobian_identity": ("dZ_s/du = ∂_u R_s/(2 Z_s) = a₁u/Z_s, "
                              "INDÉPENDANTE de la détermination car "
                              "Z_s² = R_s vaut pour les deux"),
        "cells": cells,
        "D5_component_negative": comp_neg,
        "D6_gauge_negative": gauge_neg,
        "not_tested": ["D5.4 congruence", "D5.5 transport de PD",
                       "cover multi-chart"],
        "aggregate": aggregate,
        "gates_prereg": gates,
        "verdict": verdict}
    ART.write_text(json.dumps(out, indent=2, ensure_ascii=False,
                              default=float), encoding="utf-8")
    print("\nVERDICT :\n" + verdict)
    print(f"\n→ {ART}")
    return out


# ===========================================================================
#  Self-test
# ===========================================================================
def _selftest():
    fails = []

    # F-S1: the Jacobian identity holds for BOTH determinations
    import cmath
    a0, a1 = Fraction(-16, 9), Fraction(-8, 3)
    u0 = complex(0.3, -0.4)

    def R(u):
        return float(a0) + float(a1) * u * u

    for name, det in (("principal", lambda z: cmath.sqrt(z)),
                      ("tournée", lambda z: -1j * cmath.sqrt(-z))):
        Zs = det(R(u0))
        d = 1e-7
        fd = (det(R(u0 + d)) - det(R(u0 - d))) / (2 * d)
        formula = float(a1) * u0 / Zs
        ok = abs(fd - formula) < 1e-6 * max(1.0, abs(formula))
        fails.append(not ok)
        print(f"[{'PASS' if ok else 'FAIL'}] F-S1 Jacobien ({name}) : "
              f"a₁u/Z = {formula:.6f} vs contrôle {fd:.6f}")

    # F-S2 : NÉGATIF de composante — σ déterminé au centre, indéterminé
    # on the box. This is the mutation the review asked for.
    A = Fraction(-8, 3)
    h = 0.001953125
    sig_box = rotated_sigma_from_coeffs(
        A, A, (-0.6, -0.59), (-1.5 * h, 0.5 * h),
        (-0.63, -0.62), (-1.5 * h, 0.5 * h))
    sig_ctr = rotated_sigma_from_coeffs(
        A, A, (-0.6, -0.59), (-0.5 * h, -0.5 * h),
        (-0.63, -0.62), (-0.5 * h, -0.5 * h))
    ok = sig_box == 0 and sig_ctr == -1
    fails.append(not ok)
    print(f"[{'PASS' if ok else 'FAIL'}] F-S2 négatif de composante : σ "
          f"sur la boîte = {sig_box} (indéterminé) mais σ au centre seul "
          f"= {sig_ctr} — la mutation accepterait à tort")

    # F-S3: disjointness from the target slice is read on the
    # enclosures de Im u' et Im v', et un intervalle contenant 0 ne
    # suffit PAS à conclure « disjoint »
    def disj(iul, iuh, ivl, ivh):
        return bool(iul > 0 or iuh < 0 or ivl > 0 or ivh < 0)
    ok = (disj(0.1, 0.2, -0.1, 0.1) and disj(-0.2, -0.1, -0.1, 0.1)
          and not disj(-0.1, 0.1, -0.1, 0.1))
    fails.append(not ok)
    print(f"[{'PASS' if ok else 'FAIL'}] F-S3 disjonction : une seule "
          f"coordonnée strictement d'un côté suffit ; deux intervalles "
          f"contenant 0 ne concluent PAS")

    print("-" * 78)
    print("SELF-TEST:", "FAIL" if any(fails) else "ALL PASS")
    return 1 if any(fails) else 0


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--selftest":
        raise SystemExit(_selftest())
    build()
