#!/usr/bin/env python3
"""
The CHART SELECTION CRITERION, strengthened, and swept again over
the 59 charts.

The certified component made the congruence step fall: with a genuinely certified component, the
NATIVE section of the target chart is not constructible, and the sweep gave
**0 out of 15** on the first cell. This script turns that into a selection
criterion in its own right and sweeps it again over EVERY cell.

THE GAP THIS SCRIPT REPAIRS. The earlier criterion bore on:

    (i)  le DOMAINE     — jauge `|Z_{g'}|` loin de 0, `|u'| ≤ 1`, `|v'| ≤ 1`
    (ii) le JACOBIEN    — `det J` d'enclosure excluant 0
    (iii) la DISJONCTION — `Im u'` ou `Im v'` de signe constant

Those three bear on the **coordinates** `(u', v')` of the target chart.
None says anything about the **radicands** of the chart's OWN section.
A chart can therefore lie outside its real slice while
having an uncertifiable section root, and that is exactly what
happened. The strengthened criterion adds:

    (iv) NATIVE SECTION: the three rows of the TARGET chart's section
         admit a CERTIFIED determination on the whole box
         (principale si `R'` évite (−∞,0], sinon tournée avec `σ'`
         déterminé — jamais au seul centre)

Checks pré-enregistrés :
  F1 REFINEMENT: the strengthened criterion is a SUBSET of the earlier
     one, chart by chart. A UNIVERSAL check: a criterion that
     accepted a chart refused before would not be a
     strengthening but a different criterion.
  F2 EXHAUSTIVENESS: the 59 charts are examined on EVERY cell,
     with no silent skip; the published verdict is a UNIVERSAL
     count, not an existence claim.
  F3 THE CRITERION BITES: it refuses at least one chart that the earlier
     one accepted. Otherwise it would add nothing and the
     resweep would be pointless. (If the criterion did not bite,
     this check MUST fail and say so.)
  F4 CENTRE-ONLY NEGATIVE CONTROL: the mutation that decides constructibility
     on the CENTRE alone (a box of half-width 0) must accept
     STRICTLY MORE charts than the version certified on the box.
     That is what separates "certified on the cell" from "seen at
     centre ».
  F5 SOUNDNESS OF THE CONSTRUCTIBLE ROWS: for any row declared
     constructible, whatever the determination, the enclosure of
     `Z'^2 - R'` must contain 0. A UNIVERSAL check over all accepted
     rows, hence NON-VACUOUS even if no whole chart passes.
  F6 IS THE OBSTRUCTION INHERITED? The conjecture is that the row
     blocking the target chart carries the SAME ambient coordinate `Z_s` as
     the one that blocked the source. This check MEASURES the proportion; it does not
     impose it. It fails only if the measurement is missing.
  F9 SUBDIVISION DOES NOT ERASE THE OBSTRUCTION: on a lost chart,
     on subdivise la cellule source en 16^d sous-cellules et on
     recounts the refusals. Were the obstruction an ARTEFACT
     D'ENCLOSURE (surestimation de `Im R'`), les refus DISPARAÎTRAIENT
     at increasing depth. The check requires instead that the number
     of refusals GROW like `4^d`, the signature of a branch locus
     of CODIMENSION 2, hence of a REAL GEOMETRIC obstruction. The
     check échoue si les refus disparaissent OU s'ils croissent en
     `16^d` (codimension 0).
  F7 PUBLISHED VERDICT: the final count `n_strong` is published for
     EVERY cell, including when it is 0. A sweep that did not
     publierait que ses succès mentirait par omission.

This script certifies NEITHER congruence NOR positivity transport: it selects the
(cell, chart) pairs on which those questions make sense.

Sorties : results/chart_selection_criterion.json
Usage   : chart_selection_criterion.py [--selftest]
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
SHELL_JSON = RES / "shell_enumeration.json"
ART = RES / "chart_selection_criterion.json"

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
#  (iv) THE NATIVE SECTION OF THE TARGET CHART: the added criterion
# ===========================================================================
def native_section_constructible(S2, g2, up, vp):
    """The three rows of the target chart's OWN section, with their
    determination certified on the WHOLE BOX.

    Retourne (rows, ok). Chaque ligne porte sa détermination, son `σ'`,
    the enclosures of `Re R'` and `Im R'`, and, if it is constructible,
    the enclosure of `Z'^2 - R'`, which must contain 0 (check F5).

    The component is certified in TWO steps:
      (i) DIRECT enclosure of `Im R'`: on the target chart it is
          not constrained to contain 0, so it is often sharper;
     (ii) à défaut, la règle par SIGNES DE COEFFICIENTS.
    If both are undetermined, the row is REFUSED. Never by trial.
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
            # F5: soundness of the row, whatever the determination
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
    print(f"the chart criterion CRITÈRE DE CHART RENFORCÉ : TM ({TM_ORDER},"
          f"{UNARY_SERIES_DEG}), {N_CELLS} cellules × 59 charts")
    print("=" * 78)
    load_canonical_MH()
    shell = json.loads(SHELL_JSON.read_text(encoding="utf-8"))
    residual = [x for x in shell["cells"] if x["still_branch"]]
    log(f"résidu : {len(residual)} paires")

    cells = []
    step = max(1, len(residual) // max(1, N_CELLS))
    for i in range(min(N_CELLS, len(residual))):
        r = residual[i * step]
        S, g, eps = tuple(r["S"]), r["g"], tuple(r["eps"])
        c = [float.fromhex(x) for x in r["center_hex"]]
        hw = float.fromhex(r["hw_hex"])
        Z, dZ, srows = build_section(S, g, eps, c, hw)
        # the centre-only MUTATION: the SAME cell seen at one point
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
                # centre-only mutation, on the SAME chart
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
        # is the obstruction inherited?
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

    # --- F9: does subdivision erase the obstruction? -------------------------
    # This is THE question that 62 -> 0 raises. Two competing explanations:
    #   (a) ARTEFACT D'ENCLOSURE — `Im R'` est surestimé et straddle 0 par
    #       excès de largeur. Alors subdiviser doit FAIRE DISPARAÎTRE les
    #       refusals, and the resweep need only go one level down.
    #   (b) OBSTRUCTION GÉOMÉTRIQUE — le lieu `{Im R' = 0} ∩ {Re R' < 0}`
    #       really crosses the cell. Then the refusals PERSIST, and
    #       their number grows like the measure of the locus: `4^d` for
    #       codimension 2 in a real box of dimension 4, against
    #       `16^d` sous-cellules au total.
    # The two predictions are DISJOINT and the count decides.
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

    # --- cross-check INDÉPENDANT contre l'artefact the chart criterion/the Jacobian criterion ----------------
    # `n_prev_criterion` DOIT reproduire `payoff_per_cell` de
    # `d5_fullcell` : deux scripts, même sélection de cellules, même
    # earlier criterion. A disagreement would signal drift in one or
    # the other, which is what a review calls a cross-check that
    # n'existait dans aucun script.
    cross = {"checked": False}
    fc = RES / "full_cell_charts.json"
    try:
        d = json.loads(fc.read_text(encoding="utf-8"))
        ref = list(d["aggregate"]["payoff_per_cell"])[:len(cells)]
        mine = [c["n_prev_criterion"] for c in cells]
        cross = {"checked": True, "reference": ref, "measured": mine,
                 "agrees": ref == mine, "source": fc.name}
    except Exception as exc:
        cross = {"checked": False, "error": str(exc)[:120]}

    # ------------------------------------------------------------------ checks
    all_rows = [row for c in cells for x in c["charts"]
                for row in (x.get("native_rows") or [])]
    built = [row for row in all_rows if row["determination"] is not None]
    checks = {
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
        "**%d charts** passed the earlier criterion (domain, Jacobian and "
        "disjointness) on %d cells; **%d** pass the STRENGTHENED "
        "criterion, which further requires the NATIVE section of the target chart "
        "to be constructible with CERTIFIED components on the whole "
        "box. %d charts are lost. Of %d blocking rows, **%d "
        "carry an ambient coordinate that the source already had to "
        "rotate**, so the obstruction is INHERITED in %s of cases."
        % (tot["prev"], len(cells), tot["strong"], tot["lost"],
           tot["blocking_total"], tot["blocking_inherited"],
           ("%.0f%%" % (100.0 * tot["blocking_inherited"]
                        / tot["blocking_total"]))
           if tot["blocking_total"] else "n/a"))

    art = {"artifact": "chart_selection_criterion",
           "claim": "The chart selection criterion, STRENGTHENED by the "
                    "constructibility of the NATIVE section of the target "
                    "chart, and swept again over the 59 charts.",
           "totals": tot, "verdict": verdict,
           "cross_check_d5_fullcell": cross,
           "subdivision_probe": sub,
           "cells": cells, "checks": checks,
           "checks_passed": sum(1 for v in checks.values() if v),
           "checks_total": len(checks),
           "provenance": provenance([SHELL_JSON], time.time() - T0)}
    ART.parent.mkdir(parents=True, exist_ok=True)
    ART.write_text(json.dumps(art, indent=2, ensure_ascii=False),
                   encoding="utf-8")
    print("-" * 78)
    for k, v in checks.items():
        print(f"  {'PASS' if v else 'FAIL'}  {k}")
    print("-" * 78)
    print(verdict)
    print(f"→ {ART}")
    return 0 if all(checks.values()) else 1


# ===========================================================================
#  Self-test
# ===========================================================================
def _selftest():
    ok = True

    def chk(name, cond):
        nonlocal ok
        print(f"  {'ok  ' if cond else 'FAIL'} {name}")
        ok = ok and bool(cond)

    # T1: a POSITIVE real constant is constructible on the principal branch
    from .taylor_models import CONE
    Rp = TMC.const(CIV.from_complex(complex(4.0, 0.0)))
    try:
        Rp.sqrt_principal()
        chk("T1 constante positive : principale OK", True)
    except BranchCutError:
        chk("T1 constante positive : principale OK", False)

    # T2: a NEGATIVE real constant refuses both the principal and the
    #      rotated one (undetermined component), which is the heart of the matter
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

    # T3: a radicand on the POSITIVE ray: the principal branch accepts it,
    #      but the ROTATED one is the continuation; it builds
    #      as soon as the component is certified.
    Rr = TMC.const(CIV.from_complex(complex(-4.0, 0.5)))
    try:
        got = tm_sqrt_rotated(Rr, 1)
    except BranchCutError:
        got = None
    chk("T3 Im > 0 : la tournée σ=+1 construit", got is not None)

    # T4: soundness, Z^2 - R contains 0 for the rotated branch of T3
    if got is not None:
        d = got * got - Rr
        dr, di = _rng(d), _rng(d, True)
        chk("T4 soundness Z² − R ∋ 0",
            dr[0] <= 0 <= dr[1] and di[0] <= 0 <= di[1])
    else:
        chk("T4 soundness Z² − R ∋ 0", False)

    # T5 NEGATIVE CONTROL: the strengthened criterion must NEVER accept a chart
    #      that the earlier one refused. This is checked on the logic itself.
    def strong(prev, native):
        return bool(prev and native)
    chk("T5 négatif : strong ⊆ prev (les 4 combinaisons)",
        all(not strong(p, n) or p
            for p in (False, True) for n in (False, True))
        and not strong(False, True))

    # T6 NEGATIVE CONTROL: a mutation taking `strong = native` alone
    #      would accept a chart refused by the earlier criterion. Without this
    #      contre-exemple, T5 serait satisfait par n'importe quoi.
    def strong_broken(prev, native):
        return bool(native)
    chk("T6 negative control: the mutation `strong = native` VIOLATES the refinement",
        strong_broken(False, True) and not False)

    # T7 : `native_section_constructible` REFUSE quand `σ'` est
    #      indéterminé — pas d'essai silencieux
    class _FakeUV:
        pass
    chk("T7 rotated_sigma_from_coeffs returns 0 on a straddle",
        rotated_sigma_from_coeffs(Fraction(1), Fraction(1),
                                  (-1.0, 1.0), (-1.0, 1.0),
                                  (-1.0, 1.0), (-1.0, 1.0)) in (0, None))

    print("SELF-TEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_selftest() if "--selftest" in sys.argv else build())
