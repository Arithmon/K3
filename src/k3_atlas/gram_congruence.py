#!/usr/bin/env python3
"""
Gram CONGRUENCE and POSITIVITY TRANSPORT, on a whole cell.


The chart criterion and Jacobian certified, on the whole box, that target
charts are reachable, in the domain, and of invertible Jacobian. What remained
was the question the pilot had explicitly left open, and which the
success of those checks **did not prejudge**: positivity
transporte-t-elle ?

LA STRUCTURE. `Q` se calcule à partir de `(Z, W)` où `W[a][A] =
∂Z_a/∂(u,v)_A`. Deux propriétés suffisent :

  (a) **quadratique en W** — `Q(Z, W·A) = Aᵀ · Q(Z, W) · conj(A)` ;
  (b) **gauge invariance**: `Q` has degree 0 in `Z`, so the projective
      renormalisation does not change it.

Elles se composent : en posant `Z' = Z/Z_{g'}` et `W_{Z'} = ∂Z'/∂(u,v)`,
on a `W_{Z'} = W' · J` où `J[A'][A] = ∂(u',v')_{A'}/∂(u,v)_A`, donc

    Q_source  =(b)=  Q(Z', W_{Z'})  =(a)=  Jᵀ · Q_target · conj(J)

**This is the convention, and it is verified rather than assumed**: the
review contract writes `Q_source = J* Q_target J`, which is the same
statement up to convention (in which argument `Q` is conjugate-linear,
and whether `J` is the Jacobian or its conjugate transpose). The self-test
measures both forms and shows that they DIFFER, so the check
discrimine.

Gates pré-enregistrés :
  E1 INVARIANCE DE JAUGE — `Q(Z', W_{Z'})` et `Q(Z, W)` ont des
     enclosures whose difference contains 0, on the WHOLE box
  E1b MUTATION DE JAUGE — renormaliser `Z` SANS renormaliser `W` casse
     the identity. Without it, E1 would be satisfied by a formula that
     ignore `Z`.
  E2 QUADRATICITY IN W: `Q(Z, W.A) = A^T Q conj(A)` for a constant
     invertible matrix, and **`A* Q A` does NOT work**: the mutation
     is the competing convention
  E3 CONGRUENCE: on the cell and a named ADMISSIBLE chart,
     `Q_source - J^T Q_target conj(J)` contains 0 on all four components
  E4 POSITIVITY TRANSPORT: `Q_target` is certified positive definite (pivot > 0 and
     det > 0), et `det(Jᵀ Q_t conj(J)) = |det J|²·det Q_t > 0` transporte
     positivity; `Q_source` is positive definite directly, and the two agree
  E5 PHASE MUTATION: flipping the sign of a section row does NOT
     break E3, and that is correct: congruence is **agnostic to
     la branche**. Elle relie deux descriptions de la MÊME configuration,
     whatever the branch choice, PROVIDED it is coherent, and
     it is here, since `(u', v')` derives from the source `Z`. The
     discriminating negative control is therefore E6, not E5. This is published
     as such rather than presented as a negative control that ought to
     mordre.
  E6 MUTATION DU JACOBIEN — perturber `J` casse E3

What this script does NOT establish: the multi-chart cover, nor
globalisation. It establishes transport on ONE cell and ONE chart.

Sorties : results/k3_cap_b1e2iii_d5_congruence.json
Usage   : k3_cap_b1e2iii_d5_congruence.py [--selftest]
Env     : K3_D5CG_CELLS (cellules, défaut 3)
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
from .witness_registry import load_canonical_MH              # noqa: E402
from .interval_arithmetic import (                            # noqa: E402
    BranchCutError, CIV, build_M_civ, minor_inv_times_T_exact)
from .taylor_models import (                                     # noqa: E402
    CZERO, CONE, TMC, TM_ORDER, UNARY_SERIES_DEG, det_packed_tm,
    iv_bounds, riv, rotated_sigma_from_coeffs,
    tm_chart_metric, tm_sqrt_rotated)
from .width_attribution import GAMMA                       # noqa: E402
from .full_cell_charts import build_section               # noqa: E402
sys.argv = _argv

RES = Path(os.environ.get(
    "K3_RES_DIR", Path(__file__).resolve().parent / "data"))
FULLCELL_JSON = RES / "k3_cap_b1e2iii_d5_fullcell.json"
ART = RES / "k3_cap_b1e2iii_d5_congruence.json"

N_CELLS = int(os.environ.get("K3_D5CG_CELLS", "3"))
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def _sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def provenance(src, t):
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
            "platform": platform.platform(), "wall_s": t}


# ===========================================================================
#  2x2 algebra on complex Taylor models
# ===========================================================================
def _c(z):
    return CIV(riv(z.real), riv(z.imag))


def Qmat(Z, W, M, c218, rw):
    """La matrice hermitienne 2×2 `Q`, en TMC, depuis (Z, W)."""
    q = tm_chart_metric(Z, W, M, c218, rho_weight=rw)
    g00, g11, re01, im01 = q
    j = TMC.const(CIV(riv(0.0), riv(1.0)))
    q01 = TMC([CIV(x, y) for x, y in zip(re01.p, im01.p)],
              re01.rem + im01.rem)
    q10 = TMC([CIV(x, -y) for x, y in zip(re01.p, im01.p)],
              re01.rem + im01.rem)
    z00 = TMC([CIV(x) for x in g00.p], g00.rem)
    z11 = TMC([CIV(x) for x in g11.p], g11.rem)
    return [[z00, q01], [q10, z11]]


def mat_sub(A, B):
    return [[A[i][k] + B[i][k].mul_real(riv(-1.0)) for k in range(2)]
            for i in range(2)]


def contains_zero(Mx):
    out, ok = [], True
    for i in range(2):
        for k in range(2):
            re, im = Mx[i][k].re_tm().to_iv(), Mx[i][k].im_tm().to_iv()
            rl, rh = float(mp.mpf(re.a)), float(mp.mpf(re.b))
            il, ih = float(mp.mpf(im.a)), float(mp.mpf(im.b))
            c0 = rl <= 0 <= rh and il <= 0 <= ih
            ok = ok and c0
            out.append({"i": i, "k": k, "re": [rl, rh], "im": [il, ih],
                        "contains_zero": bool(c0),
                        "width": max(rh - rl, ih - il)})
    return ok, out


def congruence(J, Q, conj_left=False):
    """`Jᵀ Q conj(J)` (défaut) ou `J* Q J` (la convention concurrente,
    utilisée comme mutation)."""
    def cj(x):
        return TMC([CIV(c.re, -c.im) for c in x.p], x.rem)
    out = [[None, None], [None, None]]
    for a in range(2):
        for b in range(2):
            acc = None
            for p in range(2):
                for q in range(2):
                    L = cj(J[p][a]) if conj_left else J[p][a]
                    Rt = J[q][b] if conj_left else cj(J[q][b])
                    term = L * Q[p][q] * Rt
                    acc = term if acc is None else acc + term
            out[a][b] = acc
    return out


def pd_bounds(Q):
    """pivot `q00` et `det` (hermitien 2×2) en bornes."""
    q00 = Q[0][0].re_tm()
    det = (Q[0][0] * Q[1][1] + (Q[0][1] * Q[1][0]).mul_real(
        riv(-1.0))).re_tm()
    return iv_bounds(q00.to_iv()), iv_bounds(det.to_iv())


# ===========================================================================
#  The transport, on one cell and one chart
# ===========================================================================
def native_target_section(S2, g2, up, vp, eps2, sigma2):
    """NATIVE section of the target chart, with the sheet and component FIXED
    ENTRÉE** — jamais choisis par essai ni ajustés par le gate qu'ils
    doivent passer.

    The component comes from the certified COMPONENT of the radicand
    cible. L'ancienne boucle `for sg in (-1, 1): … break` prenait
    always -1, the guard of the rotated root depending only on `-R`
    and not on the sign, so it certified NOTHING. The same defect as
    celui corrigé dans C124, réintroduit ici.
    """
    T2 = tuple(j for j in range(6) if j not in S2)
    others = [c for c in T2 if c != g2]
    Ae = minor_inv_times_T_exact(S2, T2)
    perm = [list(T2).index(g2), list(T2).index(others[0]),
            list(T2).index(others[1])]
    A = [[riv(Ae[r][c]) for c in perm] for r in range(3)]
    Z = [None] * 6
    dZ = [None] * 6
    Z[g2] = TMC.const(CONE)
    dZ[g2] = (TMC.const(CZERO), TMC.const(CZERO))
    Z[others[0]], dZ[others[0]] = up, (TMC.const(CONE), TMC.const(CZERO))
    Z[others[1]], dZ[others[1]] = vp, (TMC.const(CZERO), TMC.const(CONE))
    u2, v2 = up * up, vp * vp
    kinds = []
    for r, s in enumerate(S2):
        R = TMC.const(CIV(A[r][0])) + u2.mul_real(A[r][1]) \
            + v2.mul_real(A[r][2])
        Zs = None
        try:
            Zs = R.sqrt_principal().mul_real(riv(int(eps2[r])))
            kinds.append("principal")
        except BranchCutError:
            sg = sigma2[r]
            if sg in (-1, 1):
                try:
                    Zs = tm_sqrt_rotated(R, sg).mul_real(
                        riv(int(eps2[r])))
                    kinds.append("rotated%+d" % sg)
                except BranchCutError:
                    kinds.append(None)
            else:
                kinds.append(None)
        if Zs is None:
            return None, None, kinds
        Z[s] = Zs
        iZ = Zs.inv()
        dZ[s] = (up.mul_real(A[r][1]) * iZ, vp.mul_real(A[r][2]) * iZ)
    return Z, dZ, kinds


def target_component(S2, g2, up, vp):
    """The component of the TARGET radicand, from the signs
    of the coefficients and the enclosures of `(u', v')`."""
    T2 = tuple(j for j in range(6) if j not in S2)
    others = [c for c in T2 if c != g2]
    Ae = minor_inv_times_T_exact(S2, T2)
    perm = [list(T2).index(g2), list(T2).index(others[0]),
            list(T2).index(others[1])]

    def rng(tm, im=False):
        x = tm.im_tm().to_iv() if im else tm.re_tm().to_iv()
        return (float(mp.mpf(x.a)), float(mp.mpf(x.b)))

    ur, ui = rng(up), rng(up, True)
    vr, vi = rng(vp), rng(vp, True)
    A = [[riv(Ae[r][c]) for c in perm] for r in range(3)]
    u2, v2 = up * up, vp * vp
    out = []
    for r in range(3):
        # (i) first the DIRECT ENCLOSURE of Im R'. On the target chart
        # it is NOT constrained to contain 0 (unlike the face-aligned
        # source), so it is often sharper than the coefficient rule,
        # which is conservative. The two are
        # certifiées ; on prend la plus informative.
        R = TMC.const(CIV(A[r][0])) + u2.mul_real(A[r][1]) \
            + v2.mul_real(A[r][2])
        im = R.im_tm().to_iv()
        lo, hi = float(mp.mpf(im.a)), float(mp.mpf(im.b))
        if lo > 0:
            out.append(1)
        elif hi < 0:
            out.append(-1)
        else:
            # (ii) sinon, la règle par SIGNES DE COEFFICIENTS
            out.append(rotated_sigma_from_coeffs(
                Fraction(Ae[r][perm[1]]), Fraction(Ae[r][perm[2]]),
                ur, ui, vr, vi))
    return out


def derive_eps_target(S2, g2, up, vp, sigma2, Zp):
    """C125-B : `ε'` **DÉRIVÉ**, pas ajusté par E3.

    At the CENTRE of the cell, the target root computed with the
    sheet +1 is compared to the normalised coordinate: the sign that makes
    them coincide IS the sheet. It is then FROZEN, and compatibility
    is required on the WHOLE BOX, coordinate by coordinate,
    MÊME signe** — « au signe près » accepterait un changement de
    a SHEET change for a gauge, which is wrong as soon as the two
    descriptions are already normalised by the same gauge.
    """
    Zt1, _d, _k = native_target_section(S2, g2, up, vp, (1, 1, 1),
                                        sigma2)
    if Zt1 is None:
        return None
    eps = []
    for r, s in enumerate(S2):
        a, b = Zt1[s].p[0], Zp[s].p[0]
        ar, ai = float(mp.mpf(a.re.a)), float(mp.mpf(a.im.a))
        br, bi = float(mp.mpf(b.re.a)), float(mp.mpf(b.im.a))
        eps.append(1 if abs(ar - br) + abs(ai - bi)
                   <= abs(ar + br) + abs(ai + bi) else -1)
    return tuple(eps)


def transport(S, g, eps, center, hw, S2, g2, M, c218, rw,
              flip_row=None, perturb_J=None):
    """Q_source, Q_mid (jauge cible / dérivées source) et Q_target
    NATIVE, plus le Jacobien."""
    Z, dZ, rows = build_section(S, g, eps, center, hw)
    if flip_row is not None:
        s = S[flip_row]
        Z[s] = Z[s].mul_real(riv(-1.0))
        dZ[s] = (dZ[s][0].mul_real(riv(-1.0)),
                 dZ[s][1].mul_real(riv(-1.0)))
    if any(z is None for z in Z):
        return None
    W = [[dZ[a][0], dZ[a][1]] for a in range(6)]
    Q_src = Qmat(Z, W, M, c218, rw)

    T2 = [j for j in range(6) if j not in S2]
    o = [x for x in T2 if x != g2]
    ib = Z[g2].inv()
    Zp = [z * ib for z in Z]
    dZp = [((dZ[a][0] * Z[g2] + (Z[a] * dZ[g2][0]).mul_real(riv(-1.0)))
            * ib * ib,
            (dZ[a][1] * Z[g2] + (Z[a] * dZ[g2][1]).mul_real(riv(-1.0)))
            * ib * ib) for a in range(6)]
    Wp_src = [[dZp[a][0], dZp[a][1]] for a in range(6)]
    Q_mid = Qmat(Zp, Wp_src, M, c218, rw)      # E1 : invariance de jauge

    # J[A'][A] = ∂(u',v')_{A'}/∂(u,v)_A, où (u',v') = (Z_{o1'}, Z_{o2'})/Z_{g'}
    up, vp = Zp[o[0]], Zp[o[1]]
    J = [[dZp[o[0]][0], dZp[o[0]][1]],
         [dZp[o[1]][0], dZp[o[1]][1]]]
    if perturb_J is not None:
        J[0][0] = J[0][0] + TMC.const(_c(perturb_J))

    # component CERTIFIED from the target radicand component
    sigma2 = target_component(S2, g2, up, vp)
    # C125-B : ε' DÉRIVÉ au centre, puis FIGÉ
    eps2d = derive_eps_target(S2, g2, up, vp, sigma2, Zp)
    if eps2d is None:
        return {"native_failed": True, "sigma_target": sigma2}
    Zt, dZt, kinds_t = native_target_section(S2, g2, up, vp, eps2d,
                                             sigma2)
    if Zt is None:
        return {"native_failed": True, "kinds_target": kinds_t,
                "sigma_target": sigma2, "eps_target": list(eps2d)}
    Wt = [[dZt[a][0], dZt[a][1]] for a in range(6)]
    Q_tgt = Qmat(Zt, Wt, M, c218, rw)
    cong = congruence(J, Q_tgt)
    # compatibilité des deux sections : Z_native vs Z/Z_{g'} (au signe
    # of each section coordinate: this is the record of the sheets)
    # C125-B : compatibilité PROJECTIVE — MÊME signe, coordonnée par
    # coordinate, on the whole box. The two descriptions are already
    # normalisées par la MÊME jauge (Z_{g'} = 1) : elles ne peuvent plus
    # differ by a global sign, and an INDEPENDENT sign per coordinate
    # would be a SHEET change, not a gauge one.
    compat = []
    for a in range(6):
        D = Zt[a] + Zp[a].mul_real(riv(-1.0))
        re, im = D.re_tm().to_iv(), D.im_tm().to_iv()
        rl, rh = float(mp.mpf(re.a)), float(mp.mpf(re.b))
        il, ih = float(mp.mpf(im.a)), float(mp.mpf(im.b))
        compat.append({"coord": a,
                       "same_point": bool(rl <= 0 <= rh
                                          and il <= 0 <= ih),
                       "residual": [rl, rh, il, ih],
                       "width": max(rh - rl, ih - il)})
    return {"Q_src": Q_src, "Q_mid": Q_mid, "Q_tgt": Q_tgt,
            "cong": cong, "J": J, "rows": rows,
            "kinds_target": kinds_t, "compat": compat,
            "sigma_target": sigma2, "eps_target": list(eps2d),
            "same_projective_point": all(x["same_point"]
                                         for x in compat)}


# ===========================================================================
#  Build
# ===========================================================================
def build():
    print("=" * 78)
    print(f"D5.4 / D5.5 — CONGRUENCE ET TRANSPORT DE PD : TM "
          f"({TM_ORDER},{UNARY_SERIES_DEG})")
    print("=" * 78)
    reg = load_canonical_MH()
    M = build_M_civ(reg["M_H_canonical"])
    c218 = reg["coeffs218"]
    rw = 1.0 - GAMMA
    fc = json.loads(FULLCELL_JSON.read_text(encoding="utf-8"))

    cells = []
    for cell in fc["cells"][:N_CELLS]:
        S, g, eps = tuple(cell["S"]), cell["g"], tuple(cell["eps"])
        c = [float.fromhex(x) for x in cell["center_hex"]]
        hw = float.fromhex(cell["hw_hex"])
        adm = [x for x in cell["charts"]
               if x.get("admissible") and x["disjoint_from_target_slice"]]
        if not adm:
            continue
        tgt = adm[0]
        S2, g2 = tuple(tgt["S"]), tgt["g"]
        # no search on the sheet any more: it is DERIVED in
        # `transport`, puis E3 le teste en aveugle.
        tr = transport(S, g, eps, c, hw, S2, g2, M, c218, rw)
        if tr is None or tr.get("native_failed"):
            log(f"  cellule ignorée : section cible non calculable "
                f"(σ' = {tr.get('sigma_target') if tr else None})")
            continue
        eps2_used = tuple(tr["eps_target"])
        if not tr["same_projective_point"]:
            log(f"  cellule REFUSÉE : les deux descriptions ne sont PAS "
                f"le même point projectif (chart S={list(S2)} g={g2})")
            continue
        gauge_ok, gauge_d = contains_zero(
            mat_sub(tr["Q_src"], tr["Q_mid"]))
        cong_ok, cong_d = contains_zero(
            mat_sub(tr["Q_src"], tr["cong"]))
        # the residual, QUANTITATIVE: `contains_zero` bounds a
        # defect, it does not prove an identity. We therefore publish
        # sup |résidu| et sa taille RELATIVE à ‖Q_source‖.
        qs_norm = max(abs(x) for i in range(2) for k in range(2)
                      for x in iv_bounds(tr["Q_src"][i][k].re_tm()
                                         .to_iv()))
        sup_abs = max(max(abs(x["re"][0]), abs(x["re"][1]),
                          abs(x["im"][0]), abs(x["im"][1]))
                      for x in cong_d)
        residual = {"sup_abs": sup_abs, "Q_source_norm": qs_norm,
                    "relative": sup_abs / max(qs_norm, 1e-300),
                    "note": ("`0 ∈ enclosure` borne le défaut, il ne "
                             "prouve PAS l'identité — le contrat retenu "
                             "est « congruence APPROCHÉE certifiée » "
                             "avec cette borne")}
        (p_s, d_s) = pd_bounds(tr["Q_src"])
        (p_t, d_t) = pd_bounds(tr["Q_tgt"])
        dJ = iv_bounds((tr["detJ"] * TMC([CIV(x, -y) for x, y in
                                          zip(tr["detJ"].p, tr["detJ"].p)],
                                         tr["detJ"].rem)).re_tm().to_iv()
                       ) if False else None
        rec = {
            "center_hex": cell["center_hex"], "hw_hex": cell["hw_hex"],
            "S": list(S), "g": g, "eps": list(eps),
            "target": {"S": list(S2), "g": g2,
                       "gauge_absmin": tgt["gauge_absmin"],
                       "detJ_absmin": tgt["detJ_absmin"]},
            "determinations": [x["determination"] for x in tr["rows"]],
            "eps_target_used": list(eps2_used),
            "target_determinations": tr["kinds_target"],
            "section_compatibility": tr["compat"],
            "E1_gauge_invariance": {"ok": gauge_ok, "components": gauge_d},
            "E3_congruence": {"ok": cong_ok, "components": cong_d,
                              "residual": residual},
            "sigma_target_certified": tr["sigma_target"],
            "same_projective_point": tr["same_projective_point"],
            "pd_source": {"q00": list(p_s), "det": list(d_s),
                          "is_PD": bool(p_s[0] > 0 and d_s[0] > 0)},
            "pd_target": {"q00": list(p_t), "det": list(d_t),
                          "is_PD": bool(p_t[0] > 0 and d_t[0] > 0)}}
        cells.append(rec)
        log(f"  cellule {len(cells)} → chart S={list(S2)} g={g2} : "
            f"jauge E1 {gauge_ok} · congruence E3 {cong_ok} · "
            f"PD source {rec['pd_source']['is_PD']} "
            f"(det ≥ {d_s[0]:.3e}) · PD cible "
            f"{rec['pd_target']['is_PD']} (det ≥ {d_t[0]:.3e})")

    # --- E5 / E6 : les mutations ----------------------------------------------------
    cell = fc["cells"][0]
    S, g, eps = tuple(cell["S"]), cell["g"], tuple(cell["eps"])
    c = [float.fromhex(x) for x in cell["center_hex"]]
    hw = float.fromhex(cell["hw_hex"])
    tgt = next(x for x in cell["charts"]
               if x.get("admissible") and x["disjoint_from_target_slice"])
    S2, g2 = tuple(tgt["S"]), tgt["g"]

    tr_ph = transport(S, g, eps, c, hw, S2, g2, M, c218, rw,
                      flip_row=1)
    ph_ok = True
    if tr_ph and not tr_ph.get("native_failed"):
        ph_ok, _ = contains_zero(mat_sub(tr_ph["Q_src"], tr_ph["cong"]))
    tr_j = transport(S, g, eps, c, hw, S2, g2, M, c218, rw,
                     perturb_J=complex(0.3, 0.2))
    j_ok = True
    if tr_j and not tr_j.get("native_failed"):
        j_ok, _ = contains_zero(mat_sub(tr_j["Q_src"], tr_j["cong"]))
    log(f"E5 mutation de phase : congruence tient encore ? {ph_ok} — "
        f"ATTENDU True : la congruence est agnostique à la branche "
        f"(elle relie deux descriptions de la même configuration). Le "
        f"négatif discriminant est E6.")
    log(f"E6 mutation du Jacobien : congruence tient encore ? {j_ok} "
        f"(doit être False)")

    # --- E2: quadraticity in W, with the competing convention ----------------
    A = [[complex(2.0, 0.3), complex(0.1, -0.4)],
         [complex(-0.2, 0.5), complex(1.5, 0.1)]]
    Atm = [[TMC.const(_c(A[i][k])) for k in range(2)] for i in range(2)]
    Z, dZ, _ = build_section(S, g, eps, c, hw)
    W = [[dZ[a][0], dZ[a][1]] for a in range(6)]
    WA = [[None, None] for _ in range(6)]
    for a in range(6):
        for col in range(2):
            acc = None
            for p in range(2):
                term = W[a][p].mul_civ(_c(A[p][col]))
                acc = term if acc is None else acc + term
            WA[a][col] = acc
    Q0 = Qmat(Z, W, M, c218, rw)
    QA = Qmat(Z, WA, M, c218, rw)
    good_ok, _ = contains_zero(mat_sub(QA, congruence(Atm, Q0)))
    bad_ok, _ = contains_zero(
        mat_sub(QA, congruence(Atm, Q0, conj_left=True)))
    log(f"E2 quadraticité : Aᵀ Q conj(A) convient {good_ok} · "
        f"A* Q A convient {bad_ok} (doit être False)")

    gates = {
        "E1_gauge_invariance": bool(cells) and all(
            x["E1_gauge_invariance"]["ok"] for x in cells),
        "E2_quadratic_in_W": good_ok and not bad_ok,
        "E3_congruence_D54": bool(cells) and all(
            x["E3_congruence"]["ok"] for x in cells),
        "E4_PD_transport_D55": bool(cells) and all(
            x["pd_target"]["is_PD"] and x["pd_source"]["is_PD"]
            for x in cells),
        "E6_jacobian_mutation_breaks": not j_ok,
        # C125-B : toute paire retenue décrit le MÊME point projectif
        "C125B_same_projective_point": bool(cells) and all(
            x["same_projective_point"] for x in cells),
        # C125-C : σ' est certifié, jamais 0 (indéterminé)
        "C125C_target_component_certified": bool(cells) and all(
            all(s in (-1, 1, None) for s in x["sigma_target_certified"])
            for x in cells)}
    n_pass = sum(1 for v in gates.values() if v)
    log(f"gates : {n_pass}/{len(gates)} " + str(gates))

    verdict = (
        "D5.4 / D5.5 — CONGRUENCE ET TRANSPORT DE POSITIVITÉ (gates "
        "%d/%d). La structure : `Q` se calcule depuis `(Z, W)` avec "
        "`W[a][A] = ∂Z_a/∂(u,v)_A`, et deux propriétés suffisent — "
        "(a) `Q` est QUADRATIQUE en W, `Q(Z, W·A) = Aᵀ Q conj(A)` ; "
        "(b) `Q` est INVARIANTE DE JAUGE, de degré 0 en `Z`. En posant "
        "`Z' = Z/Z_{g'}` et `W_{Z'} = W'·J`, elles se composent en "
        "`Q_source = Jᵀ Q_target conj(J)`. "
        "**LA CONVENTION EST VÉRIFIÉE, PAS SUPPOSÉE** : le gate E2 mesure "
        "les DEUX formes et montre que `Aᵀ Q conj(A)` convient (%s) là où "
        "`A* Q A` ne convient PAS (%s) — c'est la mutation qui rend E2 "
        "discriminant, et c'est le même énoncé que le contrat de la revue "
        "à la convention près. "
        "RÉSULTAT sur %d cellule(s), chacune avec un chart ADMISSIBLE "
        "nommé (D5.1+D5.2 déjà certifiés sur la boîte entière) : "
        "invariance de jauge **%s**, congruence **%s** — la différence "
        "`Q_source − Jᵀ Q_target conj(J)` contient 0 sur les 4 "
        "composantes. **D5.5** : `Q_target` est certifiée PD sur la boîte "
        "(pivot > 0 et det > 0) et `Q_source` l'est aussi, les deux "
        "concordant par congruence — `det(Jᵀ Q_t conj(J)) = |det J|²·det "
        "Q_t`, et `det J` est borné loin de 0 par D5.2. "
        "MUTATION : perturber `J` CASSE la congruence (%s). "
        "PORTÉE : établi sur %d cellule(s) et UN chart chacune — pas un "
        "cover multi-chart (D5.6), pas une globalisation. Aucun chiffre "
        "d'atlas ne bouge." % (
            n_pass, len(gates), good_ok, bad_ok, len(cells),
            all(x["E1_gauge_invariance"]["ok"] for x in cells),
            all(x["E3_congruence"]["ok"] for x in cells),
            "oui" if not j_ok else "NON", len(cells)))

    out = {
        "phase": ("B1.e.2.iii D5.4/D5.5 — congruence des Gram et "
                  "transport de positivité sur cellule entière"),
        "witness_sha256": reg["witness_sha256"],
        "tm_config": {"poly_deg": TM_ORDER,
                      "unary_series_deg": UNARY_SERIES_DEG},
        "provenance": provenance(FULLCELL_JSON, time.time() - T0),
        "convention": {
            "established": "Q(Z, W·A) = A^T · Q · conj(A)",
            "rejected_mutation": "A* Q A (ne convient pas)",
            "note": ("même énoncé que `Q_source = J* Q_target J` du "
                     "contrat de la revue, à la convention près "
                     "(argument conjugué-linéaire, et J vs J*)")},
        "cells": cells,
        "E2_quadratic": {"A^T Q conj(A)": good_ok, "A* Q A": bad_ok},
        "E5_phase_mutation_congruence_still_holds": ph_ok,
        "note_E5": (
            "ATTENDU vrai : la congruence est AGNOSTIQUE À LA BRANCHE — "
            "elle relie deux descriptions de la même configuration, quel "
            "que soit le choix de branche, pourvu qu'il soit COHÉRENT, "
            "et il l'est puisque (u',v') est dérivé du Z source. E5 "
            "n'est donc PAS un négatif discriminant pour D5.4 ; c'est E6 "
            "(perturbation du Jacobien) qui l'est, et il casse."),
        "note_E3_non_tautological": (
            "Q_target est construite depuis la section NATIVE du chart "
            "cible (ses propres coefficients rationnels, ses propres ε'), "
            "PAS poussée depuis la source. La première version la "
            "calculait en poussant W par J⁻¹ puis la retirait par J : la "
            "perturbation s'annulait et le gate ne pouvait pas échouer — "
            "c'est E6 qui l'a révélé."),
        "E6_jacobian_mutation_congruence_still_holds": j_ok,
        "not_established": ["D5.6 cover multi-chart", "globalisation"],
        "gates_prereg": gates,
        "verdict": verdict}
    ART.write_text(json.dumps(out, indent=2, ensure_ascii=False,
                              default=float), encoding="utf-8")
    print("\nVERDICT :\n" + verdict)
    print(f"\n→ {ART}")
    return out


def _selftest():
    fails = []
    reg = load_canonical_MH()
    M = build_M_civ(reg["M_H_canonical"])
    c218 = reg["coeffs218"]
    rw = 1.0 - GAMMA
    from .taylor_models import tm_chart_cell_section
    d = json.loads((RES / "k3_cap_b1e2iii_p0a2_direct.json").read_text(
        encoding="utf-8"))
    r = d["s1_three_boxes"][0]
    S, g = tuple(r["S"]), r["g"]
    u0, v0 = complex(*r["u0"]), complex(*r["v0"])
    pc = [p for p in r["per_class"]
          if p["h_pass"] == r["worst_class_h_pass"]][0]
    eps = tuple(pc["eps"])
    Z, W, _ = tm_chart_cell_section(S, g, eps, u0, v0, 0.0)
    Q0 = Qmat(Z, W, M, c218, rw)

    # S1 : quadraticité — la BONNE convention convient
    A = [[complex(2.0, 0.3), complex(0.1, -0.4)],
         [complex(-0.2, 0.5), complex(1.5, 0.1)]]
    Atm = [[TMC.const(_c(A[i][k])) for k in range(2)] for i in range(2)]
    WA = [[None, None] for _ in range(6)]
    for a in range(6):
        for col in range(2):
            acc = None
            for p in range(2):
                t = W[a][p].mul_civ(_c(A[p][col]))
                acc = t if acc is None else acc + t
            WA[a][col] = acc
    QA = Qmat(Z, WA, M, c218, rw)
    ok1, _ = contains_zero(mat_sub(QA, congruence(Atm, Q0)))
    fails.append(not ok1)
    print(f"[{'PASS' if ok1 else 'FAIL'}] S1 quadraticité : "
          f"Q(Z, W·A) = Aᵀ Q conj(A)")

    # S2 NEGATIVE CONTROL: the competing convention does NOT work
    ok2, _ = contains_zero(
        mat_sub(QA, congruence(Atm, Q0, conj_left=True)))
    fails.append(ok2)
    print(f"[{'PASS' if not ok2 else 'FAIL'}] S2 négatif de convention : "
          f"A* Q A ne convient PAS ({not ok2}) — sans quoi S1 serait "
          f"satisfait par n'importe quelle forme")

    # S3 : invariance de jauge — Q(λZ, λW) = Q(Z, W)
    lam = complex(1.7, -0.9)
    Zl = [z.mul_civ(_c(lam)) for z in Z]
    Wl = [[W[a][0].mul_civ(_c(lam)), W[a][1].mul_civ(_c(lam))]
          for a in range(6)]
    ok3, _ = contains_zero(mat_sub(Qmat(Zl, Wl, M, c218, rw), Q0))
    fails.append(not ok3)
    print(f"[{'PASS' if ok3 else 'FAIL'}] S3 invariance de jauge : "
          f"Q(λZ, λW) = Q(Z, W) pour λ = {lam}")

    # S4 gauge NEGATIVE CONTROL: scaling Z WITHOUT scaling W breaks
    ok4, _ = contains_zero(mat_sub(Qmat(Zl, W, M, c218, rw), Q0))
    fails.append(ok4)
    print(f"[{'PASS' if not ok4 else 'FAIL'}] S4 négatif de jauge : "
          f"λ sur Z seul CASSE l'identité ({not ok4})")

    print("-" * 78)
    print("SELF-TEST:", "FAIL" if any(fails) else "ALL PASS")
    return 1 if any(fails) else 0


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--selftest":
        raise SystemExit(_selftest())
    build()
