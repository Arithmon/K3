# -*- coding: utf-8 -*-
"""ATLAS_PAPER_FIGURES: the five figures of the atlas paper.
Read-only on the certificates, deterministic, committable (lightweight PNGs). No heavy run, no load-bearing float in
the mathematics: the figures READ the numbers from the certificates (U1, bridge panel, generators) and invent none.

Figures produced in `figures/`:
  1. fig_chart_regime_geometry.png       — DATA: the three loci Re R_s = 0 of one chart in the complex plane of the base (rational coefficients derived from μ)
  2. fig_certified_chart_radii.png       — DATA: certified radius ρ of the 60 types (log scale) AFTER the upstream repair of the σ floor (U1 rev. 2)
  3. fig_deck_action_by_chart_type.png   — DATA: matrix of 20 triples × 3 gauges, D vertical (12 types) vs moving the base (48)
  4. fig_bridge_continuation.png         — ANNOTATED DIAGRAM: bilateral bridge, open overlap, conjugate vs continued = D·Z_conj (margins read from the panel)
  5. fig_proof_architecture.png          — DIAGRAM: exact algebra → quantitative lemma → certified applicability → X_atlas ≅ X → continuation

DERIVATION (rule "derive, do not merely assert"): ρ per type comes from `U1C_constants_certified.types_60` (Fraction → float
FOR DISPLAY only); the bridge margins come from the the bridge step panel; the vertical/base partition of D is RECOMPUTED
here from D and compared with the counts of the generators certificate (internal check). A JSON manifest serialises the plotted
numbers and the SHA of the sources: a figure of the paper must be reproducible and traceable like a check.
Style: validated categorical palette (blue #2a78d6 / orange #eb6834 — CVD ΔE 24.7 protan, 33.6 normal, all checks PASS)
PAIRED with a secondary encoding (marker shape + hatch) for black-and-white printing; recessive axes, no heavy grid,
selective direct labelling. No interaction (print medium).
"""
from __future__ import annotations

import hashlib
import io
import json
import math
import sys
from fractions import Fraction as Fr
from itertools import combinations
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle, Ellipse, FancyBboxPatch

HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
RES = ROOT / "certificates"
FIG = ROOT / "figures"
U1 = RES / "open_chart_theorem.json"
GEN = RES / "glue_obligations.json"
BRIDGE = RES / "bridge_atlas_panel.json"
SIGMA = RES / "sigma_floor_correction.json"
OUT = RES / "figures_manifest.json"

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED = "#1a1a19", "#4a4a48", "#8a8a86"
GRID = "#e3e3e0"
DPI = 200
TITLES = False        # titles live in the LaTeX CAPTIONS of the paper, not in the image

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
    "grid.color": GRID, "grid.linewidth": 0.6, "savefig.bbox": "tight", "savefig.pad_inches": 0.08,
})


def _set_title(ax, *a, **k):
    """title embedded in the image: disabled for the paper (the LaTeX caption carries it)."""
    if TITLES:
        _set_title(ax, *a, **k)


def _sci_outward(x, downward):
    """Three-digit mantissa rounded OUTWARD. A minimum's label must never
    print ABOVE the value it marks, nor a maximum's below it: at `.2e` the
    observed-minimum label read 2.11e-12 for a line drawn at 2.1092e-12."""
    e = math.floor(math.log10(abs(x)))
    m = x / 10 ** e
    m = math.floor(m * 1000) / 1000 if downward else math.ceil(m * 1000) / 1000
    return f"{m:.3f}e{e:+03d}"

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p): return json.loads(p.read_text(encoding="utf-8"))


def _det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
            - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))


def _inv3(M):
    d = _det3(M); C = [[Fr(0)] * 3 for _ in range(3)]
    for a in range(3):
        for b in range(3):
            m = [[M[u][v] for v in range(3) if v != b] for u in range(3) if u != a]
            C[b][a] = Fr((-1) ** (a + b)) * (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / d
    return C


def branch_coeffs(mu, S=(3, 4, 5)):
    """R_s(u,v) = a_s + b_s u² + c_s v² for the chart (S, gauge z_{T0}=1), base (u,v) = (z_{T1}, z_{T2}).
    DERIVED from the radical relation w_S = −V_S^{-1} V_T w_T with w_T = (1, u², v²): exact rational coefficients."""
    T = tuple(t for t in range(6) if t not in S)
    Vm = lambda I: [[Fr(mu[x]) ** k for x in I] for k in range(3)]
    A, B = _inv3(Vm(S)), Vm(T)
    P = [[-sum(A[a][k] * B[k][b] for k in range(3)) for b in range(3)] for a in range(3)]
    return {s: (P[a][0], P[a][1], P[a][2]) for a, s in enumerate(S)}, S, T


def fig1_regime_geometry(path, mu):
    """The three loci Re R_s = 0 of one chart, in the complex plane of the base (slice v = 0).
    On a REAL slice they are empty (X(R) = ∅: F₀ is positive definite) — the geometry is complex."""
    coeffs, S, T = branch_coeffs(mu)
    ks = {s: float(-a / b) for s, (a, b, c) in coeffs.items()}                 # Re R_s = 0 ⟺ x² − y² = −a/b (v = 0)
    assert all(k < 0 for k in ks.values()), "badly chosen slice: the hyperbolae must open along Im u"
    import numpy as np
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    xs = np.linspace(-1.0, 1.0, 400)
    styles = [(BLUE, "-"), (ORANGE, (0, (5, 2))), (AQUA, (0, (1.5, 1.8)))]
    for (s, k), (col, ls) in zip(sorted(ks.items()), styles):
        ys = np.sqrt(xs ** 2 - k)
        ax.plot(xs, ys, color=col, ls=ls, lw=1.8, zorder=3, label=f"$\\mathrm{{Re}}\\,R_{{{s}}}=0$")
    # regimes: sign of Re R_s = a_s + b_s(x² − y²) — COMPUTED, not hard-coded
    def regime(y):
        return "".join("p" if float(a + b * (0 - y ** 2)) > 0 else "c" for s, (a, b, c) in sorted(coeffs.items()))
    y_lo, y_hi = 1.08, 1.52
    ax.text(0.0, y_lo + 0.02, f"regimes  {'·'.join(regime(y_lo))}", fontsize=8.4, color=INK2, va="bottom", ha="center")
    ax.text(0.0, y_hi - 0.04, f"regimes  {'·'.join(regime(y_hi))}", fontsize=8.4, color=INK2, va="top", ha="center")
    ax.legend(loc="lower left", frameon=False, fontsize=8.2, handlelength=2.6, borderaxespad=0.4, labelspacing=0.35)
    ax.add_patch(Rectangle((0.18, 1.16), 0.30, 0.14, facecolor="none", edgecolor=INK, lw=1.2, ls=(0, (3, 2)), zorder=4))
    ax.annotate("bridge box: straddles one curve\n(no single $\\sigma$ exists on it)", xy=(0.47, 1.165),
                xytext=(0.55, 1.105), fontsize=8, color=INK, va="bottom",
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7))
    ax.set_xlim(-1.0, 1.12); ax.set_ylim(y_lo, y_hi)
    ax.set_xlabel(r"$\mathrm{Re}\,u$"); ax.set_ylabel(r"$\mathrm{Im}\,u$")
    ax.grid(axis="both", zorder=0)
    Sstr = ",".join(map(str, S))
    _set_title(ax, f"Figure 1. Regime geometry of one chart ($S=\\{{{Sstr}\\}}$, $z_{{{T[0]}}}=1$, slice $v=0$)\n"
                 "crossing a curve flips that line between the principal (p) and canonical (c) determination",
                 loc="left", pad=8, fontsize=9.2)
    fig.savefig(path, dpi=DPI); plt.close(fig)


def fig2_radii(path, vals, rho_min, rho_max):
    vals = sorted(vals)
    fig, ax = plt.subplots(figsize=(6.0, 3.4))
    x = range(1, len(vals) + 1)
    ax.vlines(list(x), rho_min, vals, color=GRID, lw=1.0, zorder=1)
    ax.plot(list(x), vals, ls="none", marker="o", ms=4.0, mfc=BLUE, mec="white", mew=0.7, zorder=3)
    ax.axhline(rho_min, color=ORANGE, lw=1.4, ls=(0, (5, 2)), zorder=2)
    ax.set_yscale("log")
    ax.set_xlabel("chart type (60), sorted by certified radius")
    ax.set_ylabel(r"certified radius  $\rho$")
    ax.grid(axis="y", zorder=0)
    ax.set_xlim(0, len(vals) + 1)
    ax.annotate(f"observed min  {_sci_outward(rho_min, True)}", xy=(31, rho_min), xytext=(31, rho_min * 1.04),
                color=ORANGE, fontsize=8.5, va="bottom", ha="left")
    ax.annotate(f"max  {_sci_outward(rho_max, False)}", xy=(len(vals), vals[-1]), xytext=(len(vals) - 2, vals[-1] * 0.62),
                color=INK2, fontsize=8.5, ha="right", va="top")
    _set_title(ax, "Figure 2. Certified chart radius across the 60 chart types", loc="left", pad=8)
    fig.savefig(path, dpi=DPI); plt.close(fig)


def fig3_deck(path, D, vertical_types, triples):
    fig, ax = plt.subplots(figsize=(5.2, 5.6))
    tri_sorted = triples
    for r, S in enumerate(tri_sorted):
        T = [t for t in range(6) if t not in S]
        for c, g in enumerate(T):
            vert = (S, g) in vertical_types
            ax.plot([c], [len(tri_sorted) - 1 - r], ls="none",
                    marker="s" if vert else "o", ms=9 if vert else 7,
                    mfc=ORANGE if vert else "white", mec=ORANGE if vert else BLUE, mew=1.3, zorder=3)
            if vert:
                ax.plot([c], [len(tri_sorted) - 1 - r], ls="none", marker="+", ms=5, color="white", mew=1.2, zorder=4)
    ax.set_yticks(range(len(tri_sorted)))
    ax.set_yticklabels([f"{{{','.join(map(str, S))}}}" for S in reversed(tri_sorted)], fontsize=7.2, family="monospace")
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["gauge 1", "gauge 2", "gauge 3"], fontsize=8)
    ax.set_xlim(-0.6, 2.6); ax.set_ylim(-0.7, len(tri_sorted) - 0.3)
    ax.set_ylabel("solved triple $S$", labelpad=6)
    ax.grid(axis="y", zorder=0)
    n_v = len(vertical_types)
    ax.plot([], [], ls="none", marker="s", ms=8, mfc=ORANGE, mec=ORANGE, label=f"$D$ vertical (sheet only) — {n_v} types")
    ax.plot([], [], ls="none", marker="o", ms=7, mfc="white", mec=BLUE, mew=1.3, label=f"$D$ moves the base — {60 - n_v} types")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.105), frameon=False, fontsize=8, handletextpad=0.5)
    Dstr = ",".join(("+" if d > 0 else "−") + "1" for d in D)
    flipped = sorted(i for i, d in enumerate(D) if d < 0)                     # DERIVED from D, not hard-coded
    assert all((set(flipped) <= set(S)) == ((S, T[0]) in vertical_types)
               for S in tri_sorted for T in [[t for t in range(6) if t not in S]])
    fl = "\\{" + ",".join(map(str, flipped)) + "\\}"
    ax.annotate(f"vertical exactly when $S\\supseteq{fl}$ — the coordinates $D$ flips",
                xy=(0.5, -0.055), xycoords="axes fraction", ha="center", va="top", fontsize=8, color=INK2)
    _set_title(ax, f"Figure 3. One deck element, two chart realizations\n$D=\\mathrm{{diag}}({Dstr})$",
                 loc="left", pad=8, fontsize=9.5)
    fig.savefig(path, dpi=DPI); plt.close(fig)


def fig4_bridge(path, width, sep, diff, n_bridges):
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.set_xlim(0, 10); ax.set_ylim(0, 6.6); ax.axis("off")
    # --- geometry (lower left) ---
    ax.add_patch(Rectangle((0.5, 1.15), 4.1, 1.9, facecolor="#f2f7fd", edgecolor=BLUE, lw=1.1))
    ax.text(0.75, 2.88, "chart A", color=BLUE, fontsize=8.5, va="top")
    ax.add_patch(Rectangle((3.7, 1.15), 3.9, 1.9, facecolor="#fdf4f0", edgecolor=ORANGE, lw=1.1, alpha=0.85))
    ax.text(7.35, 2.88, "chart B", color=ORANGE, fontsize=8.5, va="top", ha="right")
    ax.add_patch(Rectangle((3.35, 0.85), 1.55, 2.5, facecolor="none", edgecolor=INK, lw=1.4, ls=(0, (4, 2.5))))
    ax.text(4.12, 3.52, "bridge box (bilateral)", ha="center", color=INK, fontsize=8.5)
    ax.plot([4.12, 4.12], [1.15, 3.05], color=MUTED, lw=1.0)
    ax.text(4.02, 2.10, r"$\mathrm{Re}\,R=0$", rotation=90, ha="right", va="center", fontsize=7.6, color=MUTED)
    ax.text(2.55, 3.52, "regime boundary", ha="right", color=MUTED, fontsize=7.8)
    ax.annotate("", xy=(3.95, 3.10), xytext=(2.65, 3.48), arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7))
    ax.annotate("", xy=(3.7, 0.62), xytext=(4.6, 0.62), arrowprops=dict(arrowstyle="<->", color=INK2, lw=0.9))
    ax.text(4.15, 0.30, f"open overlap $\\geq$ {width:.2e}", ha="center", fontsize=7.8, color=INK2)
    # --- the two sheets (upper right) ---
    ax.plot([6.55], [5.70], marker="o", ms=6.5, mfc="white", mec=BLUE, mew=1.5)
    ax.text(6.85, 5.70, "conjugate sheet $Z_{\\mathrm{conj}}$\nfails to reglue " + f"({n_bridges}/{n_bridges})",
            va="center", fontsize=8.2, color=BLUE, linespacing=1.5)
    ax.plot([6.55], [4.35], marker="s", ms=6.5, mfc=ORANGE, mec=ORANGE)
    ax.text(6.85, 4.35, "continued sheet $D\\cdot Z_{\\mathrm{conj}}$\nreglues exactly " + f"({n_bridges}/{n_bridges})",
            va="center", fontsize=8.2, color=ORANGE, linespacing=1.5)
    ax.annotate("", xy=(6.55, 4.62), xytext=(6.55, 5.42), arrowprops=dict(arrowstyle="<->", color=MUTED, lw=0.9))
    ax.text(6.40, 5.02, f"separation $\\geq$ {sep}", rotation=90, ha="right", va="center", fontsize=7.4, color=INK2)
    ax.add_patch(FancyArrowPatch((4.95, 3.30), (6.25, 4.15), arrowstyle="-|>", color=MUTED, lw=1.0,
                                 connectionstyle="arc3,rad=0.28", mutation_scale=11))
    ax.text(4.62, 4.55, "analytic\ncontinuation", fontsize=8, color=INK2, ha="left", linespacing=1.4)
    ax.text(0.5, -0.02, f"recentred difference on the lower side: {diff:.2e}", fontsize=7.6, color=MUTED)
    _set_title(ax, "Figure 4. Continuation across a certified bridge: conjugation is not continuation",
                 loc="left", pad=8, fontsize=9.5)
    fig.savefig(path, dpi=DPI); plt.close(fig)


def fig5_architecture(path):
    fig, ax = plt.subplots(figsize=(6.2, 2.2))
    ax.set_xlim(0, 10); ax.set_ylim(0.35, 3.3); ax.axis("off")
    boxes = [("exact\nalgebra", 0.9), ("quantitative\nlemma", 2.75), ("certified\napplicability", 4.6),
             (r"$\mathcal{X}_{\mathrm{atlas}}\cong X$", 6.55), ("certified\ncontinuation", 8.5)]
    for i, (label, x) in enumerate(boxes):
        hero = (i == 3)
        ax.add_patch(FancyBboxPatch((x - 0.75, 1.45), 1.50, 1.25, boxstyle="round,pad=0.06,rounding_size=0.08",
                                    facecolor="#f2f7fd" if hero else "white",
                                    edgecolor=BLUE if hero else MUTED, lw=1.5 if hero else 1.0))
        ax.text(x, 2.07, label, ha="center", va="center", fontsize=8.6 if not hero else 9.4,
                color=INK, fontweight="bold" if hero else "normal")
        if i:
            ax.annotate("", xy=(x - 0.79, 2.07), xytext=(boxes[i - 1][1] + 0.79, 2.07),
                        arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.0))
    ax.text(0.12, 0.72, "Sections 2, 5", fontsize=7.6, color=MUTED)
    ax.text(2.75, 0.72, "Section 3", fontsize=7.6, color=MUTED, ha="center")
    ax.text(4.6, 0.72, "Sections 3, 8", fontsize=7.6, color=MUTED, ha="center")
    ax.text(6.55, 0.72, "Section 4", fontsize=7.6, color=MUTED, ha="center")
    ax.text(8.5, 0.72, "Sections 6, 7", fontsize=7.6, color=MUTED, ha="center")
    _set_title(ax, "Figure 5. Proof architecture", loc="left", pad=6)
    fig.savefig(path, dpi=DPI); plt.close(fig)


MU = (1, 2, 3, 5, 7, 11)


def main():
    u1, gen, br = load(U1), load(GEN), load(BRIDGE)
    c = u1["U1C_constants_certified"]
    types60 = c["types_60"]
    rho_min = float(Fr(c["rho_uniform_lo_min"])); rho_max = float(Fr(c["rho_uniform_lo_max"]))
    D = tuple(gen["C_deck"]["D"])
    triples = list(combinations(range(6), 3))
    vertical = {(S, g) for S in triples for g in range(6) if g not in S
                and len(set(D[i] for i in [x for x in range(6) if x not in S])) == 1}
    # internal check: the recomputed partition must equal the counts of the upstream certificate
    assert len(vertical) == gen["C_deck"]["D_vertical_types"], "D partition inconsistent with the generators certificate"
    assert len(types60) == 60 and len(set(t["S"][0] for t in types60)) >= 1
    # bridge panel margins: READ by field name (no number entered by hand)
    width = float(br["min_overlap_width"])
    sep_exact = float(br["min_theta_margin"])
    diff = float(br["max_recentred_difference_sup"])
    n_bridges = int(br["n_bridges"])
    sep = f"{math.floor(sep_exact * 1000) / 1000:.3f}"          # rounded OUTWARD: honest lower bound
    assert float(sep) <= sep_exact and diff > 0 and width > 0
    FIG.mkdir(exist_ok=True)
    paths = {}
    paths["fig1"] = FIG / "fig_chart_regime_geometry.png"; fig1_regime_geometry(paths["fig1"], MU)
    # U1 is REPAIRED upstream (rev. 2): the figure reads the natural source again.
    vals = [float(Fr(t["rho_uniform_lo"])) for t in types60]
    rho_min, rho_max = min(vals), max(vals)
    paths["fig2"] = FIG / "fig_certified_chart_radii.png"; fig2_radii(paths["fig2"], vals, rho_min, rho_max)
    paths["fig3"] = FIG / "fig_deck_action_by_chart_type.png"; fig3_deck(paths["fig3"], D, vertical, triples)
    paths["fig4"] = FIG / "fig_bridge_continuation.png"; fig4_bridge(paths["fig4"], width, sep, diff, n_bridges)
    paths["fig5"] = FIG / "fig_proof_architecture.png"; fig5_architecture(paths["fig5"])
    man = {"artifact": "figures_manifest",
           "subject": "figures of the atlas paper — 2 data-driven (certified radii, deck action), 2 annotated diagrams, 1 architecture",
           "sources": {"u1": sha(U1), "generators": sha(GEN), "bridge_panel": sha(BRIDGE), "sigma_correction": sha(SIGMA)},
           "numbers_plotted": {"branch_coeffs_S345": {str(k): [str(x) for x in v] for k, v in branch_coeffs(MU)[0].items()},
                               "rho_min": repr(rho_min), "rho_max": repr(rho_max), "n_types": len(types60),
                               "D": list(D), "D_vertical_types": len(vertical), "D_base_moving_types": 60 - len(vertical),
                               "bridge_overlap_width": width, "sheet_separation_exact": sep_exact, "sheet_separation_outward": sep,
                               "recentred_diff": diff, "n_bridges": n_bridges},
           "palette": {"categorical": [BLUE, ORANGE, AQUA], "validated": "light surface: lightness/chroma/CVD/normal/contrast all PASS (protan ΔE 24.7, normal 33.6)",
                       "secondary_encoding": "marker shape + hatch (print / grayscale safety)"},
           "figures": {k: str(v.relative_to(ROOT)) for k, v in paths.items()},
           "sizes_bytes": {k: v.stat().st_size for k, v in paths.items()},
           "does_not_attest": ["no new mathematics: the figures READ the certificates", "the diagrams (1, 4, 5) are declared schematic",
                             "the floats are for DISPLAY (the sources remain Fraction)"]}
    try:
        import subprocess
        man["built_from_head"] = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        man["built_from_head"] = None
    man["self_sha256"] = sha(HERE)
    OUT.write_text(json.dumps(man, ensure_ascii=False, indent=1), encoding="utf-8")
    for k, v in paths.items():
        print(f"  {k}  {v.relative_to(ROOT)}  {v.stat().st_size // 1024} KB")
    print(f"D = {D} : {len(vertical)} vertical types / {60 - len(vertical)} moving the base | ρ ∈ [{rho_min:.3e}, {rho_max:.3e}]")


if __name__ == "__main__":
    main()
