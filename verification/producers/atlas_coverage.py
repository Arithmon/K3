#!/usr/bin/env python3
"""
Coverage certificate for the chart atlas: uniform pivot floor.

Certifies   inf_{X}  max_{i<j<k} |Z_i Z_j Z_k| * |V_ijk|  >=  TAU,
where V_ijk = (mu_j-mu_i)(mu_k-mu_i)(mu_k-mu_j) and mu = (1,2,3,5,7,11).

For the Vandermonde CI(2,2,2), the 3x3 minor on columns (i,j,k) of the
complex Jacobian M = [2 mu_j^{k-1} Z_j] factors EXACTLY:
|det M_S| = 8 |Z_i Z_j Z_k| |V_ijk|. The certificate therefore gives a
uniform floor on the largest minor: every point of the surface admits a pivot
above the threshold. This is the quantitative coverage criterion of the atlas
("min over X of the max minor >= tau").

Method: an AMBIENT exclusion branch-and-bound — no Krawczyk operator, no
sheet. Per gauge g in {0..5}: the slice { Z_g = t real, t in [1/sqrt6, 1],
|Z_j| <= t } (11 real dimensions), reduced by the sign group
H_g = { +-eps : eps in Z_2^3 } x {conj} (16 elements, 4 half-spaces; the
over-cover is valid, since covering too much is not wrong). Each box is:
  - EXCLUDED if one interval residual excludes 0:
      Re Q_k, Im Q_k (k=1..3), the sphere Sum|Z|^2 - 1, or |Z_j|^2 > t^2;
  - CERTIFIED if the interval lower bound on P^2 (through squared
    mignitudes) exceeds TAU^2;
  - otherwise bisected along its widest dimension.
Arithmetic: float64 intervals, outward-rounded by one ulp after every
operation (np.nextafter; IEEE round-to-nearest, hence a rigorous enclosure),
vectorised over batches of boxes.

Usage: atlas_coverage.py [TAU=0.5] [N_CORES=4] [W_MIN=1e-4] [BOX_BUDGET_M=200]

THE BOX BUDGET IS NOT OPTIONAL. At the default (200 M) gauge 1 stops early on
`budget_abort` and the verdict falls to FAILED on 16 unresolved boxes, while
the other five gauges reproduce identically. The published run uses 1000, and
the one-command verifier passes it explicitly.

Output: certificates/atlas_coverage.json
Verdict CERTIFIED iff no box is left unresolved on any of the six gauges.
"""
from __future__ import annotations
import io, json, sys, time
from fractions import Fraction
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path
import numpy as np

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "certificates"

TAU = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
N_CORES = int(sys.argv[2]) if len(sys.argv) > 2 else 4
W_MIN = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-4
BOX_BUDGET = int(float(sys.argv[4]) * 1e6) if len(sys.argv) > 4 else 200_000_000

MU = [1, 2, 3, 5, 7, 11]
LAMBDA_INT = [[m ** k for m in MU] for k in range(3)]        # lignes k=0,1,2
COORD_CHARS = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 1, 0], [0, 0, 1], [0, 0, 1]]
TRIPLES = list(combinations(range(6), 3))
V2 = [float(((MU[j] - MU[i]) * (MU[k] - MU[i]) * (MU[k] - MU[j])) ** 2)
      for i, j, k in TRIPLES]                                # exact en float64
CHUNK = 120_000          # boxes per vectorised batch
PRESPLIT_DIMS = 5        # 2^5 = 32 sous-racines par jauge

NEG, POS = -np.inf, np.inf


def dn(a):
    return np.nextafter(a, NEG)


def up(a):
    return np.nextafter(a, POS)


# ---------------------------------------------------------------------------
#  Sign group on the gauge-g slice: half-spaces (over-cover, 16 elements)
# ---------------------------------------------------------------------------
def sign_group_cuts(g):
    """Real coordinates c to restrict to x_c >= 0.

    H_g = { s·eps : eps in Z_2^3, s in {+-1} t.q. (s·eps)_g = +1 } x {id, conj},
    sign patterns on the 12 real coordinates (X_j, Y_j are coordinates 2j and
    2j+1; conjugation flips every Y). Greedy reduction: pick a flipped
    coordinate, cut x_c >= 0, keep the stabiliser; every orbit retains a
    representative in the intersection of the half-spaces (induction on the chain
    de sous-groupes) => over-cover valide.
    """
    pats = []
    for c0 in range(2):
        for c1 in range(2):
            for c2 in range(2):
                s6 = [(-1) ** (COORD_CHARS[j][0] * c0
                               + COORD_CHARS[j][1] * c1
                               + COORD_CHARS[j][2] * c2) for j in range(6)]
                if s6[g] < 0:
                    s6 = [-v for v in s6]
                for conj in (1, -1):
                    p12 = []
                    for j in range(6):
                        p12.append(s6[j])
                        p12.append(s6[j] * conj)
                    pats.append(tuple(p12))
    pats = sorted(set(pats))
    pinned = 2 * g + 1                       # Y_g = 0 on the slice
    cuts = []
    group = pats
    while True:
        cand = None
        for p in group:
            for c in range(12):
                if c != pinned and p[c] < 0:
                    cand = c
                    break
            if cand is not None:
                break
        if cand is None:
            break
        cuts.append(cand)
        group = [p for p in group if p[cand] > 0]
    return cuts


# ---------------------------------------------------------------------------
#  Vectorised kernel: (L, H) of shape (K, 12) -> excluded / certified masks
# ---------------------------------------------------------------------------
def process_batch(L, H, g, tau2_up):
    """Evaluate a batch of boxes. Returns (excl, cert, excl_counts)."""
    K = L.shape[0]
    XL, XH = L[:, 0::2], H[:, 0::2]
    YL, YH = L[:, 1::2], H[:, 1::2]

    # interval squares (K,6)
    posX = XL > 0
    negX = XH < 0
    sqXL = np.where(posX, dn(XL * XL), np.where(negX, dn(XH * XH), 0.0))
    sqXH = up(np.maximum(XL * XL, XH * XH))
    posY = YL > 0
    negY = YH < 0
    sqYL = np.where(posY, dn(YL * YL), np.where(negY, dn(YH * YH), 0.0))
    sqYH = up(np.maximum(YL * YL, YH * YH))

    # squared moduli m_j = X_j^2 + Y_j^2 (K,6)
    mL = dn(sqXL + sqYL)
    mH = up(sqXH + sqYH)

    excl = np.zeros(K, dtype=bool)
    counts = {}

    # sphere: Sum m_j - 1
    sL = mL[:, 0].copy()
    sH = mH[:, 0].copy()
    for j in range(1, 6):
        sL = dn(sL + mL[:, j])
        sH = up(sH + mH[:, j])
    e = (sL > 1.0) | (sH < 1.0)
    counts['sphere'] = int(e.sum())
    excl |= e

    # modulus : |Z_j|^2 > t^2  (t = X_g > 0)
    t2H = sqXH[:, g]
    e = np.zeros(K, dtype=bool)
    for j in range(6):
        if j == g:
            continue
        e |= mL[:, j] > t2H
    counts['modulus'] = int((e & ~excl).sum())
    excl |= e

    # quadriques : Re Q_k = Sum lam (X^2 - Y^2), Im Q_k = 2 Sum lam X Y
    dL = dn(sqXL - sqYH)
    dH = up(sqXH - sqYL)
    p1, p2, p3, p4 = XL * YL, XL * YH, XH * YL, XH * YH
    pL = dn(np.minimum(np.minimum(p1, p2), np.minimum(p3, p4)))
    pH = up(np.maximum(np.maximum(p1, p2), np.maximum(p3, p4)))
    nq_re = nq_im = 0
    for k in range(3):
        lam = LAMBDA_INT[k]
        rL = dn(lam[0] * dL[:, 0])
        rH = up(lam[0] * dH[:, 0])
        iL = dn(lam[0] * pL[:, 0])
        iH = up(lam[0] * pH[:, 0])
        for j in range(1, 6):
            rL = dn(rL + dn(lam[j] * dL[:, j]))
            rH = up(rH + up(lam[j] * dH[:, j]))
            iL = dn(iL + dn(lam[j] * pL[:, j]))
            iH = up(iH + up(lam[j] * pH[:, j]))
        e = (rL > 0.0) | (rH < 0.0)
        nq_re += int((e & ~excl).sum())
        excl |= e
        e = (iL > 0.0) | (iH < 0.0)          # the factor 2 does not affect the sign
        nq_im += int((e & ~excl).sum())
        excl |= e
    counts['quad_re'] = nq_re
    counts['quad_im'] = nq_im

    # certification : P^2_lo = max_triples V^2 * mig2_i * mig2_j * mig2_k
    mig2 = np.maximum(mL, 0.0)
    p2max = np.zeros(K)
    for t, (i, j, k) in enumerate(TRIPLES):
        q = dn(dn(dn(mig2[:, i] * mig2[:, j]) * mig2[:, k]) * V2[t])
        np.maximum(p2max, q, out=p2max)
    cert = p2max >= tau2_up
    return excl, cert & ~excl, counts


def bnb_task(args):
    """Depth-first branch and bound on one subroot. Returns stats and unresolved boxes."""
    g, L0, H0, tau, w_min, budget = args
    tau2_up = up(np.float64(tau) * np.float64(tau))
    pinned = 2 * g + 1
    stack = [(L0[None, :].copy(), H0[None, :].copy())]
    n_proc = n_excl = n_cert = 0
    counts_tot = {}
    unresolved = []
    t0 = time.time()
    while stack:
        L, H = stack.pop()
        if L.shape[0] > CHUNK:
            stack.append((L[CHUNK:], H[CHUNK:]))
            L, H = L[:CHUNK], H[:CHUNK]
        n_proc += L.shape[0]
        excl, cert, counts = process_batch(L, H, g, tau2_up)
        for key, v in counts.items():
            counts_tot[key] = counts_tot.get(key, 0) + v
        n_excl += int(excl.sum())
        n_cert += int(cert.sum())
        undec = ~(excl | cert)
        if not undec.any():
            continue
        L, H = L[undec], H[undec]
        w = H - L
        w[:, pinned] = -1.0
        too_small = w.max(axis=1) < w_min
        if too_small.any():
            for idx in np.where(too_small)[0][:64]:
                unresolved.append((L[idx].tolist(), H[idx].tolist()))
            keep = ~too_small
            L, H = L[keep], H[keep]
            w = w[keep]
            if L.shape[0] == 0:
                continue
        if n_proc > budget:
            unresolved.append(('budget_abort', int(L.shape[0])))
            break
        dims = w.argmax(axis=1)
        mid = 0.5 * (L[np.arange(L.shape[0]), dims]
                     + H[np.arange(L.shape[0]), dims])
        L1, H1 = L.copy(), H.copy()
        H1[np.arange(L.shape[0]), dims] = mid
        L2, H2 = L.copy(), H.copy()
        L2[np.arange(L.shape[0]), dims] = mid
        stack.append((np.concatenate([L1, L2]), np.concatenate([H1, H2])))
    return {'g': g, 'n_proc': n_proc, 'n_excl': n_excl, 'n_cert': n_cert,
            'counts': counts_tot, 'unresolved': unresolved,
            'seconds': time.time() - t0}


# ---------------------------------------------------------------------------
#  Roots per gauge, plus the presplit
# ---------------------------------------------------------------------------
def gauge_roots(g):
    # A PROVEN bound, not an assumed one. One `nextafter` is NOT enough: sqrt
    # and division compose two roundings, and the result stays ABOVE 1/sqrt(6)
    # by about 8.6e-19 (checked in exact rationals: t^2 - 1/6 = +7.05e-19).
    # The slice 1/sqrt(6) <= t < t_lo was therefore formally outside the
    # enumerated domain, which is enough to cost the certificate the word
    # "exhaustive". Step down until t_lo^2 <= 1/6 holds EXACTLY in rational
    # arithmetic — two steps here, but the loop relies on no platform property.
    t_lo = 1.0 / np.sqrt(6.0)
    while Fraction(float(t_lo)) ** 2 > Fraction(1, 6):
        t_lo = np.nextafter(t_lo, NEG)               # t_lo^2 <= 1/6, proven
    L = -np.ones(12)
    H = np.ones(12)
    L[2 * g], H[2 * g] = t_lo, 1.0                   # X_g = t
    L[2 * g + 1] = H[2 * g + 1] = 0.0                # Y_g = 0
    for c in sign_group_cuts(g):
        L[c] = 0.0
    roots = [(L, H)]
    for _ in range(PRESPLIT_DIMS):
        new = []
        for (l, h) in roots:
            w = h - l
            w[2 * g + 1] = -1.0
            d = int(np.argmax(w))
            m = 0.5 * (l[d] + h[d])
            l1, h1 = l.copy(), h.copy()
            h1[d] = m
            l2, h2 = l.copy(), h.copy()
            l2[d] = m
            new += [(l1, h1), (l2, h2)]
        roots = new
    return roots


def main():
    print("=" * 72)
    print(f"Coverage certificate — uniform pivot floor  TAU={TAU}  "
          f"W_MIN={W_MIN}  cores={N_CORES}")
    print("=" * 72)
    tasks = []
    for g in range(6):
        cuts = sign_group_cuts(g)
        roots = gauge_roots(g)
        print(f"[g={g}] half-space cuts x_c>=0: {cuts}  "
              f"subroots: {len(roots)}")
        budget = BOX_BUDGET // (6 * len(roots))
        for (l, h) in roots:
            tasks.append((g, l, h, TAU, W_MIN, budget))
    t0 = time.time()
    results = []
    with Pool(N_CORES) as pool:
        for i, r in enumerate(pool.imap_unordered(bnb_task, tasks)):
            results.append(r)
            if (i + 1) % 24 == 0 or (i + 1) == len(tasks):
                np_ = sum(x['n_proc'] for x in results)
                nu = sum(len(x['unresolved']) for x in results)
                print(f"  [{i+1}/{len(tasks)}] boxes processed {np_:.3e}  "
                      f"non-résolues {nu}  ({time.time()-t0:.0f}s)")
    per_gauge = {}
    for g in range(6):
        rs = [r for r in results if r['g'] == g]
        per_gauge[g] = {
            'n_proc': sum(r['n_proc'] for r in rs),
            'n_excl': sum(r['n_excl'] for r in rs),
            'n_cert': sum(r['n_cert'] for r in rs),
            'n_unresolved': sum(len(r['unresolved']) for r in rs),
            'counts': {k: sum(r['counts'].get(k, 0) for r in rs)
                       for k in ['sphere', 'modulus', 'quad_re', 'quad_im']},
        }
    n_unres = sum(v['n_unresolved'] for v in per_gauge.values())
    verdict = 'CERTIFIED' if n_unres == 0 else 'FAILED'
    elapsed = time.time() - t0
    out = {
        'phase': 'K3 CAP D3-atlas — coverage certificate (pivot floor)',
        'date': '2026-07-10',
        'statement': ('inf over gauge-fixed CI(2,2,2) sphere variety of '
                      'max_{i<j<k} |Z_i Z_j Z_k| * |V_ijk|  >=  TAU ; '
                      'equivalently max 3x3 complex Jacobian minor >= 8*TAU'),
        'tau': TAU,
        'minor_floor_8tau': 8 * TAU,
        'verdict': verdict,
        'w_min': W_MIN,
        'arithmetic': 'float64 intervals, 1-ulp outward rounding (np.nextafter) per op',
        'domain': ('per gauge g: slice Z_g = t in [1/sqrt6, 1] real, '
                   'Y_g = 0, |Z_j| <= t ; sign-group reduction x16 '
                   '(Z_2^3 x {+-1} x conj), 4 half-space cuts, over-cover'),
        'mu': MU,
        'design_min_ref': 'k3_cap_d3_atlas_pivot_min_design.json',
        'per_gauge': per_gauge,
        'total_boxes': sum(v['n_proc'] for v in per_gauge.values()),
        'seconds': elapsed,
        'n_cores': N_CORES,
        'unresolved_sample': [u for r in results for u in r['unresolved']][:16],
    }
    RES.mkdir(parents=True, exist_ok=True)   # without `parents`, the producer dies AFTER the computation
    path = RES / "atlas_coverage.json"
    path.write_text(json.dumps(out, indent=1))
    print("-" * 72)
    print(f"VERDICT: {verdict}  (unresolved {n_unres})")
    print(f"total boxes {out['total_boxes']:.3e}  in {elapsed:.0f}s")
    print(f"-> {path}")
    return 0 if verdict == 'CERTIFIED' else 1


if __name__ == '__main__':
    sys.exit(main())
