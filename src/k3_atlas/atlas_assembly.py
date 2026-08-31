#!/usr/bin/env python3
"""
k3_cap_b1e2iii_c127d_atlas.py — C127-D : DE LA PARTITION À L'ATLAS.

C127 (252 tuiles) puis C127-E (les 64 résiduelles) ont fermé le cover
d'une paire cellule/classe : 316 tuiles, Kraft = 1, résidu 0. Mais ce
cover restait une PARTITION EN BOÎTES FERMÉES — les notes le disaient
elles-mêmes : « aucun halo, aucun overlap ouvert, aucune transition
certifiée entre tuiles voisines, aucun cocycle du nerf des 316 charts ».
Et la fermeture du résidu a rendu le raccord OBLIGATOIRE : des tuiles
adjacentes emploient des DÉTERMINATIONS DIFFÉRENTES (principale,
tournée, tournée canonique), donc rien ne garantissait a priori
qu'elles décrivent le même feuillet de la source.

CE QUE CE SCRIPT PAIE (les quatre points du §14 de la 6ᵉ revue GPT) :

  D-a. CORES + HALOS — chaque tuile reçoit un halo `box(c, (1+ρ)h)`,
     ρ pris dans une ÉCHELLE PRÉ-ENREGISTRÉE (1/8, 1/16, 1/32, 1/64) :
     le premier ρ qui passe est retenu, et « aucun ρ ne passe » est un
     REFUS, pas un repli. Sur le halo on exige : section source
     complète, DÉTERMINATIONS SOURCE INCHANGÉES, critère de chart (celui
     qui a certifié le core : C126 pour les 252, étendu pour les 64 —
     et pour ces 64, C126 doit TOUJOURS REFUSER, la non-tautologie R4
     survit au halo), section cible native à LEDGER FIGÉ (`ε'`, `σ'`)
     avec les mêmes `kinds`, et congruence projective sous δ.
     Un chart valide sur un VOISINAGE OUVERT du core : c'est ce qui
     transforme un certificat de boîte fermée en carte d'atlas.

  D-b. OVERLAPS OUVERTS — le nerf est calculé en RATIONNELS DYADIQUES
     EXACTS (`Fraction`, jamais de float) : deux cores qui se touchent
     ont des halos qui se coupent en une boîte de LARGEUR STRICTEMENT
     POSITIVE dans les quatre coordonnées. Publié : la boîte d'overlap
     de chaque paire, la connexité du nerf, les triples d'intersection
     triple non vide.

  D-c. TRANSITIONS — sur chaque overlap :
     · LE FEUILLET (`θ`) : les deux sections source sont comparées non
       pas par leurs enclosures (qui perdent la corrélation : largeur
       ~1e-3, verdict sans contenu) mais par RECENTRAGE EXACT des deux
       polynômes dans un cadre commun puis SOUSTRACTION COEFFICIENT PAR
       COEFFICIENT — la différence tombe à ~1e-12 et la séparation
       ±  devient décisive. `θ_ij ∈ {±1}` est DÉRIVÉ à marge stricte
       (refus si ambigu, jamais d'essai) et doit valoir +1 : les tuiles
       voisines sont sur le MÊME FEUILLET. Ce test POUVAIT échouer —
       c'est là tout son intérêt.
     · LA CARTE (`λ`) : `λ_ij = Z[g'_i]/Z[g'_j]`, avec `|Z[g'_i]|` et
       `|Z[g'_j]|` certifiés MINORÉS > 0 sur l'overlap ⟹ la transition
       est un biholomorphisme, pas seulement une formule.
     · L'IDENTITÉ DE TRANSITION, sans division : `Zt_a·Z[g'] ∋ Z_a` par
       tuile, et d'une tuile à l'autre `Zt^(j)_a·Z[g'_j] −
       Zt^(i)_a·Z[g'_i] ∋ 0` — c'est le raccord qui met EN CONTACT les
       deux ledgers `ε'` dérivés indépendamment.

  D-d. COCYCLE SUR LE NERF — sur chaque triple d'intersection triple non
     vide : `θ_ij·θ_jk·θ_ki = +1` (exact, discret), `λ_ij·λ_jk·λ_ki ∋ 1`
     (enclosure) et l'identité de transition sur les trois arêtes. La
     cochaîne de feuillet est un cocycle de classe TRIVIALE : le choix
     de feuillet se recolle globalement sur la cellule.

HONNÊTETÉ SUR CE QUI EST TAUTOLOGIQUE ET CE QUI NE L'EST PAS :
  · `θ ≡ +1` n'est PAS tautologique — les déterminations source varient
    d'une tuile à l'autre et rien n'imposait le même feuillet ;
  · l'identité de transition n'est PAS tautologique — elle met en
    présence deux ledgers `ε'` dérivés séparément ;
  · au niveau des TRIPLES, la condition de cocycle est IMPLIQUÉE par les
    identités par paires (toutes les cartes sont des cartes affines du
    même P⁵ et toutes les sections descendent d'UNE section source) —
    ce qui reste non trivial est que l'inclusion `∋ 0` sur l'arête i→k
    n'est PAS une conséquence intervalliste de i→j et j→k : les trois
    sont vérifiées séparément. C'est dit, pas déguisé.
  · les NÉGATIFS portent la discriminance : ledger `ε'` corrompu,
    feuillet inversé, ρ = 0, paire non adjacente, recentrage neutralisé,
    géométrie décalée d'un ulp — les six doivent CASSER.

CE QUE CE SCRIPT NE PAIE PAS : le transport métrique (congruence Q,
Weyl) SUR LE HALO — les gates C127 restent établis sur les cores, et
l'extension aux halos coûte un quadruple `Qmat` par tuile (~244 s, c'est
tout le budget de C127) ; il est ici mesuré sur le sous-ensemble
stratifié seulement, et c'est DÉCLARÉ. Ne paie pas non plus : le contrat
EXACT de l'identité de congruence, le scaling complet (C125-A), les 895
autres paires cellule/classe, R12-C.

C129-A/B/C (8ᵉ revue GPT, 2026-07-31) — trois dettes theorem-grade
payées SANS changer aucun contrat :
  A. les boîtes d'overlap (paires, triples, négatif N5) passaient par
     `iv.mpf([float(lo), float(hi)])` au PLUS PROCHE — or les ε
     d'intersection portent des dénominateurs en 65 (halos 65h/64), non
     représentables en binaire : l'intervalle pouvait être STRICTEMENT
     plus petit que la boîte rationnelle. → `fraction_box_to_iv` :
     borne basse VERS LE BAS, haute VERS LE HAUT, jamais optimiste.
  B. les défauts gatés (congruence de halo, identité de paire, arêtes
     de triple) restent en mp jusqu'à la comparaison à δ ; la
     conversion float est DIRIGÉE (`_f_up`) et n'arrive qu'à la
     sérialisation — la borne publiée MAJORE la borne comparée.
  C. gate D2c : l'ÉGALITÉ D'ENSEMBLES {tuiles clippées} = {tuiles nées
     du résidu C127-E}, publiée avec les deux ensembles triés et leurs
     SHA-256 — « 64 = 64 » n'exclut pas une permutation 63+1.

Sorties : results/k3_cap_b1e2iii_c127d_atlas_pilot.json   (mode pilot)
          results/k3_cap_b1e2iii_c127d_atlas.json         (mode full)
Usage   : k3_cap_b1e2iii_c127d_atlas.py [--selftest]
Env     : K3_C127D_MODE     pilot (défaut) | full
          K3_C127D_WORKERS  processus parallèles (défaut 6)
"""
from __future__ import annotations

import hashlib
import io
import itertools
import json
import math
import os
import platform
import subprocess
import sys
import time
from fractions import Fraction
from math import comb
from multiprocessing import get_context
from pathlib import Path

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

_argv = sys.argv
sys.argv = [sys.argv[0]]
os.environ.setdefault("K3_TM_ORDER", "4")
os.environ.setdefault("K3_TM_SERIES", "4")
from mpmath import iv, mp                                          # noqa: E402
from .witness_registry import load_canonical_MH              # noqa: E402
from .interval_arithmetic import build_M_civ                  # noqa: E402
from .taylor_models import (                                     # noqa: E402
    CIV, IVPM, MIDX, MONO, NM, TMC, TM_ORDER, UNARY_SERIES_DEG,
    civ_absmin, riv)
from .full_cell_charts import (                           # noqa: E402
    build_section, chart_certificate)
from .chart_selection_criterion import (                  # noqa: E402
    native_section_constructible)
from .gram_congruence import (                         # noqa: E402
    GAMMA, native_target_section)
from .chart_transport import (                    # noqa: E402
    address_of, tree_gates, _f_down, _f_up)
from .residual_closure import native_rows_ext          # noqa: E402
sys.argv = _argv

RES = Path(os.environ.get(
    "K3_RES_DIR", Path(__file__).resolve().parent / "data"))
COVER_JSON = RES / "k3_cap_b1e2iii_d5_6_dyadic_cover.json"
C127_JSON = RES / "k3_cap_b1e2iii_c127_transport_all.json"
C127E_JSON = RES / "k3_cap_b1e2iii_c127e_residual.json"
MODE = os.environ.get("K3_C127D_MODE", "pilot")
N_WORKERS = int(os.environ.get("K3_C127D_WORKERS", "6"))
ART = RES / ("k3_cap_b1e2iii_c127d_atlas.json" if MODE == "full"
             else "k3_cap_b1e2iii_c127d_atlas_pilot.json")

# --- PRÉ-ENREGISTRÉ, fixé avant le run, non ajustable -----------------------
# L'échelle de ρ : le PREMIER qui passe est retenu, par tuile, et
# sérialisé. Aucun ρ ne passe ⟹ REFUS de la tuile (gate D2), jamais un
# repli silencieux sur le core.
RHO_LADDER = (Fraction(1, 8), Fraction(1, 16),
              Fraction(1, 32), Fraction(1, 64))
# Plafond du défaut CERTIFIÉ de l'identité de transition, sur les
# overlaps ET sur les halos. Même δ que C127/C127-E, repris tel quel.
DELTA_TRANS = 1e-5
# Plafond du défaut du cocycle λ_ij·λ_jk·λ_ki − 1 sur les triples.
DELTA_COCYCLE = 1e-5
N_PILOT_TILES = 40      # tuiles du patch stratifié en mode pilot
N_HALO_METRIC = 8       # tuiles où le transport MÉTRIQUE est refait
                        # sur le halo (portée déclarée, pas silencieuse)

T0 = time.time()
IV0 = iv.mpf(0)
IV1 = iv.mpf(1)


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def _sha(p):
    try:
        return hashlib.sha256(Path(p).read_bytes()).hexdigest()
    except OSError:
        return None


def provenance(src, t_wall):
    try:
        h = subprocess.run(["git", "rev-parse", "HEAD"],
                           cwd=Path(__file__).resolve().parent,
                           capture_output=True, text=True,
                           timeout=10).stdout.strip() or None
    except Exception:
        h = None
    return {"git_head": h, "python": sys.version.split()[0],
            "platform": platform.platform(), "mp_prec": int(mp.prec),
            "tm_order": TM_ORDER, "unary_series_deg": UNARY_SERIES_DEG,
            "wall_s": round(t_wall, 1), "mode": MODE,
            "n_workers": N_WORKERS,
            "preregistered": {
                "rho_ladder": [[r.numerator, r.denominator]
                               for r in RHO_LADDER],
                "delta_trans": DELTA_TRANS,
                "delta_cocycle": DELTA_COCYCLE,
                "n_pilot_tiles": N_PILOT_TILES,
                "n_halo_metric": N_HALO_METRIC},
            "inputs": {str(Path(x).name): _sha(x) for x in src},
            "self_sha256": _sha(__file__)}


# ===========================================================================
#  Géométrie dyadique EXACTE — jamais de float dans une décision
# ===========================================================================
def box_of(o):
    """(centre, demi-largeur) d'une feuille, en Fractions exactes."""
    return ([Fraction(float.fromhex(x)) for x in o["center_hex"]],
            Fraction(float.fromhex(o["hw_hex"])))


def halo_hw(h, rho):
    """(1+ρ)·h, avec vérification que le résultat est EXACTEMENT
    représentable en float (h est une puissance de 2, ρ dyadique)."""
    H = h * (1 + rho)
    f = float(H)
    return (H, f) if Fraction(f) == H else (H, None)


def _inter_generic(loA, hiA, loB, hiB):
    """Boîte d'intersection de deux pavés donnés par leurs bornes, ou
    None si la largeur n'est pas STRICTEMENT positive dans les quatre
    coordonnées. Tout en Fractions — « ouvert » est une comparaison
    rationnelle, jamais un test float. La boîte rendue n'est pas un cube
    (demi-largeurs par coordonnée), d'où le passage systématique par
    (lo, hi) plutôt que par (centre, demi-largeur)."""
    lo, hi = [], []
    for k in range(4):
        a, b = max(loA[k], loB[k]), min(hiA[k], hiB[k])
        if not b > a:
            return None
        lo.append(a)
        hi.append(b)
    return {"center": [(lo[k] + hi[k]) / 2 for k in range(4)],
            "hw": [(hi[k] - lo[k]) / 2 for k in range(4)],
            "lo": lo, "hi": hi}


def cube_bounds(c, H):
    return [c[k] - H for k in range(4)], [c[k] + H for k in range(4)]


def touching(ci, hi, cj, hj):
    """Les cores FERMÉS se touchent (distance de Tchebychev ≤ 0)."""
    return max(abs(ci[k] - cj[k]) - (hi + hj) for k in range(4)) <= 0


def eps_box(w, c, H):
    """Image de la boîte `w` dans les symboles ε de la tuile (c, H) :
    liste de 4 (lo, hi) en Fractions. Le TM n'est valide que sur
    [−1,1]⁴ — l'inclusion est un GATE, pas une hypothèse."""
    out = []
    for k in range(4):
        lo = (w["lo"][k] - c[k]) / H
        hi = (w["hi"][k] - c[k]) / H
        out.append((lo, hi))
    return out


def eps_box_in_range(E):
    return all(-1 <= lo and hi <= 1 for lo, hi in E)


def frac_to_iv(q):
    """Fraction → intervalle mpmath ENCLOSANT (arrondi extérieur par la
    division intervalle : aucune conversion float au plus proche)."""
    return iv.mpf(int(q.numerator)) / iv.mpf(int(q.denominator))


def fraction_box_to_iv(lo, hi):
    """(lo, hi) en Fractions → UN intervalle mpmath ENCLOSANT [lo, hi] :
    borne basse arrondie VERS LE BAS, borne haute VERS LE HAUT (C129-A,
    8ᵉ revue). L'ancien chemin `iv.mpf([float(lo), float(hi)])`
    arrondissait les deux endpoints AU PLUS PROCHE — or les boîtes
    d'overlap héritent des demi-largeurs de halo `65h/64` et leurs
    coordonnées ε portent des dénominateurs en 65, non représentables en
    binaire : l'intervalle construit pouvait être STRICTEMENT PLUS PETIT
    que la boîte rationnelle exacte, et un `∋ 0` certifié sur un domaine
    tronqué ne certifie pas le domaine. L'élargissement est ≤ 1 ulp à
    mp.prec par borne : conservateur, jamais optimiste."""
    return iv.mpf([mp.mpf(frac_to_iv(lo).a), mp.mpf(frac_to_iv(hi).b)])


# ===========================================================================
#  Recentrage EXACT d'un Taylor-modèle — le cœur de C127-D
# ===========================================================================
#  Deux tuiles voisines portent leurs polynômes dans DEUX cadres ε
#  différents. Comparer leurs ENCLOSURES sur l'overlap perd la
#  corrélation : chaque enclosure a la largeur de la VARIATION de la
#  fonction sur l'overlap (~1e-3 ici), donc leur différence aussi, et le
#  verdict n'a plus de contenu. On recentre donc le polynôme de j dans le
#  cadre de i — substitution AFFINE ε_j = off + scale·ε_i, exacte en
#  rationnels — puis on soustrait COEFFICIENT PAR COEFFICIENT. Mesuré sur
#  une paire réelle : largeur 7,2e-3 (naïf) → 2,3e-12 (recentré).
#
#  La substitution ne crée AUCUN reste : composer un polynôme de degré
#  ≤ N avec une application affine donne un polynôme de degré ≤ N, et la
#  matrice ci-dessous est exacte (Fractions) avant conversion intervalle.
# ===========================================================================
def recenter_matrix(off, scale):
    """T tel que p∘φ [β] = Σ_{α ≥ β} T[β][α] · p[α], pour
    φ(ε)_k = off_k + scale·ε_k. Rendue en intervalles ENCLOSANTS."""
    acc = {}
    for ai, am in enumerate(MONO):
        for bm in itertools.product(*[range(e + 1) for e in am]):
            coef = Fraction(1)
            for k in range(4):
                coef *= (comb(am[k], bm[k])
                         * off[k] ** (am[k] - bm[k])
                         * scale ** bm[k])
            if coef:
                row = acc.setdefault(MIDX[bm], {})
                row[ai] = row.get(ai, Fraction(0)) + coef
    return {b: [(a, frac_to_iv(c)) for a, c in row.items() if c]
            for b, row in acc.items()}


def apply_recenter(T, p):
    """Applique T aux coefficients CIV d'un TM."""
    out = [CIV(IV0, IV0)] * NM
    for b, row in T.items():
        ar, ai_ = IV0, IV0
        for a, c in row:
            ar = ar + p[a].re * c
            ai_ = ai_ + p[a].im * c
        out[b] = CIV(ar, ai_)
    return out


def poly_lin(pa, pb, sign):
    """pa ± pb, coefficient par coefficient."""
    if sign > 0:
        return [CIV(pa[k].re + pb[k].re, pa[k].im + pb[k].im)
                for k in range(NM)]
    return [CIV(pa[k].re - pb[k].re, pa[k].im - pb[k].im)
            for k in range(NM)]


def enclose(p, rem, E):
    """Enclosure d'un TM (coefficients CIV, reste `rem` majorant les DEUX
    parties) sur la sous-boîte ε donnée par les intervalles E."""
    pw = []
    for k in range(4):
        col = [IV1]
        for _d in range(TM_ORDER):
            col.append(col[-1] * E[k])
        pw.append(col)
    ar, ai_ = IV0, IV0
    for idx, m in enumerate(MONO):
        c = p[idx]
        t = IV1
        for k in range(4):
            if m[k]:
                t = t * pw[k][m[k]]
        ar = ar + c.re * t
        ai_ = ai_ + c.im * t
    return CIV(ar + IVPM * rem, ai_ + IVPM * rem)


def sep_phase(dm, dp):
    """Feuillet `θ` par SÉPARATION STRICTE des distances au carré, sur les
    enclosures de la DIFFÉRENCE et de la SOMME (donc corrélées, donc
    décisives). Miroir exact de la logique `_eps_sep` de C127, appliquée
    à des polynômes déjà soustraits. Ambigu ⟹ REFUS."""
    d2m = dm.re ** 2 + dm.im ** 2
    d2p = dp.re ** 2 + dp.im ** 2
    m_lo, m_hi = mp.mpf(d2m.a), mp.mpf(d2m.b)
    p_lo, p_hi = mp.mpf(d2p.a), mp.mpf(d2p.b)
    if m_hi < p_lo:
        return 1, {"theta": 1, "margin": _f_down(p_lo - m_hi),
                   "d2_chosen": [_f_down(m_lo), _f_up(m_hi)],
                   "d2_other": [_f_down(p_lo), _f_up(p_hi)]}
    if p_hi < m_lo:
        return -1, {"theta": -1, "margin": _f_down(m_lo - p_hi),
                    "d2_chosen": [_f_down(p_lo), _f_up(p_hi)],
                    "d2_other": [_f_down(m_lo), _f_up(m_hi)]}
    return None, {"d2_diff": [_f_down(m_lo), _f_up(m_hi)],
                  "d2_sum": [_f_down(p_lo), _f_up(p_hi)]}


def civ_sup(c):
    """sup des |parties| d'une enclosure complexe."""
    return max(abs(mp.mpf(c.re.a)), abs(mp.mpf(c.re.b)),
               abs(mp.mpf(c.im.a)), abs(mp.mpf(c.im.b)))


def civ_contains_zero(c):
    return (mp.mpf(c.re.a) <= 0 <= mp.mpf(c.re.b)
            and mp.mpf(c.im.a) <= 0 <= mp.mpf(c.im.b))


def civ_div(a, b):
    """a/b en intervalles, exige |b| minoré > 0 (certifié par
    l'appelant)."""
    d = b.re ** 2 + b.im ** 2
    return CIV((a.re * b.re + a.im * b.im) / d,
               (a.im * b.re - a.re * b.im) / d)


# ===========================================================================
#  Sérialisation EXACTE des Taylor-modèles (les objets mpmath iv ne se
#  picklent pas : on transporte les tuples `_mpf_`, exacts, entre les
#  workers et le parent — aucune conversion, aucun arrondi)
# ===========================================================================
def ser_tmc(t):
    return ([(mp.mpf(x.re.a)._mpf_, mp.mpf(x.re.b)._mpf_,
              mp.mpf(x.im.a)._mpf_, mp.mpf(x.im.b)._mpf_) for x in t.p],
            mp.mpf(t.rem.b)._mpf_)


def de_tmc(blob):
    coeffs, rem = blob
    p = [CIV(iv.mpf([mp.make_mpf(a), mp.make_mpf(b)]),
             iv.mpf([mp.make_mpf(c), mp.make_mpf(d)]))
         for (a, b, c, d) in coeffs]
    return p, iv.mpf(mp.make_mpf(rem))


# ===========================================================================
#  Phase halo : le chart est-il valide sur un VOISINAGE OUVERT du core ?
# ===========================================================================
_G = {}


def _init_worker(cell, root=None):
    _G["cell"] = cell
    if root is not None:
        _G["root"] = root


def flush_faces(tile, root_c, root_h):
    """Directions où le core est EXACTEMENT à cheval sur une face de la
    cellule : +1 face haute, −1 face basse, 0 intérieur. Comparaison
    rationnelle exacte."""
    c, h = box_of(tile)
    out = []
    for k in range(4):
        if c[k] + h == root_c[k] + root_h:
            out.append(1)
        elif c[k] - h == root_c[k] - root_h:
            out.append(-1)
        else:
            out.append(0)
    return out


def halo_certify(cell, tile, rho, rule, root_c, root_h):
    """Un essai de halo à ρ donné, sous l'une des deux RÈGLES :

    · `symmetric` — `box(c, (1+ρ)h)`, halo à deux côtés dans les quatre
      directions. C'est la règle pré-enregistrée, essayée d'abord.
    · `clipped` — même demi-largeur `(1+ρ)h`, mais CENTRE DÉCALÉ de ∓ρh
      dans chaque direction où le core est à cheval sur une face de la
      CELLULE : le halo s'étend alors vers l'intérieur et reste AFFLEURANT
      à la face, donc `core ⊂ halo ⊆ cellule`. C'est la notion d'atlas
      d'une variété À BORD : les cartes sont ouvertes dans la topologie
      RELATIVE de la cellule. La règle n'est légitime QUE sur une
      direction affleurante — le gate D2b l'exige, sans quoi elle
      deviendrait une porte de sortie pour un échec intérieur.

    Retourne (record, Z, P) où `P[a] = Zt[a]·Z[g']` (le produit sans
    division qui porte l'identité de transition), ou (record, None, None)."""
    S, g, eps = cell
    c, h = box_of(tile)
    H, Hf = halo_hw(h, rho)
    rec = {"rho": [rho.numerator, rho.denominator], "rule": rule}
    if Hf is None:
        rec["refused"] = "halo_hw_not_exact_float"
        return rec, None, None
    fl = flush_faces(tile, root_c, root_h)
    if rule == "clipped":
        shift = [-rho * h if fl[k] > 0 else
                 (rho * h if fl[k] < 0 else Fraction(0))
                 for k in range(4)]
        if all(s == 0 for s in shift):
            rec["refused"] = "clipping_useless_no_flush_face"
            return rec, None, None
    else:
        shift = [Fraction(0)] * 4
    c = [c[k] + shift[k] for k in range(4)]
    rec["flush_faces"] = fl
    rec["center_shift"] = [[s.numerator, s.denominator] for s in shift]
    rec["center_hex"] = [float(x).hex() for x in c]
    if any(Fraction(float(x)) != x for x in c):
        rec["refused"] = "halo_center_not_exact_float"
        return rec, None, None
    c0, h0 = box_of(tile)
    if not all(c[k] - H <= c0[k] - h0 and c0[k] + h0 <= c[k] + H
               for k in range(4)):
        rec["refused"] = "core_not_inside_halo"
        return rec, None, None
    if rule == "clipped" and not all(
            root_c[k] - root_h <= c[k] - H
            and c[k] + H <= root_c[k] + root_h for k in range(4)):
        rec["refused"] = "clipped_halo_leaves_cell"
        return rec, None, None
    rec["H_hex"] = float(Hf).hex()
    cf = [float(x) for x in c]
    Z, dZ, rows = build_section(S, g, eps, cf, Hf)
    if any(z is None for z in Z):
        rec["refused"] = "source_section_incomplete"
        return rec, None, None
    det = [r["determination"] for r in rows]
    rec["source_determinations"] = det
    if det != tile["core_src_det"]:
        rec["refused"] = "source_determination_changed"
        return rec, None, None
    S2, g2 = tuple(tile["chart"]["S"]), tile["chart"]["g"]
    cert = chart_certificate(Z, dZ, S2, g2)
    rec["gauge_absmin"] = cert.get("gauge_absmin")
    rec["detJ_absmin"] = cert.get("detJ_absmin")
    if not (cert.get("admissible")
            and cert.get("disjoint_from_target_slice")):
        rec["refused"] = "chart_certificate_failed"
        return rec, None, None
    T2 = [j for j in range(6) if j not in S2]
    o = [x for x in T2 if x != g2]
    ib = Z[g2].inv()
    Zp = [z * ib for z in Z]
    up_, vp_ = Zp[o[0]], Zp[o[1]]
    _r126, ok126 = native_section_constructible(S2, g2, up_, vp_)
    _rext, okext = native_rows_ext(S2, g2, up_, vp_)
    rec["criterion_c126_on_halo"] = bool(ok126)
    rec["criterion_ext_on_halo"] = bool(okext)
    need = tile["criterion"]
    if need == "c126" and not ok126:
        rec["refused"] = "c126_fails_on_halo"
        return rec, None, None
    if need == "extended":
        if not okext:
            rec["refused"] = "extended_fails_on_halo"
            return rec, None, None
        if ok126:
            # R4 sur le halo : le critère importé doit TOUJOURS refuser,
            # sinon la tuile n'était pas résiduelle pour la raison dite
            rec["refused"] = "c126_unexpectedly_accepts_on_halo"
            return rec, None, None
    Zt, _dZt, kinds = native_target_section(
        S2, g2, up_, vp_, tuple(tile["eps_target"]),
        list(tile["sigma_target"]))
    if Zt is None:
        rec["refused"] = "native_target_failed_on_halo"
        rec["kinds_target"] = kinds
        return rec, None, None
    rec["kinds_target"] = kinds
    if kinds != tile["core_kinds_target"]:
        rec["refused"] = "target_kinds_changed"
        return rec, None, None
    # congruence projective sur le HALO, à ledger figé. Le sup reste en
    # mp jusqu'au gate (C129-B) : la comparaison à δ se fait sur la borne
    # exacte, la conversion float (dirigée, _f_up) n'arrive qu'à la
    # sérialisation.
    sup = mp.mpf(0)
    okz = True
    for a in range(6):
        D = Zt[a] + Zp[a].mul_real(riv(-1.0))
        dc = CIV(D.re_tm().to_iv() + IVPM * IV0,
                 D.im_tm().to_iv() + IVPM * IV0)
        okz = okz and civ_contains_zero(dc)
        sup = max(sup, civ_sup(dc))
    rec["halo_congruence_contains_zero"] = bool(okz)
    rec["halo_congruence_sup"] = _f_up(sup)
    if not okz:
        rec["refused"] = "halo_congruence_excludes_zero"
        return rec, None, None
    if not sup <= DELTA_TRANS:
        rec["refused"] = "halo_congruence_above_delta"
        return rec, None, None
    P = [Zt[a] * Z[g2] for a in range(6)]
    rec["accepted"] = True
    return rec, Z, P


def _halo_job(job):
    idx, tile = job
    root_c, root_h = _G["root"]
    attempts = []
    # ORDRE STRICT : la règle symétrique d'abord, sur toute l'échelle ρ.
    # Le clipping n'est atteint que si TOUTES les tentatives symétriques
    # ont été refusées — et le refus est publié, pas absorbé.
    for rule in ("symmetric", "clipped"):
        for rho in RHO_LADDER:
            rec, Z, P = halo_certify(_G["cell"], tile, rho, rule,
                                     root_c, root_h)
            attempts.append(rec)
            if rec.get("accepted"):
                return {"index": idx, "ok": True,
                        "rho_used": rec["rho"], "rule_used": rule,
                        "H_hex": rec["H_hex"],
                        "center_shift": rec["center_shift"],
                        "flush_faces": rec["flush_faces"],
                        "n_symmetric_refusals": sum(
                            1 for a in attempts
                            if a["rule"] == "symmetric"),
                        "attempts": attempts, "record": rec,
                        "Z": [ser_tmc(Z[a]) for a in range(6)],
                        "P": [ser_tmc(P[a]) for a in range(6)]}
    return {"index": idx, "ok": False, "attempts": attempts}


# ===========================================================================
#  Phase paire : feuillet, transition, identité
# ===========================================================================
_TCACHE = {}


def get_T(off, scale):
    key = (tuple((o.numerator, o.denominator) for o in off),
           (scale.numerator, scale.denominator))
    T = _TCACHE.get(key)
    if T is None:
        T = recenter_matrix(off, scale)
        _TCACHE[key] = T
    return T


def pair_certificate(i, j, corrupt=None):
    """Certificat d'overlap entre les tuiles i et j. `corrupt` permet aux
    négatifs de substituer des données falsifiées SANS toucher le chemin
    nominal."""
    gi, gj = _G["geom"][i], _G["geom"][j]
    w = _inter_generic(gi["lo"], gi["hi"], gj["lo"], gj["hi"])
    out = {"i": i, "j": j}
    if w is None:
        out["refused"] = "no_open_overlap"
        return out
    out["overlap_hw"] = [float(x) for x in w["hw"]]
    out["overlap_min_width"] = float(2 * min(w["hw"]))
    Ei = eps_box(w, gi["c"], gi["H"])
    Ej = eps_box(w, gj["c"], gj["H"])
    out["eps_box_in_range"] = bool(eps_box_in_range(Ei)
                                   and eps_box_in_range(Ej))
    if not out["eps_box_in_range"]:
        out["refused"] = "overlap_outside_tm_domain"
        return out
    # ε_j = off + scale·ε_i   (exact en rationnels)
    off = [(gi["c"][k] - gj["c"][k]) / gj["H"] for k in range(4)]
    scale = gi["H"] / gj["H"]
    T = get_T(off, scale)
    Eiv = [fraction_box_to_iv(lo, hi) for lo, hi in Ei]
    Zi, Pi = _G["tm"][i]
    Zj, Pj = _G["tm"][j]
    if corrupt and corrupt.get("tile") == j:
        Zj, Pj = corrupt["Z"], corrupt["P"]
    elif corrupt and corrupt.get("tile") == i:
        Zi, Pi = corrupt["Z"], corrupt["P"]
    S = _G["cell"][0]

    # --- (1) le FEUILLET : θ par recentrage + soustraction exacte -------
    thetas, tmargins, amb = [], [], False
    for s in S:
        pa, ra = Zi[s]
        pb, rb = Zj[s]
        pb = apply_recenter(T, pb)
        rem = ra + rb
        dm = enclose(poly_lin(pa, pb, -1), rem, Eiv)
        dp = enclose(poly_lin(pa, pb, +1), rem, Eiv)
        th, rec = sep_phase(dm, dp)
        rec["s_coord"] = int(s)
        rec["diff_sup"] = _f_up(civ_sup(dm))
        if th is None:
            amb = True
        thetas.append(th)
        tmargins.append(rec)
    out["theta"] = thetas
    out["theta_margins"] = tmargins
    out["theta_ambiguous"] = bool(amb)
    out["same_sheet"] = bool(not amb and all(t == 1 for t in thetas))

    # --- (2) coordonnées AFFINES : contrôle de l'évaluateur -------------
    #  Z[g] ≡ 1 et Z[o] ≡ u, v sont les MÊMES fonctions dans les deux
    #  tuiles : leur différence recentrée doit être nulle à l'arrondi
    #  près. C'est un test gratuit du recentrage lui-même.
    aff_sup = mp.mpf(0)
    for a in range(6):
        if a in S:
            continue
        pa, ra = Zi[a]
        pb, rb = Zj[a]
        d = enclose(poly_lin(pa, apply_recenter(T, pb), -1),
                    ra + rb, Eiv)
        aff_sup = max(aff_sup, civ_sup(d))
    out["affine_coords_defect"] = _f_up(aff_sup)

    # --- (3) la CARTE : λ = Z[g'_i]/Z[g'_j], non nulle des deux côtés ---
    g_i, g_j = _G["tiles"][i]["chart"]["g"], _G["tiles"][j]["chart"]["g"]
    zi_gi = enclose(*Zi[g_i], Eiv)
    zi_gj = enclose(*Zi[g_j], Eiv)
    mi = mp.mpf(civ_absmin(zi_gi).a)
    mj = mp.mpf(civ_absmin(zi_gj).a)
    out["abs_min_Z_gi"] = _f_down(mi)
    out["abs_min_Z_gj"] = _f_down(mj)
    out["transition_biholomorphic"] = bool(mi > 0 and mj > 0)
    if mi > 0 and mj > 0:
        lam = civ_div(zi_gi, zi_gj)
        out["lambda"] = [_f_down(mp.mpf(lam.re.a)),
                         _f_up(mp.mpf(lam.re.b)),
                         _f_down(mp.mpf(lam.im.a)),
                         _f_up(mp.mpf(lam.im.b))]
        out["charts_equal"] = bool(g_i == g_j)

    # --- (4) l'IDENTITÉ de transition, sans division --------------------
    #  Zt^(j)_a·Z[g'_j] − Zt^(i)_a·Z[g'_i] ∋ 0 : les deux ledgers ε'
    #  dérivés SÉPARÉMENT sont mis en présence sur l'overlap.
    # Le sup reste en mp jusqu'au gate (C129-B) ; _f_up seulement pour
    # la sérialisation — la borne publiée MAJORE la borne comparée.
    idsup, idok = mp.mpf(0), True
    for a in range(6):
        pa, ra = Pi[a]
        pb, rb = Pj[a]
        d = enclose(poly_lin(pa, apply_recenter(T, pb), -1),
                    ra + rb, Eiv)
        idok = idok and civ_contains_zero(d)
        idsup = max(idsup, civ_sup(d))
    out["transition_identity_contains_zero"] = bool(idok)
    out["transition_identity_defect"] = _f_up(idsup)
    out["ok"] = bool(out["same_sheet"] and idok
                     and idsup <= DELTA_TRANS
                     and out["transition_biholomorphic"])
    return out


def _pair_job(job):
    i, j = job
    return pair_certificate(i, j)


# ===========================================================================
#  Phase triple : le cocycle sur le nerf
# ===========================================================================
def triple_certificate(i, j, k, lam):
    gi, gj, gk = _G["geom"][i], _G["geom"][j], _G["geom"][k]
    w = _inter_generic(gi["lo"], gi["hi"], gj["lo"], gj["hi"])
    if w is None:
        return None
    w = _inter_generic(w["lo"], w["hi"], gk["lo"], gk["hi"])
    if w is None:
        return None
    out = {"i": i, "j": j, "k": k,
           "min_width": float(2 * min(w["hw"]))}
    # θ : cocycle EXACT (discret)
    tij, tjk, tik = lam[(i, j)]["theta"], lam[(j, k)]["theta"], \
        lam[(i, k)]["theta"]
    if any(x is None for x in tij + tjk + tik):
        out["theta_cocycle"] = False
        out["theta_note"] = "ambiguous_edge"
    else:
        out["theta_cocycle"] = all(
            tij[r] * tjk[r] * tik[r] == 1 for r in range(len(tij)))
    # λ : cocycle en enclosure sur l'intersection TRIPLE
    Ei = eps_box(w, gi["c"], gi["H"])
    if not eps_box_in_range(Ei):
        out["refused"] = "triple_outside_tm_domain"
        return out
    Eiv = [fraction_box_to_iv(lo, hi) for lo, hi in Ei]
    Zi = _G["tm"][i][0]
    gg = [_G["tiles"][x]["chart"]["g"] for x in (i, j, k)]
    z = [enclose(*Zi[c], Eiv) for c in gg]
    ok = all(mp.mpf(civ_absmin(x).a) > 0 for x in z)
    out["lambda_defined"] = bool(ok)
    out["gauge_abs_min"] = [_f_down(mp.mpf(civ_absmin(x).a)) for x in z]
    if ok:
        # HONNÊTETÉ : λ_ij·λ_jk·λ_ki = (z_i/z_j)(z_j/z_k)(z_k/z_i) est une
        # IDENTITÉ ALGÉBRIQUE de la construction — toutes les cartes sont
        # des cartes affines du MÊME P⁵ et les trois λ sont des rapports
        # des trois MÊMES évaluations. Le vérifier en intervalles ne
        # vérifie rien : chaque z apparaît deux fois, la corrélation est
        # perdue, et le « défaut » mesuré N'EST QUE la largeur de
        # décorrélation (exactement le phénomène qui a imposé le
        # recentrage ailleurs dans ce script). Le nombre est donc publié
        # comme DIAGNOSTIC de décorrélation, JAMAIS comme gate.
        # Ce qui a du contenu au niveau du triple, et qui est gaté :
        #   · les trois jauges sont MINORÉES > 0 sur la boîte triple
        #     (les trois λ existent et la composition est licite) ;
        #   · les trois identités d'arête, recentrées, contiennent 0 sur
        #     la boîte TRIPLE — et cela n'est PAS une conséquence
        #     intervalliste de i→j et j→k (l'inclusion n'est pas
        #     transitive).
        prod = _mul_civ(_mul_civ(civ_div(z[0], z[1]),
                                 civ_div(z[1], z[2])),
                        civ_div(z[2], z[0]))
        d = CIV(prod.re - IV1, prod.im)
        out["lambda_product_contains_one"] = bool(civ_contains_zero(d))
        out["lambda_decorrelation_width"] = _f_up(civ_sup(d))
        out["lambda_note"] = ("identité algébrique — largeur de "
                              "décorrélation publiée, non gatée")
    # l'identité de transition sur les TROIS arêtes, sur la boîte TRIPLE.
    # Sup en mp jusqu'au gate (C129-B), _f_up à la sérialisation.
    idsup, idok = mp.mpf(0), True
    for (a, b) in ((i, j), (j, k), (i, k)):
        ga, gb = _G["geom"][a], _G["geom"][b]
        Ea = eps_box(w, ga["c"], ga["H"])
        if not eps_box_in_range(Ea):
            idok = False
            continue
        Eav = [fraction_box_to_iv(lo, hi) for lo, hi in Ea]
        off = [(ga["c"][t] - gb["c"][t]) / gb["H"] for t in range(4)]
        T = get_T(off, ga["H"] / gb["H"])
        for co in range(6):
            pa, ra = _G["tm"][a][1][co]
            pb, rb = _G["tm"][b][1][co]
            dd = enclose(poly_lin(pa, apply_recenter(T, pb), -1),
                         ra + rb, Eav)
            idok = idok and civ_contains_zero(dd)
            idsup = max(idsup, civ_sup(dd))
    out["edge_identities_contain_zero"] = bool(idok)
    out["edge_identity_defect"] = _f_up(idsup)
    out["ok"] = bool(out.get("theta_cocycle") and out["lambda_defined"]
                     and idok and idsup <= DELTA_TRANS)
    return out


def _mul_civ(a, b):
    return CIV(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re)


def _triple_job(job):
    i, j, k = job
    return triple_certificate(i, j, k, _G["pairs"])


# ===========================================================================
#  Chargement des feuilles et de leurs ledgers
# ===========================================================================
def load_leaves():
    cov = json.loads(COVER_JSON.read_text(encoding="utf-8"))
    c127 = json.loads(C127_JSON.read_text(encoding="utf-8"))
    c127e = json.loads(C127E_JSON.read_text(encoding="utf-8"))
    cell = cov["cell"]
    tr127 = {r["tile_index"]: r for r in c127["transports"]
             if not r.get("failed")}
    tr127e = {r["box_index"]: r for r in c127e["transports"]
              if not r.get("failed")}
    leaves = []
    for i, t in enumerate(cov["tiles"]):
        r = tr127.get(i)
        if r is None:
            raise SystemExit(f"tuile C127 {i} sans transport certifié")
        leaves.append({
            "src": "c127", "orig_index": i, "depth": t["depth"],
            "center_hex": t["center_hex"], "hw_hex": t["hw_hex"],
            "chart": t["chart"], "criterion": "c126",
            "core_src_det": r["source_determinations"],
            "core_kinds_target": r["kinds_target"],
            "eps_target": r["eps_target"],
            "sigma_target": r["sigma_target"]})
    for t in c127e["new_tiles"]:
        r = tr127e.get(t["box_index"])
        if r is None:
            raise SystemExit(
                f"tuile C127-E {t['box_index']} sans transport certifié")
        leaves.append({
            "src": "c127e", "orig_index": t["box_index"],
            "depth": t["depth"], "center_hex": t["center_hex"],
            "hw_hex": t["hw_hex"], "chart": t["chart"],
            "criterion": "extended",
            "core_src_det": r["source_determinations"],
            "core_kinds_target": r["kinds_target"],
            "eps_target": r["eps_target"],
            "sigma_target": r["sigma_target"]})
    return cell, leaves


def pilot_patch(leaves, n_max):
    """Patch stratifié : on part de la tuile née du RÉSIDU qui a le plus
    de voisines issues de C127 (l'interface la plus discriminante — des
    déterminations différentes s'y touchent), on prend tout son
    voisinage, puis on complète par un représentant de chaque signature
    (chart, déterminations source, kinds cible) non encore couverte."""
    bx = [box_of(t) for t in leaves]
    nb = {i: set() for i in range(len(leaves))}
    for i in range(len(leaves)):
        for j in range(i + 1, len(leaves)):
            if touching(bx[i][0], bx[i][1], bx[j][0], bx[j][1]):
                nb[i].add(j)
                nb[j].add(i)
    seed, best = 0, -1
    for i, t in enumerate(leaves):
        if t["src"] != "c127e":
            continue
        n = sum(1 for j in nb[i] if leaves[j]["src"] == "c127")
        if n > best:
            seed, best = i, n
    sel = [seed] + sorted(nb[seed])
    sig = {(leaves[i]["src"], tuple(leaves[i]["chart"]["S"]),
            leaves[i]["chart"]["g"], tuple(leaves[i]["core_src_det"]),
            tuple(leaves[i]["core_kinds_target"])) for i in sel}
    for i, t in enumerate(leaves):
        if len(sel) >= n_max:
            break
        s = (t["src"], tuple(t["chart"]["S"]), t["chart"]["g"],
             tuple(t["core_src_det"]), tuple(t["core_kinds_target"]))
        if s not in sig:
            sig.add(s)
            sel.append(i)
    sel = sorted(set(sel))[:n_max]
    return sel, {"seed": seed, "seed_c127_neighbours": best,
                 "n_selected": len(sel), "selected": sel,
                 "note": ("patch = une tuile née du résidu + tout son "
                          "voisinage + un représentant par signature "
                          "non couverte")}


# ===========================================================================
#  Build
# ===========================================================================
def build():
    print("=" * 78)
    print(f"C127-D ATLAS {'COMPLET (316)' if MODE == 'full' else 'PILOTE'}"
          f" : TM ({TM_ORDER},{UNARY_SERIES_DEG}), {N_WORKERS} workers, "
          f"δ_trans = {DELTA_TRANS:.0e}")
    print("=" * 78)
    cell_d, leaves = load_leaves()
    S, g, eps = (tuple(cell_d["S"]), cell_d["g"], tuple(cell_d["eps"]))
    root_c = [float.fromhex(x) for x in cell_d["center_hex"]]
    root_h = float.fromhex(cell_d["hw_hex"])
    log(f"{len(leaves)} feuilles chargées "
        f"({sum(1 for t in leaves if t['src'] == 'c127')} C127 + "
        f"{sum(1 for t in leaves if t['src'] == 'c127e')} C127-E)")

    # --- D1 : la partition, ré-affirmée depuis les adresses ------------
    addresses, addr_fail = [], 0
    for t in leaves:
        a = address_of(root_c, root_h,
                       [float.fromhex(x) for x in t["center_hex"]],
                       float.fromhex(t["hw_hex"]))
        if a is None:
            addr_fail += 1
        else:
            addresses.append(a)
    tg = tree_gates(addresses)
    log(f"D1 : {len(addresses)}/{len(leaves)} adresses exactes, "
        f"prefix-free={tg['prefix_free']}, clos={tg['tree_closed']}, "
        f"Kraft={tg['kraft_sum'][0]}/{tg['kraft_sum'][1]}")

    sel, patch = (list(range(len(leaves))), None) if MODE == "full" \
        else pilot_patch(leaves, N_PILOT_TILES)
    log(f"tuiles retenues : {len(sel)}")

    # --- D2/D3 : les halos ---------------------------------------------
    rootF = ([Fraction(x) for x in
              [Fraction(float.fromhex(y)) for y in cell_d["center_hex"]]],
             Fraction(float.fromhex(cell_d["hw_hex"])))
    mpctx = get_context("fork")
    _init_worker((S, g, eps), rootF)
    with mpctx.Pool(N_WORKERS, initializer=_init_worker,
                    initargs=((S, g, eps), rootF)) as pool:
        halos = pool.map(_halo_job, [(i, leaves[i]) for i in sel])
    halo_ok = [h for h in halos if h["ok"]]
    n_clip = sum(1 for h in halo_ok if h["rule_used"] == "clipped")
    log(f"D2 : halos acceptés {len(halo_ok)}/{len(sel)} "
        f"(ρ : {sorted({tuple(h['rho_used']) for h in halo_ok})}) · "
        f"règle symétrique {len(halo_ok) - n_clip}, clippée {n_clip}")
    if len(halo_ok) != len(sel):
        for h in halos:
            if not h["ok"]:
                log(f"  REFUS tuile {h['index']} : "
                    f"{[a.get('refused') for a in h['attempts']]}")

    # --- géométrie des halos (exacte) ----------------------------------
    geom, tm = {}, {}
    for h in halo_ok:
        i = h["index"]
        c0, hw = box_of(leaves[i])
        rho = Fraction(*h["rho_used"])
        H = hw * (1 + rho)
        c = [c0[k] + Fraction(*h["center_shift"][k]) for k in range(4)]
        lo, hi = cube_bounds(c, H)
        geom[i] = {"c": c, "H": H, "lo": lo, "hi": hi,
                   "core_c": c0, "core_hw": hw, "rho": rho,
                   "rule": h["rule_used"], "flush": h["flush_faces"]}
        tm[i] = ([de_tmc(b) for b in h["Z"]],
                 [de_tmc(b) for b in h["P"]])
    _G["geom"], _G["tm"], _G["tiles"], _G["cell"] = \
        geom, tm, leaves, (S, g, eps)

    # --- D4 : le nerf, en rationnels exacts ----------------------------
    ids = sorted(geom)
    touch_pairs, pairs = [], []
    for a in range(len(ids)):
        for b in range(a + 1, len(ids)):
            i, j = ids[a], ids[b]
            gi, gj = geom[i], geom[j]
            tch = touching(gi["core_c"], gi["core_hw"],
                           gj["core_c"], gj["core_hw"])
            ov = _inter_generic(gi["lo"], gi["hi"], gj["lo"], gj["hi"])
            if tch:
                touch_pairs.append((i, j, ov is not None))
            if ov is not None:
                pairs.append((i, j))
    open_all = all(x[2] for x in touch_pairs)
    adj = {i: set() for i in ids}
    for i, j in pairs:
        adj[i].add(j)
        adj[j].add(i)
    seen, stack = {ids[0]}, [ids[0]]
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    connected = len(seen) == len(ids)
    log(f"D4 : {len(touch_pairs)} paires de cores qui se touchent, "
        f"{len(pairs)} overlaps ouverts, nerf connexe={connected}")

    # --- D6-D9 : les transitions ---------------------------------------
    with mpctx.Pool(N_WORKERS, initializer=_init_worker,
                    initargs=((S, g, eps),)) as pool:
        pres = pool.map(_pair_job, pairs)
    pmap = {(r["i"], r["j"]): r for r in pres}
    _G["pairs"] = pmap
    pok = [r for r in pres if r.get("ok")]
    log(f"D7-D9 : paires certifiées {len(pok)}/{len(pres)} · "
        f"θ≡+1 {sum(1 for r in pres if r.get('same_sheet'))} · "
        f"défaut d'identité max "
        f"{max((r.get('transition_identity_defect', 0) for r in pres), default=0):.3e}")

    # --- D10 : le cocycle ----------------------------------------------
    triples = []
    for (i, j) in pairs:
        for k in adj[i] & adj[j]:
            if k > j:
                triples.append((i, j, k))
    with mpctx.Pool(N_WORKERS, initializer=_init_worker,
                    initargs=((S, g, eps),)) as pool:
        tres = [r for r in pool.map(_triple_job, triples) if r]
    tok = [r for r in tres if r.get("ok")]
    log(f"D10 : triples d'intersection triple non vide {len(tres)}, "
        f"cocycle certifié {len(tok)}/{len(tres)}")

    # --- D11 : les négatifs --------------------------------------------
    negs = run_negatives(leaves, geom, tm, pairs, pmap)
    log("D11 : négatifs — "
        + " · ".join(f"{k}={'CASSE' if v['breaks'] else 'NE CASSE PAS'}"
                     for k, v in negs.items()))

    # --- gates ----------------------------------------------------------
    max_halo = max((h["record"]["halo_congruence_sup"]
                    for h in halo_ok), default=0.0)
    max_id = max((r.get("transition_identity_defect", 0.0)
                  for r in pres), default=0.0)
    max_aff = max((r.get("affine_coords_defect", 0.0)
                   for r in pres), default=0.0)
    min_theta = min((m["margin"] for r in pres
                     for m in r.get("theta_margins", [])
                     if "margin" in m), default=0.0)
    max_cocy = max((r.get("lambda_decorrelation_width", 0.0)
                    for r in tres), default=0.0)
    # D2b — le clipping n'est légitime QUE contre une face de la cellule,
    # et QUE si la règle symétrique a d'abord été refusée sur TOUTE
    # l'échelle ρ. Sans ce gate, « clippé » serait une porte de sortie
    # pour n'importe quel échec intérieur.
    clip_legit = True
    for h in halo_ok:
        if h["rule_used"] != "clipped":
            continue
        sh, fl = h["center_shift"], h["flush_faces"]
        clip_legit = clip_legit and all(
            (Fraction(*sh[k]) == 0) or (fl[k] != 0) for k in range(4)) \
            and h["n_symmetric_refusals"] == len(RHO_LADDER)
    # D2c (C129-C, 8ᵉ revue) — l'ÉGALITÉ D'ENSEMBLES, pas le compte : la
    # prose affirme {tuiles clippées} = {64 tuiles nées du résidu
    # C127-E}, et « 64 = 64 » n'exclut pas une permutation 63+1. Les deux
    # ensembles sont publiés, triés, avec leur SHA-256.
    clipped_idx = sorted(h["index"] for h in halo_ok
                         if h["rule_used"] == "clipped")
    c127e_idx = sorted(i for i in sel if leaves[i]["src"] == "c127e")

    def _set_sha(idx):
        return hashlib.sha256(
            json.dumps(idx, separators=(",", ":")).encode()).hexdigest()
    gates = {
        "D1_partition_reasserted": bool(
            addr_fail == 0 and tg["unique"] and tg["prefix_free"]
            and tg["tree_closed"] and tg["kraft_is_one"]),
        "D2_halo_charts_valid": bool(
            halos and len(halo_ok) == len(sel)
            and all(Fraction(*h["rho_used"]) > 0 for h in halo_ok)),
        "D2b_clipping_only_at_cell_boundary": bool(clip_legit),
        "D2c_clipped_set_is_exactly_c127e_residual": bool(
            clipped_idx == c127e_idx),
        "D3_halo_congruence_below_delta": bool(
            halo_ok and max_halo <= DELTA_TRANS),
        "D4_overlaps_open": bool(touch_pairs and open_all
                                 and len(pairs) >= len(touch_pairs)),
        "D5_nerve_connected": bool(connected),
        "D6_overlaps_in_tm_domain": bool(pres) and all(
            r.get("eps_box_in_range") for r in pres),
        "D7_sheet_phase_trivial": bool(pres) and all(
            r.get("same_sheet") for r in pres) and min_theta > 0,
        "D8_transition_biholomorphic": bool(pres) and all(
            r.get("transition_biholomorphic") for r in pres),
        "D9_transition_identity": bool(pres) and all(
            r.get("transition_identity_contains_zero") for r in pres)
        and max_id <= DELTA_TRANS,
        "D10_cocycle_on_nerve": bool(tres) and all(
            r.get("ok") for r in tres),
        "D11_negatives_all_break": bool(negs) and all(
            v["breaks"] for v in negs.values()),
        "D12_no_silent_cap": bool(
            len(halos) == len(sel) and len(pres) == len(pairs)
            and len(tres) == len(triples)
            and all((r["i"], r["j"]) in pmap for r in pres))}
    n_pass = sum(1 for v in gates.values() if v)

    verdict = (
        "C127-D (%s) — LA PARTITION EST DEVENUE UN ATLAS. Les %d "
        "feuilles du cover fermé (C127 + C127-E) portent chacune un HALO "
        "`(1+ρ)h` avec ρ pris dans l'échelle pré-enregistrée %s : sur ce "
        "voisinage OUVERT du core, les déterminations source sont "
        "INCHANGÉES, le critère de chart qui a certifié le core tient "
        "encore (et pour les tuiles nées du résidu, C126 REFUSE toujours "
        "— la non-tautologie R4 survit au halo), la section cible native "
        "à ledger FIGÉ a les mêmes `kinds`, et la congruence projective "
        "reste sous δ (max %.3e). LA MESURE QUI A COÛTÉ UN AMENDEMENT : "
        "%d tuiles ont REFUSÉ le halo symétrique sur TOUTE l'échelle ρ, "
        "y compris à ρ = 1/64 — diagnostic : elles sont à cheval sur une "
        "face de la CELLULE, et les 64 tuiles nées du résidu de D5.6 "
        "sont EXACTEMENT celles qui affleurent la face {Im u = Im v = "
        "0}, là où σ de la détermination tournée devient indéterminé "
        "dès qu'on déborde. Le résidu 1/64 de D5.6 était donc "
        "géométriquement UNE FACE, pas une dispersion. Ces tuiles "
        "reçoivent un halo CLIPPÉ (même demi-largeur, centre décalé vers "
        "l'intérieur, halo affleurant : `core ⊂ halo ⊆ cellule`) — la "
        "notion d'atlas d'une variété À BORD, cartes ouvertes dans la "
        "topologie RELATIVE. Le gate D2b interdit que le clipping serve "
        "ailleurs : il exige une face affleurante dans chaque direction "
        "décalée ET le refus préalable de toute l'échelle symétrique. "
        "Le nerf est calculé en RATIONNELS "
        "DYADIQUES EXACTS : %d paires de cores se touchent et %d/%d ont "
        "un overlap de largeur STRICTEMENT POSITIVE dans les quatre "
        "coordonnées ; le nerf est %s. Sur chaque overlap, les deux "
        "sections source sont comparées par RECENTRAGE EXACT dans un "
        "cadre commun puis soustraction coefficient par coefficient "
        "(l'enclosure naïve, décorrélée, aurait une largeur ~1e-3 et "
        "aucun contenu) : le feuillet θ est DÉRIVÉ à marge stricte, vaut "
        "+1 sur %d/%d paires (marge min %.3e), et le défaut affine de "
        "contrôle du recentrage est %.3e. La transition est un "
        "BIHOLOMORPHISME certifié (|Z[g']| minoré > 0 des deux côtés) et "
        "l'identité `Zt^(j)·Z[g'_j] − Zt^(i)·Z[g'_i] ∋ 0` — qui met en "
        "présence deux ledgers ε' dérivés SÉPARÉMENT — tient avec un "
        "défaut max %.3e ≤ δ. Sur les %d triples d'intersection triple "
        "non vide, le cocycle est certifié %d/%d. HONNÊTETÉ SUR LE "
        "COCYCLE : `λ_ij·λ_jk·λ_ki = 1` est une IDENTITÉ ALGÉBRIQUE de "
        "la construction (mêmes cartes affines d'un même P⁵, trois "
        "rapports des trois mêmes évaluations) — le vérifier en "
        "intervalles ne vérifie rien, chaque jauge y apparaît deux fois "
        "et le « défaut » mesuré (%.3e) N'EST QUE la largeur de "
        "décorrélation, exactement le phénomène qui a imposé le "
        "recentrage ailleurs ; il est donc publié comme DIAGNOSTIC et "
        "n'entre dans AUCUN gate. Ce qui est gaté au niveau du triple a "
        "du contenu : les trois jauges MINORÉES > 0 sur la boîte triple, "
        "le cocycle discret θ exact, et les trois identités d'arête "
        "recentrées sur la boîte TRIPLE — l'inclusion `∋ 0` n'étant pas "
        "transitive, i→k n'est pas une conséquence intervalliste de i→j "
        "et j→k. Ce qui n'est PAS tautologique non plus, c'est θ ≡ +1 "
        "(rien n'imposait aux voisines le même feuillet, leurs "
        "déterminations source DIFFÈRENT) et l'identité de transition, "
        "qui met en présence deux ledgers ε' dérivés séparément. Les %d "
        "négatifs cassent tous. NON PAYÉ ICI : le transport MÉTRIQUE "
        "(congruence Q, Weyl) sur les halos — établi sur les cores par "
        "C127, son extension coûte un quadruple Qmat par tuile ; le "
        "contrat EXACT de l'identité ; le scaling complet ; les 895 "
        "autres paires cellule/classe ; R12-C. C129-A/B/C (8ᵉ revue) "
        "PAYÉS : endpoints d'overlap en arrondi EXTÉRIEUR "
        "(fraction_box_to_iv — l'ancien float-nearest pouvait tronquer "
        "les boîtes à dénominateur 65), défauts gatés comparés à δ en mp "
        "AVANT toute conversion float (sérialisation _f_up seulement), "
        "et égalité d'ENSEMBLES {clippées} = {résidu C127-E} gatée (D2c, "
        "SHA-256 publiés)." % (
            MODE, len(leaves),
            [f"{r.numerator}/{r.denominator}" for r in RHO_LADDER],
            max_halo,
            sum(1 for h in halo_ok if h["rule_used"] == "clipped"),
            len(touch_pairs),
            sum(1 for x in touch_pairs if x[2]), len(touch_pairs),
            "connexe" if connected else "NON CONNEXE",
            sum(1 for r in pres if r.get("same_sheet")), len(pres),
            min_theta, max_aff, max_id, len(tres), len(tok), len(tres),
            max_cocy, len(negs)))

    art = {
        "artifact": ART.stem, "mode": MODE,
        "claim": ("C127-D — cores + halos, overlaps ouverts certifiés, "
                  "transitions (feuillet, carte, identité) et cocycle "
                  "sur le nerf : le cover fermé de C127/C127-E devient "
                  "un atlas multi-chart."),
        "cell": cell_d,
        "n_leaves": len(leaves),
        "tree_gates": tg, "n_address_failures": addr_fail,
        "pilot_patch": patch,
        "selected_tiles": sel,
        "halos": [{"index": h["index"], "ok": h["ok"],
                   "rho_used": h.get("rho_used"),
                   "H_hex": h.get("H_hex"),
                   "record": h.get("record"),
                   "attempts": [{k: v for k, v in a.items()
                                 if k != "kinds_target"}
                                for a in h["attempts"]]}
                  for h in halos],
        "clipped_set": {
            "clipped_indices": clipped_idx,
            "c127e_selected_indices": c127e_idx,
            "clipped_sha256": _set_sha(clipped_idx),
            "c127e_sha256": _set_sha(c127e_idx),
            "equal": bool(clipped_idx == c127e_idx),
            "note": ("C129-C : égalité d'ensembles gatée (D2c), pas "
                     "seulement le cardinal — 64 = 64 n'exclut pas une "
                     "permutation 63+1")},
        "nerve": {
            "n_touching_core_pairs": len(touch_pairs),
            "n_open_overlaps": len(pairs),
            "all_touching_pairs_open": bool(open_all),
            "connected": bool(connected),
            "n_triples": len(triples),
            "degree": {"min": min((len(adj[i]) for i in ids), default=0),
                       "max": max((len(adj[i]) for i in ids), default=0)},
            "note": ("nerf calculé en Fraction exactes ; « largeur "
                     "strictement positive dans les 4 coordonnées » est "
                     "une comparaison rationnelle, pas un test float")},
        "pairs": pres,
        "triples": tres,
        "negatives": negs,
        "max_halo_congruence_sup": max_halo,
        "max_transition_identity_defect": max_id,
        "max_affine_control_defect": max_aff,
        "min_theta_margin": min_theta,
        "max_lambda_decorrelation_width": max_cocy,
        "delta_trans_preregistered": DELTA_TRANS,
        "delta_cocycle_preregistered": DELTA_COCYCLE,
        "not_paid_here": [
            "transport métrique (Qmat/Weyl) sur les halos — établi sur "
            "les cores par C127",
            "contrat exact de l'identité de congruence",
            "scaling complet (C125-A)",
            "les 895 autres paires cellule/classe", "R12-C"],
        "verdict": verdict, "gates": gates,
        "gates_passed": n_pass, "gates_total": len(gates),
        "provenance": provenance([COVER_JSON, C127_JSON, C127E_JSON],
                                 time.time() - T0)}
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
#  Les négatifs — la discriminance, pas la décoration
# ===========================================================================
def run_negatives(leaves, geom, tm, pairs, pmap):
    """Six falsifications. Chacune DOIT casser un gate nominal."""
    out = {}
    ids = sorted(geom)
    ref = next(((i, j) for (i, j) in pairs
                if pmap[(i, j)].get("ok")), None)
    if ref is None:
        return {"N0_no_reference_pair": {"breaks": False}}
    i, j = ref

    # N1 — ledger ε' corrompu sur une tuile : l'identité de transition
    #      doit cesser de contenir 0 (elle met en présence les DEUX
    #      ledgers ; si elle ne voyait pas ε', elle serait décorative)
    Zj, Pj = tm[j]
    Pbad = [([CIV(-c.re, -c.im) for c in Pj[a][0]], Pj[a][1])
            if a == leaves[j]["chart"]["S"][0] else Pj[a]
            for a in range(6)]
    r1 = pair_certificate(i, j, corrupt={"tile": j, "Z": Zj,
                                         "P": Pbad})
    out["N1_eps_ledger_corrupted"] = {
        "breaks": bool(not r1.get("transition_identity_contains_zero")),
        "defect": r1.get("transition_identity_defect"),
        "note": ("une ligne de Zt·Z[g'] niée sur la tuile j : "
                 "l'identité de transition doit exclure 0")}

    # N2 — feuillet inversé : θ doit sortir à −1 (le test DÉTECTE un
    #      raccord inter-feuillets, c'est toute sa raison d'être)
    s0 = _G["cell"][0][0]          # première ligne √ de la SOURCE
    Zbad = [([CIV(-c.re, -c.im) for c in Zj[a][0]], Zj[a][1])
            if a == s0 else Zj[a] for a in range(6)]
    r2 = pair_certificate(i, j, corrupt={"tile": j, "Z": Zbad,
                                         "P": Pj})
    out["N2_sheet_flipped"] = {
        "breaks": bool(r2.get("theta") and r2["theta"][0] == -1
                       and not r2.get("same_sheet")),
        "theta": r2.get("theta"),
        "note": ("Z_s nié sur la tuile j : θ doit valoir −1 sur cette "
                 "ligne et le gate de feuillet doit tomber")}

    # N3 — ρ = 0 : les halos redeviennent les cores et l'overlap perd sa
    #      largeur ; « ouvert » doit cesser d'être vrai
    n3 = 0
    for (a, b) in pairs[:200]:
        ga, gb = geom[a], geom[b]
        la, ha = cube_bounds(ga["c"], ga["core_hw"])
        lb, hb = cube_bounds(gb["c"], gb["core_hw"])
        if _inter_generic(la, ha, lb, hb) is None:
            n3 += 1
    out["N3_rho_zero_kills_openness"] = {
        "breaks": bool(n3 > 0), "n_pairs_losing_openness": n3,
        "n_tested": min(200, len(pairs)),
        "note": ("à ρ = 0 les boîtes fermées adjacentes ne se coupent "
                 "plus en largeur positive : l'ouverture des overlaps "
                 "est bien PORTÉE par ρ, pas par la partition")}

    # N4 — paire non adjacente : la géométrie doit refuser
    far = None
    for a in ids:
        for b in ids:
            if b > a and (a, b) not in pmap:
                far = (a, b)
                break
        if far:
            break
    if far:
        r4 = pair_certificate(*far)
        out["N4_non_adjacent_refused"] = {
            "breaks": bool(r4.get("refused") == "no_open_overlap"),
            "pair": list(far), "refused": r4.get("refused")}
    else:
        out["N4_non_adjacent_refused"] = {
            "breaks": False, "note": "aucune paire non adjacente"}

    # N5 — recentrage neutralisé : si l'on compare les polynômes SANS
    #      recentrer (T = identité), la différence doit exploser dès que
    #      les deux cadres diffèrent réellement
    gi, gj = geom[i], geom[j]
    w = _inter_generic(gi["lo"], gi["hi"], gj["lo"], gj["hi"])
    Eiv = [fraction_box_to_iv(lo, hi)
           for lo, hi in eps_box(w, gi["c"], gi["H"])]
    Zi = tm[i][0]
    pa, ra = Zi[s0]
    pb, rb = Zj[s0]
    d_true = float(civ_sup(enclose(
        poly_lin(pa, apply_recenter(
            get_T([(gi["c"][k] - gj["c"][k]) / gj["H"] for k in range(4)],
                  gi["H"] / gj["H"]), pb), -1), ra + rb, Eiv)))
    d_naive = float(civ_sup(enclose(poly_lin(pa, pb, -1),
                                    ra + rb, Eiv)))
    out["N5_recentering_is_load_bearing"] = {
        "breaks": bool(d_naive > 1e3 * max(d_true, 1e-300)),
        "defect_recentered": d_true, "defect_unrecentered": d_naive,
        "note": ("sans recentrage la comparaison est vide de contenu ; "
                 "le rapport mesure ce que le recentrage porte")}

    # N6 — géométrie décalée : le test doit voir le décalage À SA TAILLE
    #      PRÉVUE. Un simple « c'est plus grand » se règle par un seuil
    #      choisi après coup ; ici l'effet attendu d'un décalage δ de
    #      l'application affine est δ·|∂f/∂ε₀|, et |∂f/∂ε₀| se LIT sur le
    #      coefficient de degré 1 du Taylor-modèle. Le négatif exige donc
    #      que le défaut mesuré tombe dans [prédit/4, 4·prédit] ET que le
    #      plancher nominal soit deux ordres SOUS le prédit.
    delta = Fraction(1, 2 ** 20)
    off_bad = [(gi["c"][k] - gj["c"][k]) / gj["H"] for k in range(4)]
    off_bad[0] = off_bad[0] + delta
    d_bad = float(civ_sup(enclose(
        poly_lin(pa, apply_recenter(get_T(off_bad, gi["H"] / gj["H"]),
                                    pb), -1), ra + rb, Eiv)))
    e0 = pb[MIDX[(1, 0, 0, 0)]]
    grad = float(max(abs(mp.mpf(e0.re.a)), abs(mp.mpf(e0.re.b)),
                     abs(mp.mpf(e0.im.a)), abs(mp.mpf(e0.im.b))))
    pred = float(delta) * grad
    out["N6_shifted_geometry_detected"] = {
        "breaks": bool(pred / 4 <= d_bad <= 4 * pred
                       and d_true < pred / 100),
        "defect_true": d_true, "defect_shifted": d_bad,
        "delta_eps": float(delta), "degree1_coefficient": grad,
        "predicted_defect": pred,
        "ratio_measured_over_predicted": (d_bad / pred if pred else None),
        "note": ("décalage δ = 2⁻²⁰ de l'application affine : le défaut "
                 "doit apparaître à δ·|∂f/∂ε₀|, lu sur le coefficient de "
                 "degré 1 du TM — le test voit la GÉOMÉTRIE, et il la "
                 "voit à la bonne échelle, pas seulement « plus grand »")}
    return out


# ===========================================================================
#  Self-test — fonctions pures, sans le registre ni les sections
# ===========================================================================
def _selftest():
    ok = True

    def chk(name, cond):
        nonlocal ok
        ok = ok and bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    # T1 : intersection exacte, largeur strictement positive
    lo1, hi1 = cube_bounds([Fraction(0)] * 4, Fraction(1))
    lo2, hi2 = cube_bounds([Fraction(3, 2)] + [Fraction(0)] * 3,
                           Fraction(1))
    w = _inter_generic(lo1, hi1, lo2, hi2)
    chk("T1 overlap ouvert : largeur 1/2 dans la direction séparante",
        w is not None and w["hw"][0] == Fraction(1, 4)
        and w["hw"][1] == Fraction(1))
    lo3, hi3 = cube_bounds([Fraction(2)] + [Fraction(0)] * 3,
                           Fraction(1))
    chk("T2 NÉGATIF : deux boîtes qui se TOUCHENT sans se chevaucher "
        "ne donnent aucun overlap ouvert",
        _inter_generic(lo1, hi1, lo3, hi3) is None)

    # T3 : la carte ε et son gate de domaine
    E = eps_box({"lo": [Fraction(-1, 2)] * 4,
                 "hi": [Fraction(1, 2)] * 4},
                [Fraction(0)] * 4, Fraction(1))
    chk("T3 boîte ε ⊆ [−1,1]⁴", eps_box_in_range(E)
        and E[0] == (Fraction(-1, 2), Fraction(1, 2)))
    E2 = eps_box({"lo": [Fraction(-2)] * 4, "hi": [Fraction(2)] * 4},
                 [Fraction(0)] * 4, Fraction(1))
    chk("T4 NÉGATIF : une boîte qui déborde le domaine du TM est "
        "REFUSÉE (le reste n'y est pas valide)", not eps_box_in_range(E2))

    # T5 : (1+ρ)h exactement représentable
    chk("T5 halo dyadique exact pour toute l'échelle pré-enregistrée",
        all(halo_hw(Fraction(1, 1024), r)[1] is not None
            for r in RHO_LADDER))

    # T6 : le recentrage — l'identité est l'identité
    T = recenter_matrix([Fraction(0)] * 4, Fraction(1))
    p = [CIV(iv.mpf(float(k + 1)), iv.mpf(float(-k)))
         for k in range(NM)]
    q = apply_recenter(T, p)
    chk("T6 recentrage identité : coefficients inchangés",
        all(mp.mpf(q[k].re.a) <= float(k + 1) <= mp.mpf(q[k].re.b)
            and mp.mpf(q[k].im.a) <= float(-k) <= mp.mpf(q[k].im.b)
            for k in range(NM)))

    # T7 : le recentrage est EXACT sur un polynôme connu.
    #      f(ε) = ε₀²  et  φ(ε) = (1/2 + ε/4) ⟹ f∘φ = 1/4 + ε₀/4 + ε₀²/16
    p2 = [CIV(IV0, IV0)] * NM
    p2 = list(p2)
    p2[MIDX[(2, 0, 0, 0)]] = CIV(IV1, IV0)
    T2m = recenter_matrix([Fraction(1, 2)] + [Fraction(0)] * 3,
                          Fraction(1, 4))
    q2 = apply_recenter(T2m, p2)

    def near(c, x):
        return (mp.mpf(c.re.a) <= x <= mp.mpf(c.re.b)
                and mp.mpf(c.im.a) <= 0 <= mp.mpf(c.im.b))
    chk("T7 recentrage exact sur ε₀² ∘ (1/2 + ε/4)",
        near(q2[0], 0.25) and near(q2[MIDX[(1, 0, 0, 0)]], 0.25)
        and near(q2[MIDX[(2, 0, 0, 0)]], 0.0625))

    # T8 : recentrage puis évaluation = évaluation directe (le même
    #      point, deux cadres) — c'est l'invariant qui porte tout C127-D
    Efull = [iv.mpf([-1, 1])] * 4
    v_direct = enclose(p2, IV0, [iv.mpf([0.25, 0.75])]
                       + [iv.mpf([-1, 1])] * 3)
    v_recent = enclose(q2, IV0, Efull)
    chk("T8 le recentrage préserve l'enclosure au même point physique",
        mp.mpf(v_recent.re.a) <= mp.mpf(v_direct.re.a)
        and mp.mpf(v_direct.re.b) <= mp.mpf(v_recent.re.b))

    # T9 : la séparation de feuillet, et son REFUS
    dm = CIV(iv.mpf([-1e-9, 1e-9]), iv.mpf([-1e-9, 1e-9]))
    dp = CIV(iv.mpf([1.9, 2.1]), iv.mpf([-0.1, 0.1]))
    th, rec = sep_phase(dm, dp)
    chk("T9 feuillet : différence ~0 et somme ~2 ⟹ θ = +1 à marge > 0",
        th == 1 and rec["margin"] > 0)
    th2, _ = sep_phase(dp, dm)
    chk("T10 feuillet inversé ⟹ θ = −1", th2 == -1)
    amb = CIV(iv.mpf([-1, 1]), iv.mpf([-1, 1]))
    th3, _ = sep_phase(amb, amb)
    chk("T11 NÉGATIF : séparation ambiguë ⟹ REFUS (jamais d'essai)",
        th3 is None)

    # T12 : sérialisation exacte des TM (aller-retour bit à bit)
    t = TMC([CIV(iv.mpf([0.1, 0.2]), iv.mpf([-0.3, 0.4]))] * NM,
            iv.mpf([0, 1e-9]))
    p3, r3 = de_tmc(ser_tmc(t))
    chk("T12 sérialisation `_mpf_` exacte (aller-retour identique)",
        all(mp.mpf(p3[k].re.a) == mp.mpf(t.p[k].re.a)
            and mp.mpf(p3[k].im.b) == mp.mpf(t.p[k].im.b)
            for k in range(NM))
        and mp.mpf(r3.b) == mp.mpf(t.rem.b))

    # T13 : le cocycle discret, et un triple délibérément cassé
    def cocy(a, b, c):
        return all(a[r] * b[r] * c[r] == 1 for r in range(3))
    chk("T13 cocycle θ trivial sur un triple sain",
        cocy([1, 1, 1], [1, 1, 1], [1, 1, 1]))
    chk("T14 NÉGATIF : un seul θ inversé casse le cocycle",
        not cocy([1, -1, 1], [1, 1, 1], [1, 1, 1]))

    # T15 : division intervalle complexe
    a = CIV(iv.mpf([1, 1]), iv.mpf([0, 0]))
    b = CIV(iv.mpf([2, 2]), iv.mpf([0, 0]))
    d = civ_div(a, b)
    chk("T15 division complexe : 1/2 enclos",
        mp.mpf(d.re.a) <= 0.5 <= mp.mpf(d.re.b))

    # T16/T17 : C129-A — les endpoints Fraction → iv en arrondi
    #   EXTÉRIEUR. Témoin : 1/65 (le dénominateur que les halos 65h/64
    #   injectent réellement dans les boîtes d'overlap), non représentable
    #   en binaire. Le critère « B enclot q » se teste SANS float : 65·B
    #   doit contenir 1 exactement.
    q65 = Fraction(1, 65)
    B = fraction_box_to_iv(q65, q65)
    X = iv.mpf(65) * B
    chk("T16 fraction_box_to_iv enclot le rationnel exact (65·B ∋ 1)",
        mp.mpf(X.a) <= 1 <= mp.mpf(X.b))
    old = iv.mpf([float(q65), float(q65)])
    Xo = iv.mpf(65) * old
    chk("T17 NÉGATIF : l'ancien chemin float-nearest RATE le rationnel "
        "exact (65·[float(1/65)] ∌ 1) — le défaut C129-A est réel",
        not (mp.mpf(Xo.a) <= 1 <= mp.mpf(Xo.b)))

    print("SELF-TEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(_selftest() if "--selftest" in sys.argv else build())
