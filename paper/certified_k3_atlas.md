# Certified Analytic Geometry on an Explicit K3 Surface

### Quantitative Charts, Branch Transitions, and Computer-Assisted Continuation

**Brieuc de La Fournière**

Independent researcher, Beaune, France · ORCID 0009-0000-0641-9740

**DOI**: [10.5281/zenodo.22047469](https://doi.org/10.5281/zenodo.22047469) · **Repository**: [github.com/Arithmon/K3](https://github.com/Arithmon/K3)

*Preprint. All quantitative hypotheses of this paper are certified by machine-checkable
artifacts; the verification entry point is a single command, described in Appendix F.*

---

## Abstract

We equip an explicit algebraic K3 surface with a finite holomorphic atlas
whose chart domains, transition maps and branch continuations are certified
rather than asserted. The surface is the smooth complete intersection of three
diagonal quadrics in $\mathbf{P}^5$ with integer coefficients. Exploiting the
diagonal quadratic structure -- for which the Taylor remainder is a single
explicit term and the relevant Hessian is constant -- we prove a uniform
quantitative chart lemma with an explicit guaranteed radius, by a
self-contained Newton contraction whose uniqueness step is an exact identity
rather than an estimate. Sixty explicitly specified chart types cover the
surface; on every nonempty overlap the transitions are restrictions of exact
algebraic coordinate changes satisfying the cocycle law identically, and they
carry the branch bookkeeping with them. We identify the chart pivots as the
Plücker coordinates of the row space of the Jacobian, which places the
combinatorics of chart changes inside $\mathrm{Gr}(3,6)$, and we determine how
a global sign automorphism acts chart by chart -- vertically on some types,
moving the base point on others. A validated computation layer, exact over
$\mathbf{Q}$ and outward-rounded elsewhere, certifies the quantitative
hypotheses together with analytic continuations across chart boundaries; there
the continued sheet and the conjugated one are shown to differ by an explicit
deck transformation. No Ricci-flat or hyperkähler metric is claimed. The
result is a reproducible bridge between an explicit projective model of a K3
surface and theorem-grade local analytic data, suitable for subsequent
validated geometric analysis.

---

## 1. Introduction

### 1.1 The gap between explicit algebraic models and usable analytic geometry

A K3 surface given by explicit projective equations is "explicit" only
algebraically. The equations alone do not provide: chart domains with
controlled radius; a constructive chart rule valid at every point; transition
maps that are effectively verified rather than asserted; robust tracking of
root branches; or objects directly consumable by rigorous analytic
computation. This paper closes that gap for one concrete surface.

None of this is automatic. A qualitative application of the implicit
function theorem yields charts but no certified numerical radius for this
instance; computer algebra yields exact local data with no analytic control;
a numerical atlas yields neither exactness nor guarantees. What is missing in
each case is the *quantitative* statement: a domain of definite size on which
a named holomorphic map is proved to exist, transitions that are verified
rather than assumed, and a branch of every square root that is followed
rather than hoped for. No single ingredient below is new; what is new is
assembling them into one globally compatible object whose hypotheses are
machine-checkable and whose verification a third party can re-run.

### 1.2 Main result

> **Theorem A (Certified finite holomorphic atlas).** Let
> $X \subset \mathbf{P}^5$ be the explicit smooth complete intersection of
> three diagonal quadrics defined in Section 2. Then $X$ carries a finite
> holomorphic atlas assembled from $60$ explicitly specified chart types,
> every local chart of every type admitting one and the same certified lower
> bound on its radius $\rho_{\mathrm{unif}} > 0$
> ($\rho_{\mathrm{unif}} \ge 2.10 \times 10^{-12}$ in the normalization of
> Section 3.5; not sharp, strictly positive, reproducible). On every nonempty
> overlap, the transition maps are restrictions of explicit algebraic
> coordinate transformations and satisfy the cocycle identities exactly. The
> resulting glued analytic space is canonically biholomorphic to $X$.

The atlas is *finite by compactness* and *quantitative by certificate*: the
theorem asserts quantitative existence with a guaranteed radius, not an
enumerated list of charts. (The distinction matters and is kept throughout;
see Section 4.2.)

### 1.3 Contributions

**(A)** A completely explicit K3 model and its branched-cover structure
(Section 2). **(B)** A uniform quantitative chart lemma giving true chart
domains with certified radius, whose proof is self-contained and exact for
quadrics (Section 3). **(C)** Global atlas closure with exact transition and
cocycle laws and explicit branch bookkeeping (Sections 4-6). **(D)** A
reproducible computer-assisted architecture certifying the quantitative
hypotheses without turning the theorem into brute-force enumeration
(Sections 7-8).

**Scope.** We do not claim an explicit Ricci-flat K3 metric, nor a global
hyperkähler family. Metric computations appear only as downstream validated
probes (Section 9) and are clearly separated from the atlas theorem.

---

## 2. The Explicit K3 Surface

### 2.1 A diagonal CI(2,2,2) model

Let $\mu = (\mu_0, \dots, \mu_5) = (1, 2, 3, 5, 7, 11)$ and, for
$k = 0, 1, 2$,

$$F_k(z) \;=\; \sum_{j=0}^{5} \mu_j^{\,k}\, z_j^2 .$$

Explicitly:

$$
\begin{aligned}
F_0 &= z_0^2 + z_1^2 + z_2^2 + z_3^2 + z_4^2 + z_5^2,\\
F_1 &= z_0^2 + 2 z_1^2 + 3 z_2^2 + 5 z_3^2 + 7 z_4^2 + 11 z_5^2,\\
F_2 &= z_0^2 + 4 z_1^2 + 9 z_2^2 + 25 z_3^2 + 49 z_4^2 + 121 z_5^2,
\end{aligned}
$$

and $X = \{F_0 = F_1 = F_2 = 0\} \subset \mathbf{P}^5$. The reader can copy
the three equations onto a sheet of paper; every constant in the paper is
derived from these six integers.

### 2.2 Smoothness and the K3 property

The Jacobian criterion reduces to Vandermonde nonvanishing: for distinct
$\mu_j$ the $3 \times 3$ minors $V_S = \det(\mu_s^k)_{k, s \in S}$ are
nonzero for every triple $S$, and smoothness follows.

> **Lemma 2.1 (at least three nonzero coordinates).** At every point of $X$,
> at least three of the six homogeneous coordinates are nonzero.

*Proof.* Suppose at most two are nonzero, say $z_i$ and $z_j$. Writing
$w = z^2$, the three equations read $W \cdot (w_i, w_j)^{\mathsf T} = 0$ with
$W = (\mu_i^k, \mu_j^k)_{k=0,1,2}$. Every $2 \times 2$ minor of $W$ is a
Vandermonde determinant in two distinct $\mu$'s, hence nonzero (checked
exactly on all $15$ pairs), so $W$ has rank $2$ and $w_i = w_j = 0$. Then all
coordinates vanish, which is not a point of $\mathbf{P}^5$. $\square$

Smoothness follows at once: for a triple $S$ of nonzero coordinates the
corresponding Jacobian minor is the algebraic pivot
$\det M_S^{\mathrm{alg}} = 8 V_S z_{s_1} z_{s_2} z_{s_3}$ (Section 3.5), which is nonzero, so the Jacobian has rank $3$
everywhere and $X$ is a smooth surface. Adjunction gives
$K_X \simeq \mathcal{O}_X$, and the Lefschetz hyperplane theorem gives
$H^1(X, \mathbf{Z}) = 0$, hence $h^{1,0} = h^{0,1} = 0$; so $X$ is K3. This
part is entirely classical and involves no computation. Lemma 2.1 will be used
a second time, for a different purpose, in Section 4.1.

### 2.3 Projection and branched-cover structure

Choosing three base coordinates and solving the three quadratic equations for
the squares of the other three exhibits, on each elimination patch,

$$z_{s_i}^2 = R_i(t_1, t_2, t_3), \qquad i = 1, 2, 3,$$

with $R_i$ explicit rational quadratic expressions (Cramer over the
Vandermonde pivot). This gives branched projections
$\pi_S : X \dashrightarrow \mathbf{P}^2$ with generic fibre the $2^3$ sign
choices, a natural sheet action of $(\mathbf{Z}/2)^3$ per patch, and the
global diagonal sign action on $X$ itself. Branch divisors are cut by
$z_s = 0$.

For the distinguished patch $S = \{3,4,5\}$, $T = \{0,1,2\}$, gauge
$z_0 = 1$ and base $(u, v) = (z_1, z_2)$, Cramer over the Vandermonde pivot
($|V_S| = 48$) gives the three radicands exactly:

$$
R_3 = -5 - \tfrac{15}{4} u^2 - \tfrac{8}{3} v^2, \qquad
R_4 = 5 + \tfrac{27}{8} u^2 + 2 v^2, \qquad
R_5 = -1 - \tfrac{5}{8} u^2 - \tfrac{1}{3} v^2 .
$$

A remark that shapes the whole paper: since $F_0 = \sum_j z_j^2$ is positive
definite, $X(\mathbf{R}) = \varnothing$ -- the surface has no real points at
all for this $\mu$. Accordingly the three radicands above have *constant sign*
on the real $(u,v)$ slice, and the branch loci are invisible there. The
regime structure that the continuation machinery of Section 7 navigates lives
in the complex directions: in the slice $v = 0$ of the complex $u$-plane the
loci $\mathrm{Re}\, R_s = 0$ are three nested rectangular hyperbolas
$x^2 - y^2 = -a_s/b_s$, and crossing one flips exactly one solved line between
the principal and the canonical determination (Figure 1).

The deck symmetries are the diagonal sign matrices
$\sigma_j = \mathrm{diag}(1, \dots, -1, \dots, 1)$ (the $-1$ in slot $j$):
each preserves every $F_k$, since the coordinates enter only through their
squares. Modulo the projective scalar $-1$ they generate a group of order
$32$, and $\sigma_1, \dots, \sigma_5$ suffice, because
$\sigma_0 = -\sigma_1\sigma_2\sigma_3\sigma_4\sigma_5$ in $\mathbf{P}^5$.

The branch locus of $\pi_S$ is $\{R_{s_1} R_{s_2} R_{s_3} = 0\}$. For the
distinguished patch this discriminant is the explicit degree-six polynomial

$$
\tfrac{2025}{256}u^6 + \tfrac{465}{32}u^4v^2 + \tfrac{2235}{64}u^4
+ \tfrac{53}{6}u^2v^4 + \tfrac{1031}{24}u^2v^2 + \tfrac{205}{4}u^2
+ \tfrac{16}{9}v^6 + \tfrac{118}{9}v^4 + \tfrac{95}{3}v^2 + 25 .
$$

This explains *before any computation* why sign patterns exist and why every
chart carries branch bookkeeping.

### 2.4 Chart notation

A chart type is a pair $(S, g)$: a triple $S \subset \{0, \dots, 5\}$ of
solved coordinates ($20$ triples) and a gauge coordinate $g \notin S$
normalized to $z_g = 1$ ($3$ choices per triple after the selector; $60$
types in total). Solved variables $w = z_S \in \mathbf{C}^3$; base variables
$b = (u, v) \in \mathbf{C}^2$; sheet sign pattern $\varepsilon \in \{\pm 1\}^3$.
Two layers must be kept apart, and the paper keeps them apart throughout.
The **algebraic layer** uses the equations as written: $V_S \in M_3(\mathbf{Z})$
and every identity of Sections 2, 4 and 5 lives over $\mathbf{Q}$. The
**normalized layer** rescales rows, $\tilde F_k = F_k / c_k$ with
$c_k^2 = \sum_j \mu_j^{2k}$; this is a conditioning convention for the
certificate, not part of the geometry. The *squares* $c_k^2$ are integers; the
$c_k$ themselves are irrational, and the certificate handles them as algebraic
numbers by outward-rounded rational bracketing.
*Worked example.* Take $S = \{3,4,5\}$ and $g = 0$: the solved variables are
$w = (z_3, z_4, z_5)$, the gauge is $z_0 = 1$, the base is
$(u, v) = (z_1, z_2)$, and a point of the chart is
$z_s = \varepsilon_s \sqrt{R_s(u,v)}$ for the radicands of Section 2.3 and a
sign pattern $\varepsilon \in \{\pm 1\}^3$. The normalization constants have exact integer
squares, $c_0^2 = 6$, $c_1^2 = 209$, $c_2^2 = 17765$ (the $c_k$ themselves
are irrational -- see the two layers above), and the Vandermonde
pivot of this triple is $|V_S| = 48$. Every constant appearing later in the
certificate for this chart is a rational function of these.

---

## 3. Quantitative Local Charts

### 3.1 Local elimination

For an admissible type $(S, g)$ write $G(w; b) := \tilde F(w, b, 1)$. Two
Jacobians must be kept apart, and confusing them is easy because they differ
only by the row normalization of Section 2.4:

$$M_S^{\mathrm{alg}} := \partial_w F, \qquad
\tilde M_S := \partial_w \tilde F = \partial_w G,
\qquad
\det M_S^{\mathrm{alg}} = c_0 c_1 c_2 \cdot \det \tilde M_S .$$

**Admissibility is the algebraic pivot condition**

$$q_S(Z) \;:=\; \frac{|\det M_S^{\mathrm{alg}}(Z)|}{\lVert Z \rVert^3}
\;>\; m_{\mathrm{alg}}, \qquad m_{\mathrm{alg}} = 4$$

($m_{\mathrm{alg}}$ preregistered, below the observed floor $4.8$), stated
once and for all on the **projective invariant** $q_S$ -- equivalently, on the
unit-norm representative $\lVert Z \rVert = 1$, which is the slice the
coverage certificate enumerates. Since $\det M_S^{\mathrm{alg}}$ is homogeneous
of degree $3$ (Section 3.5), the bare inequality
$|\det M_S^{\mathrm{alg}}| > 4$ would depend on the choice of representative,
whereas the vanishing condition $\det M_S^{\mathrm{alg}} \neq 0$ of Lemma 4.1
does not, and is used as such. Throughout, $U_S$ denotes the open
$\{q_S > m_{\mathrm{alg}}\}$. It is the
algebraic pivot that enters the coverage argument of Section 4.1 and the
coordinate floor of Lemma 3.1, both of which use the integer factorization
$\det M_S^{\mathrm{alg}} = 8 V_S z_{s_1} z_{s_2} z_{s_3}$ with
$V_S \in \mathbf{Z}$. The **Newton estimates of Section 3.4 use the
row-normalized $\tilde M_S$**, whose smallest singular value is the $\sigma$
of Proposition 3.2. The same symbol is never used for both.



### 3.2 Why the present model is unusually tractable

For diagonal quadrics the following identities are *exact*, not estimates:

* (E1) $\partial_w G(w; b) = \tilde M(w) = 2 \tilde V_S \,\mathrm{diag}(w)$,
  independent of $b$;
* (E2) $G(w; b) - G(w_0; b) - \tilde M(w_0)(w - w_0)
  = \tilde V_S \cdot ((w - w_0)^2)$ (componentwise square): the Taylor
  remainder is a single explicit quadratic term, the "Hessian" is constant;
* (E3) $G(w_0; b) - G(w_0; b_0) = \tilde V_{uv} \cdot (b^2 - b_0^2)$.

No numerical approximation of any higher derivative is ever needed. This is
the structural reason the algebraic layer closes over $\mathbf{Q}$, the row
normalization being the only place where algebraic numbers enter -- and there
they are bracketed, never approximated.

### 3.3 The pivot-to-coordinate floor

The condition bounds a determinant; the Newton estimates need a bound on the
*coordinates*. Bridging the two is a short argument that will be used twice --
here for $\sigma$, and in Section 4.3 for the distinctness of sheet sign patterns.
Sheets of $\pi_S$ do merge, exactly where a solved coordinate vanishes, so a
uniform bound cannot hold on all of $X$ and has to be earned on the domains
actually used.

> **Lemma 3.1 (pivot-to-coordinate floor).** Let $a_s, b_s, c_s$ be the
> rational coefficients of the radicand $R_s = a_s + b_s u^2 + c_s v^2$ and
> set $B_s := (|a_s| + |b_s| + |c_s|)^{1/2}$. On a chart domain
> $U_S \cap C_{S,g}$ every solved coordinate is bounded below:
> $$|z_{s_i}| \;>\; \max\!\left\{\;
> \frac{m_{\mathrm{alg}}}{8\,|V_S| \prod_{j \neq i} B_{s_j}}, \;\;
> \big(|a_{s_i}| - |b_{s_i}| - |c_{s_i}|\big)^{1/2} \;\right\} \;>\; 0,$$
> the second entry being used only when the quantity under the root is
> positive. Over the $60$ chart types the worst such bound is
> $9.5711 \times 10^{-4}$ (the certified value $9.571170\ldots \times 10^{-4}$
> rounded outward), attained at $S = \{3,4,5\}$, gauge $z_0$, solved
> coordinate $z_5$.
> Consequently two sign patterns
> differing in slot $s$ give points at distance at least twice that in the
> $s$-th coordinate: the eight sheets are uniformly separated on certified
> domains.

*Proof.* First the threshold has to be turned into a statement about the
representative actually used. In the gauge $z_g = 1$ we have
$\lVert Z \rVert \ge 1$, so
$q_S(Z) > m_{\mathrm{alg}}$ gives
$|\det M_S^{\mathrm{alg}}(Z)| = q_S(Z)\,\lVert Z \rVert^3 > m_{\mathrm{alg}}$:
the projective condition implies the plain inequality on the gauge representative,
which is the form the factorization consumes. From there, the pivot threshold
bounds a **product**, not its factors: by the factorization of Section 3.5,
$|\det M_S^{\mathrm{alg}}| = 8|V_S| \prod_i |z_{s_i}| > m_{\mathrm{alg}}$ gives
$\prod_i |z_{s_i}| > m_{\mathrm{alg}} / (8|V_S|)$. The sector selector supplies the
missing upper bounds: $|u|, |v| \le 1$, so
$|z_s|^2 = |R_s(u,v)| \le |a_s| + |b_s| + |c_s| = B_s^2$. Dividing the product
bound by the two upper bounds $B_{s_j}$, $j \neq i$, gives the first entry.
The second is independent of the condition: the same triangle inequality read from
below gives $|R_{s_i}| \ge |a_{s_i}| - |b_{s_i}| - |c_{s_i}|$, which is the
sharper bound on the lines where the constant term dominates. All constants
are exact rationals with outward-rounded roots, and the worst case is computed
over all $60$ types. The bound does not depend on the gauge: changing $g$
permutes the columns of $w_T$, which leaves each row sum
$|a_s| + |b_s| + |c_s|$ unchanged, so the worst case is attained at the three
gauges of $T = \{0,1,2\}$ alike and the certificate records the first of
them. $\square$

So the pivot threshold is precisely the guard that keeps every chart away from
its own branch locus.


### 3.4 A quantitative inverse-function lemma

> **Proposition 3.2 (Uniform quantitative chart lemma).** Fix a type
> $(S, g)$, a lower bound $\sigma \le s_{\min}(\tilde M(w_0))$ -- supplied by
> $s_{\min}(\tilde M(w_0)) \ge 2\, s_{\min}(\tilde V_S) \cdot
> \min_{s \in S} |z_s|$ together with Lemma 3.1 -- and the
> exact quadratic remainder bound (E2) with constant
> $L = 2\,\sigma_{\max}(\tilde V_S)$. For every admissible centre the
> simplified Newton iteration for $G(\cdot\,; b) = 0$ contracts on a ball of
> explicit radius: with residual $\eta \le a_{S,g}\,\rho(2 + \rho)/\sigma$,
> where $a_{S,g}$ bounds the norm of the base-variable block $\tilde V_{uv}$
> of the normalized Vandermonde and $|b_0| \le 1$ by the sector selector
> of Section 3.5, and
> $h = L\eta/\sigma < 1/2$, there is a unique zero at distance
> $r^* = (1 - \sqrt{1 - 2h}\,)\,\sigma / L \le 2\eta$, unique in the open
> ball $\lVert w - w_0 \rVert < \sigma/L$ by the midpoint identity
> $0 = \tilde M\big((w + w')/2\big)(w - w')$, which is exact for quadrics.
> The local solution is holomorphic in $b$: with the centre held fixed and
> $G$ quadratic, every simplified-Newton iterate is a polynomial in the base
> variables $b$, and the limit is uniform. All hypotheses are verified by exact rational checks on all 60
> types.

The proof is self-contained (no external quotable constant is load-bearing)
and is given in Appendix B.

The floor of Lemma 3.1 and the displacement bound of Proposition 3.2 close,
together, the step that obligation **O5** of Section 4.3 consumes:

> **Corollary 3.3 (charts avoid their own branch locus).** On a chart domain
> the solved coordinates stay bounded away from zero *throughout the certified
> ball*, not merely at its centre. Indeed the branch locus of $\pi_S$ is
> $\{R_{s_1} R_{s_2} R_{s_3} = 0\}$, that is, the vanishing of a solved
> coordinate; Lemma 3.1 bounds each one below by $\underline{z}_{S,g}$ at the
> centre, and Proposition 3.2 confines the solution to
> $\lVert w - w_0 \rVert \le r^* \le 2\eta$, so every point of the chart
> satisfies $|z_{s_i}| \ge \underline{z}_{S,g} - 2\eta$. Over the $60$ types
> the ratio $2\eta / \underline{z}_{S,g}$ is at most
> $1.72 \times 10^{-2}$, attained at $S = \{0,3,5\}$, $g = 1$
> ($2\eta = 2.22 \times 10^{-5}$ against a floor of
> $1.29 \times 10^{-3}$): the displacement is two orders of magnitude below the
> floor it would have to cross. Hence every certified chart domain lies in
> $X^\circ$ for its own projection, the eight signed sheets are distinct on it,
> and the sign pattern is unambiguous.

The constants $\eta$ and $\underline{z}_{S,g}$ are both serialized per type in
the chart certificate, so the ratio above is read off the artifact rather than
recomputed by hand.

### 3.5 Projective normalization

The radius must not depend on a naive choice of projective representative.
The gauge is fixed by $z_g = 1$ on the affine open $A_g = \{z_g \neq 0\}$,
and centres are selected in the closed sector
$C_{S,g} = \{|z_g| = \max_t |z_t|\}$ (selector only; it guarantees
$|u|, |v| \le 1$ and is never used as an open cover).
The homogeneity is exactly right for this. Since
$\partial \tilde F_k / \partial z_s = 2(\mu_s^k/c_k) z_s$, we have
$\tilde M_S = 2 \tilde V_S \, \mathrm{diag}(z_S)$ and therefore the exact
factorization

$$\det M_S^{\mathrm{alg}} \;=\; 8\, V_S\, z_{s_1} z_{s_2} z_{s_3},
\qquad V_S \in \mathbf{Z},$$

for the unnormalized equations -- an identity entirely over $\mathbf{Z}$ --
and $\det \tilde M_S = \det M_S^{\mathrm{alg}} / (c_0 c_1 c_2)$ after row
normalization. Both were verified symbolically on all $20$ triples. It is homogeneous of degree $3$ in
$Z$, so

$$q_S(Z) \;=\; \frac{|\det M_S^{\mathrm{alg}}(Z)|}{\lVert Z \rVert^3}$$

-- the quantity on which Section 3.1 states the admissibility condition --
is invariant under $Z \mapsto \lambda Z$: the certified radius is read off a
genuine projective invariant, not off a choice of representative. (Replacing
the exponent $3$ by $2$ destroys the invariance -- this is one of the negative
controls.) The sector $|z_g| = \max_t |z_t|$ selects a representative; it is
never used as an open set of a cover.

### 3.6 The certified radius

| Quantity | Exact / certified bound |
| --- | ---: |
| algebraic pivot threshold $m_{\mathrm{alg}}$ (preregistered) | $4$ (observed floor $4.8$) |
| chart types | $60$ ($20$ triples $\times$ selector) |
| $\sigma_{\min}(\tilde V_S)$ | $\ge |\det \tilde V_S| / \lVert \tilde V_S \rVert_F^2$ (exact rational) |
| Lipschitz constant $L$ | $\le 2\,\sigma_{\max}(\tilde V_S) \le 2 \lVert \tilde V_S \rVert_F$ |
| contraction certificate | $h < 1/2$ on all $60$ types |
| per-type radius $\rho_{S,g}$ | $\in [2.10 \times 10^{-12},\ 1.71 \times 10^{-4}]$ (outward-rounded) |
| guaranteed uniform radius $\rho_{\mathrm{unif}} := \min_{S,g} \rho_{S,g}$ | $\ge 2.10 \times 10^{-12}$ |

The uniform radius is a Frobenius/determinant bound: strictly positive and
reproducible, deliberately not sharp. A uniform radius is a single number, not an interval: the *per-type* radii
$\rho_{S,g}$ range over eight decades (Figure 2) and $\rho_{\mathrm{unif}}$ is
their minimum. The per-type values are in the certificate artifact.

---

## 4. From Local Charts to a Global Atlas

### 4.1 Global admissibility

Every point of $X$ admits at least one admissible triple: the pivot opens
$U_S = \{q_S > m_{\mathrm{alg}}\}$ cover $X$ (coverage margin $\tau = 0.6$; the zero
set of all pivots is empty). Coverage is *qualitatively* a two-line consequence of Lemma 2.1. By the
factorization of Section 3.5, $\det M_S^{\mathrm{alg}} \neq 0$ if and only if the three
solved coordinates are all nonzero; Lemma 2.1 provides such a triple at every
point; hence

> **Lemma 4.1 (coverage).** The nonvanishing pivot opens
> $W_S := \{\det M_S^{\mathrm{alg}} \neq 0\}$, $S$ ranging over the $20$
> triples, cover $X$. (Vanishing is a projective condition, so $W_S$ needs no
> normalization; the symbol $U_S$ is reserved throughout for the *quantitative*
> open $\{q_S > m_{\mathrm{alg}}\} \subset W_S$.)

What the large enumeration audit adds is not the covering -- it is the
*quantitative floor*: a positive threshold $m_{\mathrm{alg}}$ (preregistered at
$m_{\mathrm{alg}} = 4$, against an observed minor floor $4.8$) and a coverage margin $\tau = 0.6$, so
that the opens $\{q_S > m_{\mathrm{alg}}\}$ still cover, with room to spare. That
audit is an independent certificate, referenced by hash from the chart
certificate, and its box enumeration belongs in the supplement.

### 4.2 Uniformity and finiteness

The chart types form a finite family (60) with uniformly controlled
constants; compactness of $X$ extracts a finite subcover of the local chart
balls. The theorem is *quantitative existence with guaranteed radius*: no
enumerated atlas is claimed, and the certification architecture explicitly
forbids the stronger wording (a firewall check turns red if "explicit
enumerated atlas" is asserted).

### 4.3 The atlas theorem, and the glued space as a corollary

The charts of Proposition 3.2 are parametrizations of open subsets of one and
the same $X \subset \mathbf{P}^5$, and on every overlap their transitions are
the ambient projective coordinate identifications of Section 5. Nothing has to
be reconstructed:

> **Theorem 4.2 (Global certified atlas).** The local parametrizations of
> Proposition 3.2 cover $X$ (Lemma 4.1), and on every nonempty overlap their
> transition maps agree with the ambient projective coordinate
> identifications and satisfy the cocycle law exactly (Section 5). They
> therefore constitute a finite holomorphic atlas on $X$, all of whose charts
> carry the uniform certified radius $\rho_{\mathrm{unif}} > 0$.

The gluing statement is then a corollary rather than the mechanism of proof:

> **Corollary 4.3.** The analytic space
> $\mathcal{X}_{\mathrm{atlas}} = \big(\bigsqcup_\alpha U_\alpha\big)/\sim$
> obtained by gluing these parametrizations along the certified transitions is
> canonically biholomorphic to $X$:
> $$\boxed{\;\mathcal{X}_{\mathrm{atlas}} \;\cong_{\mathrm{bihol}}\; X.\;}$$

Reading the corollary as a construction from data -- hand someone the chart
types, the radii and the transition formulas, and they rebuild $X$ -- is what
makes it worth stating; it is not what carries Theorem 4.2.

For the corollary the obligations are the standard ones, and it is worth
recording which are algebra and which consume certified data.

| | Obligation | Proof mode | Support |
| --- | --- | --- | --- |
| **O1** | $\sim$ is an equivalence relation | exact algebra | the cocycle of Proposition 5.1, together with $P_{SS} = \mathrm{Id}$ and $P_{S'S}^{-1} = P_{SS'}$ |
| **O2** | $\Phi$ is well defined | exact algebra | every chart section solves $F_0 = F_1 = F_2 = 0$ by construction (elimination), hence lands in $X$ |
| **O3** | $\Phi$ is a local biholomorphism | certificate | Proposition 3.2: certified radius, contraction $h < 1/2$ |
| **O4** | $\Phi$ is surjective | certificate | the pivot opens cover $X$ (Lemma 4.1) |
| **O5** | $\Phi$ is injective | exact algebra + sign pattern | on a chart domain the three solved coordinates are bounded away from zero, so the eight signed sheets are *distinct there* and the sign pattern separates them; $\sim$ therefore captures all coincidences of chart points |
| **O6** | $\mathcal{X}_{\mathrm{atlas}}$ is second countable and Hausdorff | topology | the atlas is finite by compactness (Section 4.2); a continuous bijection that is a local homeomorphism is open, hence a homeomorphism, so $\mathcal{X}_{\mathrm{atlas}}$ inherits Hausdorffness from $X$ |

Obligation **O5** uses the pivot-to-coordinate floor of Lemma 3.1: on a chart
domain the three solved coordinates are bounded away from zero, so the eight
signed sheets are distinct there. The sheets that merge over $B$ are not lost --
they are covered by charts of a different triple, which Lemma 2.1 always
provides.

Three obligations are exact algebra, two are certificates, one is
topological. Note that Hausdorffness is *not* bought with a separation
estimate: once $\Phi$ is a bijective local homeomorphism it is open, and the
quotient inherits the topology of $X$. Lemma 3.1 is therefore not needed
for O6 -- it is a quantitative lemma in its own right, used for the sign pattern
distinctness of **O5** and for the lower bound on $\sigma$ in Section 3.

### 4.4 Topology as an audit, not the proof

Homology of the nerve is used only as an independent consistency control.
Homology matching is *not* used to identify the glued analytic space with
$X$; the identification is the canonical map $\Phi$ of Theorem 4.2.

---

## 5. Exact Transition Maps

### 5.1 Algebraic chart changes

Between two types $(S, g)$ and $(T, h)$ the transition is the restriction of
an exact algebraic coordinate change: a permutation of coordinates, a
projective rescaling by $z_h/z_g$, and a re-selection of solved roots.
Two closed formulas make this concrete. An **elementary swap**
$S = \{3,4,5\} \to S' = \{2,4,5\}$ acts on the squares by

$$w_2 = -6 w_3 - 15 w_4 - 45 w_5, \qquad w_4 = w_4, \qquad w_5 = w_5,$$

that is, by a matrix whose last two rows are those of the identity: exactly
one row changes, the one carrying the newly solved coordinate. On the
coordinates themselves the transition is
$z_{s'} = \varepsilon_{s'} \sqrt{(P_{S'S} w_S)_{s'}}$, so the sheet sign pattern
$\varepsilon$ is part of the data of the map -- this is the branch-crossing
content, and it is why transitions are recorded together with their sign patterns.
A **gauge change** $(S, g) \to (S, g')$ is the projective rescaling by
$z_g / z_{g'}$ and touches no sign pattern.

### 5.2 Exact cocycle law

Underlying the chart transitions there is a purely algebraic cocycle, at the
level of the squares. Since the quadrics are linear in $w_j = z_j^2$, the
constraint $V_{\mathrm{full}}\, w = 0$ determines all six squares from any
admissible triple; hence for every ordered pair $(S, S')$ there is an exact
rational $3 \times 3$ matrix $P_{S'S}$ with $w_{S'} = P_{S'S}\, w_S$ on $X$.
These satisfy $P_{SS} = \mathrm{Id}$ (20/20), $P_{SS'} P_{S'S} = \mathrm{Id}$
(400/400) and

$$P_{S''S'}\, P_{S'S} \;=\; P_{S''S} \qquad\text{on all } 20^3 = 8000
\text{ ordered triples, in exact rational arithmetic.}$$

On every triple overlap of charts,
$\Phi_{\gamma\beta} \circ \Phi_{\beta\alpha} = \Phi_{\gamma\alpha}$
*exactly* -- these are algebraic identities in exact arithmetic, not
numerical statements with tolerance. The certification layer verified the
full transfer/cocycle panel ($3540/3540$ identities). The paper carefully
distinguishes exact algebraic identity from certified numerical
applicability on a given domain.

### 5.3 Finite transition generators

The twenty solved triples are not a combinatorial device imported from
outside: they are the standard coordinate charts of a Grassmannian. The
Jacobian $dF = 2\,V \mathrm{diag}(z)$ is a $3 \times 6$ matrix, so its row
space is a point of $\mathrm{Gr}(3,6)$, and its $S$-minor is exactly the pivot

$$\det M_S^{\mathrm{alg}} \;=\; 8\, V_S\, z_{s_1} z_{s_2} z_{s_3}$$

(verified on all $20$ minors). **The pivots are the Plücker coordinates of the
row space of $dF$.** The solved triples are therefore the standard coordinate
charts of $\mathrm{Gr}(3,6)$, pulled back along
$p \mapsto \operatorname{rowspace} dF(p)$, and an elementary swap
$|S \cap S'| = 2$ is precisely the standard chart change there; the Johnson
graph $J(6,3)$ is the adjacency graph of those elementary chart changes. (It
is not the nerve of the cover: any two standard charts of $\mathrm{Gr}(3,6)$
meet, so that nerve would be the full simplex on $20$ vertices.)

With that reading, we exhibit a small generating set. Two elementary moves act
on chart types $(S, g)$:

* an **elementary swap** $S \to S'$ with $|S \cap S'| = 2$ (one solved
  coordinate is exchanged against one base coordinate);
* a **gauge change** $(S, g) \to (S, g')$ with $g' \in T \setminus \{g\}$
  (projective rescaling by $z_g/z_{g'}$: the source gauge has $z_g = 1$,
  so reaching $z_{g'} = 1$ divides the representative by $z_{g'}$).

Together with the sheet flips of Section 6.2 these generate everything, with
an explicit bound on word length:

> **Proposition 5.1 (finite generation with length bound).** The $20$ solved
> triples with adjacency $|S \cap S'| = 2$ form the Johnson graph $J(6,3)$:
> $20$ nodes, $90$ edges, connected, of diameter $3$. The $60$ chart types
> under $\{$elementary swap, gauge change$\}$ -- the swap keeping $g$ fixed,
> the gauge change keeping $S$ fixed -- form a connected graph of
> diameter $4$; the extra step over $J(6,3)$ is real, a witness being
> $(\{0,1,2\}, 3) \to (\{0,2,4\}, 3) \to (\{0,2,4\}, 1) \to
> (\{2,3,4\}, 1) \to (\{3,4,5\}, 1)$, where the gauge $3$ must be moved
> before the triple can reach $\{3,4,5\}$. Consequently, by the cocycle of
> Section 5.2, **every chart transition is a word of length at most $4$ in
> elementary swaps and gauge changes, composed with a single sheet-group
> element.**

*Remark.* Both moves are necessary. Removing gauge changes breaks the type
graph into **six connected components of ten types each** -- freezing the
gauge $g$ confines $S$ to $J(5,3)$ on the five remaining coordinates -- so no
composition of swaps alone realizes a gauge change.

---

## 6. Branch Geometry

### 6.1 Sheet labels and the radical description

The three quadrics are *linear in the squares* $Z_j^2$: solving the
Vandermonde system gives, for any two complementary triples $S, T$,

$$w_S \;=\; -V_S^{-1} V_T\, w_T \qquad (w_S = (z_s^2)_{s \in S}),$$

with exact rational coefficients (the largest coefficient over all triples is
$C = \max_S \lVert V_S^{-1} V_T \rVert_\infty = 112$, exactly). Each patch
therefore carries the eight signed sheets
$z_s = \varepsilon_s \sqrt{R_s}$, $\varepsilon \in \{\pm 1\}^3$, and the
signed radical descriptions are exhaustive. They are *distinct* exactly off
the branch locus: over $X^\circ$ a point determines one sign pattern, whereas along
a branch stratum the sign patterns are identified precisely in the slots whose
solved coordinate vanishes (there $+\sqrt{0} = -\sqrt{0}$). Every certified
chart domain lies in $X^\circ$ for its own projection, by Corollary 3.3, so on
chart data the sign pattern is unambiguous. The sheet sign pattern
$\varepsilon$ is part of the chart data.

A discipline that the certification imposed and the paper keeps: **sign patterns
are derived, never defaulted.** In the bridge construction of Section 7 the
sign pattern left at the default $(1,1,1)$ was *refused* by the exact regluing
test, and the derived sign pattern $(1,-1,-1)$ closed it; the refusal was correct
behaviour, not an error.

### 6.2 Deck transformations

Diagonal sign matrices $Z_j \mapsto \pm Z_j$ preserve each quadric
$Q_m(Z) = \sum_a \mu_a^m Z_a^2$ trivially (coordinates enter only through
their squares), so every such matrix induces a holomorphic automorphism of
$X$. Modulo the projective scalar $-1$ they form the **global diagonal sign
automorphism group**
$$G_{\mathrm{sign}} \simeq (\mathbf{Z}/2)^5, \qquad |G_{\mathrm{sign}}| = 32,$$
which is *not* the deck group of any one projection. For a fixed $\pi_S$ the
deck group is the vertical subgroup
$G_{\mathrm{deck}}(S) \simeq (\mathbf{Z}/2)^3$ of order $8$ and index $4$
(Section 6.3), acting as the sheet permutations over $X^\circ$. The certification layer works with one distinguished non-scalar
representative,

$$D \;=\; \mathrm{diag}(+1, -1, +1, +1, +1, -1),$$

an involution which is *not* a projective scalar, and certifies its exact
reuse of chart data: the deck-transported charts remain valid with a
certified radius $\rho_D \ge 1/4096$ (conservative), with the transported
rows split into $36$ ambient and $28$ relative-open cases. What is *not*
trivial about $D$ -- and is measured, not asserted -- is Section 7.3: $D$ is
the exact discrepancy between conjugation and continuation.

### 6.3 One abstract transformation, different chart realizations

The same abstract deck element acts differently depending on the chart type.
In a chart where a flipped coordinate is a *solved* variable, the element is
a sheet permutation over a fixed base point (a vertical arrow in the sign pattern);
in a chart where that coordinate is a *base* or *gauge* variable, the same
element moves the base point. The certification history contains a sharp
instance of the distinction: an audit revealed that $316$ atlas nodes had
been implicitly typed as $316$ vertical-$T$ arrows, which is false -- the
correct statement types the *carrier* of the transformation chart by chart,
and one derived predicate is not universal across chart types. In the paper's
language:

$$\text{abstract deck action} \;\neq\; \text{chartwise vertical action},$$

and the atlas bookkeeping keeps the two straight by typing the carrier
$(S, g)$ of every arrow. The split is quantitative: the deck group has order
$32$ (diagonal signs $\{\pm 1\}^6$ modulo the projective scalar), generated
by five flips; for each chart type the *vertical* subgroup -- the elements
fixing the base point, i.e. those whose sign vector is constant on the
complement $T$ -- has order $8$ and **index $4$**, the quotient
$(\mathbf{Z}/2)^2$ moving the base. For the distinguished element
$D = \mathrm{diag}(+1,-1,+1,+1,+1,-1)$ of Section 6.2 this is completely
explicit: $D$ is vertical for exactly $4$ of the $20$ triples, i.e. $12$ of
the $60$ chart types, and moves the base point on the remaining $48$
(Figure 3, where the pattern is visible: vertical exactly on the triples
containing both coordinates that $D$ flips). Deck statements about *sheets* are made relative to
a chart's sign pattern convention (source and target sheets qualified), never as
absolute class labels.

### 6.4 Away from the branch locus

Let $X^\circ = X \setminus \pi^{-1}(B)$ with $B$ the branch divisor of the
chosen projection. Only on $X^\circ$ is the projection an etale covering
with a genuine local system of sheets; on $B$ the correct language is
ramification. All local-system statements in this paper are restricted to
$X^\circ$.

---

## 7. Computer-Assisted Analytic Continuation

### 7.1 Why continuation needs certification

A root formula chosen at a centre does not, by itself, stay on the correct
branch over a whole domain: the radicand $R$ can cross the branch cut. We
fix the convention throughout: $\sqrt{\cdot}$ denotes the principal square
root, with cut $(-\infty, 0]$, so the principal determination of $\sqrt{R}$ is
available exactly where $R$ avoids $(-\infty, 0]$, and a convexity lemma for
the certified boxes turns this into a checkable sign condition. The rotated determination

$$\sqrt{R}_{\mathrm{rot}} \;=\; \sigma\, i\, \sqrt{-R}_{\mathrm{principal}}$$

costs no new guard: it is available exactly where $R$ avoids
$[0, +\infty)$, since that holds iff $-R$ avoids $(-\infty, 0]$, so the
existing guard applies verbatim to $-R$. Together the two determinations
cover $\mathbf{C} \setminus \{0\}$. The sign $\sigma$ is not free -- it is
fixed by the argument half-plane: $\sigma = +1$ for $\arg R \in (0, \pi]$,
$\sigma = -1$ for $\arg R \in (-\pi, 0)$.

### 7.2 Certified bridge domains

A *bridge* is a box straddling a regime boundary on which the section is
built bilaterally. In the certified panel ($64$ bridge charts over a witness
arc): the analytic regime of each solved line is assigned *by the certified
sign of* $\mathrm{Re}\, R$, without trial (pattern
(principal, canonical, canonical) on all $64$); the section is bilateral by
construction -- on a bilateral box $\mathrm{Im}\, R$ changes sign, so a
one-sided constructor with a fixed $\sigma$ cannot exist, and the negative
control that forces one is refused $64/64$; the target chart gauge is
bounded away from zero; recentring is anisotropic ($2H$ in the imaginary
direction, $H$ in the real one) by an exact affine substitution in rational
arithmetic, adding no remainder; and the bridge sign pattern is derived
($(1,-1,-1)$ on all $64$), as in Section 6.1.

Two senses of "exact" occur in this section and are kept apart. An
*algebraic* identity -- the cocycle law, the deck relation of Proposition 7.1
-- is derived symbolically and holds identically. A *certified residual* is a
number with an enclosure: when a regluing is reported at
$4.48 \times 10^{-15}$, that is an outward-rounded bound on a difference, not
a proof that the difference vanishes. Both are certified; only the first is an
equality.

Each bridge recloses on the lower side with a certified residual: open overlap
of certified minimal width $5.49 \times 10^{-4}$ in all four coordinates, anchor strictly
interior, sheet separation at least $1.042$, recentred difference
$4.48 \times 10^{-15}$.

### 7.3 Continuation versus conjugation

> **Proposition 7.1 (the deck discrepancy).** On the upper side of the
> bridge panel, the analytically continued section does *not* reglue to the
> conjugated atlas: the conjugate sections fail the exact regluing test on
> all $64$ bridges. It reglues exactly to the *deck-translated* conjugate
> $D \cdot Z_{\mathrm{conj}}$, with
> $D = \mathrm{diag}(+1, -1, +1, +1, +1, -1)$, on all $64$ bridges, with
> $O(1)$ margins. The identity is discriminating: changing any single sign
> of $D$ breaks it (negative control).

In slogan form: the *derived mirror chart* and the *analytically continued
sheet* are different objects, and their exact discrepancy is the deck
transformation $D$ (Figure 4). Conjugation of the chart data is a legitimate construction -- but it is
antiholomorphic in nature and it lands on the other sheet, so the conjugated
description is a *diagnostic*, related to the continued one by $D$.

The sign pattern has a predictable geometric reason. On the real corner
$\{\mathrm{Im}\, u = \mathrm{Im}\, v = 0\}$ the radicand $R$ is real: where
$\mathrm{Re}\, R > 0$ the principal root is real and the antiholomorphic
involution *fixes* it; where $\mathrm{Re}\, R < 0$ the canonical root
$i\sqrt{-R}$ is purely imaginary and the involution *negates* it. The
mixed sign pattern of $D$ is exactly the record of which solved lines sit in
which regime.

### 7.4 Exact overlap closure and the certified nerve

Beyond the bridges, the panel closes into a nerve in which *an edge is a
certified transition*, never a mere geometric intersection: $380$ nodes
($316$ lower charts $+\ 64$ bridges); $5396$ lower-lower edges imported
from the prior panel and re-verified; $64$ bridge-lower edges certified;
$210/210$ bridge-bridge edges certified (minimal margin $1.043$, sup
difference $7.94 \times 10^{-15}$); $588$ new triples; the nerve is
connected. The purely geometric intersection graph ($338$ edges) is
published separately under its true name and carries no analytic claim.

Continuation also crosses codimension-one real faces: on a witness tile the
continued section traverses the $\mathrm{Re}$-face into the neighbouring
canonical cell, and the sheet reached is *identified by a control
experiment, not assumed*: neither recorded label closes naively; running the
same experiment with the home label on the home side reproduces the exact
predicted mismatch pattern $\theta = \varepsilon_{\mathrm{derived}} \odot
\text{label}$, elucidating the discrepancy as a measured convention
difference between the register convention and the sign-certified regime
convention. The identification is therefore *relative* and exact.

The body exhibits three worked instances --
one bridge with its certified margins, one face crossing with its control
experiment, and the deck identity of Proposition 7.1 -- each with the exact
numbers quoted above. The full $64$-bridge panel, the $380$-node nerve with
its $5{,}670$ certified edges, and every box enumeration stay in the
supplement and the artifact record.

---

## 8. Certification Architecture and Reproducibility

### 8.1 Proof layers

The layering is the architecture of the paper (Figure 5): exact algebra
supplies the identities, the quantitative lemma turns them into a chart with a
radius, certification discharges the hypotheses, and only then does the atlas
close.

| Claim | Proof mode |
| --- | --- |
| $X$ is K3 | exact algebra |
| branched-cover structure | exact algebra |
| transition formulas | exact symbolic identities |
| chart applicability | interval / directed arithmetic |
| quantitative radius | exact bound + certified inputs |
| branch regime on a box | interval sign certificate |
| overlap continuation | exact identity after certified applicability |
| exploratory metric probes | validated numerics, not theorem input |

### 8.2 Negative controls

Each principal certificate is accompanied by at least one perturbation known
to invalidate the relevant hypothesis (wrong pivot, mutated sign pattern, wrong
branch, wrong transition, box crossing a branch locus), verifying that the
test is discriminating rather than vacuous.

### 8.3 Reproducibility package

Exact equations; certificate JSONs; scripts; environment pins; SHA-256
hashes. Verification is a single command,

```
python3 verification/verify.py
```

which distinguishes three levels, and the distinction is the point. The cheap
certificates (each under a second) are **re-executed** and checked for all
checks green, all negative controls green, and an unchanged outcome; a
regenerated artifact is deliberately *not* hash-compared, since its provenance
field changes with every commit and comparing it would test the repository
rather than the mathematics. The coverage enumeration is **recomputed** and
its counters compared one by one with the shipped file, gauge by gauge: it is
the one certificate whose verdict would otherwise merely be read back. The
expensive certificates -- the bridge panel, the metric path, the face
crossing -- are **hash-verified** against recorded SHA-256 values; reproducing
them needs the compute machine, not this script.
The verifier fails, with a nonzero exit code, on any altered hash or any red
condition, and it is read-only on the working tree: replaying rewrites provenance
fields, so the original bytes are snapshotted and restored. Large box-enumeration artifacts never enter the PDF: dataset/supplement
only.

---

## 9. Validated Metric Compatibility on Selected Continuation Domains

As a downstream demonstration that the certified analytic atlas can support
validated tensorial computations, we report certified relative residuals of
order $10^{-9}$ on bridge/edge domains and a positive Weyl margin, under a
preregistered threshold $\delta = 10^{-5}$.
The panel runs on one certified cell ($S = \{0,1,5\}$, gauge $z_2$, sign pattern
$\varepsilon = (1,1,1)$) against a threshold **fixed before the run**,
$\delta = 10^{-5}$. On the $64$ bridge charts the relative residual of the
metric path is at most $5.93 \times 10^{-9}$; on the $274$ certified nerve
edges it is at most $3.19 \times 10^{-9}$; the Weyl slack stays above
$2.98 \times 10^{-2}$ throughout, and the frozen target-kind pattern
(principal, rotated, principal) holds on $64/64$.

Two honest qualifications. First, metric congruence is certified *below*
$\delta$, not exactly: these are enclosures of a residual, not an identity.
Second, the certificate runs the metric path from one representative on each
certified domain; **no independent two-sided metric comparison is claimed**,
and the artifact records that fact explicitly.

**This section does not establish a Ricci-flat or hyperkähler metric on
$X$.**

---

## 10. Relation to Existing Explicit and Numerical K3 Constructions

### 10.1 Algebraic explicitness

Quartic, CI(2,2,2), Kummer and elliptic K3 models are classical, and the
general theory of these surfaces is standard [1, 2]. The novelty here is not
the existence of the equations.

### 10.2 Analytic and asymptotic constructions

Yau's theorem [3] gives existence of the Ricci-flat Kähler metric without any
explicit description; Gross-Wilson [4] and the Kummer gluing tradition, up to
recent work [5], prove existence and asymptotics by analytic means. Our object
is not a new abstract existence proof, and this paper claims no metric at all.

### 10.3 Numerical K3 geometry

Donaldson's balanced-metric algorithm [6], the Headrick-Wiseman K3
computations [7], the energy-functional approach of Headrick and Nassar [8] --
applied by them to the Fermat quartic K3 -- the projective-embedding methods
of Douglas-Karp-Lukic-Reinbacher [9] and the machine-learning approximations
of the `cymetric` line [10] produce global numerical fits of Calabi-Yau
metrics. The cited implementations do not provide interval-enclosed chart
domains or exact transition certificates of the kind assembled here, and their
reported residuals are diagnostics rather than enclosures. The two lines of
work answer different questions.

### 10.4 Computer-assisted proofs in complex geometry

The closest methodological neighbour is recent: Ishige and Takayasu [11]
compute the monodromy of Picard-Fuchs equations for a family of K3 toric
hypersurfaces by rigorous analytic continuation along contours, controlling
truncation and rounding with interval arithmetic. Their certified continuation
lives in the *base parameter of a family*; ours lives in *charts on a fixed
surface*, and the two are complementary. The general validated-numerics
toolkit we rely on is standard [12-15], as is the Newton-Kantorovich/Krawczyk
lineage behind Proposition 3.2 [16, 17].

Against that background the claim of this paper is deliberately narrow:

> To our knowledge, the combination of an explicit CI(2,2,2) model,
> quantitative certified charts with a guaranteed radius, exact branch-aware
> transition data, and reproducible computer-assisted continuation between
> charts has not previously been provided for a concrete K3 surface.

We do *not* claim priority for certified continuation on K3 data in general
[11], and we make no claim about computer-assisted proofs of Calabi-Yau
metrics, which are the subject of separate ongoing work.

## 11. Discussion

### 11.1 What has actually been achieved

$$\text{exact algebraic model} \rightarrow \text{quantitative local charts}
\rightarrow \text{certified transitions} \rightarrow \text{global analytic
closure}.$$

### 11.2 What remains open

No certified global Ricci-flat metric; no exact hyperkähler identity; no K3
family with controlled constants; no seven-dimensional gluing result. (The
internal sign pattern separates the closed fixed-K3 atlas theorem from the open
metric and family programmes.)

### 11.3 Why this infrastructure matters

The atlas is designed as an analytic substrate for validated PDE, metric and
deformation computations on the same explicit surface.

---

## Use of large language models

This work was carried out in sustained iteration with two large language models:
Claude (Anthropic) and GPT (OpenAI). Their use went well beyond assisted copy
editing, and it is therefore documented here rather than left undeclared. The
certificate scripts, the negative controls, the derivations reported in the
appendices and this manuscript were developed over repeated cycles of drafting,
criticism and revision between the author and the two models; the models also
assisted in checking references and code outputs against their sources. Two
defects in the certificates were found in exactly this way, while writing
the argument out in continuous prose for an outside reader, and both are
recorded in the paper: an invariance checked under one group generator instead of
the whole group, and a bound on a product used where a bound on a minimum was
required. Neither model is an author: authorship carries accountability, and
every condition, every run and every sentence of this paper was reviewed, verified
and decided by the author, who takes sole responsibility for its accuracy,
integrity and originality.

## Acknowledgements

The author thanks Anthropic for Claude and OpenAI for GPT. The nature and extent
of the two models' contribution to this work is documented in the section "Use
of large language models" above.

## Statements and Declarations

**Competing interests.** The author has no financial interests, and no
institutional or funding relationship, bearing on this work.

**Funding.** No funding was received for this work.

**Data and code availability.** Every quantitative claim in this paper is backed
by a machine-checkable certificate. The certificate scripts, their JSON outputs,
the figure sources and the verification entry point of Appendix F are publicly
available in the repository accompanying this paper; a timestamped archival
deposit identifier will be supplied at submission. Verification is a single
command, and it distinguishes three regimes: the certificates it re-executes,
the coverage enumeration it recomputes and compares counter by counter, and
those it verifies by hash (Appendix F).

**Ethics approval.** Not applicable. The work involves no human participants, no
animal subjects and no personal data.

**Author contributions.** Sole author: the author designed the certificates,
wrote the analysis code, ran the verifications and wrote the manuscript. Large
language models were used throughout, in the manner and to the extent documented
above; they meet no authorship criterion, and accountability for the work is the
author's alone.

---

## References

Entries were checked against publisher, society or catalogue metadata on
2026-08-20.

1. W. P. Barth, K. Hulek, C. A. M. Peters, A. Van de Ven, *Compact Complex
   Surfaces*, 2nd ed., Ergebnisse der Mathematik und ihrer Grenzgebiete,
   3. Folge / A Series of Modern Surveys in Mathematics **4**, Springer,
   2004. doi:10.1007/978-3-642-57739-0.
2. D. Huybrechts, *Lectures on K3 Surfaces*, Cambridge Studies in Advanced
   Mathematics **158**, Cambridge University Press, 2016.
   doi:10.1017/CBO9781316594193.
3. S.-T. Yau, On the Ricci curvature of a compact Kähler manifold and the
   complex Monge-Ampère equation, I, *Comm. Pure Appl. Math.* **31** (1978),
   339-411. doi:10.1002/cpa.3160310304.
4. M. Gross, P. M. H. Wilson, Large complex structure limits of K3 surfaces,
   *J. Differential Geom.* **55** (2000), no. 3, 475-546.
   doi:10.4310/jdg/1090341262. arXiv:math/0008018.
5. B. Shackleton, A Calabi-Yau metric on the Kummer surface,
   arXiv:2605.02046 (2026).
6. S. K. Donaldson, Some numerical results in complex differential geometry,
   *Pure Appl. Math. Q.* **5** (2009), no. 2, 571-618.
   doi:10.4310/pamq.2009.v5.n2.a2.
7. M. Headrick, T. Wiseman, Numerical Ricci-flat metrics on K3, *Class.
   Quantum Grav.* **22** (2005), no. 23, 4931-4960.
   doi:10.1088/0264-9381/22/23/002. arXiv:hep-th/0506129.
8. M. Headrick, A. Nassar, Energy functionals for Calabi-Yau metrics,
   *Adv. Theor. Math. Phys.* **17** (2013), no. 5, 867-902.
   doi:10.4310/ATMP.2013.v17.n5.a1. arXiv:0908.2635.
9. M. R. Douglas, R. L. Karp, S. Lukic, R. Reinbacher, Numerical Calabi-Yau
   metrics, *J. Math. Phys.* **49** (2008), 032302. doi:10.1063/1.2888403.
   arXiv:hep-th/0612075.
10. M. Larfors, A. Lukas, F. Ruehle, R. Schneider, Numerical metrics for
    complete intersection and Kreuzer-Skarke Calabi-Yau manifolds, *Mach.
    Learn.: Sci. Technol.* **3** (2022), 035014.
    doi:10.1088/2632-2153/ac8e4e. arXiv:2205.13408.
11. T. Ishige, A. Takayasu, Computer-assisted proofs for finding the monodromy
    of Picard-Fuchs differential equations for a family of K3 toric
    hypersurfaces, *Commun. Nonlinear Sci. Numer. Simul.* **152** (2026),
    109408. doi:10.1016/j.cnsns.2025.109408. arXiv:2501.03792.
12. W. Tucker, *Validated Numerics: A Short Introduction to Rigorous
    Computations*, Princeton University Press, 2011.
    ISBN 978-0-691-14781-9.
13. R. E. Moore, R. B. Kearfott, M. J. Cloud, *Introduction to Interval
    Analysis*, SIAM, 2009, 223 pp. doi:10.1137/1.9780898717716.
14. S. M. Rump, Verification methods: rigorous results using floating-point
    arithmetic, *Acta Numerica* **19** (2010), 287-449.
    doi:10.1017/S096249291000005X.
15. M. T. Nakao, M. Plum, Y. Watanabe, *Numerical Verification Methods and
    Computer-Assisted Proofs for Partial Differential Equations*, Springer
    Series in Computational Mathematics **53**, Springer, 2019.
    doi:10.1007/978-981-13-7669-6.
16. R. Krawczyk, Newton-Algorithmen zur Bestimmung von Nullstellen mit
    Fehlerschranken, *Computing* **4** (1969), 187-201.
    doi:10.1007/BF02234767.
17. J. M. Ortega, The Newton-Kantorovich theorem, *Amer. Math. Monthly*
    **75** (1968), no. 6, 658-660. doi:10.2307/2313800.

---

## Appendix A -- Exact defining equations

**A.1 The data.** $\mu = (\mu_0, \dots, \mu_5) = (1, 2, 3, 5, 7, 11)$ and
$F_k(z) = \sum_{j=0}^{5} \mu_j^{\,k} z_j^2$ for $k = 0, 1, 2$;
$X = \{F_0 = F_1 = F_2 = 0\} \subset \mathbf{P}^5$. Everything else in the
paper is a rational function of these six integers.

**A.2 Vandermonde determinants.** For a triple $S$,
$V_S = \det(\mu_s^k)_{k=0,1,2;\, s \in S} \in \mathbf{Z}$. All twenty are
nonzero (this is Lemma 2.1's ingredient); their absolute values are

| $S$ | $|V_S|$ | $S$ | $|V_S|$ | $S$ | $|V_S|$ | $S$ | $|V_S|$ |
| --- | ---: | --- | ---: | --- | ---: | --- | ---: |
| 012 | 2 | 013 | 12 | 014 | 30 | 015 | 90 |
| 023 | 16 | 024 | 48 | 025 | 160 | 034 | 48 |
| 035 | 240 | 045 | 240 | 123 | 6 | 124 | 20 |
| 125 | 72 | 134 | 30 | 135 | 162 | 145 | 180 |
| 234 | 16 | 235 | 96 | 245 | 128 | 345 | 48 |

with extremes $\min |V_S| = 2$ (at $S = \{0,1,2\}$) and $\max |V_S| = 240$.
The largest transition coefficient over all triples is
$C = \max_S \lVert V_S^{-1} V_T \rVert_\infty = 112$, attained at
$S = \{0,1,2\}$; the largest single entry is $80$.

**A.3 Normalization.** $c_k^2 = \sum_j \mu_j^{2k}$ gives
$c_0^2 = 6$, $c_1^2 = 209$, $c_2^2 = 17765$. The squares are integers, the
$c_k$ are not. The **algebraic layer** (Sections 2, 4, 5 and Appendices A, C)
lives over $\mathbf{Q}$ with $V_S \in M_3(\mathbf{Z})$; the **normalized
layer** (Section 3 and Appendix B) is a conditioning convention, and the
certificate brackets the $c_k$ by outward-rounded rationals rather than
approximating them.

**A.4 Distinguished patch.** $S = \{3,4,5\}$, gauge $z_0 = 1$, base
$(u,v) = (z_1, z_2)$: the radicands and the degree-six discriminant are
displayed in Section 2.3, and the corresponding regime geometry is Figure 1.

---

## Appendix B -- Proof of Proposition 3.2

**B.1 Setup.** Fix a chart type $(S, g)$. Write $w = z_S \in \mathbf{C}^3$ for
the solved variables, $b = (u,v) \in \mathbf{C}^2$ for the base, $z_g = 1$ for
the gauge, and $G(w; b) := \tilde F(w, b, 1) : \mathbf{C}^3 \times
\mathbf{C}^2 \to \mathbf{C}^3$ for the normalized equations. Let $w_0$ solve
$G(w_0; b_0) = 0$ at a centre $b_0$ in the sector, so $|b_0| \le 1$.

**B.2 Three exact identities.** Because the $F_k$ are diagonal quadrics, the
following hold identically -- they are not estimates:

* **(E1)** $\partial_w G(w; b) = \tilde M(w) = 2\, \tilde V_S \,
  \mathrm{diag}(w)$, independent of $b$;
* **(E2)** $G(w; b) - G(w_0; b) - \tilde M(w_0)(w - w_0) =
  \tilde V_S \cdot \big((w - w_0)^2\big)$, the square taken componentwise:
  the Taylor remainder is one explicit quadratic term and the "Hessian" is
  constant;
* **(E3)** $G(w_0; b) - G(w_0; b_0) = \tilde V_{uv} \cdot (b^2 - b_0^2)$.

**B.3 Residual.** Let $\sigma \le s_{\min}(\tilde M(w_0))$ and let $a_{S,g}$
bound $\lVert \tilde V_{uv} \rVert$. For $|b - b_0| \le \rho$, (E3) and
$G(w_0;b_0) = 0$ give $\lVert G(w_0; b) \rVert \le a_{S,g}\,\rho\,(2|b_0| +
\rho) \le a_{S,g}\,\rho\,(2 + \rho)$, hence the simplified-Newton residual

$$\eta \;:=\; \lVert \tilde M(w_0)^{-1} G(w_0; b) \rVert
\;\le\; \frac{a_{S,g}\,\rho\,(2 + \rho)}{\sigma}.$$

**B.4 Contraction.** Let $T(w) = w - \tilde M(w_0)^{-1} G(w; b)$. By (E1),
$$DT(w) = I - \tilde M(w_0)^{-1}\tilde M(w)
= 2\,\tilde M(w_0)^{-1} \tilde V_S \,\mathrm{diag}(w_0 - w),$$
so $\lVert DT(w) \rVert \le (L/\sigma)\lVert w - w_0 \rVert$ with
$L = 2\,\sigma_{\max}(\tilde V_S)$. With
$$h \;=\; \frac{L\,\eta}{\sigma} \;<\; \frac{1}{2},$$
$T$ maps the closed ball of radius $r^* = (1 - \sqrt{1 - 2h}\,)\,\sigma/L$
around $w_0$ into itself and contracts there; $r^* \le 2\eta$. Hence a zero
$w(b)$ exists in that ball, and it is unique in it.

**B.5 Uniqueness in the large ball.** For diagonal quadrics the midpoint
identity is exact: for any $w, w'$,
$$G(w; b) - G(w'; b) \;=\; \tilde M\!\left(\frac{w + w'}{2}\right)(w - w'),$$
since $w_s^2 - w_s'^2 = 2 \cdot \frac{w_s + w_s'}{2} \cdot (w_s - w_s')$
term by term. If $w, w'$ are two zeros with
$\lVert w - w_0 \rVert < \sigma/L$ and $\lVert w' - w_0 \rVert < \sigma/L$
(open ball, so the perturbation estimate stays strictly below the
invertibility threshold), then $\tilde M((w+w')/2)$ is invertible there and $w = w'$. Uniqueness
therefore holds on a ball much larger than $r^*$.

**B.6 Holomorphy.** Start the iteration at the constant $w_0$. Since $G$ is
quadratic and $\tilde M(w_0)$ is a fixed matrix, each iterate $T^n(w_0)$ is a
*polynomial* in the base variables $b$. The convergence of B.4 is uniform on
the ball, so the limit $w(b)$ is holomorphic in $b$, and the chart map
$b \mapsto (w(b), b)$ is a holomorphic section.

**B.7 Certified constants.** All quantities above are bounded by exact
rationals: $\sigma_{\max}(\tilde V_S) \le \lVert \tilde V_S \rVert_F$
(upper), $\sigma_{\min}(\tilde V_S) \ge |\det \tilde V_S| /
\lVert \tilde V_S \rVert_F^2$ (lower), square roots bracketed by integer
square roots to $10^{-40}$. The checks verify $h < 1/2$, the contraction
factor, $r^* \le 2\eta$ and $r^* < \sigma/L$ on all $60$ types, and a
negative control multiplying $\rho$ by $100$ drives $h$ above $1/2$. The
resulting uniform radius is
$\rho_{\mathrm{unif}} = \min_{S,g} \rho_{S,g} \ge 2.10 \times 10^{-12}$, the
per-type values spanning $[2.10 \times 10^{-12},\, 1.71 \times 10^{-4}]$
(Figure 2); these are Frobenius/determinant bounds, deliberately not sharp.

The lower bound on $\sigma$ deserves its own line, because it is where an
earlier version of this certificate was optimistic. Since
$\tilde M(w_0) = 2\,\tilde V_S \,\mathrm{diag}(z_S)$,
$$s_{\min}(\tilde M(w_0)) \;\ge\; 2\, s_{\min}(\tilde V_S)\,
\min_{s \in S}|z_s|,$$
so what is needed is the **minimum** of the solved coordinates. The pivot
threshold bounds only their **product**; substituting one for the other is
invalid whenever a solved coordinate exceeds $1$, which happens on certified
domains, and an explicit point of $X$ on a certified chart domain witnesses
the failure. The valid bound on $\min_s |z_s|$ is Lemma 3.1, and the radius
quoted above is the one it yields.

---

## Appendix C -- Transition formulas

**C.1 The square-level cocycle.** The quadrics are linear in $w_j = z_j^2$,
so $V_{\mathrm{full}}\, w = 0$ determines all six squares from any admissible
triple. For each ordered pair $(S, S')$ this yields an exact rational matrix
$P_{S'S} \in M_3(\mathbf{Q})$ with $w_{S'} = P_{S'S} w_S$ on $X$, obtained by
composing $w_T = -V_T^{-1} V_S w_S$ with the selection of the $S'$ entries.
Verified in exact arithmetic: $P_{SS} = \mathrm{Id}$ (20/20),
$P_{SS'}P_{S'S} = \mathrm{Id}$ (400/400), and
$P_{S''S'}P_{S'S} = P_{S''S}$ on all $20^3 = 8000$ ordered triples.

**C.2 Elementary swap.** For $S = \{3,4,5\} \to S' = \{2,4,5\}$,

$$P_{S'S} = \begin{pmatrix} -6 & -15 & -45 \\ 0 & 1 & 0 \\ 0 & 0 & 1
\end{pmatrix}, \qquad\text{i.e.}\qquad w_2 = -6w_3 - 15w_4 - 45w_5 .$$

Exactly one row differs from the identity -- the row of the newly solved
coordinate. On coordinates the map is
$z_{s'} = \varepsilon_{s'}\sqrt{(P_{S'S} w_S)_{s'}}$, so the sign pattern travels
with the transition.

**C.3 Gauge change.** $(S, g) \to (S, g')$ is the projective rescaling by
$z_g/z_{g'}$; it changes no sign pattern and no solved triple.

**C.4 Generation.** The twenty triples with adjacency $|S \cap S'| = 2$ form
$J(6,3)$: $20$ nodes, $90$ edges, connected, diameter $3$. The $60$ chart
types under $\{$swap, gauge change$\}$ form a connected graph of diameter
$4$ -- the swap holds $g$ fixed and the gauge change holds $S$ fixed, so a
transition whose target triple contains the source gauge needs a gauge move
of its own. With C.1, every transition is a word of length at most $4$ in
these two moves composed with a single sheet-group element
(Proposition 5.1). Removing
gauge changes disconnects the type graph into six components of ten.

**C.5 Scope.** C.1--C.4 are algebraic identities over $\mathbf{Q}$. That a
given transition *applies* on a given domain -- gauge nonvanishing, branch
regime constant, overlap open -- is a separate, certified statement
(Sections 7 and 8).

---

## Appendix D -- Branch conventions

**D.1 Chart data.** A chart carries: a solved triple $S$; a gauge $g \notin S$
with $z_g = 1$; base variables $b = (u,v)$; and a sign pattern
$\varepsilon \in \{\pm 1\}^3$ with $z_s = \varepsilon_s \sqrt{R_s(b)}$. The
sector $C_{S,g} = \{|z_g| = \max_{t \in T} |z_t|\}$ is a **selector** of
representatives (it gives $|u|, |v| \le 1$ and feeds Lemma 3.1); it is never
used as an open set of a cover.

**D.2 Determinations.** $\sqrt{\cdot}$ is the principal square root, cut
$(-\infty, 0]$. The **principal** determination of $\sqrt{R}$ is available
where $R$ avoids $(-\infty,0]$; the **rotated** (canonical) determination
$\sqrt{R} = \sigma\, i\, \sqrt{-R}$ is available where $R$ avoids
$[0,+\infty)$, the two conditions being exchanged by $R \mapsto -R$. The
sign is not free:

| $\arg R$ | $\arg(-R)$ | $i\sqrt{-R}$ | $\sigma$ |
| --- | --- | --- | ---: |
| $(0, \pi]$ | $\arg R - \pi$ | $+\sqrt{R}$ | $+1$ |
| $(-\pi, 0)$ | $\arg R + \pi$ | $-\sqrt{R}$ | $-1$ |

On a bilateral box $\mathrm{Im}\,R$ changes sign, so no single $\sigma$
exists there -- which is why bridge sections are built bilaterally
(Section 7.2), and why forcing a component $\sigma$ on a bridge is refused
$64/64$.

**D.3 Sign groups.** Diagonal sign matrices preserve every $F_k$. Modulo the
projective scalar $-1$ they form the global group
$G_{\mathrm{sign}} \simeq (\mathbf{Z}/2)^5$ of order $32$, generated by
$\sigma_1, \dots, \sigma_5$ with
$\sigma_0 = -\sigma_1\sigma_2\sigma_3\sigma_4\sigma_5$. For a fixed
projection $\pi_S$ the deck group is the vertical subgroup
$G_{\mathrm{deck}}(S) \simeq (\mathbf{Z}/2)^3$, of order $8$ and index $4$:
an element is vertical for $(S,g)$ iff its sign vector is constant on the
complement $T$. For $D = \mathrm{diag}(+1,-1,+1,+1,+1,-1)$ this holds for
exactly $4$ triples, i.e. $12$ of the $60$ types (Figure 3).

**D.4 Continuation versus conjugation.** On the real corner
$\{\mathrm{Im}\,u = \mathrm{Im}\,v = 0\}$ the radicand is real, and the
antiholomorphic involution acts on each solved line according to its regime:

| regime on the line | root | conjugation acts by |
| --- | --- | --- |
| $\mathrm{Re}\,R > 0$ (principal) | real | fixes it |
| $\mathrm{Re}\,R < 0$ (canonical) | purely imaginary $i\sqrt{-R}$ | negates it |

The mixed sign pattern of $D$ is exactly the record of which lines sit in
which regime; this is the geometric content of Proposition 7.1.

**D.5 Discipline.** Ledgers are **derived, never defaulted**: in the bridge
construction the default $(1,1,1)$ was refused by the exact regluing test and
the derived $(1,-1,-1)$ closed it. Sheet statements are made relative to a
chart's sign pattern convention, never as absolute class labels; identifications
across conventions are measured, as in the face-crossing control experiment
$\theta = \varepsilon_{\mathrm{derived}} \odot \text{label}$ (Section 7.4).

---

## Appendix E -- Certificate specification

**E.1 Anatomy.** A certificate is one JSON document. The vocabulary is
shared, but it is not a schema every file satisfies, and stating otherwise
would be a claim the reader can falsify in one command. Counted as TOP-LEVEL
fields over the fourteen shipped files -- the level at which a skeleton is a
skeleton, and not "the key appears somewhere in the document":

| field | certificates carrying it |
| --- | ---: |
| `artifact` | 13 |
| `checks`, `checks_passed`, `checks_total` | 12 |
| `built_from_head`, `self_sha256` | 10 |
| `kind`, `seconds` | 9 |
| `upstream`, `perturbation_tests` | 8 |
| `does_not_attest`, `subject` | 6 |
| `outcome`, `provenance` | 5 |
| `n_atteste_pas` | 4 |

Three of the files nest `self_sha256` inside a `provenance` block instead of
carrying it at the top level, which is why that row reads 10 and not 13.

The fields mean: an `artifact` name; a `subject` line stating the claim; a
`kind` (exact algebra, interval, closed form); `upstream`, holding the
SHA-256 of every certificate consumed; the payload blocks; `checks`, a
dictionary of **named booleans**; `checks_passed` / `checks_total`;
`perturbation_tests`, the negative controls; an `outcome` string; the explicit
list of what the certificate does *not* establish -- written
`does_not_attest` in the six most recent files and `n_atteste_pas` in the four
oldest, a rename this repository did not finish; and the provenance triple
`seconds`, `built_from_head`, `self_sha256`.

The four design records are the least uniform: they predate the convention
and carry the payload without the counters. The nine files the verification
command checks all carry `checks` and its two counters.

**E.2 Green and red.** A certificate is green iff every condition and every
negative control is true. The `outcome` string encodes the claim when green
and collapses to a `*_checks_red` marker otherwise, so a red certificate
cannot be quoted as if it were green.

**E.3 Negative controls.** Each principal condition ships with at least one
deliberate mutation known to invalidate its hypothesis -- a repeated $\mu$, a
wrong pivot threshold, a mutated sign pattern, a wrong branch, a falsified
transition entry, a box crossing a branch locus. A condition that no mutation can
redden is vacuous and is rejected rather than reported.

**E.4 Derivation discipline.** A separate lint enforces three rules on the
certificates: constants that can be derived must be derived rather than
named (D1); a certificate must not silently consume another block's output
(D2, blocking); and no condition may be constant-true by construction (D3).

**E.5 Claim boundary.** `n_atteste_pas` is not a disclaimer but part of the
specification: it is where the metric probe of Section 9 records that no
independent two-sided comparison was run, and where the atlas certificates
record that the quantitative floor, not the covering, is what the enumeration
audit supplies.

---

## Appendix F -- Reproduction instructions

**F.1 One command.**

```
python3 verification/verify.py
```

It prints one line per *checked* certificate and exits nonzero on any
failure; a `--quick` flag checks hashes only, skipping both replay and
recomputation. It is read-only on the working tree, and takes about two
minutes, almost all of it in the recomputed coverage certificate.

The repository ships fourteen certificate files and the command checks nine
of them: five replayed, one recomputed, three hash-verified (F.3). The other
five are records rather than claims -- four design documents (the uniform
chart lemma, the exact transitions, the regional gluing contract, the closure
skeleton) and the figure manifest. They do have checks of their own, and an
earlier revision of this verifier replayed them; that was reverted on
purpose. Shipping their producers would pull twenty-three further artifacts
into the repository through what those producers read -- contract amendments,
preregistrations, an internal chain none of the theorems above rest on -- and
a repository that carries what it does not need is harder to audit, not
easier.

**F.2 Environment.** Python $\ge 3.10$, NumPy, SymPy 1.14.0, mpmath pinned at
1.3.0 -- the requirements listed in the repository README; the certificates
were produced under Python 3.14 and have also been run under 3.12. The pin is
enforced by `verify.py` itself, which checks the installed mpmath before
anything else runs and refuses to continue against a different version,
because that version is serialized in the provenance of the whole chain.

**F.3 Three tiers.** Five cheap certificates (each under a second) are
**re-executed** and checked for all checks green, all negative controls green
and an unchanged outcome: the chart theorem, the atlas closeout, the
generator/obligation certificate, the smoothness/coverage/transition
certificate, and the sigma-floor correction whose repair is reported in
Appendix B.7 -- that last one is replayed precisely so that the correction is
checked by the reader rather than asserted by the author.

The coverage certificate of Section 4.1 is **recomputed**: its branch-and-bound
enumeration is re-run and its counters are compared one by one with the shipped
file, gauge by gauge. This is the exhaustive enumeration behind the observed
minor floor $4.8$ -- $71{,}807{,}792$ boxes in float64 interval arithmetic with
one-ulp outward rounding, about seventy seconds on four cores -- and it is the
computation on which the preregistered threshold $m_{\mathrm{alg}} = 4$ has its
margin. Reading its verdict would test nothing; re-running it reddens if a single
counter moves. The published run is reproduced exactly, box for box, on all six
gauges.

The three expensive ones are **hash-verified**:

| artifact | SHA-256 (first 16) |
| --- | --- |
| bridge panel (64 charts, nerve) | `58c06da48412c8b1` |
| metric path, full (64 bridges, 274 edges) | `397c19b105819f79` |
| face crossing (control experiment) | `9088e59844c7e0dd` |

**F.4 What re-running does not test.** A regenerated artifact is
deliberately not hash-compared: its provenance field changes with every
commit, so such a comparison would test the repository rather than the
mathematics. Reproducing the expensive certificates from scratch requires the
compute machine, not this script; the figures are regenerated by their own
deterministic script, whose manifest records every plotted number and the
hashes of its sources.

---

## Figures

Generated by `verification/producers/figures.py` (deterministic;
manifest with source hashes and every plotted number in
`figures_manifest.json`). Figures 1, 2 and 3 are data-driven
(computed from $\mu$, $D$ and the certificates); 4 and 5 are schematic, and say so.

1. `fig_chart_regime_geometry.png` -- the three loci $\mathrm{Re}\,R_s = 0$
   of the chart $S = \{3,4,5\}$, $z_0 = 1$ in the complex base slice $v = 0$:
   nested rectangular hyperbolas whose rational coefficients are derived from
   $\mu$, the regime pattern on either side, and a bridge box straddling one
   curve. *(data)*
2. `fig_certified_chart_radii.png` -- certified radius across the 60 chart
   types on a log scale, spanning $[2.10 \times 10^{-12},\, 1.71 \times 10^{-4}]$
   (endpoints rounded outward).
   The dashed line is the observed minimum over the types,
   $2.109 \times 10^{-12}$; the guaranteed floor quoted in the text,
   $\rho_{\mathrm{unif}} \ge 2.10 \times 10^{-12}$, is that value rounded
   outward. Regenerated from the repaired chart certificate. *(data)*
3. `fig_deck_action_by_chart_type.png` -- the 20 triples $\times$ 3 gauges,
   marking where $D$ is vertical (12 types) and where it moves the base (48),
   exhibiting the pattern: vertical exactly when $S \supseteq \{1, 5\}$, the
   coordinates $D$ flips. *(data)*
4. `fig_bridge_continuation.png` -- continuation across a certified bridge:
   open overlap, regime boundary, and the conjugate sheet versus the continued
   sheet $D \cdot Z_{\mathrm{conj}}$. *(schematic, annotated with certified
   margins)*
5. `fig_proof_architecture.png` -- exact algebra $\to$ quantitative lemma
   $\to$ certified applicability $\to$ $\mathcal{X}_{\mathrm{atlas}} \cong X$
   $\to$ certified continuation. *(schematic)*
