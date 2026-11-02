# Why vlink's answers are correct

This note states what vlink relies on from the literature and proves every
other step its trusted kernel (`vlink/kernel.py`, `vlink/field.py`) takes.
Search code (`beads.py`, `necklace.py`, `certify.py`) is never trusted: its
outputs are certificates that the kernel checks.

Throughout, M is a closed 3-manifold triangulated by T, with a squared length
assigned to each edge; every tetrahedron is given the Euclidean metric with
those edge lengths, and M the induced piecewise Euclidean (PE) metric.

## 0. Cited results

- **(G) Gromov's link condition** [Bridson–Haefliger, Ch. II.5]. A PE
  complex with finitely many shapes is locally CAT(0) iff the link of every
  cell is CAT(1). For a 3-manifold the links to check are: edge links
  (circles) and vertex links (piecewise spherical 2-spheres); the links of
  faces and of interior points are automatically fine.
- **(C) Circles and short geodesics** [Bridson–Haefliger, Ch. II.1 and II.5]. A circle of length
  ℓ is CAT(1) iff ℓ ≥ 2π. A compact piecewise spherical 2-complex that is
  locally CAT(1) is CAT(1) iff it contains no closed geodesic of length < 2π.
  A CAT(1) space contains no closed local geodesic of length < 2π.
- **(B) Bowditch** [Bowditch; Elder–McCammond–Meier (EMM) Lemma 4.5]. In a
  locally CAT(1) space that is not CAT(1), the shortest closed geodesic has
  the same length as the shortest *unshrinkable* closed geodesic.
- **(E) Gallery reduction** [Elder–McCammond 2002 §3; EMM Lemmas 4.7–4.9 and
  §4]. In a 2-dimensional piecewise spherical complex, a closed geodesic
  either passes through a vertex — then it is a *necklace* of *beads*
  (vertex-to-vertex geodesic segments through strips of triangles) — or it
  determines a disc, annular or Möbius gallery. Disc galleries reduce to edge
  links; annular galleries are shrinkable (4.7); Möbius strips do not immerse
  in a 2-sphere (EMM §4), so they cannot occur in a vertex link of a
  3-manifold. A geodesic with a bead longer than π is shrinkable; one with a
  bead of length exactly π is equivalent to one whose beads are all shorter
  than π (4.9).
- **(L) Local geodesics at a cone point** [Bridson–Haefliger, Ch. I.5 and II.5]. A path through a
  point W of a piecewise spherical surface, geodesic on either side of W, is a
  local geodesic at W iff both angles it makes at W (measured around the
  circle link of W) are ≥ π.

Combining (G), (C), (B), (E): **M is locally CAT(0) iff (i) every tetrahedron
is non-degenerate, (ii) every edge has total dihedral angle ≥ 2π, and (iii) no
vertex link contains a necklace of length < 2π whose beads are shorter than π,
or have length exactly π.** With `pi_beads=False` the last clause is dropped,
relying on (E, Lemma 4.9); certificates then record that assumption.

## 1. Exact arithmetic (field.py, kernel.AngleSum)

**Field.** Elements are Σ c_s √s over distinct square-free integers s with
rational c_s. Square roots of distinct square-free integers are linearly
independent over ℚ, so the representation is unique and `is_zero` (no stored
terms) is exact.

**Sign.** Pick a square-free p > 1 such that every key is divisible by p or
coprime to it (`_splitter`: start from a key ≠ 1 and replace p by gcd(p, k)
while some gcd lies strictly between 1 and p; no factoring). Write
x = A + B√p, A collecting the keys coprime to p and B the others divided by
p; A and B are free of √p. If B = 0 then sign x = sign A; if A = 0 then
sign x = sign B; if sign A = sign B that is the sign. Otherwise A − B√p has
the sign of A, and x(A − B√p) = A² − pB², so sign x = sign A · sign(A² − pB²).
Products of keys keep the divisible-or-coprime property, so A² − pB² is free
of √p and the recursion terminates. The same p gives the conjugation
√p ↦ −√p used by `inverse`. The float pre-check returns a sign only if
|est| > 10⁻⁹ · Σ|terms|; the float error is at most (k+2)·2.3·10⁻¹⁶·Σ|terms|
for k terms, which is smaller whenever k < 4·10⁶.

**Square-free parts.** `squarefree_split(n)` trial-divides only by primes p
with p³ ≤ n; the remaining cofactor m has no prime factor below n^{1/3}, so
it has at most two prime factors: it is p² (detected exactly by isqrt) or
square-free.

**Towers.** For irrational link fields, √(b²) with b² ∈ F not rational is a
new symbol √r (r > 0 in F). Tower elements are F-combinations of products of
such symbols, reduced by √r² = r. Their sign is decided by rigorous interval
arithmetic (each float error bounded from the magnitudes of the exact parts)
or, when the interval contains 0, by splitting on one symbol exactly as above
(it needs only r > 0). Equality is decided through the sign, so it is exact
even when symbols are algebraically related.

**AngleSum.** An angle θ ∈ [0, π] is given by exact a, b² with b ≥ 0 and
(a, b) ≠ 0, so a + ib is a positive multiple of e^{iθ}. Invariant: z is a
positive multiple of e^{iΘ} with Θ ∈ [hπ, (h+1)π). After adding θ,
Θ' ∈ [hπ, (h+2)π). With s = (−1)^h: s·sin Θ' > 0 iff Θ' ∈ (hπ, (h+1)π),
< 0 iff Θ' ∈ ((h+1)π, (h+2)π), and when sin Θ' = 0, Θ' = hπ iff s·cos Θ' > 0.
So h is updated exactly; "total ≥ kπ" is h ≥ k. Two sums with equal h differ
by less than π, so the sign of Im(z̄₁z₂) decides their order (`angle_leq`).
b = √(b²) is exact: in F when b² is rational, as a tower symbol otherwise.

## 2. Link geometry (tri.py; kernel.Geometry, Link, step)

**Combinatorics.** `tri.py` decodes the isomorphism signature into a gluing
table exactly as Regina's `fromIsoSig` does (Regina 7.4,
`triangulation/detail/isosig-impl.h`), and numbers edges and vertices in order
of first appearance over (tetrahedron, local index), as Regina does; both are
checked against Regina on all 61,848 triangulations of the closed and
Hodgson–Weeks censuses and on 3,000 random closed triangulations
(`tests/test_tri.py`). It rejects anything that is not a closed 3-manifold
triangulation: boundary facets, an edge identified with itself in reverse,
or a vertex link that is not a 2-sphere. (With every edge valid and no
boundary, a vertex link is a closed connected surface — connected because a
vertex class is a single orbit of corners under face gluings — so Euler
characteristic 2 makes it a sphere.) This is what §0 (E) needs: Möbius
galleries do not occur because vertex links are spheres.

The link triangles of vertex v are the corners (τ, w) of tetrahedra with
vertex w in the class v, in sorted order; vertex j of the corner is the j-th
smallest vertex of τ other than w. Crossing the link edge opposite vertex j
is crossing the face of τ opposite that vertex, so the neighbour is the
corner (τ', g(w)) for the tetrahedron τ' glued there by g, and the link
gluing is g restricted to the face. A link vertex is an end of an edge of M
at v (two ends of one edge are distinct link vertices).

Link triangle t of the link of vertex v is the corner of a tetrahedron τ at v;
its vertices are the edge vectors f₀, f₁, f₂ of τ from v, with Gram matrix
G_t[j][k] = (d(v,j) + d(v,k) − d(j,k))/2, rational. Validity of τ (positive
Gram determinants of τ and its faces) makes G_t positive definite.

**Development.** Let t' be glued to t across the link edge opposite vertex d,
i.e. τ' is glued to τ along the face spanned by f_a, f_b. That face is shared,
so G_t and G_t' agree on (a, b) (checked). In τ, f_d = p_a f_a + p_b f_b + n
with n ⊥ span(f_a, f_b), |n| = h; in τ', f'_{d'} = q_a f_a + q_b f_b + n' with
|n'| = h'. In the development of τ ∪ τ' across their common face the two
tetrahedra lie on opposite sides of it, so n' = −(h'/h) n. `step` computes
exactly this, with p, q rational (from the Gram data) and μ = h'/h = √(h'²/h²),
h², h'² rational. Hence developed strips are isometric copies of the union of
the tetrahedra's corners along the strip.

## 3. Bead enumeration (kernel.enumerate_beads): completeness

**Lemma (visible interval).** Let S, A, B be points of S² with A, B ≠ ±S and
S not on the great circle of A, B. A geodesic from S of length < π meets the
open minor arc AB (transversally) iff its initial direction lies strictly
inside the angle at S of the spherical triangle SAB; that angle is < π.

*Proof.* SAB is a non-degenerate triangle with sides < π, so it is convex and
lies in an open hemisphere not containing −S. A geodesic leaving S strictly
inside the angle enters the triangle and, before reaching −S, leaves it
through the open side AB (it cannot leave through SA or SB, which it meets only
at S). Conversely a point X of the open arc AB at distance < π is joined to S
by a minor arc inside the triangle, leaving S strictly inside the angle. ∎

If A ∥ S, then S lies on the great circle of AB, which a geodesic from S meets
only at ±S or runs along; it cannot cross the open arc AB transversally. So
abandoning such strips loses no bead.

**Theorem (completeness).** Every bead of length < π (and every lune of
geodesics of length exactly π) in the link is found by `enumerate_beads`.

*Proof.* Let γ go from S (corner d₀ of t₀) to E through t₀, …, t_m, crossing
the shared edges e₁, …, e_m in their relative interiors. If m = 0, γ is an edge
of t₀: a thin bead, listed. Otherwise γ leaves S strictly inside t₀'s corner
and crosses e₁, the edge opposite d₀. In the development (§2) γ is a minor arc,
and by the lemma its initial direction w lies in the open interval I_k of
directions crossing e₁, …, e_k, for every k. The search begins with I₁ and,
from triangle t_k entered through e_k, tries both other edges of t_k — both
contain t_k's new vertex, and γ leaves through one of them — keeping the strip
while I_k ≠ ∅; it contains w, so the branch following γ is never cut. At t_m,
E is t_m's new vertex and γ reaches it from inside, so E's direction lies
strictly inside I_m (by the lemma applied to the edges crossed): the bead is
recorded. If E = −S every direction in I_m gives a geodesic of length π: the
lune is recorded (with `pi_beads`). ∎

Termination: a branch is cut when I_k = ∅ or a vertex is parallel to S; a
strip longer than 4× the number of link triangles raises an error rather than
being silently dropped (never observed). Evidence beyond the proof:
`tests/brute_beads.py` (exhaustive strips up to 7 triangles agree exactly) and
Elder–McCammond's `cat.g` bead catalogue for the disphenoid and cube shapes.

## 4. Checking a bead (kernel.check_bead)

Given a strip with gluing steps and endpoints S, E (not parallel), the kernel
develops the strip (§2) and requires, for each shared edge AB, that the minor
arc SE crosses the open arc AB: with D = det, the intersection line of the
planes (S,E) and (A,B) is X = S·D(A,B,E) − E·D(A,B,S) (up to a scalar);
X lies on the open arc AB iff its coefficients in A, B have equal non-zero
sign, and on the open arc SE iff its coefficients in S, E do. Crossing every
shared edge in its interior, in a strip of convex triangles, puts γ inside the
strip, so γ is a bead. The kernel also requires that the strip never goes straight
back across the edge it just crossed. Then the crossings come in strip
order: by induction, S lies on triangle k's side of edge k, the segment
meets edge k once (at parameter τ_k), the next edge's interior lies strictly
on the far side (the next triangle is developed across edge k), so
τ_{k+1} > τ_k, and S lies on triangle k+1's side of edge k+1 (Lean:
`checkBead_order`). For a π-bead the half great circle through the given
point Q is tested instead (Q on S's side: Y·(Q − (Q·S/S·S)S) > 0).

## 5. Junction angles (kernel.junction_side)

At a link vertex W (an edge of T) the triangles around W form a circle of
total angle Σ corner angles = the edge's dihedral sum (≥ 2π by (ii)). The
kernel sums, exactly, the angle from the reversed incoming direction u to the
first edge, the corner angles of the triangles passed, and the angle from the
last edge to the outgoing direction x — on each of the two sides. Each corner
angle is intrinsic to one triangle (its Gram matrix), and each partial angle
is computed in the frame where both rays are known. When incoming and outgoing
corners coincide, the kernel first checks whether x lies between u and the
first edge (then that side's angle is ∠(u, x) < π). By (L), the junction is a
local geodesic iff both sums are ≥ π.

## 6. Certificates

**False (vertex).** The kernel checks each bead (§4), that consecutive beads
meet, each junction (§5) and that the total length is < 2π (§1). The beads
then form a closed local geodesic of length < 2π, so by (C) the vertex link
is not CAT(1) and M is not locally CAT(0) by (G). **False (edge, invalid)**
are direct computations (§1, §2).

**True.** The kernel checks validity and every edge itself, enumerates every
bead itself (§3), and for each vertex link verifies the shortest-path trees:
for each bead s, values D_s(x) (exact sums of bead lengths along tree chains)
with D_s(s) = len(s), such that for every pair x → y with a passing junction
(junctions involving a π-bead count as passing):

1. y < s, or D_s(x) + len(y) ≥ 2π, or y is in the tree with
   D_s(y) ≤ D_s(x) + len(y);
2. if y = s: D_s(x) ≥ 2π.

*Soundness.* Suppose a necklace b₀ → b₁ → … → b_k → b₀ of length < 2π exists;
let s = b₀ be its lowest-numbered bead, so every b_i ≥ s. By induction on i,
b_i is in the tree with D_s(b_i) ≤ len(b₀) + … + len(b_i): true for i = 0, and
condition 1 applied to b_i → b_{i+1} (its partial length is < 2π) gives it for
i+1. Then condition 2 for b_k → s says D_s(b_k) ≥ 2π, while
D_s(b_k) ≤ total length < 2π — a contradiction. With conditions (i), (ii) and
the completeness theorem, the reduction in §0 shows M is locally CAT(0).
Counting π-bead junctions as passing only adds edges to the graph, so it can
only make the certificate harder to satisfy, never unsound.

## 7. What remains trusted

- The cited results in §0 (and, with `pi_beads=False`, EMM Lemma 4.9).
  In the Lean formalisation (`lean/VLink/Soundness.lean`) these are five
  axioms about the opaque predicate `LocallyCAT0`; the soundness of both
  certificate checks is proved from them (`checkFalse_sound`,
  `checkTrue_sound`), together with most of this note — §1 (arithmetic, angle
  sums), §2 (development), §4 (bead checks), §5 (junction sums), §6 (trees) —
  as Lean theorems. Still argued only here: bead
  enumeration completeness (§3); positive definiteness (§2), the fan walk
  (§5) and the crossing order (§4) are now Lean theorems too.
- That `tri.py` decodes isosigs as Regina does (tested, §2: the statement
  "isosig X with squared lengths L" refers to Regina's triangulation and edge
  numbering). Regina itself is no longer used by the kernel.
- Python itself and the ~1,350 lines of `kernel.py`, `field.py` and `tri.py`.

_Theorem numbers in Bridson–Haefliger to be checked against the book before citing._
