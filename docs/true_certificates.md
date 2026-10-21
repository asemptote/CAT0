# Design: certificates for True answers

A True answer says the metric is locally CAT(0). It rests on three claims:

1. every tetrahedron is valid and every edge's dihedral angles sum to ≥ 2π;
2. **(A)** the list of beads in each vertex link is complete;
3. **(B)** no cyclic chain of beads with junction angles ≥ π on both sides
   (a necklace) has total length < 2π.

Claim 1 is a finite exact computation the kernel already performs
(`Geometry.tet_valid`, `Geometry.edge_angle`); it needs no certificate.

## (A) Bead completeness: enumeration moves into the kernel

Proving completeness means re-enumerating, so a certificate would not make
checking cheaper. Instead bead enumeration becomes kernel code (today it is
`beads.py`). Its correctness rests on one fact, to be documented next to it:
seen from S, the relative interior of an edge AB (S not on its great circle)
is hit exactly by the rays from S strictly inside the angle at S of the
spherical triangle SAB, which is < π; a strip is extended only while the
intersection of these open intervals over its crossed edges is non-empty.

Evidence so far:

- `tests/brute_beads.py`: on five links, every bead of ≤ 5 triangles found by
  exhaustive strip enumeration with the kernel's exact bead test is exactly
  the set vlink finds (deeper check running).
- Elder–McCammond's `cat.g` (GAP) per-shape catalogue: disphenoid — one bead
  type (2 tetrahedra, length π/2), matching every thick bead vlink finds in
  disphenoid links; cube orthoscheme — 2-tetrahedron types of length π/3 and
  2π/3 (matching vlink) and 3-tetrahedron types of length π/2 that do not
  occur in the links examined (confirmed by the brute force).

Later optimisation: verify the bead catalogue once per set of shapes (as
`cat.g` does) and match occurrences per triangulation.

## (B) No short necklace: shortest-path certificates

Let H be the directed graph on beads with an edge x → y when x ends at the
link vertex where y starts and the junction x → y is ≥ π on both sides;
y has weight len(y). A necklace of length < 2π is a directed cycle of total
weight < 2π. Each cycle is examined from its lowest-numbered bead s.

The certificate gives, for each bead s, values D_s(x) for the beads x > s
(each D_s(x) represented as an explicit chain of beads whose lengths it sums;
or ∞). The kernel checks:

- D_s(s) ≤ len(s);
- for every pair x → y with matching link vertices and y > s: the junction
  fails (kernel), or D_s(y) ≤ D_s(x) + len(y);
- for every x with x → s possible: the junction fails, or D_s(x) ≥ 2π.

By induction along any path, D_s is a lower bound on path lengths from s,
so every cycle through s as its lowest bead has length ≥ 2π. Checking costs
O(beads × junction pairs) comparisons, independent of how long the search
took; the search (Dijkstra) produces the D_s.

## Status

Implemented: `kernel.enumerate_beads` (A); `kernel.check_true_certificate`,
`kernel.check_vertex_trees` and `vlink/certify.py` (B);
`linkcond(..., certify=True)` attaches and kernel-checks a certificate for True
answers; flags `pi_beads` (D1) and `exact_only` (D2). Junctions involving a
π-length bead are counted as passing (conservative) when `pi_beads=True`.

First results (flat cube, H²×ℝ cube, hyperbolic disphenoid): every True answer
certified in both π settings, also with `exact_only=True`; certificates 7–187 KB,
generation and checking 3–38 s each. With rational links (cube, disphenoid)
exact-only checking is barely slower than the float path.

## Decisions (adopted behind flags)

**D1. Drop beads of length exactly π, citing Elder–McCammond–Meier Lemma 4.9**
(a necklace with a bead of length exactly π is equivalent to one whose beads
are all shorter than π). This removes the most intricate part of the search
and of the certificates, and vlink's last `UndecidedError` case. Check before
adopting: no answer changes on the 146 metrics checked so far.

**D2. Rigorous floating-point intervals inside the kernel.** The comparisons in
(B) number up to ~10⁶ per link; doing all of them in exact arithmetic is slow.
A small kernel primitive — a float enclosure of each bead length with a proven
error bound, exact fallback when enclosures overlap — makes checking fast.
The alternative keeps the kernel float-free but makes certificate checking
much slower.

Recommendation: adopt both (D1 after the no-change check).
