# Negative controls: why there are none, and what replaces them

A negative control would be a metric whose correct answer is known to be
False for topological reasons: a cube or disphenoid metric satisfying the
edge condition on a closed manifold that admits no locally CAT(0) metric at
all (spherical, S²×S¹, Nil, Sol or SL₂~ geometry; `controls/topology.py`).
The checker would have to answer False on it, independently of its own
geometry.

## None exist for these shapes

- **Census.** The full census scan (`census_full.md`) enumerates every cube
  and disphenoid metric satisfying the edge condition on all 61,848 closed
  and Hodgson–Weeks census triangulations, exhaustively. All 38 lie on flat,
  hyperbolic, H²×R or graph manifolds; none on the classes above (6,407
  spherical, 1 S²×S¹, 337 Nil, 268 Sol and 15,846 SL₂~ triangulations in the
  closed census).
- **Pachner walks** from census triangulations (`census_walks*`) reach new
  triangulations carrying shape metrics, all on flat, H²×R, hyperbolic or
  other Seifert manifolds (e.g. SFS [RP²: (3,1) (3,−1)], orbifold Euler
  characteristic −1/3, Euler number 0: H²×R), where the answer is not forced.
- **Generic lengths.** A search over 390 triangulations of these manifolds
  with generic squared lengths near the two shapes
  (`controls/generic_controls.py`) found no metric satisfying the edge
  condition; the best margins approach 0 from below (a question for C9).

## A counting constraint

Both shapes have dihedral angles summing to 7π/3 per tetrahedron (disphenoid:
2·π/2 + 4·π/3; cube orthoscheme: 3·π/2 + π/4 + π/3 + π/4). A closed
triangulated 3-manifold with n tetrahedra and v vertices has n + v edges
(Euler characteristic 0, 2n faces). Summing the edge condition over all
edges,

  (7π/3) n ≥ 2π (n + v),  so  n ≥ 6v,

with equality exactly when every edge has dihedral sum 2π, i.e. the metric
is flat. Every edge also needs degree ≥ 4 (no dihedral angle exceeds π/2).
Minimal census triangulations rarely meet these conditions — most fail at
once on an edge of degree ≤ 3 — which is why shape metrics are so scarce.

## What validates False answers instead

1. Every False answer comes with a necklace certificate, checked by two
   independent exact kernels (Python `vlink/kernel.py`; Lean `vlinkcheck`,
   whose arithmetic, geometry and bead checks are proved) and by the
   independent float checker `verify_cert.py`.
2. Wrong certificates are rejected: the 2020 code's witnesses, the
   junction-bug witnesses, and corrupted certificates (80 in all,
   `tests/cert_corpus.py`).
3. Bead enumeration is complete against brute force (`tests/brute_beads.py`,
   strips up to 7 triangles) and agrees with Elder–McCammond's `cat.g`.
4. Every True answer has a certificate too, so neither kind of answer
   depends on the search being right.

Decision (B8): no negative controls of this kind are used; the checks above
replace them, and the scarcity is recorded here and in `census_full.md`.
