# The Elder–McCammond–Meier degree metric on 5/6 triangulations (2026-11-02)

Elder–McCammond–Meier (*Combinatorial conditions that imply
word-hyperbolicity for 3-manifolds*, Topology 42 (2003), Def. 5.1) give a
5/6-triangulation (every edge of degree 5 or 6) the metric with squared
length 4 on degree-5 edges and 3 on degree-6 edges, and prove it locally
CAT(0), with no flat planes in the universal cover, when no triangle has two
degree-5 edges (Thurston's 5/6\* condition). Their tetrahedra are then the
regular one, a "mixed" one and the disphenoid.

**Census.** Of the 61,848 census triangulations (closed census and
Hodgson–Weeks), 136 are 5/6-triangulations (115 distinct isosigs) and none
is 5/6\*: every one has a triangle with two degree-5 edges, outside the
theorem. Inputs: `emm56_inputs.jsonl` (with `triangles_two_deg5`).

**Results** (`linkcond(certify=True)`, `emm56_linkcond.jsonl`; certificates on
`/import/taiji1/wardak/cat0-vertex-link/certificates/emm56`). Every answer is
certified and accepted by the Python kernel and by `vlinkcheck` (Lean); no
disagreement.

| | locally CAT(0) | not CAT(0) |
|---|---|---|
| closed census (58) | 19 | 39 |
| Hodgson–Weeks (78) | 31 | 47 |
| distinct isosigs (115) | 39 | 76 |

- Every "not CAT(0)" answer fails the edge condition (an edge with dihedral
  sum < 2π), and every failing edge lies on a triangle with two degree-5
  edges (182 of degree 5, 5 of degree 6): exactly where the 5/6\* condition
  fails.
- All 39 locally CAT(0) metrics are non-flat (some edge has dihedral sum
  > 2π). They are on hyperbolic manifolds, as they must be: closed census
  Hyp_1.75712603 (Z_7), Hyp_1.82434432, Hyp_1.83193119, Hyp_1.88541473,
  Hyp_2.45616037, Hyp_2.46823220 (four homology types), Hyp_2.56897060 (five),
  Hyp_2.57516207; Hodgson–Weeks #20, 21, 32, 113–115, 133, 135, 138, 142, 145,
  331, 344, 356, 368, 369, 378, 954, 1025, 1027, 1211, 2263, 2753, 2790, 2914,
  2939, 2984, 3163, 3175, 5649, 7002.

So the degree metric is locally CAT(0) on a third of the census
5/6-triangulations even though the 5/6\* hypothesis fails; where it fails, it
fails only at edges next to a triangle with two degree-5 edges.
