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

## Where the edge condition fails

So on these 115 triangulations the metric is locally CAT(0) exactly when the
edge condition holds: no vertex link ever fails. The degree metric has eleven
tetrahedron types (by the set of degree-5 edges, up to symmetry); their
dihedral angles, in degrees:

| degree-5 edges | 5/6\* allowed | at degree-5 edges | at degree-6 edges |
|---|---|---|---|
| none (regular) | yes | | 70.53 ×6 |
| one (mixed) | yes | 75.52 | 65.91 ×4, 83.62 |
| two opposite (disphenoid) | yes | 90 ×2 | 60 ×4 |
| two adjacent | no | 70.14 ×2 | 62.47, 64.06, 78.22 ×2 |
| three, a triangle | no | 65.91 ×3 | 75.52 ×3 |
| three, at a vertex | no | 67.38 ×3 | 73.90 ×3 |
| three, a path | no | 64.16, 83.50 ×2 | 59.37 ×2, 71.79 |
| four, a triangle plus one | no | 62.47 ×2, 78.22, 79.66 | 70.14 ×2 |
| four, a 4-cycle | no | 76.66 ×4 | 57.42 ×2 |
| five | no | 60, 73.90 ×4 | 67.38 |
| six | no | 70.53 ×6 | |

A degree-5 edge needs its five angles to average at least 72°, a degree-6
edge its six to average at least 60°. With the three 5/6\* types this always
holds (Elder–McCammond–Meier Lemma 5.3); the other types bring angles below
72° at degree-5 edges and below 60° at degree-6 edges, which is where the
failures occur.

## More 5/6-triangulations from Pachner walks

`validation/controls/five_six_walks.py`: simulated annealing over 2-3, 3-2
and 4-4 moves from the 115 triangulations above (two walks of 200,000 steps
each, cost = distance of the edge degrees from {5, 6}), 1,073 s on 4 cores of
orr. It reached only 10 new 5/6-triangulations (12–16 tetrahedra; at least 4
triangles with two degree-5 edges each, so none 5/6\*): 5/6-triangulations
are isolated among their Pachner neighbours (a 2-3 move creates a degree-3
edge). All 10 fail the edge condition (certified, Lean agrees), so the
observation "locally CAT(0) ⇔ edge condition" holds on all 125.

No 5/6\* triangulation was found. Each vertex link of one is a triangulated
sphere with twelve pairwise non-adjacent degree-5 vertices, so at least 32
vertices (the dual of the truncated icosahedron), so a one-vertex 5/6\*
triangulation has at least 15 tetrahedra and very constrained vertex links;
neither the censuses nor these walks (which never got below 4 bad triangles)
contain one. A positive control would need a construction, not a search.
