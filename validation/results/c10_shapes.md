# C10: Coxeter tetrahedra and orthoschemes on the censuses (2026-11-02)

Shapes (`validation/controls/shape_metrics.py`):

- the three Euclidean Coxeter tetrahedra: C~3 (the cube orthoscheme, `cube`),
  A~3 (the disphenoid, `disphenoid`), both scanned before (`census_full.md`),
  and B~3 (`b3`: the alcove (0,0,0), (1,0,0), (1/2,1/2,0), (1/2,1/2,1/2),
  dihedral angles π/4, π/3, π/2);
- orthoschemes `orth-a-b-c`: a path with orthogonal steps of squared lengths
  a, b, c, for a, b, c ∈ {1, …, 4} with gcd 1, up to reversal, other than the
  cube (32 shapes).

None of these shapes has a dihedral angle above π/2, so an edge of degree
≤ 3 rules each out. The scan (`validation/shape_scan.py`) therefore ran on
the 36,150 triangulations whose edges all have degree ≥ 4: 28,417 of the
closed census (`closed-or-census-11.dat`, 50,817) and 7,733 of Hodgson–Weeks
(`hodgson-weeks.dat`, 11,031); list on
`/import/taiji1/wardak/cat0-vertex-link/work/data/c10/candidates.json`.

Result: **no metric** satisfying the edge condition, for any of the 33 shapes,
on any triangulation; every search finished (no limit or timeout reached).
149 s on 4 cores of orr. Summary lines: `c10_shape_scan.jsonl`.

Checks of the search: it reproduces every cube and disphenoid metric of
`census_full_metrics.jsonl` exactly, `orth-1-1-1` gives the same metrics as
`cube`, and `orth-2-2-2` (a non-exact shape, so the floating-point path) gives
the cube metrics scaled by 2.

Why so few: an orthoscheme with unequal steps is not the cell of a cube-like
subdivision of one shape (a box a×b×c splits into six orthoschemes, one per
order of the steps, so of up to six different shapes), and the B~3 alcove has
volume 1/24 of a unit cube, so a closed flat manifold tiled by it needs many
more tetrahedra than the closed census allows.
