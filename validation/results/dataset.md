# Dataset (validation/results/dataset.jsonl)

Every metric examined: 109 distinct (isosig, squared edge lengths), from the 2020
candidate files, the full census scan (`census_full.md`), Pachner walks and the
controls. Every answer is certified in both settings (`pi_beads` True and False):
the certificate is accepted by the Python kernel and by `vlinkcheck` (Lean);
certificates are stored on `/import/taiji1/wardak/certificates/allm2`.

| answer | reason | flat | metrics |
|---|---|---|---|
| not CAT(0) | vertex | no | 21 |
| locally CAT(0) |  | no | 70 |
| locally CAT(0) |  | yes | 18 |

Geometry of the manifold (where recorded): hyperbolic 74, flat 18, other 8, graph manifold 6, H2xR 2, unknown 1.

Tetrahedra per triangulation: 6 to 18.
