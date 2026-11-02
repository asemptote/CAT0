# Cube and disphenoid metrics on the whole census (2026-10-29)

Every triangulation of the closed census (`closed-or-census-11.dat`, 50,817)
and of the Hodgson–Weeks census (`hodgson-weeks.dat`, 11,031) was searched
for all metrics built from one tetrahedron shape — the cube orthoscheme or
the Sommerville disphenoid (`controls/shape_metrics.py`) — that satisfy the
edge condition. Every search finished (no limit or timeout reached); most end
at once, since an edge of degree ≤ 3 already rules a shape out. Run through
`validation/workqueue.py` (`tasks:shape_metrics_batch`, 1,240 tasks), then
each metric through `linkcond(certify=True)` with both kernels
(`tasks:run_check`, impl `linkcond_cert`).

Files: `census_full_metrics.jsonl` (the metrics), `census_full_linkcond.jsonl`
(answers, Python and Lean verdicts; certificates stored on
`/import/taiji1/wardak/cat0-vertex-link/certificates/census`).

| | metrics | locally CAT(0) | of which flat | not CAT(0) |
|---|---|---|---|---|
| closed census, cube | 13 | 13 | 11 | 0 |
| closed census, disphenoid | 20 | 19 | 7 | 1 (vertex: `lLLPALzQcaceedfhikjjknkxknnnxkkkk`) |
| Hodgson–Weeks, disphenoid | 5 | 5 | 0 | 0 |
| Hodgson–Weeks, cube | 0 | | | |

38 metrics on 25 triangulations. Every answer is certified and its certificate
accepted by the Python kernel and by `vlinkcheck` (Lean).

The search adds nothing new: the 33 closed-census metrics are those found
earlier (`census_shape_metrics.jsonl`), and the 5 Hodgson–Weeks ones are
`examples_3_disphenoid` #20–24 of the 2020 candidate files. So those
candidate files already list every cube or disphenoid metric on these
censuses. Locally CAT(0), non-flat metrics occur on hyperbolic manifolds
`Hyp_1.01494161`, `Hyp_1.82434432`, `Hyp_2.02988321` (two triangulations) in
the closed census and on HW #174, #2027, #2980, #6705, #7720.
