# vlink in Lean 4

A Lean 4 / Mathlib formalisation of vlink's trusted kernel (`vlink/kernel.py`,
`field.py`, `tri.py`), and `vlinkcheck`, an executable checker for vlink's
certificates built on it.

| Layer | File | Status |
|---|---|---|
| L1 exact arithmetic | `VLink/MQ.lean` (Σ c √s), `VLink/Tower.lean` (nested √) | proved: every operation computes the real value; `sign` is sound (when it answers, it is the sign of the value) |
| L2 angle sums | `VLink/Angle.lean` | proved: `AngleSum.add` tracks Σθ exactly with its half-turn count; `Θ ≥ kπ ↔ h ≥ k`; `leq` decides `Θ₁ ≤ Θ₂` |
| L3 combinatorics | `VLink/Tri.lean` | executable (isosig decoding, skeleton, manifold check) |
| L3 vector predicates | `VLink/Vec.lean` | proved: the exact operations compute the real inner product, determinant and cross product; `parallel` decides parallelism and the inner product's sign (`parallel_ok`); an accepted `crosses` means the geodesic meets the open edge arc at an interior point (`crosses_ok_minor`, `crosses_ok_half`); the angles added for a bead length and for a corner are the angle between the endpoints and between the tangent directions (`length_angle`, `tangent_angle`) |
| L3 beads | `VLink/Kernel.lean` | proved: an accepted bead crosses every interior edge of its developed strip (`checkBead_ok`); exact angle sums represent the sum of their angles (`angleSum_ok`) |
| L4 True certificates | `VLink/Necklace.lean`, `VLink/Kernel.lean` | proved: if every tree's values satisfy the checked inequalities, every necklace (a `Cycle` of beads, each passing to the next) has length `≥ 2π` (`no_short_necklace`); the kernel's tree check establishes exactly those inequalities for the exact angle sums, so when it passes no necklace at that vertex is shorter than `2π` (`checkVertexTrees_ok`); an exact angle sum represents a single angle (`AngleSum.rep_unique`) |
| L3 development, junctions, L4 False certificates | `VLink/Kernel.lean`, `VLink/Cert.lean` | executable; built on the proved layers |
| bridge to CAT(0) | `docs/proofs.md` | on paper: cited results (Gromov, Bowditch, Elder–McCammond–Meier) and the kernel's geometric lemmas |

No proof uses linear independence of square roots or any normal form: zero
tests go through the sign, and the sign algorithm is sound for any
representation. Every check in `vlinkcheck` is exact (Python's
`exact_only=True`).

## Build (orr)

The toolchain and Mathlib cache live on `/import/taiji1` (not in `~`):

```
export ELAN_HOME=/import/taiji1/wardak/elan PATH=/import/taiji1/wardak/elan/bin:$PATH
export MATHLIB_CACHE_DIR=/import/taiji1/wardak/mathlib-cache
cd /import/taiji1/wardak/cat0-vertex-link/CAT0/lean
lake exe cache get    # prebuilt Mathlib
lake build            # proofs + vlinkcheck
```

## Use

```
.lake/build/bin/vlinkcheck certificates.jsonl
```

One certificate per line (or a record with a `certificate` field, as
`validation/run_check.py` writes); prints `ACCEPT`/`REJECT`, the claim
(`True`/`False`) and the kernel's message.

Cross-check against the Python kernel:

```
python tests/cert_corpus.py corpus.jsonl --workers 8
lean/.lake/build/bin/vlinkcheck corpus.jsonl > lean.txt
python tests/cert_corpus.py --compare corpus.jsonl lean.txt
```

(or per certificate through `validation/workqueue.py` with
`tasks:lean_verdict`, as done on 2026-10-22). Result: on all 169
certificates — 15 genuine False and 74 genuine True (the 37 regression
metrics answered True, both π settings), 79 wrong ones (the 2020 code's
witnesses, the junction-bug witnesses, corrupted False and True
certificates, a True certificate for a metric that is not CAT(0)) —
`vlinkcheck` and the Python kernel (`exact_only=True`) reach the same
verdict. Lean took 373 CPU-seconds in all (at most 26 s for one True
certificate); the Python kernel, exact throughout, 78 CPU-minutes.

`vlinkcheck` is also vlink's reference checker: `vlink.check.linkcond` runs
it on every certificate when it is built (`vlink/lean.py`) and raises an
error if it disagrees with the Python kernel. On the 109 validation metrics
(both π settings, 218 certificates) and the 38 census metrics it agrees
everywhere.

## Next

- Development: `step` places the next link triangle isometrically (the
  developed vectors have the triangle's Gram matrix), from the transition
  data; and the Gram matrix of a valid tetrahedron is positive definite
  (so `parallel_ok` reads as "same / opposite direction", and the symmetric
  hypotheses hold).
- Junctions: `junctionSide` sums the corner angles around the link vertex
  (`tangent_angle` for each), so `≥ π` is the local-geodesic condition; then
  a passing entry of `junctionRow` (`junctionRow_mem`) means the necklace is
  a local geodesic at that junction.
- The certificate theorems: `checkFalse c = ok → ¬ LinkCondition` and
  `checkTrue c = ok → LinkCondition`, with the cited results (Gromov,
  Bowditch, Elder–McCammond–Meier, bead enumeration completeness) as named
  axioms.
