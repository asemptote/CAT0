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
| L3 link data | `VLink/Kernel.lean` | proved: every link triangle's Gram matrix is positive definite (`link_gram_posdef`, by Sylvester's criterion from the leading minors the kernel checks); the link's gluing is consistent (`glueOk`, checked) |
| L3 development | `VLink/Develop.lean`, `VLink/Kernel.lean` | proved: each step places the neighbouring link triangle isometrically, on the far side of the shared edge (`step_iso`, `stepWith_iso`), so in the frame of a strip's first triangle every developed triangle has exactly its own Gram matrix (`develop_iso`) |
| L3 junctions | `VLink/Kernel.lean` | proved: the exact angle at each side of a junction is the sum of the corner angles listed by the fan walk (`junctionFan_ok`), each the true tangent angle (`tangent_angle`); the walk turns around the junction's link vertex one glued triangle at a time (`fan_step`), and every angle it adds is a corner at that vertex (`fanWalk_corners`) |
| L3 beads in order | `VLink/Develop.lean`, `VLink/Kernel.lean` | proved: an accepted bead's segment meets its strip's interior edges at strictly increasing parameters, so it traverses the strip in order (`order_step`, `checkBead_order`); this needs the strip never to turn straight back across the edge it just crossed, which both kernels now check |
| L3 bead enumeration | `VLink/Complete.lean` | proved: every bead the kernel accepts is found by its own enumeration (`enumerateBeads_complete`; a π-bead as one with the same strip and ends): the bead's direction lies strictly inside the interval of directions toward each edge it crosses (`visible_crossing`) and hence inside their intersections (`intersect_within`), so the search never abandons its strip, and the strip ends at the new vertex of its last triangle (`end_off_edge`) |
| L4 soundness | `VLink/Soundness.lean` | proved from the cited theorems: an accepted False certificate implies `¬ LocallyCAT0` (`checkFalse_sound`), an accepted True certificate implies `LocallyCAT0` (`checkTrue_sound`) |

`LocallyCAT0 isosig sqlens` is opaque (Mathlib has no piecewise Euclidean
complexes); five axioms, stated in terms of what the kernel computes,
characterise it (`VLink/Soundness.lean`): a degenerate tetrahedron or a short
edge rules it out; a short necklace (`ShortNecklace`: accepted beads,
junctions `≥ π`, length `< 2π`) rules it out; valid tetrahedra, edges `≥ 2π`
and no short necklace (with or without π-beads) imply it. Both directions use
the same `ShortNecklace`. A True certificate's trees show that no necklace of
the *enumerated* beads is short (`NoShortNecklace`); by completeness that
means no short necklace at all (`not_short_of_enum`). `#print axioms`
confirms the two soundness theorems use exactly these and Lean's standard
axioms (`propext`, `Classical.choice`, `Quot.sound`), no `sorry`.

The axioms are justified by the cited theorems (Gromov's link condition,
Bridson–Haefliger; Bowditch; Elder–McCammond; Elder–McCammond–Meier) together
with facts proved here — exact arithmetic and angle sums, positive definite
Gram matrices (`link_gram_posdef`), isometric development (`develop_iso`),
geodesics crossing their strips' edges in order (`checkBead_ok`,
`checkBead_order`), lengths and corner angles (`length_angle`,
`tangent_angle`), junction sums over the fan around the link vertex
(`junctionFan_ok`, `fanWalk_corners`), complete bead enumeration
(`enumerateBeads_complete`), shortest-path trees (`checkVertexTrees_ok`).
Nothing in `docs/proofs.md` is left to paper alone. What the axioms add to
the literature is that the kernel's accepted beads and junctions are the
beads and junctions of the cited theorems, which is what the Lean
descriptions of `checkBead` and `junctionFan` above say.

No proof uses linear independence of square roots or any normal form: zero
tests go through the sign, and the sign algorithm is sound for any
representation. Every check in `vlinkcheck` is exact (Python's
`exact_only=True`).

## Build (orr)

The toolchain and Mathlib cache live on `/import/taiji1` (not in `~`):

```
export ELAN_HOME=/import/taiji1/wardak/cat0-vertex-link/tools/elan PATH=/import/taiji1/wardak/cat0-vertex-link/tools/elan/bin:$PATH
export MATHLIB_CACHE_DIR=/import/taiji1/wardak/cat0-vertex-link/tools/mathlib-cache
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

- Nothing outstanding in the kernel proofs. Still trusted outside Lean: that
  `Tri.lean` (and `tri.py`) decode isosigs as Regina does (tested).
