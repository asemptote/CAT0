# vlink in Lean 4

A Lean 4 / Mathlib formalisation of vlink's trusted kernel (`vlink/kernel.py`,
`field.py`, `tri.py`), and `vlinkcheck`, an executable checker for vlink's
certificates built on it.

| Layer | File | Status |
|---|---|---|
| L1 exact arithmetic | `VLink/MQ.lean` (Σ c √s), `VLink/Tower.lean` (nested √) | proved: every operation computes the real value; `sign` is sound (when it answers, it is the sign of the value) |
| L2 angle sums | `VLink/Angle.lean` | proved: `AngleSum.add` tracks Σθ exactly with its half-turn count; `Θ ≥ kπ ↔ h ≥ k`; `leq` decides `Θ₁ ≤ Θ₂` |
| L3 combinatorics | `VLink/Tri.lean` | executable (isosig decoding, skeleton, manifold check) |
| L3 geometry, L4 certificates | `VLink/Kernel.lean`, `VLink/Cert.lean` | executable; built on the proved L1–L2 |
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

## Next

- Specifications for the geometric predicates (`crosses`, `parallel`,
  `tangentAngle`, development) as statements about real vectors, proved
  from L1.
- The certificate theorems: `checkFalse c = ok → ¬ LinkCondition` and
  `checkTrue c = ok → LinkCondition`, with the cited results as named axioms.
