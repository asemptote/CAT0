# gallery.py

A Python implementation of Elder and McCammond's method to check Gromov's link condition for vertices

Extensive documentation in code files

## Requirements
Python 3, Regina (regina-normal) 7
    
## Dependencies:
- SAGE, or
- [`constructible.py`](https://github.com/leovt/constructible): To use a different underlying constructible number object, simply
        change what `constrsphtrig.sqrt` refers to after importing. (If
        `constructible.py` is not available create a file with stub variables
        `sqrt` and `Constructible`.)

# vlink (rewrite, 2026)

`vlink/` re-implements the checker with exact arithmetic in multiquadratic
fields instead of nested constructible numbers (profiling showed 99.6% of
`gallery.py`'s time in that arithmetic). Every geometric decision is an exact
sign test; ties at length 2π and angle π are decided exactly.

```python
from vlink.check import linkcond, vertex_linkcond
linkcond('gLLPQaceedffnkxknnn', [1, 2, 3, 1, 1, 2, 2])   # validity + edges + vertices
# -> {'result': True, 'reason': None, 'flat': True}
```

- `field.py` exact Q(√q₁, √q₂, …) numbers; `link.py` exact development of
  link triangles; `beads.py` bead enumeration; `necklace.py` necklace search;
  `edges.py` exact tetrahedron validity and edge condition; `check.py` entry points.
- With `allow_repeats=False` the bead sets are identical to `gallery.py`'s and
  results agree on all validation inputs. The default `allow_repeats=True`
  drops `gallery.py`'s unproven rule that a necklace never reuses a link
  triangle; on all validation inputs both settings give the same answers.
- Necklace search is pruned exactly (shortest bead-path bounds between link
  vertices; each loop searched once): worst validation input 9 s, against
  25–50 min for `gallery.py`.
- Requirements: Python ≥ 3.10, `regina` ≥ 7.4 (PyPI wheels: Linux/macOS).
- Tests: `python tests/test_vlink_geometry.py`, `python tests/test_vlink_edges.py`.

## Validation

`validation/` holds the tooling used to check the checkers (results and a
comparison table in `validation/results/summary.md`):

- `verify_cert.py` independent float geometry: input validity, edge condition,
  flatness, and verification of every necklace witness;
- `run_check.py` runs an implementation (`current`, `vlink`, `vlink_repeats`,
  `legacy2020`) over the candidate files with resume support;
- `controls/` labelled controls from the closed census: `topology.py` classifies
  manifolds that cannot be CAT(0); `shape_metrics.py` enumerates all cube /
  disphenoid metrics satisfying the edge condition; Pachner-walk and
  generic-length searches for negative controls;
- `legacy/port_2020.py` reproduces the September 2020 checker under Python 3.

The candidate files themselves (`CAT0_Candidates/`) are not part of this repository.
