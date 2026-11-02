# Elder–McCammond's `cat.g` (GAP): independent bead catalogue

`cat.g` (version 0.7, March 2002) is Murray Elder and Jon McCammond's GAP
program for testing the link condition of piecewise Euclidean / spherical
complexes built from finitely many shapes. It enumerates, per shape, the
galleries that can carry closed geodesics ("vertex" galleries — our beads —
and Möbius galleries). It is not included here (no licence is stated); get it
from Jon McCammond's software page:

    http://web.math.ucsb.edu/~jon.mccammond/software/cat.html
    http://web.math.ucsb.edu/~jon.mccammond/software/catFiles/cat.g

A copy is kept outside the repository (`working/reference/gap/cat.g` in the
project folder; on taiji1, `cat0-vertex-link/gap-cat/cat.g`). The only change
made to run it under GAP 4: its first line `Exec("clear");` is commented out.

Run (GAP 4 from conda-forge, environment `~/envs/gap` on orr), with `cat.g` in
this folder:

    gap -q run_cat.g      # galleries for the disphenoid and the cube orthoscheme -> run_cat.log
    gap -q lengths_cat.g  # geodesic lengths of the cube's vertex galleries

What it confirmed (`docs/true_certificates.md`): for the disphenoid, one
vertex-gallery type (2 tetrahedra, length π/2), matching every thick bead
vlink finds in disphenoid links; for the cube orthoscheme, 2-tetrahedron
types of length π/3 and 2π/3 (matching vlink) and 3-tetrahedron types of
length π/2 that do not occur in the links examined (confirmed by
`tests/brute_beads.py`).
