"""Tetrahedron validity and the edge condition, decided by the kernel.

- validity: every tetrahedron is a non-degenerate Euclidean tetrahedron
  (kernel.Geometry.tet_valid: Gram determinants of the tetrahedron and its
  faces are positive);
- edge condition: the dihedral angles around every edge sum to >= 2*pi
  (kernel.Geometry.edge_angle, an exact AngleSum);
- flatness: every edge sum is exactly 2*pi.
"""

import math

from . import kernel

TWO_PI = 2 * math.pi


def _exactly_two_pi(acc):
    # total in [2*pi, 3*pi) with Im z = 0 means total == 2*pi exactly
    return acc.h == 2 and acc.im.is_zero()


def edge_report(isosig, sqlens):
    """dict(valid, edge_condition, flat, failing_edges, sums); sums are floats
    (dihedral sum / 2*pi per edge) for information only."""
    geo = kernel.Geometry(isosig, sqlens)
    for t in range(geo.T.n):
        if not geo.tet_valid(t):
            return dict(valid=False, edge_condition=None, flat=None,
                        failing_edges=None, degenerate_tet=t)
    failing, flat, sums = [], True, []
    for e in range(len(geo.T.edges)):
        acc = geo.edge_angle(e)
        if not acc.at_least(2):
            failing.append(e)
        if not _exactly_two_pi(acc):
            flat = False
        total = 0.0
        for tet, k in geo.T.edges[e]:
            a, b2 = geo.dihedral(tet, *kernel.EDGES[k])
            total += math.atan2(math.sqrt(max(0.0, float(b2))), float(a))
        sums.append(total / TWO_PI)
    return dict(valid=True, edge_condition=not failing, flat=flat,
                failing_edges=failing, sums=sums)
