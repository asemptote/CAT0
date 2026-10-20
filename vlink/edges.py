"""Exact edge and tetrahedron checks for a PE metric with rational squared lengths.

- validity: every tetrahedron is a non-degenerate Euclidean tetrahedron
  (Cayley-Menger determinants of the tetrahedron and its four faces > 0);
- edge condition: the dihedral angles around every edge sum to >= 2*pi
  (the link of every edge has no closed geodesic shorter than 2*pi);
- flatness: every edge sum is exactly 2*pi.

All exact. A dihedral angle theta of a tetrahedron at edge (i, j) is the angle
between the other two edge vectors u, w at i projected perpendicular to the
edge: with u', w' those projections, a = u'.w' is rational (Gram matrix of
rational squared lengths) and b = |u' x w'| = sqrt(|u'|^2 |w'|^2 - a^2) is the
square root of a rational, so a + i*b = |u'||w'| e^{i theta} lies in the field
Q(i, sqrt q, ...) exactly. Sums are compared with 2*pi in floating point, and
near-ties are decided exactly from the argument of the product (as in
necklace.cmp_two_pi).
"""

import math
from fractions import Fraction

import regina

from .field import MQ

TWO_PI = 2 * math.pi
EPS = 1e-9


def _det(M):
    n = len(M)
    M = [list(r) for r in M]
    det = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return Fraction(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            det = -det
        det *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            for k in range(c, n):
                M[r][k] -= f * M[c][k]
    return det


def _gram(d, base, others):
    """Gram matrix of the edge vectors from `base` to `others`."""
    return [[(d[(base, x)] + d[(base, y)] - (d[(x, y)] if x != y else 0)) / 2
             for y in others] for x in others]


def tet_sqlens(tet, sqlens):
    d = {}
    for k in range(6):
        m = tet.edgeMapping(k)
        d[(m[0], m[1])] = d[(m[1], m[0])] = Fraction(sqlens[tet.edge(k).index()])
    return d


def tet_valid(d):
    """All four faces and the tetrahedron itself have positive Gram determinant
    (i.e. positive area / volume)."""
    if _det(_gram(d, 0, (1, 2, 3))) <= 0:
        return False
    for f in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
        if _det(_gram(d, f[0], f[1:])) <= 0:
            return False
    return True


def dihedral_complex(d, i, j):
    """(a, b2): a + i*sqrt(b2) is a positive multiple of e^{i theta}, theta the
    dihedral angle at edge (i, j); a, b2 rational, b2 >= 0."""
    k, l = [x for x in range(4) if x not in (i, j)]
    G = _gram(d, i, (j, k, l))   # rows/cols: e = j, u = k, w = l
    ee, eu, ew = G[0][0], G[0][1], G[0][2]
    uu, uw, ww = G[1][1], G[1][2], G[2][2]
    # u' = u - (u.e/e.e) e, w' likewise
    a = uw - eu * ew / ee
    nu = uu - eu * eu / ee
    nw = ww - ew * ew / ee
    return a, nu * nw - a * a


def edge_report(isosig, sqlens):
    """Exact validity, edge condition and flatness.

    Returns dict(valid, edge_condition, flat, failing_edges, sums) where sums
    are floats (dihedral sum / 2*pi per edge) for information only.
    """
    T = regina.Triangulation3.fromIsoSig(isosig)
    if len(sqlens) != T.countEdges():
        raise ValueError('expected %d squared lengths' % T.countEdges())
    ds = []
    for tet in T.tetrahedra():
        d = tet_sqlens(tet, sqlens)
        if not tet_valid(d):
            return dict(valid=False, edge_condition=None, flat=None,
                        failing_edges=None, degenerate_tet=tet.index())
        ds.append(d)
    failing, flat, sums = [], True, []
    for e in T.edges():
        total, re, im = 0.0, MQ.rational(1), MQ()
        for emb in e.embeddings():
            m = emb.vertices()
            a, b2 = dihedral_complex(ds[emb.tetrahedron().index()], m[0], m[1])
            total += math.atan2(math.sqrt(float(b2)), float(a))
            b = MQ.sqrt_rational(b2)
            re, im = re * a - im * b, re * b + im * a
        sums.append(total / TWO_PI)
        if abs(total - TWO_PI) <= EPS:
            c = im.sign()            # total = 2*pi + delta, sign(delta) = sign(Im)
        else:
            c = 1 if total > TWO_PI else -1
        if c < 0:
            failing.append(e.index())
        if c != 0:
            flat = False
    return dict(valid=True, edge_condition=not failing, flat=flat,
                failing_edges=failing, sums=sums)
