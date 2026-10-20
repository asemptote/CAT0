"""Enumerate every metric on a triangulation built from one tetrahedron shape.

Shapes (squared edge lengths; dihedral angles are multiples of pi/12):
  cube        the orthoscheme of the unit cube: for a vertex path a-b-c-d,
              |ab|^2=|bc|^2=|cd|^2=1, |ac|^2=|bd|^2=2, |ad|^2=3
              (12 placements per tetrahedron: 24 paths up to reversal)
  disphenoid  the Sommerville tetrahedron with one pair of opposite edges of
              squared length 4 (dihedral pi/2) and four of squared length 3
              (dihedral pi/3)  (3 placements per tetrahedron)

A backtracking search places a shape on each tetrahedron so that shared edges
get equal lengths, and keeps only metrics satisfying the edge condition (total
dihedral angle >= 2*pi around every edge), checked exactly in units of pi/12.
"""

import itertools
import math
import os
import sys
import time

import regina

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from verify_cert import dihedral, embed_tet  # noqa: E402

UNIT = math.pi / 12
FULL = 24  # 2*pi in units of pi/12


def _units(d2):
    """Dihedral angles (in pi/12 units) of the tetrahedron with squared lengths d2."""
    p, _ = embed_tet(d2)
    out = {}
    for i, j in itertools.combinations(range(4), 2):
        u = dihedral(p, i, j) / UNIT
        assert abs(u - round(u)) < 1e-9, u
        out[(i, j)] = int(round(u))
    return out


def _placements(shape):
    """List of {(i,j): (sqlen, units)} for every placement of the shape."""
    pats = set()
    if shape == 'cube':
        for a, b, c, d in itertools.permutations(range(4)):
            if (a, b, c, d) > (d, c, b, a):
                continue
            sq = {}
            for (x, y), s in (((a, b), 1), ((b, c), 1), ((c, d), 1),
                              ((a, c), 2), ((b, d), 2), ((a, d), 3)):
                sq[tuple(sorted((x, y)))] = s
            pats.add(tuple(sorted(sq.items())))
    elif shape == 'disphenoid':
        for pair in ((0, 1), (0, 2), (0, 3)):
            opp = tuple(x for x in range(4) if x not in pair)
            sq = {e: (4 if e in (pair, opp) else 3)
                  for e in itertools.combinations(range(4), 2)}
            pats.add(tuple(sorted(sq.items())))
    else:
        raise ValueError(shape)
    out = []
    for pat in sorted(pats):
        sq = dict(pat)
        d2 = dict(sq)
        d2.update({(j, i): v for (i, j), v in sq.items()})
        u = _units(d2)
        out.append({e: (sq[e], u[e]) for e in sq})
    return out


PLACEMENTS = {s: _placements(s) for s in ('cube', 'disphenoid')}


def enumerate_metrics(isosig, shape, limit=100, timeout=30.0):
    """Returns (metrics, complete) where metrics is a list of squared-edge-length
    lists satisfying the edge condition, and complete says whether the search
    finished (not cut off by limit/timeout)."""
    T = regina.Triangulation3.fromIsoSig(isosig)
    n, ne = T.size(), T.countEdges()
    deg = [T.edge(e).degree() for e in range(ne)]
    # local edge (i,j) of tet t -> global edge index
    gl = []
    for t in range(n):
        tet = T.tetrahedron(t)
        m = {}
        for k in range(6):
            mp = tet.edgeMapping(k)
            m[tuple(sorted((mp[0], mp[1])))] = tet.edge(k).index()
        gl.append(m)
    # visit tetrahedra in BFS order so edges complete early
    order, seen = [], {0}
    queue = [0]
    while queue:
        t = queue.pop(0)
        order.append(t)
        for f in range(4):
            a = T.tetrahedron(t).adjacentSimplex(f)
            if a is not None and a.index() not in seen:
                seen.add(a.index())
                queue.append(a.index())
    sqlen = [0] * ne
    count = [0] * ne   # how many embeddings of each edge have been placed
    total = [0] * ne   # dihedral sum so far (pi/12 units)
    found, deadline = [], time.time() + timeout
    state = dict(stop=False)

    def rec(k):
        if state['stop']:
            return
        if k == n:
            found.append(list(sqlen))
            if len(found) >= limit:
                state['stop'] = True
            return
        if time.time() > deadline:
            state['stop'] = True
            return
        t = order[k]
        for pl in PLACEMENTS[shape]:
            changed, ok = [], True
            for e_loc, (s, u) in pl.items():
                e = gl[t][e_loc]
                if sqlen[e] and sqlen[e] != s:
                    ok = False
                    break
                changed.append((e, sqlen[e] == 0))
                sqlen[e] = s
                count[e] += 1
                total[e] += u
                if count[e] == deg[e] and total[e] < FULL:
                    ok = False
                    break
            if ok:
                rec(k + 1)
            for e_loc, (s, u) in list(pl.items())[:len(changed)]:
                e = gl[t][e_loc]
                count[e] -= 1
                total[e] -= u
            for e, was_new in changed:
                if was_new:
                    sqlen[e] = 0
            if state['stop']:
                return

    rec(0)
    return found, not state['stop']
