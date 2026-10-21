"""vlink.tri (the kernel's own isosig decoder and skeleton) against Regina.

    python tests/test_tri.py [census.dat ...]

For every triangulation: the gluing table equals Regina's, and edges and
vertices are numbered as Regina numbers them. With census files (lines
'isosig NAME ...') every triangulation in them is checked; otherwise the
regression inputs.
"""

import json
import os
import sys

import regina

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from vlink import tri  # noqa: E402


def compare(sig):
    T = regina.Triangulation3.fromIsoSig(sig)
    mine = tri.Triangulation(sig)
    assert mine.n == T.size(), sig
    for t in range(T.size()):
        tet = T.tetrahedron(t)
        for f in range(4):
            g = tet.adjacentGluing(f)
            assert mine.adj[t][f] == tet.adjacentTetrahedron(f).index(), (sig, t, f)
            assert mine.glue[t][f] == tuple(g[i] for i in range(4)), (sig, t, f)
        for k in range(6):
            assert mine.edge_of[t][k] == tet.edge(k).index(), (sig, t, k)
            m = tet.edgeMapping(k)
            assert {m[0], m[1]} == set(tri.EDGES[k]), (sig, k)
        for v in range(4):
            assert mine.vertex_of[t][v] == tet.vertex(v).index(), (sig, t, v)
    assert len(mine.edges) == T.countEdges() and mine.n_vertices == T.countVertices()
    for e in range(T.countEdges()):
        assert len(mine.edges[e]) == T.edge(e).degree(), (sig, e)


def test_regression_inputs():
    data = json.load(open(os.path.join(ROOT, 'tests', 'data', 'regression_cases.json')))
    sigs = sorted({c['isosig'] for c in data['cases']})
    for sig in sigs:
        compare(sig)
    return len(sigs)


def _random_closed(rng, n):
    """A random closed (no boundary facets) triangulation with n tetrahedra:
    a random pairing of the 4n faces with random gluing maps."""
    import itertools
    perms = list(itertools.permutations(range(4)))
    while True:
        T = regina.Triangulation3()
        T.newTetrahedra(n)
        faces = [(t, f) for t in range(n) for f in range(4)]
        rng.shuffle(faces)
        ok = True
        for (t, f), (u, h) in zip(faces[::2], faces[1::2]):
            if (t, f) == (u, h):
                ok = False
                break
            p = rng.choice([p for p in perms if p[f] == h])
            T.tetrahedron(t).join(f, T.tetrahedron(u), regina.Perm4(*p))
        if ok and T.isConnected():
            return T


def test_random_triangulations(trials=3000):
    """Random closed triangulations, manifold or not: tri accepts exactly the
    valid ones whose vertex links are all spheres (Regina: valid and closed),
    and agrees with Regina on those."""
    import random
    rng = random.Random(1)
    accepted = rejected = 0
    for _ in range(trials):
        T = _random_closed(rng, rng.randint(1, 4))
        sig = T.isoSig()
        manifold = T.isValid() and T.isClosed()
        try:
            tri.Triangulation(sig)
        except tri.Invalid:
            assert not manifold, sig
            rejected += 1
            continue
        assert manifold, sig
        compare(sig)
        accepted += 1
    assert accepted and rejected
    return accepted, rejected


def census(paths):
    n = 0
    for path in paths:
        for line in open(path):
            parts = line.split()
            # closed-or-census: "isosig NAME : #k"; hodgson-weeks: "volume_k isosig"
            sig = parts[1] if "_" in parts[0] else parts[0]
            compare(sig)
            n += 1
            if n % 5000 == 0:
                print(n, flush=True)
    return n


if __name__ == '__main__':
    if len(sys.argv) > 1:
        print('ok, census triangulations agree with Regina:', census(sys.argv[1:]))
    else:
        print('ok, regression triangulations agree with Regina:', test_regression_inputs())
        print('ok, random closed triangulations (manifolds, non-manifolds):',
              test_random_triangulations())
