"""Closed triangulated 3-manifolds from isomorphism signatures (trusted).

The kernel's combinatorics, with no dependence on Regina:

  * decode: Regina's isomorphism signatures (Regina 7.4,
    engine/triangulation/detail/isosig-impl.h, fromIsoSig) into a gluing
    table: adj[t][f] is the tetrahedron glued to face f of t (the face
    opposite vertex f) and glue[t][f] the vertex map (a tuple p, vertex i of
    t -> vertex p[i] of adj[t][f]);
  * skeleton: edges and vertices numbered as Regina numbers them, in order
    of first appearance over (tetrahedron, local index), local edge k joining
    EDGES[k] (tests/test_tri.py checks both against Regina on the census);
  * manifold check: closed, every edge valid (not identified with itself
    reversed), every vertex link a 2-sphere;
  * vertex links: link triangle = corner (t, v) of a tetrahedron at the
    vertex, numbered by sorting the corners; its vertices 0, 1, 2 are the
    tetrahedron's other vertices in increasing order (see corner_perm).
"""

from itertools import permutations

_B64 = {c: i for i, c in enumerate(
    'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789+-')}
PERMS4 = sorted(permutations(range(4)))     # Regina's Perm<4>::orderedSn
IDENTITY = (0, 1, 2, 3)
EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))   # Regina's edge numbering
EDGE_INDEX = {e: k for k, e in enumerate(EDGES)}
EDGE_INDEX.update({(b, a): k for (a, b), k in list(EDGE_INDEX.items())})


class Invalid(Exception):
    """Not the isosig of a closed triangulated 3-manifold."""


def inverse(p):
    q = [0] * len(p)
    for i, x in enumerate(p):
        q[x] = i
    return tuple(q)


def decode(sig):
    """(adj, glue) of the triangulation with isomorphism signature sig."""
    s = sig.strip()
    pos = 0

    def single():
        nonlocal pos
        if pos >= len(s) or s[pos] not in _B64:
            raise Invalid('incomplete or invalid base64 encoding')
        pos += 1
        return _B64[s[pos - 1]]

    def integer(nchars):
        v = 0
        for i in range(nchars):
            v |= single() << (6 * i)
        return v

    adj, glue = [], []
    while pos < len(s):
        n = single()
        nchars = 1
        if n == 63:
            nchars = single()
            n = integer(nchars)
        if n == 0:
            continue
        actions, nfacets, njoins = [], 0, 0
        while nfacets < 4 * n:
            c = single()
            for j in range(3):
                a = (c >> (2 * j)) & 3
                if nfacets == 4 * n:
                    if a:
                        raise Invalid('extraneous facet actions')
                    continue
                if a == 3:
                    raise Invalid('invalid facet action')
                nfacets += 1 if a == 0 else 2
                njoins += a == 2
                if nfacets > 4 * n:
                    raise Invalid('facet actions do not match triangulation size')
                actions.append(a)
        dest = [integer(nchars) for _ in range(njoins)]
        perm = [integer(1) for _ in range(njoins)]
        if pos < len(s) and s[pos] == '.':
            raise Invalid('lock masks are not supported')
        base = len(adj)
        adj += [[None] * 4 for _ in range(n)]
        glue += [[None] * 4 for _ in range(n)]
        act = iter(actions)
        nxt, jp = 1, 0
        for p in range(n):
            for f in range(4):
                if adj[base + p][f] is not None:
                    continue
                a = next(act)
                if a == 0:
                    raise Invalid('boundary facet: the triangulation is not closed')
                if a == 1:
                    if nxt >= n:
                        raise Invalid('gluing to non-existent simplex')
                    q, g = nxt, IDENTITY
                    nxt += 1
                else:
                    if perm[jp] >= 24:
                        raise Invalid('invalid gluing permutation')
                    q, g = dest[jp], PERMS4[perm[jp]]
                    jp += 1
                    if q >= nxt or adj[base + q][g[f]] is not None or (q == p and g[f] == f):
                        raise Invalid('invalid gluing destination')
                adj[base + p][f], glue[base + p][f] = base + q, g
                adj[base + q][g[f]], glue[base + q][g[f]] = base + p, inverse(g)
    if not adj:
        raise Invalid('empty triangulation')
    return adj, glue


class Triangulation:
    """A closed triangulated 3-manifold given by its isosig."""

    def __init__(self, sig):
        self.isosig = sig
        self.adj, self.glue = decode(sig)
        self.n = len(self.adj)
        self._edges()
        self._vertices()

    def _edges(self):
        """edges[e]: sorted list of (t, k), local edge k of tetrahedron t;
        edge_of[t][k] = e. end[(t, v, w)] = (e, 0 or 1) tells which end of its
        edge e the vertex v of the edge {v, w} of t is (consistently over e's
        embeddings). An edge identified with itself reversed is invalid."""
        self.edge_of = [[None] * 6 for _ in range(self.n)]
        self.edges, self.end = [], {}
        for t in range(self.n):
            for k in range(6):
                if self.edge_of[t][k] is not None:
                    continue
                e = len(self.edges)
                seen = {(t, EDGES[k])}
                stack, members = [(t, EDGES[k])], []
                self.edge_of[t][k] = e
                while stack:
                    u, (a, b) = stack.pop()
                    members.append((u, EDGE_INDEX[(a, b)]))
                    self.end[(u, a, b)], self.end[(u, b, a)] = (e, 0), (e, 1)
                    for f in range(4):
                        if f in (a, b):
                            continue
                        g = self.glue[u][f]
                        w, ab = self.adj[u][f], (g[a], g[b])
                        if (w, ab[::-1]) in seen:
                            raise Invalid('edge %d is identified with itself in reverse' % e)
                        if (w, ab) not in seen:
                            seen.add((w, ab))
                            self.edge_of[w][EDGE_INDEX[ab]] = e
                            stack.append((w, ab))
                self.edges.append(sorted(members))

    def _vertices(self):
        """vertex_of[t][v]; corners[i] = sorted corners (t, v) of vertex i.
        Every vertex link must be a 2-sphere: it is a closed surface (closed
        triangulation, valid edges), connected, so it suffices that its Euler
        characteristic (link vertices = ends of edges at the vertex) is 2."""
        self.vertex_of = [[None] * 4 for _ in range(self.n)]
        self.corners = []
        for t in range(self.n):
            for v in range(4):
                if self.vertex_of[t][v] is not None:
                    continue
                i = len(self.corners)
                self.vertex_of[t][v] = i
                stack, members = [(t, v)], []
                while stack:
                    u, x = stack.pop()
                    members.append((u, x))
                    for f in range(4):
                        if f != x:
                            w, y = self.adj[u][f], self.glue[u][f][x]
                            if self.vertex_of[w][y] is None:
                                self.vertex_of[w][y] = i
                                stack.append((w, y))
                self.corners.append(sorted(members))
        self.n_vertices = len(self.corners)
        for i, cs in enumerate(self.corners):
            n_link_vertices = len({self.end[(t, v, w)] for t, v in cs for w in range(4) if w != v})
            if n_link_vertices - 3 * len(cs) // 2 + len(cs) != 2:
                raise Invalid('the link of vertex %d is not a 2-sphere' % i)


def corner_perm(v):
    """Link-triangle vertex j -> tetrahedron vertex, for the corner at v:
    (others in increasing order..., v)."""
    return tuple(x for x in range(4) if x != v) + (v,)
