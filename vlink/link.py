"""Exact geometry of a vertex link, developed in one fixed multiquadratic field.

Link triangle t of the link of vertex v sits in the corner of tetrahedron
tet(t) at v; its vertices j = 0, 1, 2 are the directions of the tetrahedron's
edges from v (Regina: vertex j lies on the edge from facePerm(t)[3] to
facePerm(t)[j]). The inner products of those three edge vectors form a
*rational* Gram matrix G_t, since squared edge lengths are rational:

    G_t[j][k] = (d(3,j) + d(3,k) - d(j,k)) / 2   (local tetrahedron indices).

A strip of link triangles is developed rigidly into the frame of its first
triangle: points are stored as coordinate vectors with respect to the first
triangle's edge vectors (so inner products are x^T G_0 y). Crossing from
triangle t to its neighbour t' over the edge (a, b) only needs

  - rational projection coefficients of the third vertex onto span(a, b),
    computed in G_t and in G_t' separately, and
  - mu = h'/h, the ratio of the third vertex's distances from that plane,
    which is the square root of the rational number h'^2 / h^2,

so every coordinate lies in Q(sqrt(q) : q among those ratios); for the cube
and disphenoid shapes mu = 1 and everything stays rational. Points are
directions: vectors are never normalised (no square roots of field elements).
"""

from fractions import Fraction

import regina

from .field import MQ, ZERO

EDGE_OF = {}


def _tet_sqlens(tet, sqlens):
    d = {}
    for k in range(6):
        m = tet.edgeMapping(k)
        i, j = m[0], m[1]
        d[(i, j)] = d[(j, i)] = Fraction(sqlens[tet.edge(k).index()])
    return d


class VertexLink:
    """Combinatorics and rational Gram data of the link of one vertex."""

    def __init__(self, T, i_v, sqlens):
        v = T.vertex(i_v)
        L = v.buildLink()
        I = v.buildLinkInclusion()
        self.n = n = L.size()
        self.T, self.i_v = T, i_v
        self.gram, self.tet, self.fp = [], [], []
        self.adj, self.glue, self.lv = [], [], []
        for t in range(n):
            tet = T.tetrahedron(I.tetImage(t))
            fp = I.facePerm(t)
            d = _tet_sqlens(tet, sqlens)
            c = fp[3]
            G = [[(d[(c, fp[j])] + d[(c, fp[k])] - (d[(fp[j], fp[k])] if j != k else 0)) / 2
                  for k in range(3)] for j in range(3)]
            self.gram.append(G)
            self.tet.append(tet.index())
            self.fp.append([fp[i] for i in range(4)])
            tri = L.triangle(t)
            self.adj.append([tri.adjacentTriangle(e).index() for e in range(3)])
            self.glue.append([[tri.adjacentGluing(e)[i] for i in range(3)] for e in range(3)])
            self.lv.append([tri.vertex(j).index() for j in range(3)])
        self.n_vertices = L.countVertices()
        self._step_cache = {}

    # ---------------------------------------------------------------- strips
    def step_data(self, t, d):
        """Rational data for crossing from triangle t over the edge opposite d.

        Returns (t2, g, alpha0, beta0, alpha, beta, mu): the neighbour, the
        vertex gluing g (t vertex -> t2 vertex), the projection coefficients of
        vertex d of t onto its edge (a, b) = the other two vertices, those of
        vertex g[d] of t2 onto (g[a], g[b]), and mu = h2/h (an MQ).
        """
        key = (t, d)
        r = self._step_cache.get(key)
        if r is not None:
            return r
        a, b = [x for x in range(3) if x != d]
        t2, g = self.adj[t][d], self.glue[t][d]
        G, G2 = self.gram[t], self.gram[t2]
        a2, b2, d2 = g[a], g[b], g[d]
        # the shared edge must have the same Gram data in both triangles
        assert (G[a][a], G[a][b], G[b][b]) == (G2[a2][a2], G2[a2][b2], G2[b2][b2]), \
            'inconsistent edge lengths across link edge'
        det = G[a][a] * G[b][b] - G[a][b] ** 2
        al0 = (G[d][a] * G[b][b] - G[d][b] * G[a][b]) / det
        be0 = (G[d][b] * G[a][a] - G[d][a] * G[a][b]) / det
        h2 = G[d][d] - al0 * G[d][a] - be0 * G[d][b]
        al = (G2[d2][a2] * G[b][b] - G2[d2][b2] * G[a][b]) / det
        be = (G2[d2][b2] * G[a][a] - G2[d2][a2] * G[a][b]) / det
        h22 = G2[d2][d2] - al * G2[d2][a2] - be * G2[d2][b2]
        if h2 <= 0 or h22 <= 0:
            raise ValueError('degenerate tetrahedron in link triangle %d or %d' % (t, t2))
        r = (t2, g, al0, be0, al, be, MQ.sqrt_rational(h22 / h2))
        self._step_cache[key] = r
        return r

    def start(self, t):
        """Developed positions of triangle t in its own frame: unit coordinate
        vectors (the edge vectors themselves)."""
        one, z = MQ.rational(1), ZERO
        return (Strip(self, t, [(one, z, z), (z, one, z), (z, z, one)]))


def _lin(c1, x, c2, y):
    return tuple(c1 * xi + c2 * yi for xi, yi in zip(x, y))


class Strip:
    """A developed sequence of link triangles (all in the first one's frame)."""

    __slots__ = ('link', 'ts', 'pos', 'frame', 'steps')

    def __init__(self, link, t, pos, ts=None, frame=None):
        self.link = link
        self.ts = ts or (t,)
        self.pos = (tuple(pos),)
        self.frame = frame if frame is not None else t
        self.steps = ()  # vertex d crossed opposite at each extension (for certificates)

    def extend(self, d):
        """New strip with the neighbour across the edge opposite vertex d of
        the last triangle appended."""
        link = self.link
        t = self.ts[-1]
        t2, g, al0, be0, al, be, mu = link.step_data(t, d)
        a, b = [x for x in range(3) if x != d]
        P = self.pos[-1]
        N = tuple(pd - pa * al0 - pb * be0 for pd, pa, pb in zip(P[d], P[a], P[b]))
        new = tuple(pa * al + pb * be - mu * n for pa, pb, n in zip(P[a], P[b], N))
        P2 = [None, None, None]
        P2[g[a]], P2[g[b]], P2[g[d]] = P[a], P[b], new
        s = Strip.__new__(Strip)
        s.link, s.ts, s.pos, s.frame = link, self.ts + (t2,), self.pos + (tuple(P2),), self.frame
        s.steps = self.steps + (d,)
        return s

    # ----------------------------------------------------------- predicates
    def dot(self, x, y):
        G = self.link.gram[self.frame]
        return sum((x[j] * y[k] * G[j][k] for j in range(3) for k in range(3)
                    if G[j][k]), ZERO)

    @staticmethod
    def det(x, y, z):
        """Triple product in frame coordinates (times det of the frame basis)."""
        return (x[0] * (y[1] * z[2] - y[2] * z[1])
                - x[1] * (y[0] * z[2] - y[2] * z[0])
                + x[2] * (y[0] * z[1] - y[1] * z[0]))

    def cos_float(self, x, y):
        return float(self.dot(x, y)) / (float(self.dot(x, x)) * float(self.dot(y, y))) ** 0.5


def link_of(isosig, i_v, sqlens):
    return VertexLink(regina.Triangulation3.fromIsoSig(isosig), i_v, sqlens)
