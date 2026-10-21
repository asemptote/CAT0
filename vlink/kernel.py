"""The trusted kernel of vlink: exact geometric predicates and certificate checks.

Everything vlink concludes about a metric rests on the functions in this
module (and on field.py and Regina's combinatorics); the search modules
(beads.py, necklace.py) only decide what to ask the kernel, and every answer
they produce is checked here before it is returned (check.py). This is the
part to review closely, and the natural target for a formal proof.

All decisions are exact; there is no floating point. The one numerical
device is AngleSum: an angle theta in [0, pi] is given as theta = atan2(b, a)
with a and b^2 exact (b >= 0), and angles are added by multiplying a + i*b
while counting half-turns exactly, so comparisons with pi and 2*pi are exact.

Geometry: link triangle t of the link of a vertex is the corner of a
tetrahedron at that vertex; its vertices 0, 1, 2 are the tetrahedron's edge
vectors from the vertex (Regina: vertex j lies on the edge from facePerm[3]
to facePerm[j]), with a rational Gram matrix. A strip of link triangles is
developed in the basis of its first triangle's edge vectors; crossing to a
neighbour uses rational projection coefficients and the square root of a
rational ratio of heights. Points are unnormalised vectors.

Certificates ('vlink-certificate/1') prove a metric is not locally CAT(0):
  invalid  a tetrahedron (or face) with non-positive Gram determinant;
  edge     an edge whose dihedral angles sum to < 2*pi;
  vertex   a necklace of beads in a vertex link: geodesics crossing their
           strips' interior edges, junction angles >= pi on both sides, total
           length < 2*pi (a pi-bead carries pi_point fixing its geodesic).
"""

from fractions import Fraction

import regina

from .field import MQ

ONE, ZERO = MQ.rational(1), MQ()
CERT_FORMAT = 'vlink-certificate/1'
E3 = ((ONE, ZERO, ZERO), (ZERO, ONE, ZERO), (ZERO, ZERO, ONE))


class Reject(Exception):
    """A claim checked by the kernel is false (or cannot be decided exactly)."""


# ----------------------------------------------------------------- angle sums

class AngleSum:
    """Exact sum of angles in [0, pi].

    (re, im) is a positive multiple of exp(i * total), and total lies in
    [h*pi, (h+1)*pi). Adding theta in [0, pi] moves total into
    [h*pi, (h+2)*pi); which half follows from the sign of (-1)^h * Im z, or
    of Re z when Im z = 0."""

    __slots__ = ('re', 'im', 'h')

    def __init__(self):
        self.re, self.im, self.h = ONE, ZERO, 0

    def add(self, a, b2):
        if any(s != 1 for s in b2.t):
            raise Reject('unsupported: sin^2 not rational (irrational link field)')
        q = b2.t.get(1, Fraction(0))
        if q < 0:
            raise Reject('negative sin^2')
        b = MQ.sqrt_rational(q)
        if a.is_zero() and b.is_zero():
            raise Reject('angle between degenerate vectors')
        self.re, self.im = self.re * a - self.im * b, self.re * b + self.im * a
        s = 1 if self.h % 2 == 0 else -1
        si = s * self.im.sign()
        if si < 0 or (si == 0 and s * self.re.sign() < 0):
            self.h += 1
        return self

    def at_least(self, k):
        """total >= k * pi ?"""
        return self.h >= k


# ----------------------------------------------------------------- linear algebra

def gram_det(rows):
    M = [list(r) for r in rows]
    n, det = len(M), Fraction(1)
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


def dot(G, x, y):
    return sum((x[j] * y[k] * G[j][k] for j in range(3) for k in range(3) if G[j][k]), ZERO)


def det(x, y, z):
    return (x[0] * (y[1] * z[2] - y[2] * z[1]) - x[1] * (y[0] * z[2] - y[2] * z[0])
            + x[2] * (y[0] * z[1] - y[1] * z[0]))


def parallel(G, x, y):
    """+1 same direction, -1 opposite, 0 not parallel."""
    c = (x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0])
    if any(not v.is_zero() for v in c):
        return 0
    return dot(G, x, y).sign()


# ----------------------------------------------------------------- metric

class Geometry:
    """Squared edge lengths of each tetrahedron (local vertex labels)."""

    def __init__(self, T, sqlens):
        self.T = T = regina.Triangulation3.fromIsoSig(T) if isinstance(T, str) else T
        if len(sqlens) != T.countEdges():
            raise Reject('wrong number of edge lengths')
        self.d = []
        for tet in T.tetrahedra():
            d = {}
            for k in range(6):
                m = tet.edgeMapping(k)
                d[(m[0], m[1])] = d[(m[1], m[0])] = Fraction(sqlens[tet.edge(k).index()])
            self.d.append(d)

    def gram(self, tet, base, others):
        d = self.d[tet]
        return [[(d[(base, x)] + d[(base, y)] - (d[(x, y)] if x != y else 0)) / 2
                 for y in others] for x in others]

    def tet_valid(self, tet):
        """Non-degenerate Euclidean tetrahedron: positive Gram determinants of
        the tetrahedron and of its four faces."""
        if gram_det(self.gram(tet, 0, (1, 2, 3))) <= 0:
            return False
        return all(gram_det(self.gram(tet, f[0], f[1:])) > 0
                   for f in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)))

    def dihedral(self, tet, i, j):
        """(a, b^2) of the dihedral angle at edge (i, j) of tetrahedron tet."""
        k, l = [x for x in range(4) if x not in (i, j)]
        G = self.gram(tet, i, (j, k, l))
        ee, eu, ew = G[0][0], G[0][1], G[0][2]
        a = G[1][2] - eu * ew / ee
        nu, nw = G[1][1] - eu * eu / ee, G[2][2] - ew * ew / ee
        return MQ.rational(a), MQ.rational(nu * nw - a * a)

    def edge_angle(self, e):
        """AngleSum of the dihedral angles around edge e."""
        acc = AngleSum()
        for emb in self.T.edge(e).embeddings():
            m = emb.vertices()
            acc.add(*self.dihedral(emb.tetrahedron().index(), m[0], m[1]))
        return acc


class Link:
    """Link of vertex i_v: per link triangle, its rational Gram matrix,
    neighbours, gluing permutations and link-vertex labels."""

    def __init__(self, geo, i_v):
        v = geo.T.vertex(i_v)
        L, I = v.buildLink(), v.buildLinkInclusion()
        self.i_v, self.n, self.n_vertices = i_v, L.size(), L.countVertices()
        self.G, self.adj, self.glue, self.lv = [], [], [], []
        for t in range(self.n):
            tet, fp = I.tetImage(t), I.facePerm(t)
            if geo.T.tetrahedron(tet).vertex(fp[3]).index() != i_v:
                raise Reject('link inclusion does not map to the vertex')
            self.G.append(geo.gram(tet, fp[3], (fp[0], fp[1], fp[2])))
            tri = L.triangle(t)
            self.adj.append([tri.adjacentTriangle(e).index() for e in range(3)])
            self.glue.append([[tri.adjacentGluing(e)[i] for i in range(3)] for e in range(3)])
            self.lv.append([tri.vertex(j).index() for j in range(3)])
        self._trans = {}

    def transition(self, t, d):
        """Crossing from triangle t over the edge opposite its vertex d:
        (t2, g, p, q, mu) with g the vertex gluing t -> t2, p / q the
        projection coefficients of t's vertex d / t2's vertex g[d] onto the
        shared edge, and mu = h2/h the ratio of their heights (an MQ)."""
        key = (t, d)
        r = self._trans.get(key)
        if r is not None:
            return r
        t2, g = self.adj[t][d], self.glue[t][d]
        a, b = [x for x in range(3) if x != d]
        Gt, G2 = self.G[t], self.G[t2]
        a2, b2, d2 = g[a], g[b], g[d]
        if (Gt[a][a], Gt[a][b], Gt[b][b]) != (G2[a2][a2], G2[a2][b2], G2[b2][b2]):
            raise Reject('glued link edge has different lengths on its two sides')
        dt = Gt[a][a] * Gt[b][b] - Gt[a][b] ** 2
        p = ((Gt[d][a] * Gt[b][b] - Gt[d][b] * Gt[a][b]) / dt,
             (Gt[d][b] * Gt[a][a] - Gt[d][a] * Gt[a][b]) / dt)
        q = ((G2[d2][a2] * Gt[b][b] - G2[d2][b2] * Gt[a][b]) / dt,
             (G2[d2][b2] * Gt[a][a] - G2[d2][a2] * Gt[a][b]) / dt)
        h = Gt[d][d] - p[0] * Gt[d][a] - p[1] * Gt[d][b]
        h2 = G2[d2][d2] - q[0] * G2[d2][a2] - q[1] * G2[d2][b2]
        if h <= 0 or h2 <= 0:
            raise Reject('degenerate link triangle')
        r = self._trans[key] = (t2, g, p, q, MQ.sqrt_rational(h2 / h))
        return r


def step(link, P, t, d):
    """Develop across the edge opposite vertex d of triangle t, whose vertices
    are at P (developed frame). Returns (t2, P2)."""
    t2, g, p, q, mu = link.transition(t, d)
    a, b = [x for x in range(3) if x != d]
    Pa, Pb, Pd = P[a], P[b], P[d]
    normal = tuple(Pd[i] - Pa[i] * p[0] - Pb[i] * p[1] for i in range(3))
    new = tuple(Pa[i] * q[0] + Pb[i] * q[1] - mu * normal[i] for i in range(3))
    P2 = [None] * 3
    P2[g[a]], P2[g[b]], P2[g[d]] = Pa, Pb, new
    return t2, tuple(P2)


def develop(link, ts, steps):
    """Developed vertex positions of the strip ts glued by `steps`."""
    if len(steps) != len(ts) - 1:
        raise Reject('wrong number of strip steps')
    P, out = E3, [E3]
    for k, d in enumerate(steps):
        t = ts[k]
        if d not in (0, 1, 2) or link.adj[t][d] != ts[k + 1]:
            raise Reject('strip step %d does not glue triangle %d to %d' % (k, t, ts[k + 1]))
        _, P = step(link, P, t, d)
        out.append(P)
    return out


# ----------------------------------------------------------------- geodesics and angles

def _cone_signs(X, A, B, C):
    """Signs of (alpha, beta) in X = alpha A + beta B (X in the plane of A, B;
    C off that plane)."""
    den = det(A, B, C).sign()
    return det(X, B, C).sign() * den, det(A, X, C).sign() * den


def crosses(G, S, Q, A, B, half_circle):
    """Does the geodesic from S toward Q cross the open edge arc AB?
    half_circle False: the minor arc S -> Q (Q at distance < pi from S);
    half_circle True : the half great circle from S through Q's ray to -S."""
    dS = det(A, B, S)
    if dS.is_zero():
        return False
    X = tuple(S[i] * det(A, B, Q) - Q[i] * dS for i in range(3))
    if all(x.is_zero() for x in X):
        return False
    al, be = _cone_signs(X, A, B, S)
    if not (al == be and al != 0):
        return False
    Y = X if al > 0 else tuple(-x for x in X)   # the point on the open arc AB
    if half_circle:
        qs, ss = dot(G, Q, S), dot(G, S, S)
        return (dot(G, Y, Q) * ss - dot(G, Y, S) * qs).sign() > 0
    ga, de = _cone_signs(Y, S, Q, A if not det(S, Q, A).is_zero() else B)
    return ga > 0 and de > 0


def tangent_angle(G, W, X, Y):
    """(a, b^2) of the angle at W between the geodesics toward X and Y."""
    ww = dot(G, W, W)
    xw, yw = dot(G, X, W), dot(G, Y, W)
    a = dot(G, X, Y) * ww - xw * yw
    nx = dot(G, X, X) * ww - xw * xw
    ny = dot(G, Y, Y) * ww - yw * yw
    return a, nx * ny - a * a


def length(G, S, E):
    """(a, b^2) of the angle between S and E (a geodesic's length)."""
    a = dot(G, S, E)
    return a, dot(G, S, S) * dot(G, E, E) - a * a


# ----------------------------------------------------------------- beads and junctions

class Bead:
    """A bead as the kernel sees it: strip ts with developed positions pos
    (frame: ts[0]'s basis, Gram G), geodesic from S = pos[0][start] to
    E = pos[-1][end]; for a pi-bead (E antipodal to S) Q fixes the geodesic.
    in_dir / out_dir: points giving the geodesic's direction at E (reversed)
    and at S."""

    __slots__ = ('ts', 'start', 'end', 'G', 'pos', 'S', 'E', 'pi', 'Q', 'in_dir', 'out_dir')

    def __init__(self, ts, start, end, G, pos, Q=None):
        self.ts, self.start, self.end, self.G, self.pos = ts, start, end, G, pos
        self.S, self.E = pos[0][start], pos[-1][end]
        self.pi = parallel(G, self.S, self.E) < 0
        self.Q = Q
        self.out_dir = Q if self.pi else self.E
        self.in_dir = Q if self.pi else self.S


def check_bead(link, spec):
    """Kernel Bead from a certificate entry {ts, start, end, steps, pi_point},
    after checking that it is a bead: its geodesic crosses every interior
    edge of the strip in its relative interior."""
    ts, s, e = list(spec['ts']), spec['start'], spec['end']
    if not ts or not all(0 <= t < link.n for t in ts) or s not in (0, 1, 2) or e not in (0, 1, 2):
        raise Reject('malformed bead')
    steps = list(spec.get('steps') or [])
    pos = develop(link, ts, steps)
    G = link.G[ts[0]]
    S, E = pos[0][s], pos[-1][e]
    if parallel(G, S, E) > 0:
        raise Reject('bead endpoints coincide')
    Q = None
    if parallel(G, S, E) < 0:
        if not spec.get('pi_point'):
            raise Reject('pi-bead without pi_point')
        Q = tuple(MQ({int(k): Fraction(c) for k, c in x.items()}) for x in spec['pi_point'])
        if parallel(G, S, Q):
            raise Reject('pi_point on the axis')
    elif len(ts) == 1 and s == e:
        raise Reject('degenerate one-triangle bead')
    b = Bead(ts, s, e, G, pos, Q)
    for k in range(len(ts) - 1):
        A, B = [pos[k][x] for x in range(3) if x != steps[k]]
        if not crosses(G, S, Q if b.pi else E, A, B, half_circle=b.pi):
            raise Reject('geodesic of bead %s does not cross interior edge %d' % (ts, k))
    return b


def junction_side(link, bin_, bout, k):
    """AngleSum of the angle at the junction bin_ -> bout (a link vertex W)
    on side k in {1, 2}: from the reversed incoming direction, around W
    through the link triangles, to the outgoing direction. Summed from the
    intrinsic corner angles of the triangles between them."""
    tm, rm = bin_.ts[-1], bin_.end
    t0, d0 = bout.ts[0], bout.start
    P, G_in = bin_.pos[-1], bin_.G
    W, u = P[rm], bin_.in_dir
    d = (rm + k) % 3
    y = 3 - rm - d
    if (tm, rm) == (t0, d0):
        # is the outgoing direction inside this corner, between u and ray(y)?
        X = tuple(sum((bout.out_dir[j] * P[j][i] for j in range(3)), ZERO) for i in range(3))
        a1, b1 = tangent_angle(G_in, W, u, X)
        s12 = AngleSum().add(a1, b1).add(*tangent_angle(G_in, W, X, P[y]))
        s3 = AngleSum().add(*tangent_angle(G_in, W, u, P[y]))
        if (s12.h == s3.h and (s12.im * s3.re - s12.re * s3.im).is_zero()
                and (s12.re * s3.re + s12.im * s3.im).sign() > 0):
            return AngleSum().add(a1, b1)
    acc = AngleSum().add(*tangent_angle(G_in, W, u, P[y]))
    t, c, ent, ex = tm, rm, d, y
    for _ in range(4 * link.n + 2):
        if acc.at_least(1):
            return acc
        g = link.glue[t][ent]
        t, c, ent = link.adj[t][ent], g[c], g[ex]
        ex = 3 - c - ent
        if (t, c) == (t0, d0):
            return acc.add(*tangent_angle(bout.G, E3[d0], E3[ent], bout.out_dir))
        acc.add(*tangent_angle(link.G[t], E3[c], E3[ent], E3[ex]))
    raise Reject('junction fan walk did not reach the outgoing triangle')


def junction_ok(link, bin_, bout):
    """Angle >= pi on both sides of the junction bin_ -> bout."""
    return all(junction_side(link, bin_, bout, k).at_least(1) for k in (1, 2))


def total_length(beads):
    """AngleSum of the lengths of kernel beads."""
    acc = AngleSum()
    for b in beads:
        acc.add(*((MQ.rational(-1), ZERO) if b.pi else length(b.G, b.S, b.E)))
    return acc


# ----------------------------------------------------------------- certificates

def check_certificate(cert):
    """(True, message) if the certificate proves the metric is not locally
    CAT(0), else (False, reason). Never raises on bad certificates."""
    try:
        return True, _check(cert)
    except Reject as e:
        return False, str(e)
    except (KeyError, TypeError, ValueError, IndexError, AttributeError) as e:
        return False, 'malformed certificate: %s: %s' % (type(e).__name__, e)


def _check(cert):
    if cert.get('format') != CERT_FORMAT:
        raise Reject('unknown certificate format')
    geo = Geometry(cert['isosig'], cert['sqlens'])
    n_tet = geo.T.size()
    reason = cert['reason']
    if reason == 'invalid':
        t = cert['tet']
        if not 0 <= t < n_tet or geo.tet_valid(t):
            raise Reject('tetrahedron %s is valid' % t)
        return 'tetrahedron %d is degenerate' % t
    for t in range(n_tet):
        if not geo.tet_valid(t):
            raise Reject('tetrahedron %d is degenerate (use reason "invalid")' % t)
    if reason == 'edge':
        if geo.edge_angle(cert['edge']).at_least(2):
            raise Reject('edge %d has dihedral sum >= 2*pi' % cert['edge'])
        return 'edge %d has dihedral sum < 2*pi' % cert['edge']
    if reason != 'vertex':
        raise Reject('unknown reason')
    i_v = cert['vertex']
    if not 0 <= i_v < geo.T.countVertices():
        raise Reject('bad vertex')
    link = Link(geo, i_v)
    beads = [check_bead(link, b) for b in cert['necklace']]
    if not beads:
        raise Reject('empty necklace')
    if total_length(beads).at_least(2):
        raise Reject('necklace length >= 2*pi')
    m = len(beads)
    for k in range(m):
        b1, b2 = beads[k], beads[(k + 1) % m]
        if link.lv[b1.ts[-1]][b1.end] != link.lv[b2.ts[0]][b2.start]:
            raise Reject('beads %d and %d do not meet' % (k, (k + 1) % m))
        for side in (1, 2):
            if not junction_side(link, b1, b2, side).at_least(1):
                raise Reject('junction %d: angle < pi on side %d' % (k, side))
    return 'vertex %d: necklace of %d beads, length < 2*pi, all junctions >= pi' % (i_v, m)
