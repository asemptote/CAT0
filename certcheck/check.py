"""Exact checker for 'vlink-certificate/1' (certificates of non-CAT(0)-ness).

Everything is decided exactly; there is no floating point. The single
numerical device is AngleSum: an angle theta in [0, pi] is given as
theta = atan2(b, a) with a and b^2 exact (b >= 0), and is added by
multiplying z by a + i*b while counting half-turns exactly (see AngleSum.add).
Comparisons with pi and 2*pi are then exact.

Certificate kinds:
  invalid  a tetrahedron (or one of its faces) has non-positive Gram
           determinant, i.e. is not a non-degenerate Euclidean simplex;
  edge     the dihedral angles around an edge sum to < 2*pi;
  vertex   a necklace in the link of a vertex: a cyclic list of beads, each a
           strip of link triangles (with its gluing steps) carrying a geodesic
           from a link vertex S to a link vertex E, such that
             - the geodesic crosses every interior edge of the strip in its
               relative interior (exact cone tests with 3x3 determinants);
             - at each junction the angle is >= pi on both sides (sums of
               intrinsic corner angles of the link triangles around the vertex);
             - the total length is < 2*pi.
           A pi-bead (E antipodal to S) carries pi_point Q: its geodesic is
           the half great circle from S through the ray toward Q.
The geometry: link triangle t is the corner of a tetrahedron at the vertex;
its three vertices are that tetrahedron's edge vectors from the vertex, with
a rational Gram matrix. A strip is developed in the basis of its first
triangle's edge vectors; crossing to a neighbour uses rational projection
coefficients and the square root of a rational height ratio.
"""

from fractions import Fraction

import regina

from vlink.field import MQ

ONE, ZERO = MQ.rational(1), MQ()
FORMAT = 'vlink-certificate/1'


class Reject(Exception):
    pass


# ------------------------------------------------------------- angle sums

class AngleSum:
    """Exact sum of angles in [0, pi].

    z is (re, im), a positive multiple of exp(i * total); h counts half turns:
    total lies in [h*pi, (h+1)*pi). Adding theta in [0, pi] moves total into
    [h*pi, (h+2)*pi), and which half it is in follows from the sign of the
    new Im z (times (-1)^h), or of Re z when Im z = 0."""

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

    def at_least(self, k):
        """total >= k*pi ?"""
        return self.h >= k


# ------------------------------------------------------------- geometry

def _gram_det(rows):
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


class Geometry:
    def __init__(self, isosig, sqlens):
        self.T = T = regina.Triangulation3.fromIsoSig(isosig)
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
        if _gram_det(self.gram(tet, 0, (1, 2, 3))) <= 0:
            return False
        return all(_gram_det(self.gram(tet, f[0], f[1:])) > 0
                   for f in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)))

    def dihedral(self, tet, i, j):
        """(a, b^2) of the dihedral angle at edge (i, j) of tetrahedron tet."""
        k, l = [x for x in range(4) if x not in (i, j)]
        G = self.gram(tet, i, (j, k, l))
        ee, eu, ew = G[0][0], G[0][1], G[0][2]
        a = G[1][2] - eu * ew / ee
        nu, nw = G[1][1] - eu * eu / ee, G[2][2] - ew * ew / ee
        return MQ.rational(a), MQ.rational(nu * nw - a * a)


class Link:
    """Vertex link: per link triangle, its rational Gram matrix, adjacency,
    gluing permutations and link-vertex labels."""

    def __init__(self, geo, i_v):
        v = geo.T.vertex(i_v)
        L, I = v.buildLink(), v.buildLinkInclusion()
        self.n = L.size()
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


def _dot(G, x, y):
    return sum((x[j] * y[k] * G[j][k] for j in range(3) for k in range(3)), ZERO)


def _det(x, y, z):
    return (x[0] * (y[1] * z[2] - y[2] * z[1]) - x[1] * (y[0] * z[2] - y[2] * z[0])
            + x[2] * (y[0] * z[1] - y[1] * z[0]))


def _parallel(G, x, y):
    """+1 same direction, -1 opposite, 0 not parallel."""
    c = (x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0])
    if any(not v.is_zero() for v in c):
        return 0
    return _dot(G, x, y).sign()


def _develop(link, ts, steps):
    """Developed vertex positions of the strip ts (frame: ts[0]'s edge vectors)."""
    if len(steps) != len(ts) - 1:
        raise Reject('wrong number of strip steps')
    P = [(ONE, ZERO, ZERO), (ZERO, ONE, ZERO), (ZERO, ZERO, ONE)]
    out = [P]
    for k, d in enumerate(steps):
        t, t2 = ts[k], ts[k + 1]
        if d not in (0, 1, 2) or link.adj[t][d] != t2:
            raise Reject('strip step %d does not glue triangle %d to %d' % (k, t, t2))
        g = link.glue[t][d]
        a, b = [x for x in range(3) if x != d]
        Gt, G2 = link.G[t], link.G[t2]
        a2, b2, d2 = g[a], g[b], g[d]
        if (Gt[a][a], Gt[a][b], Gt[b][b]) != (G2[a2][a2], G2[a2][b2], G2[b2][b2]):
            raise Reject('glued link edge has different lengths on its two sides')
        det = Gt[a][a] * Gt[b][b] - Gt[a][b] ** 2
        # projections of the third vertex onto the shared edge, in t and in t2
        p = ((Gt[d][a] * Gt[b][b] - Gt[d][b] * Gt[a][b]) / det,
             (Gt[d][b] * Gt[a][a] - Gt[d][a] * Gt[a][b]) / det)
        q = ((G2[d2][a2] * Gt[b][b] - G2[d2][b2] * Gt[a][b]) / det,
             (G2[d2][b2] * Gt[a][a] - G2[d2][a2] * Gt[a][b]) / det)
        h = Gt[d][d] - p[0] * Gt[d][a] - p[1] * Gt[d][b]
        h2 = G2[d2][d2] - q[0] * G2[d2][a2] - q[1] * G2[d2][b2]
        if h <= 0 or h2 <= 0:
            raise Reject('degenerate link triangle')
        mu = MQ.sqrt_rational(h2 / h)
        Pa, Pb, Pd = P[a], P[b], P[d]
        normal = tuple(Pd[i] - Pa[i] * p[0] - Pb[i] * p[1] for i in range(3))
        new = tuple(Pa[i] * q[0] + Pb[i] * q[1] - mu * normal[i] for i in range(3))
        P2 = [None] * 3
        P2[a2], P2[b2], P2[d2] = Pa, Pb, new
        P = P2
        out.append(P)
    return out


def _cone_signs(X, A, B, C):
    """Signs of (alpha, beta) in X = alpha A + beta B (X in the plane of A, B;
    C off that plane)."""
    den = _det(A, B, C).sign()
    return _det(X, B, C).sign() * den, _det(A, X, C).sign() * den


def _crosses(G, S, Q, A, B, half_circle):
    """Does the geodesic from S toward Q cross the open edge arc AB?
    half_circle False: the minor arc S -> Q (Q a point at distance < pi);
    half_circle True : the half great circle from S through Q's ray to -S."""
    dS = _det(A, B, S)
    if dS.is_zero():
        return False          # S on the great circle of AB
    # intersection line of the planes (S, Q) and (A, B), as a vector
    X = tuple(S[i] * _det(A, B, Q) - Q[i] * dS for i in range(3))
    if all(x.is_zero() for x in X):
        return False
    al, be = _cone_signs(X, A, B, S)
    if al == be and al != 0:
        Y = X if al > 0 else tuple(-x for x in X)   # the point on the open arc AB
    else:
        return False
    if half_circle:
        # Y must lie on Q's side of S: Y . (Q - (Q.S / S.S) S) > 0
        qs, ss = _dot(G, Q, S), _dot(G, S, S)
        return (_dot(G, Y, Q) * ss - _dot(G, Y, S) * qs).sign() > 0
    ga, de = _cone_signs(Y, S, Q, A if not _det(S, Q, A).is_zero() else B)
    return ga > 0 and de > 0


def _tangent_angle(G, W, X, Y):
    """(a, b^2) of the angle at W between the geodesics toward X and Y."""
    ww = _dot(G, W, W)
    xw, yw = _dot(G, X, W), _dot(G, Y, W)
    a = _dot(G, X, Y) * ww - xw * yw
    nx = _dot(G, X, X) * ww - xw * xw
    ny = _dot(G, Y, Y) * ww - yw * yw
    return a, nx * ny - a * a


def _length(G, S, E):
    a = _dot(G, S, E)
    return a, _dot(G, S, S) * _dot(G, E, E) - a * a


def _vec(raw):
    return tuple(MQ({int(s): Fraction(c) for s, c in x.items()}) for x in raw)


# ------------------------------------------------------------- beads

class CheckedBead:
    pass


def _check_bead(link, b):
    ts, s, e = list(b['ts']), b['start'], b['end']
    if not ts or not all(0 <= t < link.n for t in ts) or s not in (0, 1, 2) or e not in (0, 1, 2):
        raise Reject('malformed bead')
    pos = _develop(link, ts, list(b.get('steps') or []))
    G = link.G[ts[0]]
    S, E = pos[0][s], pos[-1][e]
    cb = CheckedBead()
    cb.ts, cb.start, cb.end, cb.G, cb.pos, cb.S, cb.E = ts, s, e, G, pos, S, E
    par = _parallel(G, S, E)
    if par > 0:
        raise Reject('bead endpoints coincide')
    cb.pi = par < 0
    if cb.pi:
        if not b.get('pi_point'):
            raise Reject('pi-bead without pi_point')
        Q = _vec(b['pi_point'])
        if _parallel(G, S, Q):
            raise Reject('pi_point on the axis')
        cb.Q = Q
    elif len(ts) == 1:
        if s == e:
            raise Reject('degenerate one-triangle bead')
    for k in range(len(ts) - 1):
        d = b['steps'][k]
        A, B = [pos[k][x] for x in range(3) if x != d]
        if cb.pi:
            ok = _crosses(G, S, cb.Q, A, B, half_circle=True)
        else:
            ok = _crosses(G, S, E, A, B, half_circle=False)
        if not ok:
            raise Reject('geodesic of bead %s does not cross interior edge %d' % (ts, k))
    # direction of the geodesic leaving S, and arriving at E (reversed)
    cb.out_dir = cb.Q if cb.pi else E
    cb.in_dir = cb.Q if cb.pi else S
    return cb


def _side(link, bin_, bout, k):
    """Angle at the junction bin_ -> bout on side k (1 or 2) as an AngleSum."""
    tm, rm = bin_.ts[-1], bin_.end
    t0, d0 = bout.ts[0], bout.start
    P = bin_.pos[-1]
    G_in, G_out = bin_.G, bout.G
    W, u = P[rm], bin_.in_dir
    d = (rm + k) % 3
    y = 3 - rm - d
    acc = AngleSum()
    if (tm, rm) == (t0, d0):
        # outgoing direction inside this corner between u and ray(y)?
        X = tuple(sum((bout.out_dir[j] * P[j][i] for j in range(3)), ZERO) for i in range(3))
        a1, b1 = _tangent_angle(G_in, W, u, X)
        a2, b2 = _tangent_angle(G_in, W, X, P[y])
        a3, b3 = _tangent_angle(G_in, W, u, P[y])
        s12, s3 = AngleSum(), AngleSum()
        s12.add(a1, b1)
        s12.add(a2, b2)
        s3.add(a3, b3)
        same = (s12.h == s3.h and (s12.im * s3.re - s12.re * s3.im).is_zero()
                and (s12.re * s3.re + s12.im * s3.im).sign() > 0)
        if same:
            acc.add(a1, b1)
            return acc
    acc.add(*_tangent_angle(G_in, W, u, P[y]))
    t, c, ent, ex = tm, rm, d, y
    for _ in range(4 * link.n + 2):
        if acc.at_least(1):
            return acc
        g = link.glue[t][ent]
        t, c, ent = link.adj[t][ent], g[c], g[ex]
        ex = 3 - c - ent
        if (t, c) == (t0, d0):
            Wt = tuple(ONE if i == d0 else ZERO for i in range(3))
            Rt = tuple(ONE if i == ent else ZERO for i in range(3))
            acc.add(*_tangent_angle(G_out, Wt, Rt, bout.out_dir))
            return acc
        Gt = link.G[t]
        e3 = [tuple(ONE if i == j else ZERO for i in range(3)) for j in range(3)]
        acc.add(*_tangent_angle(Gt, e3[c], e3[ent], e3[ex]))
    raise Reject('junction fan walk did not reach the outgoing triangle')


# ------------------------------------------------------------- entry point

def check_certificate(cert):
    """Returns (True, message) if the certificate proves the metric is not
    locally CAT(0), else (False, reason). Never raises on bad certificates."""
    try:
        return True, _check(cert)
    except Reject as e:
        return False, str(e)
    except (KeyError, TypeError, ValueError, IndexError) as e:
        return False, 'malformed certificate: %s: %s' % (type(e).__name__, e)


def _check(cert):
    if cert.get('format') != FORMAT:
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
        e = geo.T.edge(cert['edge'])
        acc = AngleSum()
        for emb in e.embeddings():
            m = emb.vertices()
            acc.add(*geo.dihedral(emb.tetrahedron().index(), m[0], m[1]))
        if acc.at_least(2):
            raise Reject('edge %d has dihedral sum >= 2*pi' % cert['edge'])
        return 'edge %d has dihedral sum < 2*pi' % cert['edge']
    if reason != 'vertex':
        raise Reject('unknown reason')
    i_v = cert['vertex']
    if not 0 <= i_v < geo.T.countVertices():
        raise Reject('bad vertex')
    link = Link(geo, i_v)
    beads = [_check_bead(link, b) for b in cert['necklace']]
    if not beads:
        raise Reject('empty necklace')
    m = len(beads)
    total = AngleSum()
    for b in beads:
        total.add(*((MQ.rational(-1), ZERO) if b.pi else _length(b.G, b.S, b.E)))
    if total.at_least(2):
        raise Reject('necklace length >= 2*pi')
    for k in range(m):
        b1, b2 = beads[k], beads[(k + 1) % m]
        if link.lv[b1.ts[-1]][b1.end] != link.lv[b2.ts[0]][b2.start]:
            raise Reject('beads %d and %d do not meet' % (k, (k + 1) % m))
        for side in (1, 2):
            if not _side(link, b1, b2, side).at_least(1):
                raise Reject('junction %d: angle < pi on side %d' % (k, side))
    return 'vertex %d: necklace of %d beads, length < 2*pi, all junctions >= pi' % (i_v, m)
