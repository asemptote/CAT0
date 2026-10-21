"""The trusted kernel of vlink: exact geometric predicates and certificate checks.

Everything vlink concludes about a metric rests on the functions in this
module (and on field.py and tri.py, the kernel's own isosig decoder and
skeleton; Regina is not used); the search modules
(beads.py, necklace.py) only decide what to ask the kernel, and every answer
they produce is checked here before it is returned (check.py). This is the
part to review closely, and the natural target for a formal proof.

All decisions are exact; there is no floating point. The one numerical
device is AngleSum: an angle theta in [0, pi] is given as theta = atan2(b, a)
with a and b^2 exact (b >= 0), and angles are added by multiplying a + i*b
while counting half-turns exactly, so comparisons with pi and 2*pi are exact.

Combinatorics come from tri.py (the kernel's own isosig decoder; Regina is
not used). Geometry: link triangle t of the link of a vertex is a corner
(tet, v) of a tetrahedron at that vertex, the corners numbered in sorted
order; its vertices 0, 1, 2 are the tetrahedron's edge vectors from v to its
other vertices in increasing order (tri.corner_perm), with a rational Gram
matrix. A strip of link triangles is
developed in the basis of its first triangle's edge vectors; crossing to a
neighbour uses rational projection coefficients and the square root of a
rational ratio of heights. Points are unnormalised vectors.

Certificates ('vlink-certificate/2') prove a metric is not locally CAT(0):
  invalid  a tetrahedron (or face) with non-positive Gram determinant;
  edge     an edge whose dihedral angles sum to < 2*pi;
  vertex   a necklace of beads in a vertex link: geodesics crossing their
           strips' interior edges, junction angles >= pi on both sides, total
           length < 2*pi (a pi-bead carries pi_point fixing its geodesic).

Certificates ('vlink-certificate-true/2') prove a metric IS locally CAT(0):
the kernel checks validity and every edge itself, enumerates every bead
itself (enumerate_beads), and verifies shortest-path trees showing that no
cyclic chain of beads with passing junctions is shorter than 2*pi
(check_vertex_trees; docs/true_certificates.md). Junctions involving
beads of length exactly pi are counted as passing (conservative); with
pi_beads=False such beads are left out and the certificate must list the
assumption Elder-McCammond-Meier 2003, Lemma 4.9.

Length comparisons use rigorous float enclosures with exact fallback
(Lengths); exact_only=True makes every comparison exact.
"""

from fractions import Fraction

from .field import MQ, Tower
from .tri import EDGES, Triangulation, corner_perm

ONE, ZERO = MQ.rational(1), MQ()
CERT_FORMAT = 'vlink-certificate/2'     # /2: link triangles are sorted corners
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
        if b2.sign() < 0:
            raise Reject('negative sin^2')
        if all(s == 1 for s in b2.t):
            b = MQ.sqrt_rational(b2.t.get(1, Fraction(0)))     # stays in F
        else:
            b = Tower.sqrt_field(b2)                            # irrational link field
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
    """Squared edge lengths of each tetrahedron (local vertex labels) of the
    closed triangulated 3-manifold T (a tri.Triangulation or an isosig)."""

    def __init__(self, T, sqlens):
        self.T = T = Triangulation(T) if isinstance(T, str) else T
        if len(sqlens) != len(T.edges):
            raise Reject('wrong number of edge lengths')
        self.d = []
        for tet in range(T.n):
            d = {}
            for k, (a, b) in enumerate(EDGES):
                d[(a, b)] = d[(b, a)] = Fraction(sqlens[T.edge_of[tet][k]])
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
        for tet, k in self.T.edges[e]:
            acc.add(*self.dihedral(tet, *EDGES[k]))
        return acc


class Link:
    """Link of vertex i_v: per link triangle (corner), its rational Gram
    matrix, neighbours adj[t][d] (across the link edge opposite vertex d),
    gluing permutations glue[t][d] (vertex j of t -> vertex glue[t][d][j] of
    adj[t][d]) and link-vertex labels lv[t][j] (numbered by first appearance)."""

    def __init__(self, geo, i_v):
        T = geo.T
        corners = T.corners[i_v]
        index = {c: t for t, c in enumerate(corners)}
        self.i_v, self.n = i_v, len(corners)
        self.G, self.adj, self.glue, self.lv = [], [], [], []
        labels = {}
        for tet, v in corners:
            fp = corner_perm(v)
            self.G.append(geo.gram(tet, v, fp[:3]))
            adj, glue = [], []
            for d in range(3):
                g = T.glue[tet][fp[d]]
                v2 = g[v]
                fp2 = corner_perm(v2)
                adj.append(index[(T.adj[tet][fp[d]], v2)])
                glue.append([fp2.index(g[fp[j]]) for j in range(3)])
            self.adj.append(adj)
            self.glue.append(glue)
            self.lv.append([labels.setdefault(T.end[(tet, v, fp[j])], len(labels))
                            for j in range(3)])
        self.n_vertices = len(labels)
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
    (frame: ts[0]'s basis, Gram G), glued by steps, carrying the geodesic
    from S = pos[0][start] to E = pos[-1][end]. For a pi-bead (E antipodal to
    S) the geodesics form a lune: `lune` = (L, R), the open interval of rays
    from S (enumeration), or Q, one chosen geodesic (certificates).
    in_dir / out_dir: points giving the geodesic's direction at E (reversed)
    and at S."""

    __slots__ = ('ts', 'start', 'end', 'G', 'pos', 'steps', 'S', 'E', 'pi', 'Q', 'lune',
                 'in_dir', 'out_dir')

    def __init__(self, ts, start, end, G, pos, Q=None, steps=(), lune=None):
        self.ts, self.start, self.end, self.G, self.pos = ts, start, end, G, pos
        self.steps, self.lune = tuple(steps), lune
        self.S, self.E = pos[0][start], pos[-1][end]
        self.pi = parallel(G, self.S, self.E) < 0
        self.Q = Q
        self.out_dir = Q if self.pi else self.E
        self.in_dir = Q if self.pi else self.S


# ----------------------------------------------------------------- bead enumeration
#
# Completeness rests on one fact: seen from S, the relative interior of an
# edge AB (S not on its great circle) is hit exactly by the rays from S
# strictly inside the angle at S of the spherical triangle SAB, an open
# interval of width < pi; rays are compared by orientation, sign det(S, X, Y).
# A geodesic from S of length < pi (or = pi) passing through a strip crosses
# each interior edge in its interior, so its direction lies in all of these
# intervals; hence extending a strip only while their intersection is a
# non-empty open interval, and ending a bead at a vertex E of the last
# triangle exactly when E's ray lies strictly inside it, finds every bead.
# A strip is abandoned when a vertex is parallel to S (zero_zero / zero_pi):
# no geodesic from S of length < pi can cross an edge ending there.

def _visible(S, A, B):
    o = det(S, A, B).sign()
    if o == 0:
        return None
    return (A, B) if o > 0 else (B, A)


def _inside(S, L, R, X, strict=True):
    a, b = det(S, L, X).sign(), det(S, X, R).sign()
    return (a > 0 and b > 0) if strict else (a >= 0 and b >= 0)


def _intersect(S, I1, I2):
    (L1, R1), (L2, R2) = I1, I2
    if _inside(S, L1, R1, L2, strict=False):
        L = L2
    elif _inside(S, L2, R2, L1, strict=False):
        L = L1
    else:
        return None
    if _inside(S, L1, R1, R2, strict=False):
        R = R2
    elif _inside(S, L2, R2, R1, strict=False):
        R = R1
    else:
        return None
    return (L, R) if det(S, L, R).sign() > 0 else None


def enumerate_beads(link, allow_repeats=True, pi_beads=True, triangles=None):
    """All beads of the link (those starting in `triangles`, default all):
    thin beads (one triangle, along an edge) and thick beads (geodesic of
    length < pi through a strip), plus pi-beads (lunes) unless pi_beads is
    False. allow_repeats=False forbids a strip visiting a triangle twice."""
    out = []
    for t in (range(link.n) if triangles is None else triangles):
        G = link.G[t]
        for d0 in range(3):
            P = E3
            a, b = (d0 + 1) % 3, (d0 + 2) % 3
            S = P[d0]
            out.append(Bead((t,), d0, a, G, (P,)))
            out.append(Bead((t,), d0, b, G, (P,)))
            if parallel(G, S, P[a]) or parallel(G, S, P[b]):
                continue
            iv = _visible(S, P[a], P[b])
            if iv is not None:
                _grow(link, G, (t,), (P,), (), d0, S, d0, iv, out, allow_repeats, pi_beads)
    return out


def _grow(link, G, ts, pos, steps, start, S, d, interval, out, allow_repeats, pi_beads):
    """Cross the edge opposite vertex d of the strip's last triangle."""
    t = ts[-1]
    t2 = link.adj[t][d]
    if t2 in ts and not allow_repeats:
        return
    if len(ts) > 4 * link.n:
        raise Reject('bead strip longer than 4x the link: growth not terminating?')
    _, P2 = step(link, pos[-1], t, d)
    ts2, pos2, steps2 = ts + (t2,), pos + (P2,), steps + (d,)
    r = link.glue[t][d][d]            # the new vertex of the new triangle
    E = P2[r]
    par = parallel(G, S, E)
    if par < 0:
        if pi_beads:
            out.append(Bead(ts2, start, r, G, pos2, steps=steps2, lune=interval))
    elif par == 0 and _inside(S, interval[0], interval[1], E):
        out.append(Bead(ts2, start, r, G, pos2, steps=steps2))
    for dn in ((r + 1) % 3, (r + 2) % 3):
        y = 3 - r - dn                # the other endpoint of the edge opposite dn
        if par or parallel(G, S, P2[y]):
            continue
        iv = _visible(S, E, P2[y])
        if iv is None:
            continue
        new = _intersect(S, interval, iv)
        if new is not None:
            _grow(link, G, ts2, pos2, steps2, start, S, dn, new, out, allow_repeats, pi_beads)


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
    b = Bead(ts, s, e, G, pos, Q, steps=steps)
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


# ----------------------------------------------------------------- length comparisons
#
# Sums of bead lengths are compared with rigorous float enclosures first and
# exactly (AngleSum) only when the enclosures overlap; exact_only=True skips
# the float path. Each bead length theta = acos(c) has float value with error
# at most length_err: float(x) for a field element x is within ferr(x) of x,
# giving an error bound dc for c, and |acos x - acos y| <= 1.6 sqrt(|x - y|)
# on [-1, 1].

import math as _math


def ferr(x):
    """Bound on |float(x) - x| for a field element x = sum c_s sqrt(s)."""
    k = len(x.t)
    return (k + 2) * 2.3e-16 * sum(abs(float(c)) * _math.sqrt(s) for s, c in x.t.items())


def length_enclosure(b):
    """(lo, hi) floats with lo <= length(b) <= hi, rigorously."""
    if b.pi:
        return _math.pi - 1e-15, _math.pi + 1e-15
    se, ss, ee = dot(b.G, b.S, b.E), dot(b.G, b.S, b.S), dot(b.G, b.E, b.E)
    fss, fee = float(ss), float(ee)
    den = _math.sqrt(fss * fee)
    c = float(se) / den
    v = _math.acos(max(-1.0, min(1.0, c)))
    dc = ferr(se) / den + abs(c) * (ferr(ss) / fss + ferr(ee) / fee) + 1e-15
    err = min(_math.pi, 1.6 * _math.sqrt(dc))
    return v - err, v + err


class Lengths:
    """Bead lengths of one link, with cached enclosures and exact sums.

    A 'chain' is a tuple of bead indices; its value is the sum of their
    lengths."""

    def __init__(self, beads, exact_only=False):
        self.beads, self.exact_only = beads, exact_only
        self._enc = [length_enclosure(b) for b in beads]

    def enc(self, chain):
        lo = sum(self._enc[i][0] for i in chain) - 1e-15 * len(chain)
        hi = sum(self._enc[i][1] for i in chain) + 1e-15 * len(chain)
        return lo, hi

    def exact(self, chain):
        return total_length([self.beads[i] for i in chain])

    def leq(self, A, B):
        """sum(A) <= sum(B) exactly (A, B chains)."""
        if not self.exact_only:
            la, ha = self.enc(A)
            lb, hb = self.enc(B)
            if ha < lb:
                return True
            if la > hb:
                return False
        return angle_leq(self.exact(A), self.exact(B))

    def at_least_two_pi(self, A):
        if not self.exact_only:
            lo, hi = self.enc(A)
            if lo > 2 * _math.pi:
                return True
            if hi < 2 * _math.pi:
                return False
        return self.exact(A).at_least(2)


def angle_leq(x, y):
    """Exact comparison of two AngleSums: total(x) <= total(y)?"""
    if x.h != y.h:
        return x.h < y.h
    # same half-turn [h pi, (h+1) pi): the difference is in (-pi, pi), and its
    # sign is the sign of Im(conj(z_x) z_y)
    im = x.re * y.im - x.im * y.re
    s = im.sign()
    if s != 0:
        return s > 0
    return True   # same argument within the same half-turn: equal


# ----------------------------------------------------------------- True certificates

TRUE_FORMAT = 'vlink-certificate-true/2'
LEMMA_4_9 = 'Elder-McCammond-Meier 2003, Lemma 4.9'


def junction_table(link, beads):
    """For each pair (x, y) with x ending where y starts: True if the
    junction x -> y passes (angle >= pi on both sides). Pairs involving a
    pi-bead count as passing (its geodesic is not fixed): conservative."""
    by_start = {}
    for i, b in enumerate(beads):
        by_start.setdefault(link.lv[b.ts[0]][b.start], []).append(i)
    table = {}
    for x, bx in enumerate(beads):
        w = link.lv[bx.ts[-1]][bx.end]
        for y in by_start.get(w, ()):
            by = beads[y]
            table[(x, y)] = True if (bx.pi or by.pi) else junction_ok(link, bx, by)
    return table


def _chain(tree, s, x):
    """The chain s .. x of bead indices following predecessor pointers."""
    out = [x]
    while x != s:
        x = tree[x]
        out.append(x)
        if len(out) > len(tree) + 2:
            raise Reject('predecessor pointers contain a cycle')
    return tuple(reversed(out))


def check_vertex_trees(link, beads, trees, exact_only=False):
    """Verify the shortest-path certificate for one vertex link: no necklace
    of length < 2*pi. trees[s] maps bead x > s (reached from s) to its
    predecessor; s itself is the root. See docs/true_certificates.md."""
    n = len(beads)
    L = Lengths(beads, exact_only)
    J = junction_table(link, beads)
    succ = {}
    for (x, y), ok in J.items():
        if ok:
            succ.setdefault(x, []).append(y)
    for s in range(n):
        tree = {int(k): int(v) for k, v in trees.get(s, trees.get(str(s), {})).items()}
        if any(x <= s or not 0 <= x < n for x in tree):
            raise Reject('tree of %d: bad bead index' % s)
        chains = {s: (s,)}
        for x in tree:
            c = _chain(tree, s, x)
            if c[0] != s:
                raise Reject('tree of %d: chain does not start at the root' % s)
            for a, b in zip(c, c[1:]):
                if not J.get((a, b)):
                    raise Reject('tree of %d: %d -> %d is not a passing junction' % (s, a, b))
            chains[x] = c
        for x, cx in chains.items():
            for y in succ.get(x, ()):
                if y == s:
                    if not L.at_least_two_pi(cx):
                        raise Reject('bead %d closes a loop shorter than 2*pi through %d' % (s, x))
                    continue
                if y < s:
                    continue
                cxy = cx + (y,)
                if L.at_least_two_pi(cxy):
                    continue
                if y not in chains:
                    raise Reject('tree of %d: %d reachable via %d but missing' % (s, y, x))
                if not L.leq(chains[y], cxy):
                    raise Reject('tree of %d: distance to %d not a lower bound via %d' % (s, y, x))
    return True


def check_true_certificate(cert, exact_only=False):
    """(True, message) if the certificate proves the metric is locally CAT(0)
    (under its listed assumptions), else (False, reason)."""
    try:
        return True, _check_true(cert, exact_only)
    except Reject as e:
        return False, str(e)
    except (KeyError, TypeError, ValueError, IndexError, AttributeError) as e:
        return False, 'malformed certificate: %s: %s' % (type(e).__name__, e)


def _check_true(cert, exact_only):
    if cert.get('format') != TRUE_FORMAT:
        raise Reject('unknown certificate format')
    opts = cert.get('options', {})
    allow_repeats, pi_beads = opts.get('allow_repeats', True), opts.get('pi_beads', True)
    assumptions = list(cert.get('assumptions', []))
    if not pi_beads and LEMMA_4_9 not in assumptions:
        raise Reject('pi_beads=False requires the assumption: ' + LEMMA_4_9)
    if not allow_repeats:
        raise Reject('allow_repeats=False relies on an unproven rule')
    geo = Geometry(cert['isosig'], cert['sqlens'])
    for t in range(geo.T.n):
        if not geo.tet_valid(t):
            raise Reject('tetrahedron %d is degenerate' % t)
    for e in range(len(geo.T.edges)):
        if not geo.edge_angle(e).at_least(2):
            raise Reject('edge %d has dihedral sum < 2*pi' % e)
    verts = {int(v['vertex']): v for v in cert['vertices']}
    if sorted(verts) != list(range(geo.T.n_vertices)):
        raise Reject('certificate must cover every vertex')
    for i_v, v in sorted(verts.items()):
        link = Link(geo, i_v)
        beads = enumerate_beads(link, allow_repeats=True, pi_beads=pi_beads)
        if v.get('n_beads') != len(beads):
            raise Reject('vertex %d: bead count differs from the kernel enumeration' % i_v)
        check_vertex_trees(link, beads, {int(k): t for k, t in v['trees'].items()}, exact_only)
    note = ' (assuming %s)' % '; '.join(assumptions) if assumptions else ''
    return 'locally CAT(0): %d tetrahedra valid, all edges >= 2*pi, no short necklace at %d vertices%s' % (
        geo.T.n, len(verts), note)


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
    n_tet = geo.T.n
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
    if not 0 <= i_v < geo.T.n_vertices:
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
