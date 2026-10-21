"""Combine beads into necklaces: closed local geodesics of length < 2*pi.

A necklace is a cyclic sequence of beads, each ending where the next starts,
such that at every junction (a link vertex W) the angle between the incoming
and outgoing geodesics is >= pi on both sides, and whose total length is
< 2*pi. Finding one shows the vertex link is not CAT(1) (Elder-McCammond).

Junction angles are decided exactly: the triangles around W are developed
(exactly, link.Strip) starting from the incoming bead's last triangle, and the
sweep from the reversed incoming direction u to the outgoing direction x is
tracked with orientation signs det(W, u, .). The sweep reaches pi exactly when
a ray becomes opposite to u, an exact test, so the ties of flat metrics (angles
exactly pi) are decided exactly.

Lengths are arccos values. Sums are compared with 2*pi with a 1e-9 margin;
inside the margin the tie is settled exactly when possible: a necklace whose
junctions are all exactly straight is a great circle, of length exactly 2*pi.
Otherwise UndecidedError is raised (never a silent guess).

Beads of length exactly pi (a lune of geodesics, `Bead.pi_range`) are handled
numerically for now: the geodesic inside the lune is a free angle, and the
junction conditions at both ends become linear inequalities in it (floats,
with near-ties reported as UndecidedError).

Search order and pruning rules follow gallery.find_necklaces: start from each
bead, extend by beads starting at the current end vertex, never reuse a bead
or a link triangle, at most one pi-bead, stop extending once the necklace can
close at its start vertex.
"""

import math

from .link import Strip

TWO_PI = 2 * math.pi
EPS = 1e-9


class UndecidedError(Exception):
    pass


def _tangent_dot_sign(strip, W, a, b):
    """Sign of the dot product of the tangent directions at W toward a and b."""
    ww = strip.dot(W, W)
    return (strip.dot(a, b) * ww - strip.dot(a, W) * strip.dot(b, W)).sign()


def _apply(P, coords):
    """Point with `coords` in a triangle's own basis, given that triangle's
    developed vertex positions P."""
    return tuple(sum((coords[j] * P[j][k] for j in range(3)), coords[0] * 0)
                 for k in range(3))


class Junctions:
    """Exact junction tests between beads of one link, cached."""

    def __init__(self, link):
        self.link = link
        self.cache = {}

    def sides(self, b_in, b_out):
        """For the junction b_in -> b_out, returns a pair, one entry per side:
        'ge' (sweep > pi), 'eq' (exactly pi) or 'lt' (< pi)."""
        key = (id(b_in), id(b_out))
        r = self.cache.get(key)
        if r is None:
            r = self.cache[key] = tuple(self._side(b_in, b_out, k) for k in (1, 2))
        return r

    def ok(self, b_in, b_out):
        return all(s != 'lt' for s in self.sides(b_in, b_out))

    def straight(self, b_in, b_out):
        return 'eq' in self.sides(b_in, b_out)

    def _side(self, b_in, b_out, k):
        link = self.link
        tm, rm = b_in.ts[-1], b_in.end
        t0, d0 = b_out.ts[0], b_out.start
        P = b_in.strip.pos[-1]
        fan = Strip.__new__(Strip)
        fan.link, fan.ts, fan.pos, fan.frame = link, (tm,), (P,), b_in.strip.frame
        fan.steps = ()
        W, u = P[rm], b_in.S
        d = (rm + k) % 3
        y = 3 - rm - d
        o = Strip.det(W, P[d], P[y]).sign()

        def state(ray):
            s = o * Strip.det(W, u, ray).sign()
            if s > 0:
                return 'before'
            if s < 0:
                return 'after'
            return 'eq' if _tangent_dot_sign(fan, W, u, ray) < 0 else 'zero'

        # Walk the fan around W. In each triangle (t, corner c of W) the sweep
        # enters through the ray at vertex `ent` and leaves through the ray at
        # vertex `ex`, then crosses the edge {W, ex} (opposite `ent`).
        t, c, ent, ex = tm, rm, d, y
        if (t, c) == (t0, d0):
            X = _apply(P, b_out.E)
            # outgoing ray inside this same corner, between u and ray(ex)?
            if o * Strip.det(W, u, X).sign() > 0 and o * Strip.det(W, X, P[ex]).sign() >= 0:
                return 'lt'
            # outgoing ray equal to u: the geodesic turns back (angle 0)
            if Strip.det(W, u, X).is_zero() and _tangent_dot_sign(fan, W, u, X) > 0:
                return 'lt'
        for _ in range(4 * link.n + 2):
            st = state(fan.pos[-1][ex])
            if st in ('after', 'eq'):
                return 'ge' if st == 'after' else 'eq'
            g = link.glue[t][ent]
            fan = fan.extend(ent)
            t, c, ent = fan.ts[-1], g[c], g[ex]
            ex = 3 - c - ent
            if (t, c) == (t0, d0):
                st = state(_apply(fan.pos[-1], b_out.E))
                return {'after': 'ge', 'eq': 'eq'}.get(st, 'lt')
        raise RuntimeError('fan walk around link vertex did not reach the outgoing triangle')


def search(link, beads, findall=False, allow_repeats=False):
    """Returns a necklace (list of beads) or None; with findall, a list of all
    necklaces found (each as a list of beads). allow_repeats lets different
    beads of a necklace share link triangles (beads themselves are never
    reused: a closed geodesic traversing a bead twice covers a shorter one)."""
    J = Junctions(link)
    by_start = {}
    for i, b in enumerate(beads):
        by_start.setdefault(link.lv[b.ts[0]][b.start], []).append(i)
    for js in by_start.values():
        js.sort(key=lambda j: beads[j].length)   # short beads first
    found = []
    # Lower bounds for pruning: shortest bead-path length between link
    # vertices, ignoring junction conditions (so it never cuts a valid loop).
    nv = link.n_vertices
    INF = float('inf')
    dist = [[0.0 if x == y else INF for y in range(nv)] for x in range(nv)]
    for b in beads:
        x, y = link.lv[b.ts[0]][b.start], link.lv[b.ts[-1]][b.end]
        dist[x][y] = min(dist[x][y], b.length)
    for k in range(nv):
        dk = dist[k]
        for x in range(nv):
            dxk = dist[x][k]
            if dxk == INF:
                continue
            dx = dist[x]
            for y in range(nv):
                if dxk + dk[y] < dx[y]:
                    dx[y] = dxk + dk[y]

    def dfs(path, total, tris, has_pi):
        gm = beads[path[-1]]
        g0 = beads[path[0]]
        w = link.lv[gm.ts[-1]][gm.end]
        if w == link.lv[g0.ts[0]][g0.start]:
            neck = Necklace(beads[i] for i in path)
            if gm.pi_range or g0.pi_range or any(b.pi_range for b in neck):
                neck.pi_point = _pi_junctions(J, neck)
                ok = neck.pi_point is not None
            else:
                ok = J.ok(gm, g0)
            if ok:
                found.append(neck)
                return not findall
            return False
        s0 = link.lv[g0.ts[0]][g0.start]
        for j in by_start.get(w, ()):
            # each loop is searched only from its lowest-numbered bead
            if j <= path[0] or j in path:
                continue
            b = beads[j]
            if has_pi and b.pi_range:
                continue
            new_total = total + b.length
            # cannot get back to the start within 2*pi (with float slack)
            if new_total + dist[link.lv[b.ts[-1]][b.end]][s0] > TWO_PI + EPS:
                continue
            if abs(new_total - TWO_PI) <= EPS:
                # Near-tie: decide total < 2*pi exactly (see cmp_two_pi).
                if cmp_two_pi([beads[i] for i in path] + [b]) >= 0:
                    continue
            elif new_total >= TWO_PI:
                continue
            if not allow_repeats and tris & set(b.ts):
                continue
            if not (gm.pi_range or b.pi_range):
                if not J.ok(gm, b):
                    continue
            if dfs(path + [j], new_total, tris | set(b.ts), has_pi or bool(b.pi_range)):
                return True
        return False

    for i, b in enumerate(beads):
        if dfs([i], b.length, set(b.ts), bool(b.pi_range)):
            break
    if findall:
        return found
    return found[0] if found else None


def cmp_two_pi(beads):
    """Exact sign of (sum of bead lengths - 2*pi), for a sum known to be within
    EPS of 2*pi.

    A bead of length theta in [0, pi] contributes the complex number
    a + i*b = |S||E| e^{i theta} with a = S.E and b = sqrt(|S|^2 |E|^2 - a^2) >= 0.
    The product z of these has argument = the total length (mod 2*pi), so
    with the total within EPS of 2*pi: total == 2*pi iff Im z == 0, and
    otherwise sign(total - 2*pi) = sign(Im z). b is exact when |S|^2|E|^2 - a^2
    is rational, which holds whenever the link's field is Q (e.g. the cube and
    disphenoid shapes); otherwise UndecidedError.
    """
    re, im = MQone, MQzero
    for bd in beads:
        st = bd.strip
        a = st.dot(bd.S, bd.E)
        b2 = st.dot(bd.S, bd.S) * st.dot(bd.E, bd.E) - a * a
        if any(s != 1 for s in b2.t):
            raise UndecidedError('length tie in a non-rational link field')
        b = _MQ.sqrt_rational(b2.t.get(1, 0))
        re, im = re * a - im * b, re * b + im * a
    if (re.sign() <= 0) and not im.is_zero():
        raise UndecidedError('length tie: product not near the positive real axis')
    return im.sign()


# ------------------------------------------------------------------ pi-beads

def _inv3(M):
    """Inverse of a 3x3 matrix of MQ (rows)."""
    a = M
    cof = [[a[(i + 1) % 3][(j + 1) % 3] * a[(i + 2) % 3][(j + 2) % 3]
            - a[(i + 1) % 3][(j + 2) % 3] * a[(i + 2) % 3][(j + 1) % 3]
            for i in range(3)] for j in range(3)]  # adjugate (transposed cofactors)
    det = sum((a[0][k] * cof[k][0] for k in range(3)), a[0][0] * 0)
    inv_det = det.inverse()
    return [[cof[i][j] * inv_det for j in range(3)] for i in range(3)]


def _to_basis(P, X):
    """Coordinates of the point X (in some frame) w.r.t. a triangle's
    developed vertex vectors P (same frame)."""
    M = [[P[j][i] for j in range(3)] for i in range(3)]  # columns = P[j]
    Mi = _inv3(M)
    return tuple(sum((Mi[i][k] * X[k] for k in range(3)), X[0] * 0) for i in range(3))


def _opposite_ray(strip, W, u):
    """A point whose ray from W is opposite to the ray toward u."""
    s = strip.dot(u, W) / strip.dot(W, W)
    return tuple(2 * s * w - x for w, x in zip(W, u))


def _boundary(link, P, frame, t, c, u, k, target):
    """Walk the fan around the link vertex at corner c of t (positions P in
    `frame`) from the ray toward u, on side k, to the corner `target`.

    Returns ('all',) if every ray of the target corner is swept >= pi,
    ('none',) if none is, or ('from', B, Pt, ent, ex): rays of the target corner
    swept >= pi are those from the ray toward B (a point, same frame as Pt) on
    toward the exit ray `ex`; Pt are the target triangle's developed vertices.
    """
    if (t, c) == target:
        # the incoming and outgoing corners coincide: the corner is split by u
        # and would need two constraints; rare, not handled yet
        raise UndecidedError('pi-bead junction inside a single corner')
    fan = Strip.__new__(Strip)
    fan.link, fan.ts, fan.pos, fan.frame = link, (t,), (P,), frame
    fan.steps = ()
    W = P[c]
    d = (c + k) % 3
    y = 3 - c - d
    o = Strip.det(W, P[d], P[y]).sign()

    def state(ray):
        s = o * Strip.det(W, u, ray).sign()
        if s > 0:
            return 'before'
        if s < 0:
            return 'after'
        return 'eq' if _tangent_dot_sign(fan, W, u, ray) < 0 else 'zero'

    ent, ex = d, y
    for _ in range(4 * link.n + 2):
        st = state(fan.pos[-1][ex])
        if st in ('after', 'eq') and (fan.ts[-1], c) != target:
            # pi reached before the target corner
            return ('all',)
        g = link.glue[t][ent]
        fan = fan.extend(ent)
        t, c, ent = fan.ts[-1], g[c], g[ex]
        ex = 3 - c - ent
        if (t, c) == target:
            Pt = fan.pos[-1]
            st_in, st_out = state(Pt[ent]), state(Pt[ex])
            if st_in in ('after', 'eq'):
                return ('all',)
            if st_out == 'before':
                return ('none',)
            return ('from', _opposite_ray(fan, W, u), Pt, ent, ex)
    raise RuntimeError('fan walk did not reach the target corner')


def _pi_feasible(link, prev, p, nxt):
    """Is there a geodesic in the pi-bead p's lune such that the junctions
    prev -> p (at S) and p -> nxt (at -S) are both >= pi on both sides?"""
    S, E = p.S, p.E
    L, R = p.pi_range
    pstrip = p.strip
    cons = []  # (tau, B): need tau * det(S, B, w) >= 0, all in p's frame
    # at S: walk around S from prev's incoming direction to p's first corner
    for k in (1, 2):
        r = _boundary(link, prev.strip.pos[-1], prev.strip.frame, prev.ts[-1], prev.end,
                      prev.S, k, (p.ts[0], p.start))
        if r[0] == 'none':
            return False
        if r[0] == 'from':
            _, B, Pt, ent, ex = r
            Bp = _to_basis(Pt, B)                     # p's frame = its first triangle's basis
            e_ent = tuple(MQone if i == ent else MQzero for i in range(3))
            e_ex = tuple(MQone if i == ex else MQzero for i in range(3))
            sigma = Strip.det(S, e_ent, e_ex).sign()  # walk orientation in p's frame
            cons.append((sigma, Bp))
    # at -S: walk around -S from nxt's outgoing direction to p's last corner
    Plast = pstrip.pos[-1]
    for k in (1, 2):
        r = _boundary(link, nxt.strip.pos[0], nxt.strip.frame, nxt.ts[0], nxt.start,
                      nxt.E, k, (p.ts[-1], p.end))
        if r[0] == 'none':
            return False
        if r[0] == 'from':
            _, B, Pt, ent, ex = r
            Bp = _apply(Plast, _to_basis(Pt, B))
            sigma = Strip.det(E, Plast[ent], Plast[ex]).sign()
            # det(E, ., .) = -lambda det(S, ., .) with lambda > 0
            cons.append((-sigma, Bp))
    # clip the open range (L, R) by the closed half-planes
    lo, lo_closed, hi, hi_closed = L, False, R, False
    for tau, B in cons:
        f_lo = tau * Strip.det(S, B, lo).sign()
        f_hi = tau * Strip.det(S, B, hi).sign()
        if f_lo >= 0 and f_hi >= 0:
            continue
        if f_lo < 0 and f_hi < 0:
            return False
        # the boundary ray (B or its opposite about S) lies inside (lo, hi)
        Z = B if _inside_arc(S, lo, hi, B) else _opposite_ray(pstrip, S, B)
        if f_lo < 0:
            lo, lo_closed = Z, True
        else:
            hi, hi_closed = Z, True
    # Feasible: return a point Q (p's frame) whose ray from S gives a suitable
    # geodesic in the lune (the certificate records it); None if infeasible.
    s = Strip.det(S, lo, hi).sign()
    if s > 0:
        return tuple(x + y for x, y in zip(lo, hi))   # strictly inside (lo, hi)
    if s == 0 and lo_closed and hi_closed and _tangent_dot_sign(pstrip, S, lo, hi) > 0:
        return lo                                      # a single admissible ray
    return None


def _inside_arc(S, lo, hi, X):
    return Strip.det(S, lo, X).sign() >= 0 and Strip.det(S, X, hi).sign() >= 0


def _pi_junctions(J, beads):
    """Necklace (cyclic list of beads) containing exactly one pi-bead: is there
    a geodesic in its lune making every junction >= pi? Ordinary junctions
    are tested exactly; the two junctions at the pi-bead by _pi_feasible."""
    m = len(beads)
    i = next(k for k, b in enumerate(beads) if b.pi_range)
    for k in range(m):
        a, b = beads[k], beads[(k + 1) % m]
        if a.pi_range or b.pi_range:
            continue
        if not J.ok(a, b):
            return None
    return _pi_feasible(J.link, beads[i - 1], beads[i], beads[(i + 1) % m])


class Necklace(list):
    """A found necklace: list of beads; pi_point is the geodesic direction
    chosen in the pi-bead's lune (a point in that bead's frame), if any."""
    pi_point = None


from .field import MQ as _MQ  # noqa: E402
MQone, MQzero = _MQ.rational(1), _MQ()
