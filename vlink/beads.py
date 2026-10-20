"""Enumerate the beads (vertex-to-vertex galleries) of a vertex link exactly.

A bead from link vertex S (corner `start` of triangle ts[0]) to link vertex E
(corner `end` of triangle ts[-1]) is a strip of link triangles through which
the geodesic S -> E runs, crossing every interior edge of the strip in its
relative interior (Elder-McCammond 2002, Sec. 3). These are the candidate
pieces of short closed geodesics in the link.

Everything is decided exactly with sign tests in the link's multiquadratic
field (no arccos), using one fact: seen from S, the relative interior of an
edge AB (S not on its great circle) is hit by exactly the rays from S
strictly inside the angle at S of the spherical triangle SAB, which is < pi.
Rays are compared by orientation, sign det(S, X, Y). A strip stays
*eligible* while the intersection of these open angular intervals over its
crossed edges is non-empty; vertex E of the next triangle ends a bead iff its
ray is strictly inside that interval. If E is antipodal to S the bead has
length pi and carries the whole open interval of geodesics (`pi_range`).

Rules kept from gallery.py so the two implementations can be compared:
  - thin beads: each edge of each triangle, from both ends;
  - a bead never visits the same link triangle twice (tobead);
  - a strip is abandoned when an edge endpoint is parallel to S (zero_zero /
    zero_pi in Gallery.eligible).
"""

import math

from .link import Strip


class Bead:
    __slots__ = ('ts', 'start', 'end', 'strip', 'S', 'E', 'cos2', 'cos_sign',
                 'length', 'pi_range')

    def __init__(self, strip, start, end, S, E, pi_range=None):
        self.ts = strip.ts
        self.start, self.end = start, end
        self.strip, self.S, self.E = strip, S, E
        se, ss, ee = strip.dot(S, E), strip.dot(S, S), strip.dot(E, E)
        self.cos_sign = se.sign()
        self.cos2 = (se * se) / (ss * ee)          # exact cos^2 of the length
        self.pi_range = pi_range
        c = float(se) / math.sqrt(float(ss) * float(ee))
        self.length = math.pi if pi_range else math.acos(max(-1.0, min(1.0, c)))

    def key(self):
        return (self.ts, self.start, self.end)

    def __repr__(self):
        return 'Bead(ts=%s, %d->%d, len/pi=%.6f%s)' % (
            self.ts, self.start, self.end, self.length / math.pi,
            ', pi-range' if self.pi_range else '')


def _parallel(strip, X, Y):
    """+1 if X, Y point the same way, -1 if opposite, 0 if not parallel."""
    cx = (X[1] * Y[2] - X[2] * Y[1], X[2] * Y[0] - X[0] * Y[2], X[0] * Y[1] - X[1] * Y[0])
    if any(not c.is_zero() for c in cx):
        return 0
    return strip.dot(X, Y).sign()


def _visible(S, A, B):
    """Oriented visible interval (L, R) of edge AB from S, or None if S lies on
    the great circle through A, B."""
    o = Strip.det(S, A, B).sign()
    if o == 0:
        return None
    return (A, B) if o > 0 else (B, A)


def _inside(S, L, R, X, strict=True):
    a, b = Strip.det(S, L, X).sign(), Strip.det(S, X, R).sign()
    return (a > 0 and b > 0) if strict else (a >= 0 and b >= 0)


def _intersect(S, I1, I2):
    """Intersection of two oriented open angular intervals about S (each of
    width < pi), or None if empty."""
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
    return (L, R) if Strip.det(S, L, R).sign() > 0 else None


def beads_from(link, t, allow_repeats=False):
    """All beads whose first triangle is t (thin and thick), in a
    deterministic order. allow_repeats lifts the no-repeated-triangle rule."""
    out = []
    for d0 in range(3):
        s0 = link.start(t)
        P = s0.pos[0]
        a, b = (d0 + 1) % 3, (d0 + 2) % 3
        S = P[d0]
        # thin beads along the two edges at d0
        out.append(Bead(s0, d0, a, S, P[a]))
        out.append(Bead(s0, d0, b, S, P[b]))
        # thick beads: cross the edge opposite d0 first
        if _parallel(s0, S, P[a]) or _parallel(s0, S, P[b]):
            continue
        iv = _visible(S, P[a], P[b])
        if iv is not None:
            _grow(s0, d0, S, d0, iv, out, allow_repeats)
    return out


def _grow(strip, start, S, d, interval, out, allow_repeats=False):
    """Cross the edge opposite vertex d of the strip's last triangle."""
    t = strip.ts[-1]
    t2 = strip.link.adj[t][d]
    if t2 in strip.ts and not allow_repeats:
        return  # a bead never repeats a triangle (as in gallery.py)
    if len(strip.ts) > 4 * strip.link.n:
        raise RuntimeError('bead strip longer than 4x the link: growth not terminating?')
    s2 = strip.extend(d)
    g = strip.link.glue[t][d]
    r = g[d]                      # the new vertex of the new triangle
    P2 = s2.pos[-1]
    E = P2[r]
    par = _parallel(s2, S, E)
    if par < 0:
        out.append(Bead(s2, start, r, S, E, pi_range=interval))
    elif par == 0 and _inside(S, interval[0], interval[1], E):
        out.append(Bead(s2, start, r, S, E))
    # continue across the two edges of the new triangle that contain E
    for dn in ((r + 1) % 3, (r + 2) % 3):
        y = 3 - r - dn            # the other endpoint of the edge opposite dn
        if par or _parallel(s2, S, P2[y]):
            continue
        iv = _visible(S, E, P2[y])
        if iv is None:
            continue
        new = _intersect(S, interval, iv)
        if new is not None:
            _grow(s2, start, S, dn, new, out, allow_repeats)


def all_beads(link, allow_repeats=False):
    return [b for t in range(link.n) for b in beads_from(link, t, allow_repeats)]
