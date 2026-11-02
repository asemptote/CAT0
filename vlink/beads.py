"""Beads for the search: the kernel's bead enumeration plus float lengths.

A bead from link vertex S (corner `start` of triangle ts[0]) to link vertex E
(corner `end` of triangle ts[-1]) is a strip of link triangles through which
the geodesic S -> E runs, crossing every interior edge of the strip in its
relative interior (Elder-McCammond 2002, Sec. 3).

The enumeration itself is kernel code (kernel.enumerate_beads), since the
completeness of the bead list is part of what a True answer claims. This
module wraps each kernel bead for the search, adding its length as a float
with a rigorous error bound (Bead.length_err), used to order and prune the
necklace search; exact comparisons are left to the kernel.
"""

import math

from . import kernel
from .link import Strip


def _ferr(x):
    """Bound on |float(x) - x| for a field element x = sum c_s sqrt(s): each
    term is rounded (float(c), sqrt) with relative error <= 2.3e-16 and the
    summation of k terms adds at most k * 1.2e-16 * sum |terms|."""
    k = len(x.t)
    return (k + 2) * 2.3e-16 * sum(abs(float(c)) * math.sqrt(s) for s, c in x.t.items())


class Bead:
    __slots__ = ('ts', 'start', 'end', 'strip', 'S', 'E', 'cos2', 'cos_sign',
                 'length', 'length_err', 'pi_range')

    def __init__(self, strip, start, end, S, E, pi_range=None):
        self.ts = strip.ts
        self.start, self.end = start, end
        self.strip, self.S, self.E = strip, S, E
        se, ss, ee = strip.dot(S, E), strip.dot(S, S), strip.dot(E, E)
        self.cos_sign = se.sign()
        self.cos2 = (se * se) / (ss * ee)          # exact cos^2 of the length
        self.pi_range = pi_range
        if pi_range:
            self.length, self.length_err = math.pi, 1e-15
            return
        fss, fee = float(ss), float(ee)
        den = math.sqrt(fss * fee)
        c = float(se) / den
        self.length = math.acos(max(-1.0, min(1.0, c)))
        # Rigorous bound on the float error of `length`: dc bounds the error of
        # c (numerator error / den, relative errors of ss and ee, and the
        # rounding of sqrt and division); |acos x - acos y| <= 1.6 sqrt(|x - y|)
        # on [-1, 1] turns it into an error bound for the angle.
        dc = (_ferr(se) / den + abs(c) * (_ferr(ss) / fss + _ferr(ee) / fee)
              + 1e-15)
        self.length_err = min(math.pi, 1.6 * math.sqrt(dc))

    # the view kernel.junction_side / kernel.total_length expect
    @property
    def pos(self):
        return self.strip.pos

    @property
    def G(self):
        return self.strip.G

    @property
    def pi(self):
        return bool(self.pi_range)

    @property
    def in_dir(self):
        return self.S

    @property
    def out_dir(self):
        return self.E

    def key(self):
        return (self.ts, self.start, self.end)

    def __repr__(self):
        return 'Bead(ts=%s, %d->%d, len/pi=%.6f%s)' % (
            self.ts, self.start, self.end, self.length / math.pi,
            ', pi-range' if self.pi_range else '')


def _wrap(link, kb):
    strip = Strip.__new__(Strip)
    strip.link, strip.ts, strip.pos, strip.frame = link, kb.ts, tuple(kb.pos), kb.ts[0]
    strip.steps = kb.steps
    return Bead(strip, kb.start, kb.end, kb.S, kb.E, pi_range=kb.lune)


def beads_from(link, t, allow_repeats=False, pi_beads=True):
    """All beads whose first triangle is t, in a deterministic order."""
    return [_wrap(link, kb) for kb in
            kernel.enumerate_beads(link.k, allow_repeats, pi_beads, triangles=[t])]


def all_beads(link, allow_repeats=False, pi_beads=True):
    return [_wrap(link, kb) for kb in kernel.enumerate_beads(link.k, allow_repeats, pi_beads)]
