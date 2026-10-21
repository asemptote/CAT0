"""Exact arithmetic in multiquadratic fields Q(sqrt(q1), sqrt(q2), ...).

An element is a finite sum  sum_s c_s * sqrt(s)  over distinct square-free
positive integers s (s = 1 is the rational part) with rational c_s. Square
roots of distinct square-free integers are linearly independent over Q, so

  * equality / zero tests are exact coefficient comparisons, and
  * signs are decided exactly: a float estimate with an error bound settles
    almost every case; otherwise x = A + B*sqrt(p) is split on a square-free
    p > 1 that divides or is coprime to every key (_splitter: gcd refinement,
    no factoring), and sign(A + B sqrt p) follows from sign(A), sign(B) and
    sign(A^2 - p B^2), recursively with p eliminated.

This replaces the nested-radical Constructible numbers of constructible.py:
in the vertex-link computations every coordinate lives in one fixed field of
this kind (see link.py), so expressions never nest and never grow in depth.
"""

from fractions import Fraction
from math import gcd, isqrt, sqrt as fsqrt

_SQRT_CACHE = {}


def _fsqrt(s):
    r = _SQRT_CACHE.get(s)
    if r is None:
        r = _SQRT_CACHE[s] = fsqrt(s)
    return r


_SQF_CACHE = {}


def squarefree_split(n):
    """n = k^2 * s with s square-free; returns (k, s). n is a positive int.

    Trial division only up to the cube root of n: what remains then has at
    most two prime factors, both larger than that, so it is either a perfect
    square p^2 (detected exactly with isqrt) or square-free."""
    r = _SQF_CACHE.get(n)
    if r is not None:
        return r
    n0, k, s, p = n, 1, 1, 2
    while p * p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        k *= p ** (e // 2)
        if e % 2:
            s *= p
        p += 1 if p == 2 else 2
    q = isqrt(n)
    if q * q == n:
        k *= q
    else:
        s *= n
    r = _SQF_CACHE[n0] = (k, s)
    return r


def _splitter(keys):
    """A square-free s > 1 such that every key is divisible by s or coprime to
    it (keys: square-free integers, some > 1). Found by gcd refinement; no
    factoring needed."""
    s = next(k for k in keys if k != 1)
    changed = True
    while changed:
        changed = False
        for k in keys:
            g = gcd(s, k)
            if 1 < g < s:
                s, changed = g, True
    return s


class MQ:
    """Element of a multiquadratic field: {square-free s: Fraction c_s}."""

    __slots__ = ('t',)

    def __init__(self, terms=None):
        # zero coefficients are never stored, so is_zero()/== stay exact
        self.t = {s: c for s, c in terms.items() if c} if terms else {}

    # -- construction ----------------------------------------------------
    @staticmethod
    def rational(q):
        q = Fraction(q)
        return MQ({1: q} if q else {})

    @staticmethod
    def sqrt_rational(q):
        """Exact sqrt of a non-negative rational q."""
        q = Fraction(q)
        if q < 0:
            raise ValueError('sqrt of negative rational %s' % q)
        if q == 0:
            return MQ()
        # sqrt(a/b) = sqrt(a*b)/b
        k, s = squarefree_split(q.numerator * q.denominator)
        return MQ({s: Fraction(k, q.denominator)})

    # -- arithmetic --------------------------------------------------------
    def __add__(self, other):
        if isinstance(other, Tower):
            return NotImplemented
        if not isinstance(other, MQ):
            other = MQ.rational(other)
        t = dict(self.t)
        for s, c in other.t.items():
            v = t.get(s, 0) + c
            if v:
                t[s] = v
            else:
                t.pop(s, None)
        return MQ(t)

    __radd__ = __add__

    def __neg__(self):
        return MQ({s: -c for s, c in self.t.items()})

    def __sub__(self, other):
        if isinstance(other, Tower):
            return NotImplemented
        return self + (-other if isinstance(other, MQ) else MQ.rational(-Fraction(other)))

    def __rsub__(self, other):
        return MQ.rational(other) - self

    def __mul__(self, other):
        if isinstance(other, Tower):
            return NotImplemented
        if not isinstance(other, MQ):
            q = Fraction(other)
            return MQ({s: c * q for s, c in self.t.items()} if q else {})
        t = {}
        for s1, c1 in self.t.items():
            for s2, c2 in other.t.items():
                g = gcd(s1, s2)
                s = (s1 // g) * (s2 // g)
                v = t.get(s, 0) + c1 * c2 * g
                if v:
                    t[s] = v
                else:
                    t.pop(s, None)
        return MQ(t)

    __rmul__ = __mul__

    def conj(self, p):
        """Galois conjugate sqrt(p) -> -sqrt(p), for p as returned by
        _splitter (every key divisible by p or coprime to it)."""
        return MQ({s: (-c if s % p == 0 else c) for s, c in self.t.items()})

    def inverse(self):
        if not self.t:
            raise ZeroDivisionError('MQ zero')
        num, den = MQ({1: Fraction(1)}), self
        while True:
            if all(s == 1 for s in den.t):
                return num * (1 / den.t[1])
            p = _splitter(list(den.t))
            c = den.conj(p)
            num, den = num * c, den * c

    def __truediv__(self, other):
        if not isinstance(other, MQ):
            return self * (1 / Fraction(other))
        return self * other.inverse()

    def __rtruediv__(self, other):
        return MQ.rational(other) * self.inverse()

    # -- predicates ----------------------------------------------------------
    def is_zero(self):
        return not self.t

    def __eq__(self, other):
        if not isinstance(other, MQ):
            other = MQ.rational(other)
        return self.t == other.t

    def __hash__(self):
        return hash(frozenset(self.t.items()))

    def __float__(self):
        return sum((float(c) * _fsqrt(s) for s, c in self.t.items()), 0.0)

    def sign(self):
        """Exact sign (-1, 0, 1)."""
        t = self.t
        if not t:
            return 0
        if len(t) == 1:
            (s, c), = t.items()
            return 1 if c > 0 else -1
        # float estimate with a generous error bound
        est, mag = 0.0, 0.0
        for s, c in t.items():
            v = float(c) * _fsqrt(s)
            est += v
            mag += abs(v)
        if abs(est) > 1e-9 * mag:
            return 1 if est > 0 else -1
        return self._exact_sign()

    def _exact_sign(self):
        t = self.t
        if not t:
            return 0
        if all(s == 1 for s in t):
            return 1 if t[1] > 0 else -1
        p = _splitter(list(t))
        A = MQ({s: c for s, c in t.items() if s % p})
        B = MQ({s // p: c for s, c in t.items() if s % p == 0})
        sa, sb = A.sign(), B.sign()
        if sb == 0:
            return sa
        if sa == 0 or sa == sb:
            return sb if sa == 0 else sa
        # opposite signs: compare A^2 with p*B^2
        return sa * (A * A - B * B * p).sign()

    def __lt__(self, other):
        return (self - other).sign() < 0

    def __le__(self, other):
        return (self - other).sign() <= 0

    def __gt__(self, other):
        return (self - other).sign() > 0

    def __ge__(self, other):
        return (self - other).sign() >= 0

    def __repr__(self):
        if not self.t:
            return 'MQ(0)'
        parts = []
        for s in sorted(self.t):
            c = self.t[s]
            parts.append(str(c) if s == 1 else '%s*sqrt(%d)' % (c, s))
        return 'MQ(%s)' % ' + '.join(parts)


ZERO = MQ()
ONE = MQ.rational(1)


# --------------------------------------------------------------------- towers
#
# For generic edge lengths the sine of an angle is sqrt(b2) with b2 in the
# link's multiquadratic field F but not rational. A Tower element is
#     sum over sets K of symbols:  c_K * prod_{i in K} sqrt(r_i)
# with c_K in F and radicands r_i > 0 in F (registered globally, one symbol
# per distinct radicand). Products reduce with sqrt(r_i)^2 = r_i.
#
# The representation need not be unique (radicands may be algebraically
# related, e.g. sqrt(2+sqrt3) * sqrt(2-sqrt3) = 1), so equality and zero
# tests go through the sign, which is exact: rigorous interval arithmetic
# first, and otherwise split on one symbol, x = A + B sqrt(r), with
# sign x = sign(A) * sign(A^2 - r B^2) when A, B have opposite signs.

_RADICANDS = []          # i -> radicand r_i (an MQ, > 0)
_RADICAND_INDEX = {}     # MQ -> i
_U = 2.3e-16             # unit roundoff (with margin)


def _mq_interval(x):
    """Rigorous float interval containing the field element x."""
    v = float(x)
    e = (len(x.t) + 2) * _U * sum(abs(float(c)) * fsqrt(s) for s, c in x.t.items()) + 1e-300
    return v - e, v + e


def _imul(a, b):
    ps = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    lo, hi = min(ps), max(ps)
    m = 4 * _U * max(abs(lo), abs(hi))
    return lo - m, hi + m


class Tower:
    __slots__ = ('t',)

    def __init__(self, terms=None):
        self.t = {k: c for k, c in terms.items() if not c.is_zero()} if terms else {}

    @staticmethod
    def lift(x):
        if isinstance(x, Tower):
            return x
        if not isinstance(x, MQ):
            x = MQ.rational(x)
        return Tower({frozenset(): x})

    @staticmethod
    def sqrt_field(r):
        """sqrt(r) for r >= 0 in F: exact in F when r is rational, otherwise a
        new (or reused) symbol."""
        if r.sign() < 0:
            raise ValueError('sqrt of a negative field element')
        if r.is_zero():
            return Tower()
        if all(s == 1 for s in r.t):
            return Tower.lift(MQ.sqrt_rational(r.t[1]))
        key = frozenset(r.t.items())
        i = _RADICAND_INDEX.get(key)
        if i is None:
            i = _RADICAND_INDEX[key] = len(_RADICANDS)
            _RADICANDS.append(r)
        return Tower({frozenset((i,)): ONE})

    def __add__(self, other):
        other = Tower.lift(other)
        t = dict(self.t)
        for k, c in other.t.items():
            v = t.get(k, ZERO) + c
            if v.is_zero():
                t.pop(k, None)
            else:
                t[k] = v
        return Tower(t)

    __radd__ = __add__

    def __neg__(self):
        return Tower({k: -c for k, c in self.t.items()})

    def __sub__(self, other):
        return self + (-Tower.lift(other))

    def __rsub__(self, other):
        return Tower.lift(other) - self

    def __mul__(self, other):
        other = Tower.lift(other)
        t = {}
        for k1, c1 in self.t.items():
            for k2, c2 in other.t.items():
                c = c1 * c2
                for i in k1 & k2:
                    c = c * _RADICANDS[i]
                k = k1 ^ k2
                v = t.get(k, ZERO) + c
                if v.is_zero():
                    t.pop(k, None)
                else:
                    t[k] = v
        return Tower(t)

    __rmul__ = __mul__

    def _interval(self):
        lo = hi = 0.0
        sq = {}
        for k, c in self.t.items():
            iv = _mq_interval(c)
            for i in k:
                if i not in sq:
                    rl, rh = _mq_interval(_RADICANDS[i])
                    sq[i] = (fsqrt(max(rl, 0.0)) * (1 - _U), fsqrt(max(rh, 0.0)) * (1 + _U))
                iv = _imul(iv, sq[i])
            lo, hi = lo + iv[0], hi + iv[1]
        m = (len(self.t) + 1) * _U * max(abs(lo), abs(hi))
        return lo - m, hi + m

    def sign(self):
        if not self.t:
            return 0
        lo, hi = self._interval()
        if lo > 0:
            return 1
        if hi < 0:
            return -1
        syms = set().union(*self.t)
        if not syms:
            return self.t.get(frozenset(), ZERO).sign()
        i = max(syms)
        A = Tower({k: c for k, c in self.t.items() if i not in k})
        B = Tower({k - {i}: c for k, c in self.t.items() if i in k})
        sa, sb = A.sign(), B.sign()
        if sb == 0:
            return sa
        if sa == 0 or sa == sb:
            return sb if sa == 0 else sa
        return sa * (A * A - B * B * _RADICANDS[i]).sign()

    def is_zero(self):
        return self.sign() == 0

    def __float__(self):
        lo, hi = self._interval()
        return (lo + hi) / 2

    def __repr__(self):
        return 'Tower(%s)' % ', '.join('%s*%s' % (c, sorted(k)) for k, c in self.t.items())
