"""Exact arithmetic in multiquadratic fields Q(sqrt(q1), sqrt(q2), ...).

An element is a finite sum  sum_s c_s * sqrt(s)  over distinct square-free
positive integers s (s = 1 is the rational part) with rational c_s. Square
roots of distinct square-free integers are linearly independent over Q, so

  * equality / zero tests are exact coefficient comparisons, and
  * signs are decided exactly: a float estimate with an error bound settles
    almost every case; otherwise x = A + B*sqrt(p) is split on a prime p and
    sign(A + B sqrt p) follows from sign(A), sign(B) and sign(A^2 - p B^2),
    recursively in fewer primes.

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


def squarefree_split(n):
    """n = k^2 * s with s square-free; returns (k, s). n is a positive int."""
    k, s, p = 1, 1, 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        k *= p ** (e // 2)
        if e % 2:
            s *= p
        p += 1 if p == 2 else 2
    return k, s * n


def _smallest_prime(n):
    p = 2
    while p * p <= n:
        if n % p == 0:
            return p
        p += 1 if p == 2 else 2
    return n


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
        return self + (-other if isinstance(other, MQ) else MQ.rational(-Fraction(other)))

    def __rsub__(self, other):
        return MQ.rational(other) - self

    def __mul__(self, other):
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
        """Galois conjugate sqrt(p) -> -sqrt(p) for a prime p."""
        return MQ({s: (-c if s % p == 0 else c) for s, c in self.t.items()})

    def inverse(self):
        if not self.t:
            raise ZeroDivisionError('MQ zero')
        num, den = MQ({1: Fraction(1)}), self
        while True:
            primes = [s for s in den.t if s != 1]
            if not primes:
                return num * (1 / den.t[1])
            p = _smallest_prime(primes[0])
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
        p = _smallest_prime(next(s for s in t if s != 1))
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
