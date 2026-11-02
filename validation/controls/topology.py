"""Classify closed census manifolds (by Regina's name) as able / unable to be CAT(0).

A closed 3-manifold with a locally CAT(0) metric is aspherical (Cartan-Hadamard),
and its fundamental group contains no Nil, Sol or SL2~ subgroups (flat torus
theorem / solvable subgroup theorem). So for every metric on these manifolds
the vertex-link checker must return False (given the edge condition holds):

  spherical  S3, RP3, L(p,q), Seifert fibred over S^2 with orbifold chi > 0
  S2xS1      S2 x S1
  Nil        Seifert fibred, Euclidean base orbifold (chi = 0), Euler number != 0
  SL2~       Seifert fibred, hyperbolic base orbifold (chi < 0), Euler number != 0
  Sol        torus bundle T x I / [[a,b],[c,d]] with |trace| > 2

Manifolds that may be CAT(0) (flat, H2xR, hyperbolic, ...) get expected=None:
the answer depends on the metric. (A *flat metric* must give True; that is
checked from the metric itself, see verify_cert.metric_report.)

Only Seifert fibred spaces over S^2 written as 'SFS [S2: (a1,b1) ... (ak,bk)]'
and torus bundles 'T x I / [ a,b | c,d ]' are parsed; other names that are not
in the list above are classified 'other'.
"""

import re
from fractions import Fraction

SFS_S2 = re.compile(r'^SFS \[S2:((?: \(-?\d+,-?\d+\))+)\]$')
FIBRE = re.compile(r'\((-?\d+),(-?\d+)\)')
TXI = re.compile(r'^T x I / \[ *(-?\d+),(-?\d+) *\| *(-?\d+),(-?\d+) *\]$')


def classify(name):
    """Returns (geometry, expected) with expected False (cannot be CAT(0)) or None."""
    name = name.strip()
    if name in ('S3', 'RP3') or name.startswith('L('):
        return 'spherical', False
    if name == 'S2 x S1':
        return 'S2xS1', False
    m = SFS_S2.match(name)
    if m:
        fibres = [(int(a), int(b)) for a, b in FIBRE.findall(m.group(1))]
        chi = 2 - sum(1 - Fraction(1, abs(a)) for a, _ in fibres)
        e = -sum(Fraction(b, a) for a, b in fibres)
        if chi > 0:
            return 'spherical', False
        if chi == 0:
            return ('flat', None) if e == 0 else ('Nil', False)
        return ('H2xR', None) if e == 0 else ('SL2~', False)
    m = TXI.match(name)
    if m:
        a, b, c, d = map(int, m.groups())
        tr = abs(a + d)
        if tr > 2:
            return 'Sol', False
        if tr < 2 or (b, c) == (0, 0):
            return 'flat', None
        return 'Nil', False
    if name.startswith('Hyp_'):
        return 'hyperbolic', None
    return 'other', None


def load_census(path):
    """Yields (isosig, name, k) from closed-or-census-11.dat lines 'sig NAME : #k'."""
    for line in open(path):
        sig, _, rest = line.strip().partition(' ')
        name, _, k = rest.rpartition(' : #')
        yield sig, name, int(k)
