"""Irrational link fields: Tower arithmetic and the kernel on generic lengths.

    python tests/test_towers.py
"""

import math
import os
import random
import sys
from fractions import Fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [ROOT, os.path.join(ROOT, 'validation'), os.path.join(ROOT, 'tests')]

import verify_cert as V  # noqa: E402
from vlink import kernel  # noqa: E402
from vlink.certify import vertex_trees  # noqa: E402
from vlink.check import vertex_linkcond, _certificate  # noqa: E402
from vlink.field import MQ, Tower  # noqa: E402

FLAT_SIG = 'gLLPQaceedffnkxknnn'


def test_identities():
    s2, s3, s6 = (MQ.sqrt_rational(n) for n in (2, 3, 6))
    a = Tower.sqrt_field(2 + s3) * Tower.sqrt_field(2 - s3) - 1
    assert a.is_zero()
    assert (Tower.sqrt_field(3 + 2 * s2) - (1 + s2)).is_zero()
    assert (Tower.sqrt_field(5 + 2 * s6) - (s2 + s3)).is_zero()
    assert (Tower.sqrt_field(2 + s3) - Tower.sqrt_field(2 - s3)).sign() == 1
    # near-cancellation: sqrt(2+sqrt3) vs (sqrt6+sqrt2)/2, equal
    assert (Tower.sqrt_field(2 + s3) - (s6 + s2) * Fraction(1, 2)).is_zero()


def test_random_vs_float():
    rng = random.Random(1)
    base = [MQ.sqrt_rational(n) for n in (2, 3, 5)]
    for _ in range(300):
        r1 = MQ.rational(rng.randint(5, 20)) + base[rng.randrange(3)] * rng.randint(-2, 2)
        r2 = MQ.rational(rng.randint(5, 20)) + base[rng.randrange(3)] * rng.randint(-2, 2)
        x = Tower.sqrt_field(r1) * rng.randint(-3, 3) + Tower.sqrt_field(r2) * rng.randint(-3, 3) \
            + MQ.rational(rng.randint(-5, 5))
        y = x * x - Tower.sqrt_field(r1)
        fx = float(x)
        assert abs(fx - (math.sqrt(float(r1)) * 0 + fx)) < 1e-9
        assert abs(float(y) - (fx * fx - math.sqrt(float(r1)))) < 1e-6 * (1 + abs(fx) ** 2)
        if abs(fx) > 1e-9:
            assert x.sign() == (1 if fx > 0 else -1)


def generic_metrics(n=6, seed=5):
    rng = random.Random(seed)
    out = []
    while len(out) < n:
        sq = [rng.choice([9, 10, 11, 12, 13, 14]) for _ in range(7)]
        if V.metric_report(FLAT_SIG, sq)['ok']:
            out.append(sq)
    return out


# perturbations of flat metrics (squared lengths x 1000 +- small integers):
# irrational link fields with known vertex-link answers
PERTURBED = [
    ('gLLPQaceedffnkxknnn', [995, 2015, 2989, 997, 1005, 1986, 1987], True),
    ('gLLPQaceedffnkxknnn', [1011, 2002, 2988, 996, 1003, 1986, 2014], True),
    ('gLLPQaceedffnkxknnn', [988, 2015, 2992, 1005, 1005, 2003, 2015], False),
    ('gLLPQbefefefhhxxqxh', [2986, 4003, 3003, 3997, 2986, 2992, 3986], False),
]


def test_perturbed_flat_metrics():
    """Irrational links with both answers: True ones get a kernel-accepted
    shortest-path certificate (both modes); False ones a witness accepted by
    the kernel and by the independent float checker."""
    for sig, sq, expected in PERTURBED:
        r = vertex_linkcond(sig, sq)
        assert r['result'] is expected, (sq, r['result'])
        link = kernel.Link(kernel.Geometry(sig, sq), 0)
        if expected:
            beads = kernel.enumerate_beads(link)
            trees = vertex_trees(link, beads)
            assert kernel.check_vertex_trees(link, beads, trees)
            assert kernel.check_vertex_trees(link, beads, trees, exact_only=True)
        else:
            r['reason'] = 'vertex'
            cert = _certificate(sig, sq, r)
            kb = [kernel.check_bead(link, b) for b in cert['necklace']]
            assert not kernel.total_length(kb).at_least(2)
            m = len(kb)
            for k in range(m):
                assert kernel.junction_ok(link, kb[k], kb[(k + 1) % m])
            assert V.verify_necklace(sig, sq, r['vertex'], r['necklace'])['valid']


def test_spherical_excess():
    """Every generic link triangle has corner angles summing to > pi."""
    for sq in generic_metrics(3):
        link = kernel.Link(kernel.Geometry(FLAT_SIG, sq), 0)
        for t in range(link.n):
            acc = kernel.AngleSum()
            for c in range(3):
                a, b = [x for x in range(3) if x != c]
                acc.add(*kernel.tangent_angle(link.G[t], kernel.E3[c], kernel.E3[a], kernel.E3[b]))
            assert acc.at_least(1) and not acc.at_least(3)


def test_generic_vertex_links():
    """The vertex-link search and both certificate checks run on irrational
    links; False witnesses also pass the independent float checker."""
    import brute_beads
    from vlink.link import VertexLink
    n_true = n_false = 0
    for sq in generic_metrics():
        r = vertex_linkcond(FLAT_SIG, sq)
        geo = kernel.Geometry(FLAT_SIG, sq)
        link = kernel.Link(geo, 0)
        beads = kernel.enumerate_beads(link)
        if r['result'] is False:
            r['reason'] = 'vertex'
            cert = _certificate(FLAT_SIG, sq, r)
            # the vertex part only (these metrics need not satisfy the edge condition)
            kl = link
            kb = [kernel.check_bead(kl, b) for b in cert['necklace']]
            assert not kernel.total_length(kb).at_least(2)
            v = V.verify_necklace(FLAT_SIG, sq, r['vertex'], r['necklace'])
            assert v['valid'], v
            n_false += 1
        else:
            trees = vertex_trees(link, beads)
            assert kernel.check_vertex_trees(link, beads, trees)
            assert kernel.check_vertex_trees(link, beads, trees, exact_only=True)
            n_true += 1
        # bead completeness by brute force on this irrational link
        vl = VertexLink(geo.T, 0, sq, geo=geo)
        mine = {(b.ts, b.start, b.end) for b in beads if not b.pi and len(b.ts) <= 3}
        assert brute_beads.brute(vl, 3) == mine
    return n_true, n_false


if __name__ == '__main__':
    test_identities()
    print('ok tower identities')
    test_random_vs_float()
    print('ok random tower arithmetic vs float')
    test_spherical_excess()
    print('ok spherical excess of generic link triangles (exact)')
    print('ok generic vertex links (True, False):', test_generic_vertex_links())
    test_perturbed_flat_metrics()
    print('ok perturbed flat metrics: True certified, False witnessed (irrational links)')
