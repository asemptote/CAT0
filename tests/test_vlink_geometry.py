"""Tests for vlink.field and vlink.link against independent float geometry.

    python -m pytest tests/test_vlink_geometry.py      (or: python tests/test_vlink_geometry.py)
"""

import math
import os
import random
import sys
from fractions import Fraction

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [ROOT, os.path.join(ROOT, 'validation')]

import verify_cert as V  # noqa: E402
from vlink.field import MQ, squarefree_split  # noqa: E402
from vlink.link import link_of  # noqa: E402

FLAT_CUBE = ('gLLPQaceedffnkxknnn', [1, 2, 3, 1, 1, 2, 2])
HYPERBOLIC = ('oLLvzLLQQQcbehlmmklnknmnklxxahaqhqxhaqaqq',
              [1, 2, 1, 3, 2, 1, 3, 1, 2, 1, 2, 2, 1, 2, 2])
DISPHENOID = ('lLLPwvQQcaceedhikjkjknkxknnkakanx', None)  # lengths from examples_3


def _random_mq(rng, primes=(1, 2, 3, 5, 6, 10, 15, 30)):
    return MQ({s: Fraction(rng.randint(-20, 20), rng.randint(1, 9))
               for s in rng.sample(primes, rng.randint(1, 4))
               if rng.random() < 0.9}) + 0


def test_field_matches_float():
    rng = random.Random(0)
    for _ in range(2000):
        x, y = _random_mq(rng), _random_mq(rng)
        assert abs(float(x + y) - (float(x) + float(y))) < 1e-9
        assert abs(float(x * y) - float(x) * float(y)) < 1e-6 * (1 + abs(float(x) * float(y)))
        if not x.is_zero():
            assert abs(float(x * x.inverse()) - 1) < 1e-9
            assert x * x.inverse() == MQ.rational(1)


def test_field_exact_sign():
    s2, s3, s6 = (MQ.sqrt_rational(n) for n in (2, 3, 6))
    assert (s2 * s3 - s6).sign() == 0
    assert (s2 * s2 - 2).sign() == 0
    # near-cancellation that floats cannot decide: (sqrt2 + sqrt3)^2 vs 5 + 2 sqrt6
    assert ((s2 + s3) * (s2 + s3) - (5 + 2 * s6)).is_zero()
    # convergents of sqrt 2: 1393/985 is below it (1393^2 < 2*985^2), and
    # 665857/470832 is above it by ~1.6e-12, beyond the float estimate's
    # tolerance, so the exact path decides it
    assert (s2 - Fraction(1393, 985)).sign() == 1
    assert (s2 - Fraction(665857, 470832)).sign() == -1
    assert (s2 - Fraction(665857, 470832))._exact_sign() == -1
    assert squarefree_split(72) == (6, 2)
    assert MQ.sqrt_rational(Fraction(8, 3)) * MQ.sqrt_rational(Fraction(8, 3)) == MQ.rational(Fraction(8, 3))


def _check_link(isosig, sqlens, n_strips=60, max_len=10, seed=0):
    rng = random.Random(seed)
    T_fl = V.Link(isosig, sqlens, 0)
    L = link_of(isosig, 0, sqlens)
    keys = set()
    # 1. triangle sides and corner geometry
    for t in range(L.n):
        s = L.start(t)
        for a, b in ((0, 1), (0, 2), (1, 2)):
            c_exact = s.cos_float(s.pos[0][a], s.pos[0][b])
            assert abs(c_exact - math.cos(T_fl.side(t, a, b))) < 1e-12
    # 2. random strips: exact development vs float development
    for _ in range(n_strips):
        t0 = rng.randrange(L.n)
        s = L.start(t0)
        choices = []
        for _ in range(rng.randint(1, max_len)):
            d = rng.randrange(3)
            choices.append(d)
            s = s.extend(d)
        pos_fl = V.develop(T_fl, list(s.ts), choices)
        E = _frame_matrix(T_fl, isosig, sqlens, t0)
        for P, Q in zip(s.pos, pos_fl):
            for j in range(3):
                x = E @ np.array([float(c) for c in P[j]])
                x /= np.linalg.norm(x)
                assert np.allclose(x, Q[j], atol=1e-9), (s.ts, j, x, Q[j])
                for c in P[j]:
                    keys.update(c.t)
    return keys


def _frame_matrix(T_fl, isosig, sqlens, t):
    """Columns: the (unnormalised) edge vectors of link triangle t in the float
    embedding used by verify_cert.Link."""
    import regina
    T = regina.Triangulation3.fromIsoSig(isosig)
    v = T.vertex(0)
    I = v.buildLinkInclusion()
    tet = T.tetrahedron(I.tetImage(t))
    fp = I.facePerm(t)
    p, _ = V.embed_tet(V.tet_sqlens(T, sqlens, tet))
    return np.array([p[fp[j]] - p[fp[3]] for j in range(3)]).T


def test_flat_cube_link():
    keys = _check_link(*FLAT_CUBE)
    assert keys <= {1}, keys  # cube shapes: everything rational


def test_hyperbolic_candidate_link():
    _check_link(*HYPERBOLIC)


def test_disphenoid_link():
    import json
    data = json.load(open(os.path.join(ROOT, 'tests', 'data', 'regression_cases.json')))
    sq = next(c['sqlens'] for c in data['cases'] if c['isosig'] == DISPHENOID[0])
    _check_link(DISPHENOID[0], sq)


def test_generic_lengths_link():
    """Random near-regular integer squared lengths: mu != 1, a real field."""
    rng = random.Random(3)
    isosig = FLAT_CUBE[0]
    for trial in range(5):
        sq = [rng.choice([9, 10, 11, 12, 13]) for _ in range(7)]
        if not V.metric_report(isosig, sq)['ok']:
            continue
        keys = _check_link(isosig, sq, n_strips=30, max_len=6, seed=trial)
        assert keys - {1}, 'expected irrational coordinates for generic lengths'


if __name__ == '__main__':
    for name, f in sorted(globals().items()):
        if name.startswith('test_'):
            f()
            print('ok', name)
