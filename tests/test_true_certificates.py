"""True certificates: the kernel accepts genuine ones and rejects wrong ones.

    python tests/test_true_certificates.py
"""

import copy
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [ROOT, os.path.join(ROOT, 'validation')]

from vlink import kernel  # noqa: E402
from vlink.certify import true_certificate  # noqa: E402

FLAT = ('gLLPQaceedffnkxknnn', [1, 2, 3, 1, 1, 2, 2])            # flat: CAT(0)
HYPERBOLIC_CANDIDATE = ('oLLvzLLQQQcbehlmmklnknmnklxxahaqhqxhaqaqq',  # not CAT(0)
                        [1, 2, 1, 3, 2, 1, 3, 1, 2, 1, 2, 2, 1, 2, 2])


def _ok(c, **kw):
    return kernel.check_true_certificate(c, **kw)[0]


def test_genuine_accepted():
    for pi in (True, False):
        c = true_certificate(*FLAT, pi_beads=pi)
        assert _ok(c), kernel.check_true_certificate(c)
        assert _ok(c, exact_only=True)


def test_wrong_rejected():
    good = true_certificate(*FLAT)
    v = good['vertices'][0]
    # a reached bead removed from a tree
    s, tree = next((s, t) for s, t in v['trees'].items() if t)
    c = copy.deepcopy(good)
    del c['vertices'][0]['trees'][s][next(iter(tree))]
    assert not _ok(c)
    # a predecessor that is not a passing junction / not in the tree
    c = copy.deepcopy(good)
    t = c['vertices'][0]['trees'][s]
    x = next(iter(t))
    t[x] = int(x)     # points to itself: a cycle
    assert not _ok(c)
    # wrong bead count, missing vertex, bad options
    c = copy.deepcopy(good); c['vertices'][0]['n_beads'] += 1; assert not _ok(c)
    c = copy.deepcopy(good); c['vertices'] = []; assert not _ok(c)
    c = copy.deepcopy(good); c['options']['pi_beads'] = False; assert not _ok(c)
    c = copy.deepcopy(good); c['options']['allow_repeats'] = False; assert not _ok(c)
    c = copy.deepcopy(good); c['format'] = 'x'; assert not _ok(c)
    # the same trees for a different metric (edge condition fails)
    c = copy.deepcopy(good); c['sqlens'] = [1] * 7; assert not _ok(c)
    # a metric that is not CAT(0): its shortest-path trees cannot certify it
    bad = true_certificate(*HYPERBOLIC_CANDIDATE)
    ok, msg = kernel.check_true_certificate(bad)
    assert not ok and 'shorter than 2*pi' in msg, msg
    return msg


if __name__ == '__main__':
    test_genuine_accepted()
    print('ok genuine True certificates accepted (both pi settings, exact_only too)')
    print('ok wrong ones rejected; not-CAT(0) metric:', test_wrong_rejected())
