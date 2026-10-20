"""Regression tests for vlink.check.linkcond on inputs with known answers.

    python tests/test_regression.py          (all 39 cases, a few minutes)
    QUICK=1 python tests/test_regression.py  (one case of each kind)

tests/data/regression_cases.json holds the inputs (no other data needed):
  - flat metrics: locally Euclidean, so the answer must be True;
  - 'witness': not CAT(0); vlink must say False with a necklace that the
    independent checker validation/verify_cert.py accepts;
  - 'agreement of gallery.py and vlink': non-flat metrics where the old and
    new implementations agree (no independent ground truth);
  - 'invalid': edge labels mistaken for squared lengths (examples_3.txt):
    every tetrahedron is degenerate.
"""

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [ROOT, os.path.join(ROOT, 'validation')]

import verify_cert as V  # noqa: E402
from vlink.check import linkcond  # noqa: E402
from vlink.edges import edge_report  # noqa: E402

DATA = json.load(open(os.path.join(ROOT, 'tests', 'data', 'regression_cases.json')))


def _cases():
    cases = DATA['cases']
    if os.environ.get('QUICK'):
        picked = {}
        for c in cases:
            picked.setdefault(c['ground_truth'], c)
        cases = list(picked.values())
    return cases


def test_linkcond_known_answers():
    for c in _cases():
        r = linkcond(c['isosig'], c['sqlens'])
        assert r['result'] == c['expected'], (c['source'], r)
        if c['flat']:
            assert r['result'] is True and r.get('flat') is True, c['source']
        if r['result'] is False:
            assert r['reason'] == 'vertex', (c['source'], r)
            v = V.verify_necklace(c['isosig'], c['sqlens'], r['vertex'], r['necklace'])
            assert v['valid'], (c['source'], 'witness rejected by verify_cert')


def test_invalid_inputs_rejected():
    for c in DATA['invalid']:
        assert edge_report(c['isosig'], c['sqlens'])['valid'] is False, c['source']
        assert linkcond(c['isosig'], c['sqlens'])['reason'] == 'invalid', c['source']


if __name__ == '__main__':
    test_invalid_inputs_rejected()
    print('ok invalid inputs')
    test_linkcond_known_answers()
    print('ok %d known-answer cases' % len(_cases()))
