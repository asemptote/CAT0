"""Exact edge condition (vlink.edges) against the independent float version
(validation/verify_cert.metric_report) on every validation input and control.

    python tests/test_vlink_edges.py
"""

import glob
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [ROOT, os.path.join(ROOT, 'validation')]

import verify_cert as V  # noqa: E402
from vlink.edges import edge_report  # noqa: E402

DATA = json.load(open(os.path.join(ROOT, 'tests', 'data', 'regression_cases.json')))


def inputs():
    """Every metric with known data in the repository (no candidate files
    needed): the regression cases, all checked metrics, and the controls."""
    seen = set()
    for c in DATA['cases']:
        key = (c['isosig'], tuple(c['sqlens']))
        if key not in seen:
            seen.add(key)
            yield c['isosig'], c['sqlens']
    for path in [os.path.join(ROOT, 'validation', 'results', 'all_metrics2.jsonl')] + \
            glob.glob(os.path.join(ROOT, 'validation', 'results', '*controls*.jsonl')):
        for line in open(path):
            r = json.loads(line)
            key = (r['isosig'], tuple(r['sqlens']))
            if key not in seen:
                seen.add(key)
                yield r['isosig'], r['sqlens']


def test_edge_report_matches_float():
    n = 0
    for sig, sq in inputs():
        exact = edge_report(sig, sq)
        fl = V.metric_report(sig, sq)
        assert exact['valid'] == fl['ok'], (sig, sq)
        if exact['valid']:
            assert exact['edge_condition'] == fl['edge_condition'], (sig, sq)
            assert exact['flat'] == fl['flat'], (sig, sq)
        n += 1
    assert n >= 50
    return n


def test_examples_3_labels_are_degenerate():
    for c in DATA['invalid']:
        assert edge_report(c['isosig'], c['sqlens'])['valid'] is False


if __name__ == '__main__':
    print('ok, inputs checked:', test_edge_report_matches_float())
    test_examples_3_labels_are_degenerate()
    print('ok examples_3 labels degenerate')
