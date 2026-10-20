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

import examples  # noqa: E402
import verify_cert as V  # noqa: E402
from vlink.edges import edge_report  # noqa: E402


def inputs():
    seen = set()
    for e in examples.load_all(examples.FILES + ('examples_3_disphenoid',)):
        key = (e['isosig'], tuple(e['sqlens']))
        if key not in seen:
            seen.add(key)
            yield e['isosig'], e['sqlens']
    for path in glob.glob(os.path.join(ROOT, 'validation', 'results', '*controls*.jsonl')):
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
    for e in examples.load('examples_3'):
        assert edge_report(e['isosig'], e['sqlens'])['valid'] is False


if __name__ == '__main__':
    print('ok, inputs checked:', test_edge_report_matches_float())
    test_examples_3_labels_are_degenerate()
    print('ok examples_3 labels degenerate')
