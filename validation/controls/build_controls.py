"""Build labelled control inputs from the closed census.

    python validation/controls/build_controls.py --workers 5 [--geometries ...]

For every census triangulation whose manifold is classified by topology.py,
enumerate all cube / disphenoid metrics satisfying the edge condition
(shape_metrics.py) and write one JSON line per metric to
validation/results/controls.jsonl with:

  expected = False   the manifold cannot be CAT(0): the checker MUST say False
  expected = True    the metric is flat: the checker MUST say True
  expected = None    no ground truth (kept for statistics only)

Each metric is cross-checked with verify_cert.metric_report (independent float
geometry): valid tetrahedra, edge condition, flatness.
"""

import argparse
import json
import os
import sys
from collections import Counter
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
VALIDATION = os.path.dirname(HERE)
sys.path[:0] = [HERE, VALIDATION]

import verify_cert as V  # noqa: E402
from shape_metrics import enumerate_metrics  # noqa: E402
from topology import classify, load_census  # noqa: E402

CENSUS = os.path.join(os.path.dirname(VALIDATION), 'CAT0_Candidates', 'closed-or-census-11.dat')


def work(args):
    sig, name, k, geom, expected, limit, timeout = args
    out = []
    for shape in ('cube', 'disphenoid'):
        metrics, complete = enumerate_metrics(sig, shape, limit=limit, timeout=timeout)
        for sq in metrics:
            rep = V.metric_report(sig, sq)
            assert rep['ok'] and rep['edge_condition'], (sig, shape, sq, rep)
            exp = expected
            if rep['flat']:
                assert expected is not False, ('flat metric on non-CAT(0) manifold', sig, name)
                exp = True
            out.append(dict(isosig=sig, name=name, census_k=k, geometry=geom, shape=shape,
                            sqlens=sq, flat=rep['flat'], expected=exp,
                            search_complete=complete))
        if not metrics:
            out.append(dict(isosig=sig, name=name, census_k=k, geometry=geom, shape=shape,
                            sqlens=None, search_complete=complete))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--limit', type=int, default=50, help='max metrics per triangulation and shape')
    ap.add_argument('--timeout', type=float, default=20.0, help='seconds per triangulation and shape')
    ap.add_argument('--geometries', nargs='*',
                    default=['spherical', 'S2xS1', 'Nil', 'Sol', 'SL2~', 'flat'])
    ap.add_argument('--out', default=os.path.join(VALIDATION, 'results', 'controls.jsonl'))
    a = ap.parse_args()

    jobs = []
    for sig, name, k in load_census(CENSUS):
        geom, expected = classify(name)
        if geom in a.geometries:
            jobs.append((sig, name, k, geom, expected, a.limit, a.timeout))
    print('%d triangulations: %s' % (len(jobs), dict(Counter(j[3] for j in jobs))), flush=True)

    stats = Counter()
    with Pool(a.workers) as p, open(a.out, 'w') as f:
        for i, rows in enumerate(p.imap_unordered(work, jobs, chunksize=8)):
            for r in rows:
                if r['sqlens'] is None:
                    stats[(r['geometry'], r['shape'], 'none' if r['search_complete'] else 'cut off')] += 1
                    continue
                f.write(json.dumps(r) + '\n')
                stats[(r['geometry'], r['shape'], 'expected=%s' % r['expected'])] += 1
            if (i + 1) % 2000 == 0:
                print(i + 1, 'done', flush=True)
    for key in sorted(stats):
        print(key, stats[key])


if __name__ == '__main__':
    main()
