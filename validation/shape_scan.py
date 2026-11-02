"""Scan census triangulations for metrics built from one tetrahedron shape (C10).

    python validation/shape_scan.py CANDIDATES.json OUT.jsonl [--workers N] [--shapes S,...]

CANDIDATES.json: [{"isosig", "census", "name"}, ...], the triangulations whose
edges all have degree >= 4 (no shape here has a dihedral angle above pi/2).
Default shapes: the B~3 Coxeter tetrahedron and the orthoschemes orth-a-b-c
with a, b, c in 1..4, gcd 1, up to reversal, other than the cube (orth-1-1-1;
the cube and the disphenoid were scanned before, census_full.md).

Writes one line per (shape, triangulation) with a metric or an incomplete
search: {"shape", "isosig", "census", "name", "complete", "metrics"}, and a
final summary line per shape {"shape", "summary": {...}}.
"""

import argparse
import itertools
import json
import math
import multiprocessing
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.join(HERE, 'controls'), HERE]


def default_shapes():
    out = ['b3']
    for a, b, c in itertools.product(range(1, 5), repeat=3):
        if (a, b, c) > (c, b, a) or math.gcd(a, math.gcd(b, c)) != 1 or (a, b, c) == (1, 1, 1):
            continue
        out.append('orth-%d-%d-%d' % (a, b, c))
    return out


def _work(job):
    from shape_metrics import enumerate_metrics
    shape, rows = job
    out = []
    for r in rows:
        t = time.time()
        metrics, complete = enumerate_metrics(r['isosig'], shape, limit=100000, timeout=600.0)
        if metrics or not complete:
            out.append(dict(shape=shape, isosig=r['isosig'], census=r['census'], name=r['name'],
                            complete=complete, n=len(metrics), metrics=metrics[:1000],
                            seconds=round(time.time() - t, 2)))
    return shape, len(rows), out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('candidates')
    ap.add_argument('out')
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--shapes', default=None)
    ap.add_argument('--chunk', type=int, default=500)
    a = ap.parse_args()
    cands = json.load(open(a.candidates))
    shapes = a.shapes.split(',') if a.shapes else default_shapes()
    jobs = [(s, cands[i:i + a.chunk]) for s in shapes for i in range(0, len(cands), a.chunk)]
    summary = {s: dict(scanned=0, triangulations=0, metrics=0, incomplete=0) for s in shapes}
    t0 = time.time()
    with open(a.out, 'w') as f, multiprocessing.Pool(a.workers) as pool:
        for shape, n, rows in pool.imap_unordered(_work, jobs):
            sm = summary[shape]
            sm['scanned'] += n
            for r in rows:
                f.write(json.dumps(r) + '\n')
                sm['triangulations'] += bool(r['n'])
                sm['metrics'] += r['n']
                sm['incomplete'] += not r['complete']
            f.flush()
        for s in shapes:
            f.write(json.dumps(dict(shape=s, summary=summary[s])) + '\n')
    for s in shapes:
        print(s, summary[s])
    print('done in %.0f s' % (time.time() - t0))


if __name__ == '__main__':
    main()
