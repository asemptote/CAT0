"""Find shape metrics on larger triangulations of manifolds that cannot be CAT(0).

    python validation/controls/pachner_controls.py --per-geometry 40 --workers 5

Minimal census triangulations have edges of too low degree for any cube /
disphenoid metric to satisfy the edge condition. As in the original candidate
searches, random walks of Pachner moves (2-3, 3-2, 4-4) grow a triangulation
of the same manifold; every new triangulation with minimum edge degree >= 4
is searched exhaustively (shape_metrics.py) for metrics satisfying the edge
condition. Metrics found on manifolds classified by topology.py as unable to
be CAT(0) are negative controls (expected=False); flat metrics are positive
controls (expected=True). Results go to validation/results/pachner_controls.jsonl.
"""

import argparse
import json
import os
import random
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

import regina

HERE = os.path.dirname(os.path.abspath(__file__))
VALIDATION = os.path.dirname(HERE)
sys.path[:0] = [HERE, VALIDATION]

import verify_cert as V  # noqa: E402
from shape_metrics import enumerate_metrics  # noqa: E402
from topology import classify, load_census  # noqa: E402

CENSUS = os.path.join(os.path.dirname(VALIDATION), 'CAT0_Candidates', 'closed-or-census-11.dat')


def random_move(T, rng, max_tets, p_up=0.25):
    """One random Pachner move, biased like the 2016 search scripts: grow with
    2-3 moves, otherwise prefer 3-2 / 4-4 moves about low-degree edges."""
    if T.size() < max_tets and rng.random() < p_up:
        return T.pachner(T.triangle(rng.randrange(T.countTriangles())))
    deg3 = [e for e in range(T.countEdges()) if T.edge(e).degree() == 3]
    deg4 = [e for e in range(T.countEdges()) if T.edge(e).degree() == 4]
    if deg3 and rng.random() < 0.55:
        return T.pachner(T.edge(rng.choice(deg3)))
    if deg4:
        return T.move44(T.edge(rng.choice(deg4)), rng.randrange(2))
    if T.size() < max_tets:
        return T.pachner(T.triangle(rng.randrange(T.countTriangles())))
    return False


def walk(args):
    sig, name, geom, expected, seed, steps, max_tets, want, timeout = args
    rng = random.Random(seed)
    T = regina.Triangulation3.fromIsoSig(sig)
    seen, out, tried = set(), [], 0
    for _ in range(steps):
        random_move(T, rng, max_tets)
        if min(T.edge(e).degree() for e in range(T.countEdges())) < 4:
            continue
        iso = T.isoSig()
        if iso in seen:
            continue
        seen.add(iso)
        tried += 1
        for shape in ('cube', 'disphenoid'):
            metrics, _ = enumerate_metrics(iso, shape, limit=2, timeout=timeout)
            for sq in metrics:
                rep = V.metric_report(iso, sq)
                assert rep['ok'] and rep['edge_condition']
                exp = expected
                if rep['flat']:
                    assert expected is not False, ('flat metric on non-CAT(0) manifold', iso, name)
                    exp = True
                out.append(dict(isosig=iso, name=name, from_isosig=sig, geometry=geom,
                                shape=shape, sqlens=sq, flat=rep['flat'],
                                expected=exp, n_tet=T.size(), seed=seed))
        if len(out) >= want:
            break
    return dict(sig=sig, name=name, geometry=geom, tried=tried, found=out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--per-geometry', type=int, default=40, help='manifolds sampled per geometry')
    ap.add_argument('--geometries', nargs='*', default=['spherical', 'S2xS1', 'Nil', 'Sol', 'SL2~', 'flat'])
    ap.add_argument('--steps', type=int, default=3000)
    ap.add_argument('--max-tets', type=int, default=16)
    ap.add_argument('--want', type=int, default=4, help='stop a walk after this many metrics')
    ap.add_argument('--timeout', type=float, default=2.0)
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--census', choices=('closed', 'hw'), default='closed',
                    help="'closed': closed-or-census-11.dat (classified by topology.py); "
                         "'hw': hodgson-weeks.dat (closed hyperbolic, geometry 'hyperbolic')")
    ap.add_argument('--out', default=os.path.join(VALIDATION, 'results', 'pachner_controls.jsonl'))
    a = ap.parse_args()

    # one triangulation (the census's #1, i.e. minimal) per manifold
    by_geom = defaultdict(list)
    if a.census == 'hw':
        hw = os.path.join(os.path.dirname(CENSUS), 'hodgson-weeks.dat')
        for line in open(hw):
            vol_k, _, sig = line.strip().partition(' ')
            vol, _, k = vol_k.partition('_')
            by_geom['hyperbolic'].append((sig, 'HW_%s #%s' % (vol, k), 'hyperbolic', None))
    else:
        for sig, name, k in load_census(CENSUS):
            geom, expected = classify(name)
            if geom in a.geometries and k == 1:
                by_geom[geom].append((sig, name, geom, expected))
    rng = random.Random(a.seed)
    jobs = []
    for geom, items in sorted(by_geom.items()):
        for i, (sig, name, g, exp) in enumerate(rng.sample(items, min(a.per_geometry, len(items)))):
            jobs.append((sig, name, g, exp, a.seed * 100003 + len(jobs), a.steps, a.max_tets, a.want, a.timeout))
    print('%d walks: %s' % (len(jobs), dict(Counter(j[2] for j in jobs))), flush=True)

    stats = Counter()
    with Pool(a.workers) as p, open(a.out, 'a') as f:
        for res in p.imap_unordered(walk, jobs):
            for r in res['found']:
                f.write(json.dumps(r) + '\n')
            f.flush()
            stats[(res['geometry'], 'walks')] += 1
            stats[(res['geometry'], 'triangulations tried')] += res['tried']
            stats[(res['geometry'], 'metrics')] += len(res['found'])
            print(res['geometry'], res['name'], 'tried', res['tried'], 'found', len(res['found']), flush=True)
    for k in sorted(stats):
        print(k, stats[k])


if __name__ == '__main__':
    main()
