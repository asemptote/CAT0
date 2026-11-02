"""Negative controls with generic edge lengths on manifolds that cannot be CAT(0).

    python validation/controls/generic_controls.py --per-geometry 100 --workers 20

Why generic lengths: for a closed triangulation with n tetrahedra and v
vertices there are n + v edges (Euler characteristic 0), so the edge
condition needs a total dihedral angle of 2*pi*(n + v). Both cube and
disphenoid tetrahedra have dihedral sum 7*pi/3, which forces n >= 6v with
equality only for flat metrics: those shapes almost never fit the edge
condition, which is why the Pachner walks found nothing on non-CAT(0)
manifolds. (Regular tetrahedra are worse: average edge degree 6n/(n+v) < 6.)
A Euclidean tetrahedron's dihedral sum can approach 3*pi, so with free edge
lengths the edge condition becomes reachable.

For each census triangulation (>= --min-tets tetrahedra) of a manifold that
cannot be CAT(0) (topology.py), a stochastic hill climb on log squared edge
lengths maximises  min over edges (dihedral sum - 2*pi)  with all tetrahedra
non-degenerate. Once the margin exceeds --target, the lengths are rounded to
small integers (squared) and re-checked with verify_cert.metric_report. Any
metric found is a negative control: the vertex link condition must fail.
Output: validation/results/generic_controls.jsonl (expected=False).
"""

import argparse
import json
import math
import os
import random
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

import numpy as np
import regina

HERE = os.path.dirname(os.path.abspath(__file__))
VALIDATION = os.path.dirname(HERE)
sys.path[:0] = [HERE, VALIDATION]

import verify_cert as V  # noqa: E402
from topology import classify, load_census  # noqa: E402

CENSUS = os.path.join(os.path.dirname(VALIDATION), 'CAT0_Candidates', 'closed-or-census-11.dat')
TWO_PI = 2 * math.pi


class Margin:
    """min over edges of (dihedral angle sum - 2*pi) as a function of the
    squared edge lengths; -inf if some tetrahedron is degenerate."""

    def __init__(self, isosig):
        T = regina.Triangulation3.fromIsoSig(isosig)
        self.ne = T.countEdges()
        self.tets = []   # per tet: list of (local i, local j, global edge)
        for tet in T.tetrahedra():
            self.tets.append([(tet.edgeMapping(k)[0], tet.edgeMapping(k)[1], tet.edge(k).index())
                              for k in range(6)])

    def margins(self, sq):
        """Per-edge (dihedral sum - 2*pi), or None if a tetrahedron is degenerate."""
        sums = np.zeros(self.ne)
        for edges in self.tets:
            d2 = {}
            for i, j, e in edges:
                d2[(i, j)] = d2[(j, i)] = sq[e]
            p, vol6 = V.embed_tet(d2)
            if p is None or vol6 < 1e-9 * max(sq) ** 1.5:
                return None
            for i, j, e in edges:
                sums[e] += V.dihedral(p, i, j)
        return sums - TWO_PI

    def __call__(self, sq):
        m = self.margins(sq)
        return -math.inf if m is None else float(np.min(m))

    def soft(self, sq, beta=12.0):
        """Smooth lower bound of the minimum margin (log-sum-exp), used to
        guide the search; the true minimum decides success."""
        m = self.margins(sq)
        if m is None:
            return -math.inf, -math.inf
        lo = float(np.min(m))
        return lo - math.log(float(np.sum(np.exp(-beta * (m - lo))))) / beta, lo


def candidates(sig, rng, walk_steps, min_degree, max_tets, want):
    """Triangulations of the same manifold (Pachner random walk, as in
    pachner_controls.py) whose edges all have degree >= min_degree: an edge
    of degree <= 2 can never reach 2*pi since dihedral angles are < pi."""
    from pachner_controls import random_move
    T = regina.Triangulation3.fromIsoSig(sig)
    out, seen = [], set()
    for _ in range(walk_steps):
        random_move(T, rng, max_tets)
        if min(T.edge(e).degree() for e in range(T.countEdges())) < min_degree:
            continue
        iso = T.isoSig()
        if iso not in seen:
            seen.add(iso)
            out.append(iso)
            if len(out) >= want:
                break
    return out


def work(args):
    sig, name, geom, seed, steps, target, restarts, walk = args
    rng = random.Random(seed)
    res = []
    for iso in candidates(sig, rng, *walk):
        r = optimise((iso, name, geom, rng.randrange(10 ** 9), steps, target, restarts))
        r['from_isosig'] = sig
        res.append(r)
        if 'sqlens' in r:
            break
    return res or [dict(isosig=sig, name=name, geometry=geom, margin=-math.inf)]


def optimise(args):
    sig, name, geom, seed, steps, target, restarts = args
    rng = random.Random(seed)
    M = Margin(sig)
    best = None
    for _ in range(restarts):
        x = [rng.gauss(0, 0.3) for _ in range(M.ne)]
        s, f = M.soft([math.exp(v) for v in x])
        sigma = 0.5
        for step in range(steps):
            y = list(x)
            for e in rng.sample(range(M.ne), rng.choice((1, 1, 2, 3))):
                y[e] += rng.gauss(0, sigma)
            sy, fy = M.soft([math.exp(v) for v in y])
            if sy >= s:
                x, s, f = y, sy, fy
            sigma = max(0.02, sigma * 0.9995)
            if f > target:
                break
        if best is None or f > best[1]:
            best = (x, f)
        if f > target:
            break
    x, f = best
    out = dict(isosig=sig, name=name, geometry=geom, margin=f)
    if f <= 0:
        return out
    # round to small integer squared lengths that keep the edge condition
    base = min(math.exp(v) for v in x)
    for K in (2, 3, 4, 6, 9, 12, 16, 25, 36, 50, 100, 200, 500, 1000):
        sq = [max(1, round(math.exp(v) / base * K)) for v in x]
        rep = V.metric_report(sig, sq)
        if rep['ok'] and rep['edge_condition']:
            g = math.gcd(*sq)
            sq = [s // g for s in sq]
            out.update(sqlens=sq, expected=False, shape='generic',
                       int_margin=M([float(s) for s in sq]))
            break
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--per-geometry', type=int, default=100)
    ap.add_argument('--geometries', nargs='*', default=['spherical', 'S2xS1', 'Nil', 'Sol', 'SL2~'])
    ap.add_argument('--min-tets', type=int, default=7)
    ap.add_argument('--steps', type=int, default=6000)
    ap.add_argument('--restarts', type=int, default=3)
    ap.add_argument('--target', type=float, default=0.05)
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--min-degree', type=int, default=4)
    ap.add_argument('--walk-steps', type=int, default=20000)
    ap.add_argument('--max-tets', type=int, default=16)
    ap.add_argument('--per-manifold', type=int, default=3, help='triangulations optimised per manifold')
    ap.add_argument('--out', default=os.path.join(VALIDATION, 'results', 'generic_controls.jsonl'))
    a = ap.parse_args()

    by_geom = defaultdict(list)
    for sig, name, k in load_census(CENSUS):
        geom, expected = classify(name)
        if geom in a.geometries and k == 1:
            by_geom[geom].append((sig, name, geom))
    rng = random.Random(a.seed)
    walk = (a.walk_steps, a.min_degree, a.max_tets, a.per_manifold)
    jobs = []
    for geom, items in sorted(by_geom.items()):
        for sig, name, g in rng.sample(items, min(a.per_geometry, len(items))):
            jobs.append((sig, name, g, rng.randrange(10 ** 9), a.steps, a.target, a.restarts, walk))
    print('%d manifolds: %s' % (len(jobs), dict(Counter(j[2] for j in jobs))), flush=True)
    stats = Counter()
    with Pool(a.workers) as p, open(a.out, 'a') as f:
        for rs in p.imap_unordered(work, jobs):
            for r in rs:
                stats[(r['geometry'], 'triangulations optimised')] += 1
                if 'sqlens' in r:
                    stats[(r['geometry'], 'metrics')] += 1
                    f.write(json.dumps(r) + '\n')
                    f.flush()
                print(r['geometry'], r['name'], r['isosig'], 'margin %.4f' % r['margin'],
                      'FOUND %s' % r['sqlens'] if 'sqlens' in r else '', flush=True)
    for k in sorted(stats):
        print(k, stats[k])


if __name__ == '__main__':
    main()
