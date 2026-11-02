"""Random Pachner walks among triangulations whose edges have degree 5 or 6.

    python validation/controls/five_six_walks.py STARTS.jsonl OUT.jsonl
        [--steps N] [--walks-per-start K] [--workers W] [--star]

Simulated annealing over 2-3, 3-2 and 4-4 moves (vertex count is preserved).
The cost is the total distance of the edge degrees from {5, 6}, plus (with
--star) the number of triangles with two edges of degree 5 (Thurston's 5/6*
condition, Elder-McCammond-Meier 2003). Every distinct triangulation of cost
0 on degrees (a 5/6-triangulation) is recorded, with its isosig, the
Elder-McCammond-Meier degree metric (squared length 4 on degree-5 edges, 3 on
degree-6 edges) and whether it is 5/6*.
"""

import argparse
import json
import math
import multiprocessing
import random
import time

import regina


def degrees(T):
    return [e.degree() for e in T.edges()]


def bad_triangles(T):
    return sum(sum(f.edge(k).degree() == 5 for k in range(3)) >= 2 for f in T.triangles())


def cost(T, star):
    c = sum(0 if d in (5, 6) else (5 - d if d < 5 else d - 6) for d in degrees(T))
    if star and c == 0:
        c += bad_triangles(T)
    elif star:
        c += 0.25 * bad_triangles(T)
    return c


def neighbour(T, rng, max_tet):
    """A random 2-3, 3-2 or 4-4 move, or None."""
    for _ in range(200):
        r = rng.random()
        if r < 0.4 and T.size() < max_tet:
            f = T.triangle(rng.randrange(T.countTriangles()))
            U = T.withPachner(f)
        elif r < 0.8:
            e = [e for e in T.edges() if e.degree() == 3]
            if not e:
                continue
            U = T.withPachner(rng.choice(e))
        else:
            e = [e for e in T.edges() if e.degree() == 4]
            if not e:
                continue
            U = T.with44(rng.choice(e), rng.randrange(2))
        if U is not None:
            return U
    return None


def walk(job):
    start, seed, steps, star = job
    rng = random.Random(seed)
    T = regina.Triangulation3.fromIsoSig(start['isosig'])
    max_tet = max(2 * T.size(), T.size() + 30)
    c = cost(T, star)
    found = {}
    t_hi, t_lo = 4.0, 0.3
    for i in range(steps):
        temp = t_hi * (t_lo / t_hi) ** ((i % 20000) / 20000)   # repeated cooling cycles
        U = neighbour(T, rng, max_tet)
        if U is None:
            continue
        cu = cost(U, star)
        if cu <= c or rng.random() < math.exp((c - cu) / temp):
            T, c = U, cu
            deg = degrees(T)
            if set(deg) <= {5, 6}:
                sig = T.isoSig()
                if sig not in found:
                    found[sig] = dict(isosig=sig, from_isosig=start['isosig'], name=start.get('name'),
                                      census=start.get('census'), n_tet=T.size(),
                                      sqlens=[4 if d == 5 else 3 for d in deg],
                                      triangles_two_deg5=bad_triangles(T), seed=seed, step=i)
    return list(found.values())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('starts')
    ap.add_argument('out')
    ap.add_argument('--steps', type=int, default=200000)
    ap.add_argument('--walks-per-start', type=int, default=2)
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--star', action='store_true')
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    starts = [json.loads(l) for l in open(a.starts) if l.strip()]
    jobs = [(s, a.seed * 100000 + 1000 * i + k, a.steps, a.star)
            for i, s in enumerate(starts) for k in range(a.walks_per_start)]
    seen = set()
    t0 = time.time()
    with open(a.out, 'w') as f, multiprocessing.Pool(a.workers) as pool:
        for rows in pool.imap_unordered(walk, jobs):
            for r in rows:
                if r['isosig'] in seen:
                    continue
                seen.add(r['isosig'])
                f.write(json.dumps(r) + '\n')
            f.flush()
    print('%d distinct 5/6-triangulations, %d of them 5/6*, in %.0f s' % (
        len(seen), sum(1 for l in open(a.out) if json.loads(l)['triangles_two_deg5'] == 0),
        time.time() - t0))


if __name__ == '__main__':
    main()
