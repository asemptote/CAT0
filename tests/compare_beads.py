"""Compare bead enumeration: gallery.py (old) vs vlink.beads (new).

    python tests/compare_beads.py ISOSIG SQLENS_JSON [--vertex 0] [--skip-old]

For every link triangle, both implementations enumerate beads; they must give
the same set of (triangles, start corner, end corner) with the same lengths.
Prints timings and any differences.
"""

import argparse
import json
import math
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from vlink.beads import beads_from  # noqa: E402
from vlink.link import link_of  # noqa: E402


def old_beads(isosig, sqlens, i_v):
    import gallery
    gallery.DEBUG = False
    T = gallery.fromIsoSig(isosig)
    v = T.vertex(i_v)
    lt = gallery.LinkTriangles(T, v, sqlens)
    out = {}
    for t in range(v.buildLink().size()):
        for g in gallery.build_galleries(lt, t):
            key = (tuple(g.ts), g.ds[0], g.rs[-1])
            out[key] = float(g.geodesiclength)
    return out


def new_beads(isosig, sqlens, i_v):
    link = link_of(isosig, i_v, sqlens)
    out = {}
    for t in range(link.n):
        for b in beads_from(link, t):
            out[b.key()] = b.length
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('isosig')
    ap.add_argument('sqlens')
    ap.add_argument('--vertex', type=int, default=0)
    ap.add_argument('--skip-old', action='store_true')
    a = ap.parse_args()
    sq = json.loads(a.sqlens)
    t0 = time.time()
    new = new_beads(a.isosig, sq, a.vertex)
    t_new = time.time() - t0
    print('new: %d beads in %.2f s (%d thick)' % (len(new), t_new, sum(len(k[0]) > 1 for k in new)))
    if a.skip_old:
        return
    t0 = time.time()
    old = old_beads(a.isosig, sq, a.vertex)
    t_old = time.time() - t0
    print('old: %d beads in %.2f s (%d thick)' % (len(old), t_old, sum(len(k[0]) > 1 for k in old)))
    only_old = sorted(set(old) - set(new))
    only_new = sorted(set(new) - set(old))
    bad_len = [(k, old[k], new[k]) for k in set(old) & set(new) if abs(old[k] - new[k]) > 1e-9]
    print('only in old: %d, only in new: %d, length mismatches: %d' % (
        len(only_old), len(only_new), len(bad_len)))
    for k in only_old[:10]:
        print('  old only', k, 'len/pi=%.6f' % (old[k] / math.pi))
    for k in only_new[:10]:
        print('  new only', k, 'len/pi=%.6f' % (new[k] / math.pi))
    for k, lo, ln in bad_len[:10]:
        print('  length', k, lo / math.pi, ln / math.pi)
    print('speedup x%.0f' % (t_old / max(t_new, 1e-9)))
    print('IDENTICAL' if not (only_old or only_new or bad_len) else 'DIFFERENT')


if __name__ == '__main__':
    main()
