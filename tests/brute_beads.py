"""Completeness check of vlink's bead enumeration by brute force.

    python tests/brute_beads.py [K]

For a few links, every strip of at most K link triangles (any gluing
sequence without immediate backtracking) and every choice of start / end
corner is tested with the kernel's exact bead test (kernel.check_bead). The
set of beads of length < pi found this way must equal the set vlink's pruned
enumeration (beads.all_beads) produces for strips of at most K triangles.
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [ROOT, os.path.join(ROOT, 'validation')]

import examples  # noqa: E402
from vlink import kernel  # noqa: E402
from vlink.beads import all_beads  # noqa: E402
from vlink.link import link_of  # noqa: E402


def brute(link, K):
    found = set()
    k = link.k
    # all strips (ts, steps) of length <= K without immediate backtracking
    stack = [((t,), ()) for t in range(k.n)]
    while stack:
        ts, steps = stack.pop()
        for s in range(3):
            for e in range(3):
                spec = dict(ts=list(ts), start=s, end=e, steps=list(steps), pi_point=None)
                try:
                    b = kernel.check_bead(k, spec)
                except kernel.Reject:
                    continue
                if not b.pi:
                    found.add((ts, s, e))
        if len(ts) < K:
            for d in range(3):
                if steps and k.adj[ts[-1]][d] == ts[-2] and \
                        k.glue[ts[-1]][d][d] == steps[-1]:
                    continue  # immediately crossing back over the same edge
                stack.append((ts + (k.adj[ts[-1]][d],), steps + (d,)))
    return found


def main(K=6):
    cases = [e for e in examples.load('examples_1') if e['idx'] in (0, 5)] + \
            [e for e in examples.load('examples_3_disphenoid') if e['idx'] in (8, 14)] + \
            [e for e in examples.load('hyperbolic_candidate')]
    ok = True
    for e in cases:
        link = link_of(e['isosig'], 0, e['sqlens'])
        mine = {b.key() for b in all_beads(link, allow_repeats=True)
                if not b.pi_range and len(b.ts) <= K}
        theirs = brute(link, K)
        missing, extra = theirs - mine, mine - theirs
        print('%s #%d: brute %d, vlink %d, missed by vlink %d, not confirmed %d' % (
            e['file'], e['idx'], len(theirs), len(mine), len(missing), len(extra)))
        for m in sorted(missing)[:5]:
            print('   missed:', m)
        for m in sorted(extra)[:5]:
            print('   unconfirmed:', m)
        ok = ok and not missing and not extra
    print('COMPLETE' if ok else 'MISMATCH')


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
