"""Pi-length beads whose junction lies inside a single triangle corner.

    python tests/test_pi_beads.py

The search treats such a junction as infeasible (the angle on the side
between the two directions is at most the corner angle, < pi). Check that
it no longer raises, and that the kernel agrees for a concrete geodesic of
the lune.
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from vlink import kernel  # noqa: E402
from vlink.beads import all_beads  # noqa: E402
from vlink.link import link_of  # noqa: E402
from vlink.necklace import _pi_feasible  # noqa: E402

FLAT = ('gLLPQaceedffnkxknnn', [1, 2, 3, 1, 1, 2, 2])


def test_same_corner_pi_junction():
    link = link_of(FLAT[0], 0, FLAT[1])
    beads = all_beads(link, allow_repeats=True)
    pis = [b for b in beads if b.pi_range]
    thin = [b for b in beads if len(b.ts) == 1]
    n = 0
    for p in pis:
        prev = next((b for b in thin if b.ts[-1] == p.ts[0] and b.end == p.start), None)
        if prev is None:
            continue
        nxt = next(b for b in beads if not b.pi_range and
                   link.lv[b.ts[0]][b.start] == link.lv[p.ts[-1]][p.end])
        assert _pi_feasible(link, prev, p, nxt) is None
        # the kernel, for the geodesic through the middle of the lune
        L, R = p.pi_range
        Q = tuple(x + y for x, y in zip(L, R))
        kp = kernel.Bead(p.ts, p.start, p.end, p.G, p.pos, Q=Q, steps=p.strip.steps)
        sides = [kernel.junction_side(link.k, prev, kp, k).at_least(1) for k in (1, 2)]
        assert not all(sides), 'kernel finds both sides >= pi'
        n += 1
    assert n > 0, 'no same-corner configuration found'
    return n


if __name__ == '__main__':
    print('ok same-corner pi junctions infeasible (search and kernel agree):',
          test_same_corner_pi_junction())
