"""Certificates for True answers (search side): shortest-path trees.

For each vertex link and each bead s, compute shortest bead-chain distances
from s over the beads numbered above s, along passing junctions, up to
length 2*pi; record them as predecessor trees. kernel.check_true_certificate
verifies the trees (docs/true_certificates.md); nothing here is trusted.

Distances are compared as floats when they differ by more than 1e-7 and
exactly (kernel.Lengths.leq) otherwise, so the trees are true shortest-path
trees and the kernel's inequalities hold.
"""

import heapq
import math

from . import kernel

TWO_PI = 2 * math.pi
TIE = 1e-7


def vertex_trees(link, beads, exact_only=False):
    """{s: {x: predecessor}} for one vertex link."""
    L = kernel.Lengths(beads, exact_only)
    J = kernel.junction_table(link, beads)
    succ = {}
    for (x, y), ok in J.items():
        if ok:
            succ.setdefault(x, []).append(y)
    mid = [(lo + hi) / 2 for lo, hi in (kernel.length_enclosure(b) for b in beads)]
    trees = {}
    for s in range(len(beads)):
        dist, chain, pred = {s: mid[s]}, {s: (s,)}, {}
        heap = [(mid[s], s)]
        while heap:
            d, x = heapq.heappop(heap)
            if d > dist[x] + 1e-12:
                continue
            for y in succ.get(x, ()):
                if y <= s:
                    continue
                nd, nc = d + mid[y], chain[x] + (y,)
                if nd > TWO_PI + TIE:
                    continue
                old = dist.get(y)
                better = (old is None or nd < old - TIE or
                          (abs(nd - old) <= TIE and L.leq(nc, chain[y]) and not L.leq(chain[y], nc)))
                if better:
                    dist[y], chain[y], pred[y] = nd, nc, x
                    heapq.heappush(heap, (nd, y))
        trees[s] = pred
    return trees


def true_certificate(isosig, sqlens, pi_beads=True, exact_only=False):
    """A 'vlink-certificate-true/2' for a metric vlink answered True."""
    geo = kernel.Geometry(isosig, sqlens)
    vertices = []
    for i_v in range(geo.T.n_vertices):
        link = kernel.Link(geo, i_v)
        beads = kernel.enumerate_beads(link, allow_repeats=True, pi_beads=pi_beads)
        trees = vertex_trees(link, beads, exact_only)
        vertices.append(dict(vertex=i_v, n_beads=len(beads),
                             trees={str(s): {str(x): p for x, p in t.items()}
                                    for s, t in trees.items() if t}))
    cert = dict(format=kernel.TRUE_FORMAT, isosig=isosig, sqlens=list(sqlens),
                options=dict(allow_repeats=True, pi_beads=pi_beads), vertices=vertices)
    if not pi_beads:
        cert['assumptions'] = [kernel.LEMMA_4_9]
    return cert
