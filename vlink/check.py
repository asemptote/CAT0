"""Vertex link condition for a PE triangulated 3-manifold (vlink entry point).

    from vlink.check import vertex_linkcond
    vertex_linkcond('gLLPQaceedffnkxknnn', [1, 2, 3, 1, 1, 2, 2])
    -> {'result': True}

Same contract as gallery.vertex_linkcond: for each vertex, enumerate beads
(vlink.beads) and search for a necklace of length < 2*pi (vlink.necklace).
Returns {'result': True} or {'result': False, 'vertex': i, 'necklace': [...]}
with each bead as {'ts', 'start', 'end', 'length'} (the format checked by
validation/verify_cert.py). Raises necklace.UndecidedError instead of
guessing in the rare cases not yet decided exactly.
"""

import regina

from .beads import all_beads
from .link import VertexLink
from .necklace import search


def vertex_linkcond(isosig, sqlens, findall=False, allow_repeats=True):
    """allow_repeats (default) lets beads and necklaces reuse link triangles.
    gallery.py forbids reuse, an unproven assumption; on every validation
    input both settings give the same answers, and with the pruned search
    allowing repeats costs nothing extra, so the assumption-free setting is
    the default. allow_repeats=False reproduces gallery.py's search space."""
    T = regina.Triangulation3.fromIsoSig(isosig)
    allfound = []
    for i_v in range(T.countVertices()):
        link = VertexLink(T, i_v, sqlens)
        beads = all_beads(link, allow_repeats)
        res = search(link, beads, findall=findall, allow_repeats=allow_repeats)
        if findall:
            allfound += [(i_v, n) for n in res]
        elif res:
            return dict(result=False, vertex=i_v, necklace=_export(res),
                        n_beads=len(beads))
    if findall:
        return dict(result=not allfound,
                    necklaces=[dict(vertex=i, necklace=_export(n)) for i, n in allfound])
    return dict(result=True)


def linkcond(isosig, sqlens, allow_repeats=True):
    """Full local CAT(0) test of a PE metric (Gromov's link condition):
    every tetrahedron valid, every edge link of length >= 2*pi (exact), and
    every vertex link without a short closed geodesic.

    Returns dict(result=True/False, reason=None/'invalid'/'edge'/'vertex', ...)."""
    from .edges import edge_report
    rep = edge_report(isosig, sqlens)
    if not rep['valid']:
        return dict(result=False, reason='invalid', degenerate_tet=rep['degenerate_tet'])
    if not rep['edge_condition']:
        return dict(result=False, reason='edge', failing_edges=rep['failing_edges'])
    r = vertex_linkcond(isosig, sqlens, allow_repeats=allow_repeats)
    r.update(reason=None if r['result'] else 'vertex', flat=rep['flat'])
    return r


def _export(necklace):
    pi_point = getattr(necklace, 'pi_point', None)
    return [dict(ts=list(b.ts), start=b.start, end=b.end, length=b.length,
                 pi=bool(b.pi_range), steps=list(b.strip.steps),
                 pi_point=_mq_vec(pi_point) if (b.pi_range and pi_point is not None) else None)
            for b in necklace]


def _mq_vec(v):
    """Serialise a vector of field elements: [{squarefree s: 'p/q'}, ...]."""
    return [{str(s): str(c) for s, c in x.t.items()} for x in v]


CERT_FORMAT = 'vlink-certificate/1'


def certificate(isosig, sqlens):
    """linkcond, plus a certificate for a False answer that
    certcheck.check_certificate verifies exactly and independently:

      reason 'invalid': a tetrahedron that is not a non-degenerate Euclidean one
      reason 'edge'   : an edge whose dihedral angles sum to < 2*pi
      reason 'vertex' : a necklace of beads in a vertex link: each bead with its
                        link triangles, start/end corners and the strip's gluing
                        steps; a pi-bead with the point fixing its geodesic.

    True answers carry no certificate yet (certificate=None)."""
    r = linkcond(isosig, sqlens)
    cert = None
    if r['result'] is False:
        cert = dict(format=CERT_FORMAT, isosig=isosig, sqlens=list(sqlens), reason=r['reason'])
        if r['reason'] == 'invalid':
            cert['tet'] = r['degenerate_tet']
        elif r['reason'] == 'edge':
            cert['edge'] = r['failing_edges'][0]
        else:
            cert['vertex'] = r['vertex']
            cert['necklace'] = [dict(ts=b['ts'], start=b['start'], end=b['end'],
                                     steps=b['steps'], pi_point=b['pi_point'])
                                for b in r['necklace']]
    r['certificate'] = cert
    return r
