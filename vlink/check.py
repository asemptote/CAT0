"""Entry points: the link condition for a PE triangulated 3-manifold.

    from vlink.check import linkcond
    linkcond('gLLPQaceedffnkxknnn', [1, 2, 3, 1, 1, 2, 2])
    -> {'result': True, 'reason': None, 'flat': True, 'verified': False, ...}

linkcond runs the full test (tetrahedron validity, edge condition, vertex
links). Every False answer comes with a certificate (kernel.CERT_FORMAT) that
the kernel checks before the answer is returned ('verified': True); a False
answer whose certificate the kernel rejects is never returned (InternalError).
True answers are not yet certified ('verified': False).

vertex_linkcond is the vertex-link search on its own (same contract as
gallery.vertex_linkcond). Both raise necklace.UndecidedError instead of
guessing in the rare cases not decided exactly.
"""

from . import kernel
from .beads import all_beads
from .edges import edge_report
from .link import VertexLink
from .necklace import search

CERT_FORMAT = kernel.CERT_FORMAT


class InternalError(Exception):
    """The search produced an answer the kernel does not confirm."""


def vertex_linkcond(isosig, sqlens, findall=False, allow_repeats=True, stop_at_start=False,
                    pi_beads=True):
    """allow_repeats (default) lets beads and necklaces reuse link triangles
    (gallery.py forbids reuse, an unproven assumption); stop_at_start=True
    reproduces gallery.py's rule of never extending a path past its start.
    pi_beads=False leaves out beads of length exactly pi: a True answer then
    also relies on Elder-McCammond-Meier (2003) Lemma 4.9."""
    geo = kernel.Geometry(isosig, sqlens)
    allfound = []
    for i_v in range(geo.T.n_vertices):
        link = VertexLink(geo.T, i_v, sqlens, geo=geo)
        beads = all_beads(link, allow_repeats, pi_beads)
        res = search(link, beads, findall=findall, allow_repeats=allow_repeats,
                     stop_at_start=stop_at_start)
        if findall:
            allfound += [(i_v, n) for n in res]
        elif res:
            return dict(result=False, vertex=i_v, necklace=_export(res),
                        n_beads=len(beads))
    if findall:
        return dict(result=not allfound,
                    necklaces=[dict(vertex=i, necklace=_export(n)) for i, n in allfound])
    return dict(result=True)


def linkcond(isosig, sqlens, allow_repeats=True, stop_at_start=False, pi_beads=True,
             certify=False, exact_only=False):
    """Full local CAT(0) test (Gromov's link condition): every tetrahedron
    valid, every edge link of length >= 2*pi, every vertex link without a
    short closed geodesic.

    Returns dict(result, reason in {None, 'invalid', 'edge', 'vertex'},
    certificate (False answers), verified, ...)."""
    rep = edge_report(isosig, sqlens)
    if not rep['valid']:
        r = dict(result=False, reason='invalid', degenerate_tet=rep['degenerate_tet'])
    elif not rep['edge_condition']:
        r = dict(result=False, reason='edge', failing_edges=rep['failing_edges'])
    else:
        r = vertex_linkcond(isosig, sqlens, allow_repeats=allow_repeats,
                            stop_at_start=stop_at_start, pi_beads=pi_beads)
        if not pi_beads:
            r['assumptions'] = ['Elder-McCammond-Meier 2003, Lemma 4.9']
        r.update(reason=None if r['result'] else 'vertex', flat=rep['flat'])
    r['certificate'] = _certificate(isosig, sqlens, r)
    r['verified'] = False
    if r['result'] is False:
        ok, msg = kernel.check_certificate(r['certificate'])
        if not ok:
            raise InternalError('kernel rejects the search result: %s' % msg)
        r['verified'], r['kernel'] = True, msg
    elif certify:
        # a True answer: shortest-path-tree certificate, checked by the kernel
        from .certify import true_certificate
        cert = true_certificate(isosig, sqlens, pi_beads=pi_beads, exact_only=exact_only)
        ok, msg = kernel.check_true_certificate(cert, exact_only=exact_only)
        r['certificate'] = cert
        if ok:
            r['verified'], r['kernel'] = True, msg
        else:
            r['certificate_error'] = msg
    return r


def certificate(isosig, sqlens):
    """linkcond (kept for compatibility: certificates are now always attached)."""
    return linkcond(isosig, sqlens)


def _certificate(isosig, sqlens, r):
    if r['result'] is not False:
        return None
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
    return cert


def _export(necklace):
    pi_point = getattr(necklace, 'pi_point', None)
    return [dict(ts=list(b.ts), start=b.start, end=b.end, length=b.length,
                 pi=bool(b.pi_range), steps=list(b.strip.steps),
                 pi_point=_mq_vec(pi_point) if (b.pi_range and pi_point is not None) else None)
            for b in necklace]


def _mq_vec(v):
    """Serialise a vector of field elements: [{squarefree s: 'p/q'}, ...]."""
    return [{str(s): str(c) for s, c in x.t.items()} for x in v]
