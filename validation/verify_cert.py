"""Independent floating-point checks of inputs and necklace certificates.

Shares no code with gallery.py / constrsphtrig.py / the 2020 vLinkCond.py: each
tetrahedron is embedded in R^3 from its squared edge lengths, link triangles are
spherical triangles of actual unit vectors, and beads are developed onto S^2 by
rigid motions. Regina is used only for combinatorics (isosig decoding,
gluings, edge and vertex numbering), independently of vlink.tri.

metric_report(isosig, sqlens)
    - is the input well formed (one length per edge)?
    - is every tetrahedron a non-degenerate Euclidean tetrahedron?
    - total dihedral angle around every edge: the edge-link condition needs
      >= 2*pi everywhere; the metric is flat iff it is == 2*pi everywhere
      (and a flat metric is locally CAT(0), so the checker must say True).

verify_necklace(isosig, sqlens, vertex, beads)
    beads: [{'ts': [link triangle indices], 'start': i, 'end': j}, ...] where
    start/end are vertex numbers (0..2) of the first/last link triangle.
    Checks that the necklace really carries a closed local geodesic of length
    < 2*pi in the link of `vertex`:
      - each bead develops to a strip whose start-end great-circle arc
        (length < pi) crosses every interior edge of the strip;
      - consecutive beads meet at the same link vertex and the necklace closes;
      - at every bead junction the angle on both sides is >= pi;
      - total length < 2*pi.
    Borderline comparisons (within TOL) are reported, not silently decided.
"""

import math
import sys

import numpy as np
import regina

TOL = 1e-9
ANGLE_TOL = 1e-7  # angles are sums of several acos() values
TWO_PI = 2 * math.pi


# ---------------------------------------------------------------- geometry

def embed_tet(d2):
    """Points p0..p3 in R^3 for squared lengths d2[(i,j)], or None if degenerate.

    Returns (points, volume6) with volume6 = 6*signed volume (>0 by construction).
    """
    p = [np.zeros(3) for _ in range(4)]
    a = d2[(0, 1)]
    if a <= 0:
        return None, 0.0
    x1 = math.sqrt(a)
    p[1] = np.array([x1, 0, 0])
    # p2 in the xy-plane
    x2 = (d2[(0, 2)] - d2[(1, 2)] + a) / (2 * x1)
    y2sq = d2[(0, 2)] - x2 * x2
    if y2sq <= TOL:
        return None, 0.0
    y2 = math.sqrt(y2sq)
    p[2] = np.array([x2, y2, 0])
    x3 = (d2[(0, 3)] - d2[(1, 3)] + a) / (2 * x1)
    y3 = (d2[(0, 3)] - d2[(2, 3)] + x2 * x2 + y2 * y2 - 2 * x2 * x3) / (2 * y2)
    z3sq = d2[(0, 3)] - x3 * x3 - y3 * y3
    if z3sq <= TOL:
        return None, z3sq
    p[3] = np.array([x3, y3, math.sqrt(z3sq)])
    return p, float(np.dot(np.cross(p[1], p[2]), p[3]))


def unit(v):
    return v / np.linalg.norm(v)


def tet_sqlens(T, sqlens, tet):
    d2 = {}
    for k in range(6):
        m = tet.edgeMapping(k)
        i, j = sorted((m[0], m[1]))
        d2[(i, j)] = d2[(j, i)] = float(sqlens[tet.edge(k).index()])
    return d2


def dihedral(p, i, j):
    """Interior dihedral angle of the tetrahedron p along edge (i, j)."""
    k, l = [x for x in range(4) if x not in (i, j)]
    e = unit(p[j] - p[i])
    u = p[k] - p[i]; u = unit(u - np.dot(u, e) * e)
    w = p[l] - p[i]; w = unit(w - np.dot(w, e) * e)
    return math.acos(max(-1.0, min(1.0, float(np.dot(u, w)))))


# ---------------------------------------------------------------- inputs

def metric_report(isosig, sqlens):
    T = regina.Triangulation3.fromIsoSig(isosig)
    rep = dict(n_tet=T.size(), n_edges=T.countEdges(), n_vertices=T.countVertices(),
               n_sqlens=len(sqlens), orientable=T.isOrientable())
    if len(sqlens) != T.countEdges():
        rep['ok'] = False
        rep['problem'] = 'expected %d edge lengths, got %d' % (T.countEdges(), len(sqlens))
        return rep
    bad, embeds = [], {}
    for tet in T.tetrahedra():
        p, vol6 = embed_tet(tet_sqlens(T, sqlens, tet))
        if p is None:
            bad.append(tet.index())
        else:
            embeds[tet.index()] = p
    rep['degenerate_tets'] = bad
    if bad:
        rep['ok'] = False
        rep['problem'] = '%d of %d tetrahedra are not non-degenerate Euclidean tetrahedra' % (
            len(bad), T.size())
        return rep
    sums = []
    for e in T.edges():
        s = 0.0
        for emb in e.embeddings():
            m = emb.vertices()
            s += dihedral(embeds[emb.tetrahedron().index()], m[0], m[1])
        sums.append(s)
    rep['edge_angle_sums_over_2pi'] = [round(s / TWO_PI, 12) for s in sums]
    rep['edge_condition'] = all(s >= TWO_PI - TOL for s in sums)
    rep['flat'] = all(abs(s - TWO_PI) <= TOL for s in sums)
    rep['ok'] = True
    return rep


# ---------------------------------------------------------------- link geometry

class Link:
    """Vertex link of T at vertex i_v, with each link triangle given as three
    unit vectors (in its own tetrahedron's frame)."""

    def __init__(self, isosig, sqlens, i_v, numbering='canonical'):
        self.T = T = regina.Triangulation3.fromIsoSig(isosig)
        v = T.vertex(i_v)
        if numbering == 'regina':
            # Regina's own vertex link (the 2020 code's numbering)
            L, I = v.buildLink(), v.buildLinkInclusion()
            fps = [(I.tetImage(t), I.facePerm(t)) for t in range(L.size())]
            self._adj = [[L.triangle(t).adjacentTriangle(d).index() for d in range(3)]
                         for t in range(L.size())]
            self._glue = [[tuple(L.triangle(t).adjacentGluing(d)[j] for j in range(3))
                           for d in range(3)] for t in range(L.size())]
            self._lv = [[L.triangle(t).vertex(j).index() for j in range(3)]
                        for t in range(L.size())]
        else:
            # vlink's numbering: the corners (tet, v) at the vertex, sorted;
            # vertex j of a corner is the j-th smallest other tet vertex
            corners = sorted((e.tetrahedron().index(), e.vertex()) for e in v.embeddings())
            index = {c: t for t, c in enumerate(corners)}
            fps = [(tet, [x for x in range(4) if x != w] + [w]) for tet, w in corners]
            self._adj, self._glue, self._lv = [], [], []
            for tet, fp in fps:
                rt = T.tetrahedron(tet)
                adj, glue = [], []
                for d in range(3):
                    g = rt.adjacentGluing(fp[d])
                    w2 = g[fp[3]]
                    fp2 = [x for x in range(4) if x != w2]
                    adj.append(index[(rt.adjacentTetrahedron(fp[d]).index(), w2)])
                    glue.append(tuple(fp2.index(g[fp[j]]) for j in range(3)))
                self._adj.append(adj)
                self._glue.append(glue)
                lv = []
                for j in range(3):
                    k = [i for i in range(6) if {rt.edgeMapping(i)[0], rt.edgeMapping(i)[1]}
                         == {fp[3], fp[j]}][0]
                    lv.append((rt.edge(k).index(), 0 if rt.edgeMapping(k)[0] == fp[3] else 1))
                self._lv.append(lv)
        self.n = len(fps)
        self.vecs = []
        for t in range(self.n):
            tet = T.tetrahedron(fps[t][0])
            fp = fps[t][1]
            p, _ = embed_tet(tet_sqlens(T, sqlens, tet))
            if p is None:
                raise ValueError('degenerate tetrahedron %d' % tet.index())
            c = p[fp[3]]
            assert tet.vertex(fp[3]).index() == i_v
            self.vecs.append([unit(p[fp[j]] - c) for j in range(3)])

    def side(self, t, a, b):
        return math.acos(max(-1.0, min(1.0, float(np.dot(self.vecs[t][a], self.vecs[t][b])))))

    def corner(self, t, c):
        """Spherical angle of link triangle t at its vertex c (= a dihedral angle)."""
        a, b = [x for x in range(3) if x != c]
        o = self.vecs[t][c]
        u = self.vecs[t][a] - np.dot(self.vecs[t][a], o) * o
        w = self.vecs[t][b] - np.dot(self.vecs[t][b], o) * o
        return math.acos(max(-1.0, min(1.0, float(np.dot(unit(u), unit(w))))))

    def lv(self, t, c):
        """Label of the link vertex at corner c of link triangle t."""
        return self._lv[t][c]

    def adj(self, t, d):
        """The link triangle across the edge opposite vertex d of t."""
        return self._adj[t][d]

    def gluing(self, t, d):
        """Vertex j of t -> vertex gluing(t, d)[j] of adj(t, d)."""
        return self._glue[t][d]

    def star(self, t, c):
        """Cyclic list of (triangle, corner, entry_vertex, exit_vertex) around
        the link vertex at corner c of triangle t; entry/exit are the other two
        corners, ordered by the walk."""
        out = []
        a, b = [x for x in range(3) if x != c]
        t0, c0, entry = t, c, a
        while True:
            exitv = 3 - c - entry
            out.append((t, c, entry, exitv))
            nxt = self.adj(t, entry)  # across the edge opposite `entry`
            g = self.gluing(t, entry)
            t, c, entry = nxt, g[c], g[exitv]
            if (t, c) == (t0, c0):
                return out
            if len(out) > 4 * self.n:
                raise RuntimeError('star walk did not close')


def frame(a, b):
    e1 = unit(a)
    e2 = unit(b - np.dot(b, e1) * e1)
    return np.array([e1, e2, np.cross(e1, e2)])


def develop(link, ts, choices):
    """Positions on S^2 of every triangle of a strip, developed rigidly.

    choices[k] = vertex number d of ts[k] whose opposite edge is glued to ts[k+1].
    Returns list of 3x3 arrays (rows = developed unit vectors of vertices 0,1,2).
    """
    pos = [np.array(link.vecs[ts[0]])]
    for k in range(len(ts) - 1):
        t, d = ts[k], choices[k]
        g = link.gluing(t, d)
        t2 = ts[k + 1]
        a, b = [x for x in range(3) if x != d]
        P = pos[-1]
        Q = np.array(link.vecs[t2])
        F = frame(P[a], P[b])
        Fq = frame(Q[g[a]], Q[g[b]])
        coords = Q @ Fq.T  # coordinates of t2's vertices in its own (a,b) frame
        # third vertex of t2 must land on the opposite side of the edge from P[d]
        side_d = np.dot(P[d], F[2])
        s = -1.0 if np.sign(coords[g[d]][2]) == np.sign(side_d) else 1.0
        coords[:, 2] *= s
        pos.append(coords @ F)
    return pos


def arc_crosses_edge(ps, pe, a, b):
    """Does the great-circle arc ps->pe (length < pi) meet the edge arc a->b in
    its relative interior?  Returns (bool, slack) with slack the angular margin
    to the edge's endpoints."""
    n1 = np.cross(ps, pe)
    n2 = np.cross(a, b)
    x = np.cross(n1, n2)
    if np.linalg.norm(x) < 1e-14:
        return False, -1.0
    x = unit(x)
    best = (False, -1.0)
    for X in (x, -x):
        on_geod = (abs(math.acos(np.clip(np.dot(ps, X), -1, 1)) +
                       math.acos(np.clip(np.dot(X, pe), -1, 1)) -
                       math.acos(np.clip(np.dot(ps, pe), -1, 1))) < 1e-9)
        da = math.acos(np.clip(np.dot(a, X), -1, 1))
        db = math.acos(np.clip(np.dot(X, b), -1, 1))
        on_edge = abs(da + db - math.acos(np.clip(np.dot(a, b), -1, 1))) < 1e-9
        if on_geod and on_edge:
            best = (True, min(da, db))
    return best


def tangent(at, to):
    v = to - np.dot(to, at) * at
    return unit(v)


def angle_between(u, w):
    return math.acos(max(-1.0, min(1.0, float(np.dot(u, w)))))


def verify_bead(link, bead):
    """Develop a bead and check its geodesic. Returns dict."""
    ts, s, e = bead['ts'], bead['start'], bead['end']
    # Enumerate the gluing choices that make ts a connected strip.
    options = [[d for d in range(3)
                if link.adj(ts[k], d) == ts[k + 1]]
               for k in range(len(ts) - 1)]
    results = []

    def rec(k, choice):
        if k == len(ts) - 1:
            results.append(list(choice))
            return
        for d in options[k]:
            rec(k + 1, choice + [d])
    rec(0, [])
    attempts = []
    for choice in results:
        pos = develop(link, ts, choice)
        ps, pe = pos[0][s], pos[-1][e]
        length = angle_between(ps, pe)
        crossings = []
        ok = length < math.pi - TOL or len(ts) == 1
        for k, d in enumerate(choice):
            a, b = [x for x in range(3) if x != d]
            hit, slack = arc_crosses_edge(ps, pe, pos[k][a], pos[k][b])
            crossings.append(slack)
            ok = ok and hit and slack > TOL
        # a one-triangle bead runs along the edge (s, e) of that triangle
        if len(ts) == 1:
            ok = ok and s != e
        attempts.append(dict(choice=choice, ok=bool(ok), length=length,
                             min_crossing_slack=min(crossings) if crossings else None,
                             pos=pos))
    good = [a for a in attempts if a['ok']]
    if not attempts:
        return dict(ok=False, reason='triangles do not form a strip')
    best = good[0] if good else attempts[0]
    best['n_strip_choices'] = len(attempts)
    if not good:
        best['reason'] = 'geodesic leaves the strip or has length >= pi'
    return best


def direction_in_star(link, star, t, c, toward, pos_t):
    """Angular position (from the start of the star walk) of the direction at
    the corner c of triangle t pointing toward `toward` (developed coords
    matching pos_t, the developed vertices of t)."""
    offset = 0.0
    for (tt, cc, entry, exitv) in star:
        if (tt, cc) == (t, c):
            w = pos_t[c]
            u = tangent(w, pos_t[entry])
            x = tangent(w, toward)
            theta = angle_between(u, x)
            return offset + theta, theta, link.corner(t, c)
        offset += link.corner(tt, cc)
    raise RuntimeError('corner not in star')


def regina_to_canonical(isosig, vertex, beads):
    """One-triangle beads in Regina's link numbering (the 2020 code, and vlink
    before certificate format /2) -> vlink's numbering."""
    T = regina.Triangulation3.fromIsoSig(isosig)
    v = T.vertex(vertex)
    L, I = v.buildLink(), v.buildLinkInclusion()
    corners = sorted((e.tetrahedron().index(), e.vertex()) for e in v.embeddings())
    out = []
    for b in beads:
        if len(b['ts']) != 1:
            raise ValueError('only one-triangle beads are translated')
        t = b['ts'][0]
        fp = I.facePerm(t)
        c = (I.tetImage(t), fp[3])
        others = [x for x in range(4) if x != fp[3]]
        out.append(dict(b, ts=[corners.index(c)], start=others.index(fp[b['start']]),
                        end=others.index(fp[b['end']])))
    return out


def verify_necklace(isosig, sqlens, vertex, beads, numbering="canonical"):
    link = Link(isosig, sqlens, vertex, numbering)
    out = dict(beads=[], junctions=[])
    bead_res = [verify_bead(link, b) for b in beads]
    for b, r in zip(beads, bead_res):
        out['beads'].append(dict(ts=b['ts'], start=b['start'], end=b['end'],
                                 ok=r['ok'], length=r.get('length'),
                                 reason=r.get('reason'),
                                 min_crossing_slack=r.get('min_crossing_slack'),
                                 strip_choices=r.get('n_strip_choices')))
    total = sum(r.get('length', 0) for r in bead_res)
    out['total_length_over_2pi'] = total / TWO_PI
    ok = all(r['ok'] for r in bead_res) and total < TWO_PI - TOL
    m = len(beads)
    for k in range(m):
        b1, b2 = beads[k], beads[(k + 1) % m]
        r1, r2 = bead_res[k], bead_res[(k + 1) % m]
        t1, c1 = b1['ts'][-1], b1['end']
        t2, c2 = b2['ts'][0], b2['start']
        j = dict(link_vertex=(link.lv(t1, c1), link.lv(t2, c2)))
        if j['link_vertex'][0] != j['link_vertex'][1]:
            j['ok'] = False
            j['reason'] = 'beads do not meet'
            ok = False
            out['junctions'].append(j)
            continue
        if 'pos' not in r1 or 'pos' not in r2:
            j['ok'] = False
            ok = False
            out['junctions'].append(j)
            continue
        star = link.star(t1, c1)
        total_angle = sum(link.corner(tt, cc) for tt, cc, _, _ in star)
        # incoming geodesic, reversed: from the junction back toward bead 1's start
        P1 = r1['pos'][-1]
        back = r1['pos'][0][b1['start']]
        P2 = r2['pos'][0]
        fwd = r2['pos'][-1][b2['end']]
        a_in, _, _ = direction_in_star(link, star, t1, c1, back, P1)
        a_out, _, _ = direction_in_star(link, star, t2, c2, fwd, P2)
        side1 = (a_out - a_in) % total_angle
        side2 = total_angle - side1
        # A local geodesic needs angle >= pi on both sides; equality is allowed,
        # so compare with a float tolerance and flag the equality cases.
        j.update(cone_angle_over_pi=total_angle / math.pi,
                 sides_over_pi=(side1 / math.pi, side2 / math.pi),
                 ok=bool(min(side1, side2) >= math.pi - ANGLE_TOL),
                 borderline=bool(abs(min(side1, side2) - math.pi) <= ANGLE_TOL))
        ok = ok and j['ok']
        out['junctions'].append(j)
    out['valid'] = bool(ok)
    return out


if __name__ == '__main__':
    import json
    for line in sys.stdin:
        r = json.loads(line)
        rep = metric_report(r['isosig'], r['sqlens'])
        r['metric'] = rep
        if r.get('result') is False and rep.get('ok'):
            r['verification'] = verify_necklace(r['isosig'], r['sqlens'],
                                                r['vertex'], r['necklace'])
        print(json.dumps(r, default=float))
