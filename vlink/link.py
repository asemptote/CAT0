"""Vertex links for the search: a thin layer over the kernel's link data.

The geometry (rational Gram matrices of link triangles, gluings, strip
transitions) lives in kernel.py; this module only adds the search's view of
a developed strip (Strip), which records the gluing steps used so that the
kernel can re-develop and check any bead the search reports.
"""

from . import kernel
from .field import ZERO


class VertexLink:
    """The link of vertex i_v of T (an isosig or tri.Triangulation) with
    squared edge lengths sqlens.

    Attributes mirror kernel.Link: n (triangles), n_vertices, gram[t] (rational
    Gram matrix of triangle t), adj, glue, lv; `k` is the kernel Link itself."""

    def __init__(self, T, i_v, sqlens, geo=None):
        geo = geo or kernel.Geometry(T, sqlens)
        self.T, self.i_v, self.geo = geo.T, i_v, geo
        self.k = k = kernel.Link(geo, i_v)
        self.n, self.n_vertices = k.n, k.n_vertices
        self.gram, self.adj, self.glue, self.lv = k.G, k.adj, k.glue, k.lv

    def start(self, t):
        """Triangle t developed in its own frame (its edge vectors)."""
        return Strip(self, t, kernel.E3)


class Strip:
    """A developed sequence of link triangles (all in the first one's frame).

    steps[k] is the vertex of ts[k] opposite the edge crossed into ts[k+1]."""

    __slots__ = ('link', 'ts', 'pos', 'frame', 'steps')

    def __init__(self, link, t, pos, ts=None, frame=None):
        self.link = link
        self.ts = ts or (t,)
        self.pos = (tuple(pos),)
        self.frame = frame if frame is not None else t
        self.steps = ()

    def extend(self, d):
        """New strip with the neighbour across the edge opposite vertex d of
        the last triangle appended (kernel.step)."""
        t2, P2 = kernel.step(self.link.k, self.pos[-1], self.ts[-1], d)
        s = Strip.__new__(Strip)
        s.link, s.ts, s.pos, s.frame = self.link, self.ts + (t2,), self.pos + (P2,), self.frame
        s.steps = self.steps + (d,)
        return s

    @property
    def G(self):
        return self.link.gram[self.frame]

    def dot(self, x, y):
        return kernel.dot(self.G, x, y)

    det = staticmethod(kernel.det)

    def cos_float(self, x, y):
        return float(self.dot(x, y)) / (float(self.dot(x, x)) * float(self.dot(y, y))) ** 0.5


def link_of(isosig, i_v, sqlens):
    return VertexLink(isosig, i_v, sqlens)


__all__ = ['VertexLink', 'Strip', 'link_of', 'ZERO']
