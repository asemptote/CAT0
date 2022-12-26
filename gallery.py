"""Constructible non-Euclidean geometry with Regina.

Author: Asem Wardak

Regina focuses on computational topology, so to perform geometric calculations
one is required to use an auxiliary package such as SnapPy, which uses
arbitrary-precision rational numbers, or Sage [1], which supports algebraic
number fields. Constructible numbers comprise a strict superset of the rational
numbers and a strict subset of the algebraic numbers. This package uses an
implementation of constructible spherical trigonometry to demonstrate the
construction of galleries based on the outline given in Elder and McCammond
(2002, 2004), along with how it can be used to check Gromov's link condition
for vertices (Elder and McCammond's GAP program [2] generating candidate
galleries without constructing a loop).


Links:
    [1] https://doc.sagemath.org/html/en/reference/number_fields/index.html
    [2] http://web.math.ucsb.edu/~jon.mccammond/software/cat.html
        http://web.math.ucsb.edu/~jon.mccammond/software/catFiles/
        http://web.math.ucsb.edu/~jon.mccammond/software/doc3.pdf

References:
    `Curvature Testing in 3-Dimensional Metric Polyhedral Complexes`
    Elder and McCammond, Experimental Mathematics (2002)
    https://doi.org/10.1080/10586458.2002.10504476
    
    `CAT(0) is an Algorithmic Property`
    Elder and McCammond, Geometriae Dedicata (2004)
    https://doi.org/10.1023/B:GEOM.0000049096.63639.e3
    Open Question 7.2: `Is there a direct _geometric_ algorithm which
    determines the local curvature of a finite $M_k$-complex in dimensions
    greater than 3?` TODO: for the 4D case, perhaps edge-to-edge `beads` are
    unshrinkable iff there are only two such `beads` with orthogonal edges?
    
    `Combinatorial conditions that imply word-hyperbolicity for 3-manifolds`
    Elder, McCammond and Meier, Topology (2003)
    https://doi.org/10.1016/S0040-9383(02)00100-3
    Section 4; p.1254: `Because we are interested in triangulated 3-manifolds,
    and Mobius strips cannot be immersed into 2-spheres, the Mobius strips in
    the output can be safely ignored.` As a result, we only need to check
    vertex-to-vertex galleries for short unshrinkable loops.

Dependencies: constrsphtrig.py, Regina
"""

import math

try: from collections.abc import Sequence
except ImportError: from collections import Sequence

try: from itertools import imap
except ImportError: imap = map  # Python 3

from contextlib import contextmanager
from functools import partial
from pickle import dumps, loads

from multiprocessing import Pool
import sys

from constrsphtrig import (ConstrAngle, pi, zero, minuspi,
    ConstrAngleLength, twopi, onepi, zerolength,
    ConstrSphPoint, zero_zero, zero_pi,
    sqrt)

if 'regina' not in dir(): import regina
fromIsoSig = (regina.NTriangulation.fromIsoSig
              if 'NTriangulation' in dir(regina)
              else regina.Triangulation3.fromIsoSig)

DEBUG = True


class Gallery(object):
    """A sequence of convex polyhedral cells.
    
    For more information, see Elder and McCammond (2002), Section 3.
    `Thus, the entire focus will be on the link of a 0-cell which is a
    2-dimensional spherical complex.`
    
    We consider vertex-to-vertex galleries (beads) and vertex-to-edge galleries
    (bags).
    
    If len(ds) == len(rs), the last element of rs is the
    index of the final bottom cell (vertex); otherwise, the last element of
    ds is the index of the final bottom cell (edge): in this case rs is
    one element smaller.
    
    Edges are indexed by the opposite vertex.
    
    Thick beads: ds and rs may be viewed as input-output pairs of the vertex
    mapping between pairs of triangles. In this way, vertex-to-edge galleries
    contain enough information to be completed into a vertex-to-vertex gallery.
    
    To add a triangle to a bead, use .tobag(d) then .tobead()
    
    Some of these functions are designed for bags, others for beads.
    
    Attributes:
        T: The triangulation of the 3-manifold.
        L: The triangulation of the vertex link.
        I: The mapping between the vertex link L and the triangulation T.
        i_v: The index of the vertex of T.
        sqedgelens: Squares of edge lengths of tetrahedra using Regina's
               numbering.
        ts: A tuple of triangles of L (top cells) comprising the gallery.
            NOTE: these are represented by indices of L.triangles()
        ps: A tuple of tuples of spherical positions of triangle points when
            developed onto the sphere.
        ds: A tuple of direction vectors (vertex indices of opposite edge)
            denoting the direction of the next glued triangle. The first
            element of ds is the first bottom cell (vertex).
        rs: A tuple of reverse direction vectors denoting the direction of
            the previous glued triangle, starting from the second triangle. 
    The following attributes are only computed on-demand:
        geodesiclength: The length of the geodesic, when there is one.
        geodesicangles: ((minimum angle between (e_0+1)%3, e_0 and geodesic,
          minimum angle between (e_0+2)%3, e_0 and geodesic),
         (minimum angle between (e_m+1)%3, e_m and geodesic,
          minimum angle between (e_m+2)%3, e_m and geodesic), [supp.],
          [e_0 = ds[0], e_m=rs[-1]] # firstbottomcellindex, lastbottomcellindex
          if the geodesic exists, or False.  TODO: polish
        geodesicrange: size of angle over which geodesics exist if
                       geodesiclength==pi, or False.
    """
    
    def __init__(self, link_triangles,
                 ts=tuple(), ps=tuple(), ds=tuple(), rs=tuple(),
                 geodesiclength=None, geodesicangles=None, geodesicrange=None):
        self.link_triangles = link_triangles
        self.ts = ts
        self.ps = ps
        self.ds = ds
        self.rs = rs
        self.geodesiclength = geodesiclength
        self.geodesicangles = geodesicangles
        self.geodesicrange = geodesicrange
    
    @property
    def isbead(self):
        return len(self.ds) == len(self.rs)
    
    @property
    def start(self): return self.ps[0][self.ds[0]]
    
    @property
    def end(self): return (self.ps[-1][self.rs[-1]] if self.isbead
                           else self.ps[-1][self.ds[-1]])
    
    @property
    def eligible(self):
        """Tests a vertex-to-edge gallery for eligibility.
        
        A partial vertex-to-vertex gallery is eligible if any galleries built
        from it can contain a short geodesic. To check this, the straddling
        condition is checked for pairs of points on each triangle opposite the
        directional point.
        
        If the gallery is eligible, then by viewing bags as beads of length pi
        with the last vertex at the antipode, self.geodesicangles and
        self.geodesicrange is updated; otherwise self.geodesicangles is set to
        False.
        
        Returns: a boolean
        
        Side effects: initialises self.geodesicangles, self.geodesicrange
        """
        # We guarantee that the triangles have angles smaller than 180* and are thus oriented towards the right (the tets are convex)
        d0 = self.ds[0]
        p0 = self.ps[0]
        theta = p0[(d0+1)%3].ang().bisect(p0[(d0+2)%3].ang())
        # orient first tri perfectly to the right
        self.ps = tuple(tuple(pt.rotate(theta) for pt in p) for p in self.ps)
        # assume not self.isbead
        pairs = [(p[(d+1)%3], p[(d+2)%3]) for p, d in zip(self.ps, self.ds)]
        if (any(zero_pi in pair for pair in pairs)
            or any(zero_zero in pair for pair in pairs)
            ):
            self.geodesicangles = False
            return False
        angs = [[pt.ang() for pt in pair] for pair in pairs]
        arange = [minuspi, pi]  # min ang in range, max ang in range
        for angmin, angmax in map(sorted, angs):
            if angmax - angmin < zero:  # pair range goes through left equator
                self.geodesicangles = False
                return False
            else:
                arange[0] = max(arange[0], angmin)
                arange[1] = min(arange[1], angmax)
            # We require max_mins to be strictly smaller than min_maxs, for if
            # they are equal there is a vertex bottom cell along the gallery,
            # and so the gallery is reducible into two (not necess. thick) v2v
            # galleries.
#            print(arange)
            if arange[1] <= arange[0]:
                self.geodesicangles = False
                return False
        max_mins, min_maxs = arange
        assert max_mins <= min_maxs
        # Angles about a point are equal to (supp.) angles about the antipodal
        # point. The angles are positive if they are above the geodesic lune,
        # and vice versa.
        endangs = [angs[0], angs[-1]]
        m_ms = (max_mins, min_maxs)  # sorted
        # If the first angle in the pair is higher, the offset is zero
        offsets = [int(endang[0]<endang[1]) for endang in endangs]
        ang0m = tuple((a[0] - m_ms[1-i], a[1] - m_ms[i])
                      for a, i in zip(endangs, offsets))
        self.geodesicangles = ang0m
        self.geodesicrange = min_maxs - max_mins
        return True
    
    @property
    def contains_shortgeodesic(self):
        """Tests a vertex-to-vertex gallery for a short geodesic.
        
        Assume gallery is very short (< pi) (Elder 2003 Sec. 4).
        Assume self.isbead == True.
        
        Returns: a boolean
        
        Side effects: initialises self.geodesicangles, self.geodesicrange
        """
        assert self.isbead
        self_end = self.end
        if self_end == zero_zero:  # i.e. 2*pi
            self.geodesiclength = zero
            return False
        # Geodesic length is not pi so there is only one geodesic path
        # rotate geodesic onto equator
        self_end_ang = self_end.ang()
        ps = [[pt.rotate(self_end_ang) for pt in p] for p in self.ps]
        # ensure that the first triangle is pointing to the right after rot.
        d0 = self.ds[0]
        rm = self.rs[-1]
        theta = ps[0][(d0+1)%3].ang().bisect(ps[0][(d0+2)%3].ang())
        if theta.cos < 0: ps = [[pt.rotate(pi) for pt in p] for p in ps]
        lon = ps[-1][rm].lon  # endpoint lon. == geodesic length
        self.geodesiclength = lon
        if lon < zero: return False  # geodesic longer than pi
        p00, p01 = (ps[0][(d0+i)%3] for i in (1,2))
        pm0, pm1 = (ps[-1][(rm+i)%3] for i in (1,2))
        # Hemispheres of nonendpoints for straddling
        hs = ([[p00.hemi, p01.hemi]] +
              [[pt.hemi for pt in p] for p in ps[1:-1]] +
              [[pm0.hemi, pm1.hemi]])
        if all(hemi for h in hs for hemi in h) and all(len(set(h))==2 for h in hs):
            # Straddling condition passed; compute (signed) boundary angles
            self.geodesicangles = ((p00.ang(), p01.ang()),
                                   (pi - pm0.ang(lon), pi - pm1.ang(lon)))
            self.geodesicrange = None
            return True
        return False  # straddling failed
    
    def tobead(self):  # assume not isbead
        """Returns the bead formed from a vertex-to-edge gallery (bag).
        
        Returns False if the bead would contain a triangle twice.
        
        Assumes that the property `eligible` has been checked.
        """
        dn = self.ds[-1]
#        tri_n = self.L.triangle(self.ts[-1])
        pn = self.ps[-1]
#        t = tri_n.adjacentTriangle(dn).index()
        t = self.link_triangles.adjacentTriangles[self.ts[-1]][dn]
        if t in self.ts: return False  # new triangle already in list
        # Add new triangle to list.
        p = [zero_zero for _ in range(3)]
        # Construct the vertex map from triangle n to n+1, where gluing occurs
        # on the two vertices opposite ds_n.
#        gl = tri_n.adjacentGluing(dn)
        gl = self.link_triangles.adjacentGluings[self.ts[-1]][dn]
        # Add the two overlapping points between triangle n, n+1.
        for i in (dn+1)%3, (dn+2)%3: p[gl[i]] = pn[i]
        # To determine orientation, rotate the two intermediate points common to
        # both triangles so that they are both on the equator, then find the point
        # on the opposite hemisphere to the point on the previous triangle.
        r = gl[dn]  # complement of self.ds[-1]
        p_r0, p_r1 = p[(r+1)%3], p[(r+2)%3]
#        ltl = link_trilengths(self.T, self.I, self.sqedgelens, t)  # only need 2 out of 3 l's
        ltl = self.link_triangles.lengths[t]
        # ltl uses edge ordering, so that e.g. ltl[(r+2)%3] is the edge length
        # between the two vertices opposite to (r+2)%3, i.e. r and (r+1)%3.
        # The possible third points of triangle n+1:
#        from time import time
#        start = time()
        assert ConstrAngle(*ltl[r], exact=False) == p_r0.dist(p_r1, exact=False)# Really expensive to compute exactly, adding algebraic nunmbers is like adding logs
#        print(f'{time()-start} seconds elapsed.')
        p_r_ = p_r0.linear(p_r1, ltl[(r+2)%3], ltl[(r+1)%3])
        if len(p_r_) < 2: raise Exception('triangle points are geodesically co'
        'llinear -> some of the tetrahedra are flat')  # TODO: use something other than exceptions
        # To see which third point is on the opposite side of the edge relative
        # to the complementary point, rotate the edge onto the equator by
        # rotating the first edge vertex onto the equator about zero_zero, then
        # rotating the second edge vertex onto the equator about the first
        # edge vertex. It then suffices to check if the two points are on
        # different hemispheres. First rotation:
        if p_r0 not in (zero_zero, zero_pi):
            ang1 = p_r0.ang()
            p2_r0 = p_r0.rotate(ang1)  # first vertex on common edge
            p2_r1 = p_r1.rotate(ang1)  # second vertex on common edge
            pt2_n = pn[dn].rotate(ang1)  # point on last triangle
            p2_r_ = [p_r.rotate(ang1) for p_r in p_r_]  # compl. pts on next tri
        else:
            p2_r0 = p_r0
            p2_r1 = p_r1
            pt2_n = pn[dn]
            p2_r_ = p_r_
        # Second rotation:
        ang2 = p2_r1.ang(p2_r0.lon)
        pt2_n = pt2_n.rotate(ang2, axis=p2_r0)
        p2_r_0 = p2_r_[0].rotate(ang2, axis=p2_r0)  # first poss. compl. pt
        p[r] = p_r_[0] if p2_r_0.lat.sgn != pt2_n.lat.sgn else p_r_[1]
        return Gallery(self.link_triangles,
                       self.ts+(t,), self.ps+(tuple(p),), self.ds, self.rs+(r,),
                       pi if p[r] == zero_pi else None, self.geodesicangles,
                       self.geodesicrange)
        
    def tobag(self, d): 
        """Returns the bag formed from a vertex-to-vertex gallery."""
        return Gallery(self.link_triangles,
                       self.ts, self.ps, self.ds+(d,), self.rs)
    
    def plot(self, title=None, geodesic=False):
        """Shows the triangles of the gallery developed onto the sphere.
    
        TODO:
            test it out
            improve using spherical voronoi diagrams
        
        Links: https://stackoverflow.com/q/4622057
        """
        map = lambda f, x: [f(y) for y in x]
        if 'plt' not in dir():
            import mpl_toolkits.mplot3d as a3
            import matplotlib.colors as colors
            import matplotlib.pyplot as plt
            import numpy as np
        fig = plt.figure()
        ax = a3.Axes3D(fig)
        for p in self.ps:
            tri = a3.art3d.Poly3DCollection([np.array([map(float, pt.cart) for pt in p])])
            tri.set_color(colors.rgb2hex(np.random.rand(3)))
            tri.set_edgecolor('k')
            ax.add_collection3d(tri)
        if geodesic: # show the geodesic range too?
            g = [map(float, self.start.cart), map(float, self.end.cart)]
            ax.scatter(*(tuple(g[i][j] for i in range(2)) for j in range(3)), c='r', s=[100,100], edgecolors='black')
        ax.set_xlim3d(-1,1)
        ax.set_ylim3d(-1,1)
        ax.set_zlim3d(-1,1)
        ax.set_title(title, fontsize=7)
        plt.show()
    
    def find_thickbeads(self):  # (partial gallery) bag -> Galleries(...) (beads)
        """Returns all the thick beads which contain a short geodesic.
        
        `self` is assumed to be a bag.
        """
#        G = Galleries(*self.context)
        G = []
        #print(self)
        if not self.eligible:
#            print(f'bag of {len(self.ps)} tris not eligible: {len(G)}')
            return G
        if DEBUG and len(self.ps) > 1: self.plot(geodesic=True)
        #print(self.isbead)
        bead = self.tobead()  # False if new triangle already in gallery
        if bead == False:
            if DEBUG: print(f'{len(self.ps)+1}-th tri already in gal: {len(G)}')
            return G  # faster than `if not bead:`
        if bead.geodesiclength == pi or bead.contains_shortgeodesic:
            if DEBUG: print(f'appending gal of {len(bead.ps)} tris to {len(G)} gals')
            G.append(bead)
        for d in ((bead.rs[-1]+i)%3 for i in (1,2)):
            G.extend(bead.tobag(d).find_thickbeads())
#        if DEBUG and G: print(f'{len(G)} gals')
        return G
  
  
class LinkTriangles(object):
    """The triangles of a vertex link. Triangles are referred to by index.
    """
    def __init__(self, T, v, sqedgelens):
        L = v.buildLink()
        I = v.buildLinkInclusion()
        len_ts = L.triangles().size()
        self.lengths = [link_trilengths(T, I, sqedgelens, t) for t in range(len_ts)]
        self.adjacentGluings = [[[int(x) for x in str(L.triangle(t).adjacentGluing(d))] for d in range(3)] for t in range(len_ts)]
        self.adjacentTriangles = [[L.triangle(t).adjacentTriangle(d).index() for d in range(3)] for t in range(len_ts)]
        self.vertices = [[L.triangle(t).vertex(d).index() for d in range(3)] for t in range(len_ts)]
    
    def geodesicAngle(self, tm, t0, rm, d):
        """Computes the angle between two triangles with common vertex (tet edge).
        
                    ------->
                   /         tri_n (represents all the intermediate triangles)
       (direction)/    3-rm-d______
                 /          |\    /\
                |           | \../  \  tri_0
                |           |  \/____\
                |    tri_m  |  /rm
                |           | /       (angles dotted, added recursively)
                            |/
                            d
        
        Args:
            tri_m: The last triangle of the previous gallery.
            tri_0: The first triangle of the next gallery.
            r_m: The index of the bottom cell in the last triangle of the gallery.
            d: The direction vertex of the next triangle in the chosen path between
               the two triangles (i.e. index of common link edge)
        
        Returns:
            (angle along given direction, closest edge of t_0 reached)
        """
#        tri_n = tri_m.adjacentTriangle(d)
        tn = self.adjacentTriangles[tm][d]
#        gl = tri_m.adjacentGluing(d)
        gl = self.adjacentGluings[tm][d]
        if tn == t0: return (zerolength, gl[d])
#        if tri_n == tri_0: return (zerolength, gl[d])#gl[3-rm-d]) (bug-causing mistake)
        rn = gl[rm]
        ang, r0 = self.geodesicAngle(tn, t0, rn, gl[3-rm-d])
#        ltl = link_trilengths(T, I, sqedgelens, tri_n.index())
        ltl = self.lengths[tn]
        cos = ((ltl[rn].cos - ltl[(rn+1)%3].cos * ltl[(rn+2)%3].cos)
               / (ltl[(rn+1)%3].sin * ltl[(rn+2)%3].sin))  # spherical cosine law
        return (ang + ConstrAngleLength(ConstrAngle(cos=cos)), r0) 


def build_galleries(link_triangles, t):  # -> Galleries(...)
#    T = fromIsoSig(self.isosig)
#    i_v = self.i_v
#    L, I = (T.vertex(i_v).buildLink(), T.vertex(i_v).buildLinkInclusion())
#    G = Galleries(*self.context)
    G = []
#    sqedgelens = self.sqedgelens
#    ltl = link_trilengths(T, I, sqedgelens, t)
    ltl = link_triangles.lengths[t]
    for d0 in range(3):
        # construct the first triangle on the sphere
        p = [zero_zero for _ in range(3)]
        d00, d01 = ((d0+i)%3 for i in (1,2))
        #a01, a02, a12 = ltl[d0:] + ltl[:d0]  # list rotation
        a12, a02, a01 = ltl[d0:] + ltl[:d0]  # list rotation
        p[d00] = ConstrSphPoint(lon=a01)
        p[d01] = zero_zero.linear(p[d00], a02, a12)[0]
        # Add thin beads. Since tet edge angles are not reflex, these beads
        # have length at most pi (<pi if tet is not flat)
        G.append(Gallery(link_triangles,
                     (t,), (tuple(p),), (d0,), (d00,),
                     a01, ((zero, p[d01].ang()),
                           (pi-p[(d00+1)%3].ang(a01),
                            pi-p[(d00+2)%3].ang(a01)))))
        p = [pt.rotate(p[d01].ang()) for pt in p]
        assert float(p[d01].lon) == float(a02)
        # After rotation p[ds_01]==p[(ds_0+2)%3] is on equator
        G.append(Gallery(link_triangles,
                     (t,), (tuple(p),), (d0,), (d01,),
                     a02, ((p[d00].ang(), zero),
                           (pi-p[(d01+1)%3].ang(a02),
                            pi-p[(d01+2)%3].ang(a02)))))
        G.extend(Gallery(link_triangles, (t,), (tuple(p),), (d0,)).find_thickbeads())
    return G


    

#def find_necklaces(self, iss_G, findall, piflag, currentsum, is_G):
def find_necklaces(gs, link_triangles, iss_G, findall, piflag, currentsum, is_G):
    # def find_necklace(self, grouping, is_gs) -> Galleries(..), False, or raise Exception(Galleries(..))?
    # visualise the numbering of exprs in angs with a diagram
    # piflag: "how much extra space do you need to take from the geodesic range in order to get the minimum pi angle between the geodesics"
    # in __doc__ describe some key variables used in the algorithm instead of repeating comments
    # `top angle` is relative to the gallery in question, not preserved around a link vertex
    gm = gs[is_G[-1]]  # last bead of list
#    tri_m = gm.L.triangle(gm.ts[-1])  # final triangle of last bead of list
    tm = gm.ts[-1]
    rm = gm.rs[-1]  # final vertex index of last bead
    angm = tuple(map(ConstrAngleLength, gm.geodesicangles[-1]))
    if len(is_G) == 1:
        currentsum = ConstrAngleLength(gm.geodesiclength)
        if currentsum == onepi: piflag = (zerolength, zerolength)
#    T = gm.T
#    I = gm.I
#    sqedgelens = self.sqedgelens
    # (base case) Check if first vertex of gs = final vertex of gs (currentsum ensured <2pi)
    g0 = gs[is_G[0]]
#    tri_0 = g0.L.triangle(g0.ts[0])
    t0 = g0.ts[0]
    d0 = g0.ds[0]
#    if tri_m.vertex(rm) == tri_0.vertex(d0):
    if link_triangles.vertices[tm][rm] == link_triangles.vertices[t0][d0]:
        ang0 = tuple(map(ConstrAngleLength, g0.geodesicangles[0]))
#        geodangs = tuple(geodangle(T, I, sqedgelens, tri_m, tri_0,
#                                   rm, (rm+i)%3)+(i,) for i in (1,2))
        geodangs = tuple(link_triangles.geodesicAngle(tm, t0, rm, (rm+i)%3)+(i,) for i in (1,2))
        # 2-i instead of i-1 because g.geodesicangles uses vertex indices,
        # while geodangle(..) uses edge indices
        angs = tuple(abs(angm[2-i]) + ang + abs(ang0[int(r0==(d0+1)%3)])
                     for ang, r0, i in geodangs)
        if piflag and len(is_G) == 1:
            # Special case: there is only one gallery and it has length pi
            #top_i = int(angm[0].sgn == 1)  # index of top angle for gm in `angs` (might not be the case for g0 if twisted, but that's fine)
            #newpiflag = tuple(max(zerolength, onepi-angs[(top_i+j)%2])
            #                  for j in (0,1))
            # order of top angle no longer matters here
            newpiflag = tuple(max(zerolength, onepi-angs[j]) for j in (0,1))
            # Check if there is a twist.
            gdr = ConstrAngleLength(g0.geodesicrange)
            # if it is twisted, then there must exist some \theta such that \theta + (gdr-\theta) \geq piflag[i] \forall i.
            # if it is not twisted, it suffices to show that 2*gdr\geq piflag[0]+piflag[1] (since by selecting \theta=piflag[0]/2 one satisfies 2\theta\geq piflag[0], 2*(gdr-\theta)\geq piflag[1]).
            twisted = angm[0].sgn != ang0[int(geodangs[1][1]==(d0+1)%3)].sgn
            assert not twisted  # cannot be immersed in 2-sphere
            condition = (all(gdr >= p for p in newpiflag)
                         if twisted else
                         (gdr + gdr >= newpiflag[0] + newpiflag[1]))
            if DEBUG: print('singlepi', twisted, condition)
        elif gm.geodesiclength == pi:
            top_i = int(angm[0].sgn == 1)  # index of top angle for gm in `angs`
            newpiflag = tuple(max(piflag[j], onepi-angs[(top_i+j)%2])
                              for j in (0,1))
            condition = (ConstrAngleLength(gm.geodesicrange)
                         >= newpiflag[0] + newpiflag[1])
            if DEBUG: print('leftpi', condition)
        elif g0.geodesiclength == pi:
            top_i = int(ang0[int(geodangs[1][1]==(d0+1)%3)].sgn == 1)  # index of top angle for g0 in `angs`
            newpiflag = tuple(max(piflag[j], onepi-angs[(top_i+j)%2])
                              for j in (0,1))
            condition = (ConstrAngleLength(g0.geodesicrange)
                         >= newpiflag[0] + newpiflag[1])
            if DEBUG: print('rightpi', condition)
        else:  # both end beads have line geodesics
            condition = all(a >= onepi for a in angs)
            #if DEBUG: print('line', condition)
        return [is_G] if condition else []
    # (recursive case) Search for candidate bead to add to gallery list.
    if findall: necklaces = []  # else return a zero- or one-tuple/list
#    for i_G in iss_G[tri_m.vertex(rm).index()]:
    for i_G in iss_G[link_triangles.vertices[tm][rm]]:
        # Check that the bead is new.
        if i_G in is_G: continue
        # Check that there is at most one bead of length pi.
        g0 = gs[i_G]
        if piflag and g0.geodesiclength == pi: continue
#        tri_0 = g0.L.triangle(g0.ts[0])
        t0 = g0.ts[0]
        d0 = g0.ds[0]
        # Check that the total length is less than 2pi.
        newsum = currentsum + ConstrAngleLength(g0.geodesiclength)
        if newsum >= twopi: continue
        # Check that no candidate triangles appear in the current list.
        ts = tuple(t for i in is_G for t in gs[i].ts)
        if any(t in ts for t in g0.ts): continue
        # Compute upper and lower geodesic angles at the vertex.
        ang0 = tuple(map(ConstrAngleLength, g0.geodesicangles[0]))
#        geodangs = tuple(geodangle(T, I, sqedgelens, tri_m, tri_0,
#                                   rm, (rm+i)%3)+(i,) for i in (1,2))
        geodangs = tuple(link_triangles.geodesicAngle(tm, t0, rm, (rm+i)%3)+(i,) for i in (1,2))
        angs = tuple(abs(angm[2-i]) + ang + abs(ang0[int(r0==(d0+1)%3)])
                     for ang, r0, i in geodangs)
        if gm.geodesiclength == pi:  # pi geodesic on the left
            top_i = int(angm[0].sgn == 1)  # index of top angle in `angs`
            newpiflag = tuple(max(piflag[j], onepi-angs[(top_i+j)%2])
                              for j in (0,1))
            condition = (ConstrAngleLength(gm.geodesicrange)
                         >= newpiflag[0] + newpiflag[1])
        elif g0.geodesiclength == pi:  # pi geodesic on the right
            top_i = int(ang0[int(geodangs[1][1]==(d0+1)%3)].sgn == 1)  # index of top angle for g0 in `angs`
            newpiflag = tuple(max(zerolength, onepi-angs[(top_i+j)%2])
                              for j in (0,1))
            condition = (ConstrAngleLength(g0.geodesicrange)
                         >= newpiflag[0] + newpiflag[1])
        else:  # both beads have line geodesics
            newpiflag = piflag
            condition = all(a >= onepi for a in angs)
        if condition:
            loop = find_necklaces(gs, link_triangles, iss_G, findall, newpiflag, newsum, is_G+(i_G,))
            if findall: necklaces.extend(loop)
            elif loop: return loop
    return necklaces if findall else []


def link_trilengths(T, I, sqedgelens, t):
    """Returns the length between two vertices in the vertex link.
    
    This corresponds to the angle between two tetrahedra edges coming out of a
    vertex.
    
    https://regina-normal.github.io/engine-docs/classregina_1_1Face_3_013_00_010_01_4.html
    'Moreover, p->facePerm(i) will indicate exactly where the ith triangle sits
    within tet: it will send 3 to the vertex of t(et?) that the triangle links,
    and it will send 0,1,2 to the vertices of tet that are parallel to vertices
    0,1,2 of this triangle.'
    
    Args:
        triIndex: The index of the triangle in the vertex link, representing a
                  tetrahedron in the triangulation.
        a, b (symmetric): The indices of the two vertices in the vertex link.
        data: Triangulation (T), vLink-triangulation map (I), squares of edge
              lengths of tetrahedra using Regina's numbering (edgelengths2)
    """
    tet = T.tetrahedron(I.tetImage(t))
    fp = I.facePerm(t)
    fpd2 = {}
    for i in range(3):
        for j in range(i+1, 4):
            fpd2[(i,j)] = sqedgelens[tet.edge(triangularindex(fp[i], fp[j]))
                                     .index()]
    # tet eAngles assumed to be not reflex
    return tuple(ConstrAngle(1, (fpd2[(i,3)]+fpd2[(j,3)]-fpd2[(i,j)]) / 
                                (2*sqrt(fpd2[(i,3)]*fpd2[(j,3)]))
                             ) for i,j in ((1,2), (0,2), (0,1))
                 )

from functools import partial

def triangularindex(i, j, N=3):
    # Edge ordering of tetrahedron (confirmed)
    # EDGES = [[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]]  (for N=3)
    if not i < j: i, j = j, i
    return sum(N-k for k in range(i)) + j-i-1

# TODO: separate into gallery generating and loop-checking methods
# (to have a way to retain the gallery data)
#print(len(set(id(g.T) for g in G.gs)))
def vertex_linkcond(isosig, sqedgelens,
                    parallel=False, progress=True, mode=''):
    # mode = '', 'cert' or 'findall'
    tqdm_ = tqdm if progress else yield_iter
    T = fromIsoSig(isosig)
    with Pool() if parallel else yieldNone() as p:
        imap_ = p.imap_unordered if parallel else imap
        for i_v, v in enumerate(T.vertices()):
            # Enumerate candidate beads.
            L, I = (v.buildLink(), v.buildLinkInclusion())
#            len_link_tris = len(L.triangles())
            len_link_tris = L.triangles().size()
#            G = Galleries(isosig, i_v, sqedgelens)
            G = []
            link_triangles = LinkTriangles(T, v, sqedgelens)
            print(link_triangles.__dict__)
            for i_tris, Gpart in tqdm_(enumerate(imap_(partial(build_galleries, link_triangles),
                                                       range(len_link_tris))),
                                       total=len_link_tris, desc='Triangle'):
                G.extend(Gpart)
            print(len(G))
            # Ensure that all galleries refer to the same Regina objects.
#            G = loads(dumps(G))
            # Group galleries by the starting link vertex.
#            iss_G = tuple([] for _ in range(len(L.vertices())))
            iss_G = tuple([] for _ in range(L.vertices().size()))
            for i_G, g in enumerate(G):
                iss_G[L.triangle(g.ts[0]).vertex(g.ds[0]).index()].append(i_G)
            # Combine beads into a necklace containing a short geodesic loop.
            len_G = len(G)
            # G.find_necklaces always returns a (possibly empty) list of Galleries
            # TODO: incorporate 'findall' in G.find_necklaces as boolean arg findall=True
            if mode == 'findall': certs = []
            for i, necklaces in tqdm_(enumerate(imap_(partial(find_necklaces,
                                                              G,
                                                              link_triangles,
                                                              iss_G,
                                                              mode=='findall',
                                                              None,
                                                              None),
                                                      ((i,) for i in range(len_G)))),
                                      total=len_G, desc='Gallery'):
                #print(necklaces)
                if necklaces:
                    if mode == 'findall': certs.extend(necklaces)
                    elif mode == 'cert':
                        print()
                        return necklaces[0]
                    else: return False  # vertex link cond. failed (not CAT(0))
    if mode == 'findall': return certs
    return True  # vertex link condition passed


@contextmanager
def yieldNone(): yield None

def yield_iter(iterable, **kwargs):
    for e in iterable:
        yield e

def diy_tqdm(iterable, total=0, desc=''):
    if not total: total = len(iterable)
    progressbar(desc, -1, total)
    for i, e in enumerate(iterable):
        progressbar(desc, i, total)
        yield e
    print()

try: from tqdm import tqdm
except ImportError: tqdm = diy_tqdm

def progressbar(name, i, N, barlength=20):
    # use print() when done
    # i ranges from 0 to N-1
    sys.stdout.write(('\r[%-'+str(barlength)+'s] %s %s/%s ')
                     % ('#'*int(barlength*(i+1)/N), name, (i+1), N))
    sys.stdout.flush()    












