"""Constructible spherical trigonometry in Python.

Author: Asem Wardak

While there exist a range of spherical trigonometry libraries for areas such as
geodesy, these libraries use floating point, or at best rational numbers.
Constructible numbers, on the other hand, have historical and geometric
importance, being those numbers which can be 'constructed' using a ruler and
compass. This implementation thus generalises the exact arithmetic used by
arbitrary-precision rational-number-based implementations.

The mathematical derivations underlying the algorithms in
ConstrSphPoint.ang(..) and ConstrSphPoint.linear(..) are presented in the
file constrsphtrig.ipynb.

Current state:
    Write unit tests
    Some usage examples

Related links:
    https://github.com/bth2008/python_spheric
    https://pypi.org/project/algebraics/#description
    http://abstract.ups.edu/aata/section-constructions.html
    
Dependencies:
    constructible.py (https://github.com/leovt/constructible)
        To use a different underlying constructible number object, simply
        change what constrsphtrig.sqrt refers to after importing. (If
        constructible.py is not available create a file with stub variables
        `sqrt` and `Constructible`.)
"""

import math

from ast import literal_eval
from numbers import Real

# Constructible numbers
from constructible import sqrt as sqrt2
from constructible import Constructible

counter = 0

def simplify(self):
    """Simplifies the internal representation of Constructible(..) numbers.
    
    The algorithm is adapted from Constructible.__str__.
    """
    if not isinstance(self, Constructible):
        return self
    elif not self.b:
        return simplify(self.a)
    elif self.b == 1:
        if not self.a:
            return sqrt2(simplify(self.r))
        else:
            return simplify(self.a) + sqrt2(simplify(self.r))
    else:
        if not self.a:
            return simplify(self.b) * sqrt2(simplify(self.r))
        return simplify(self.a) + simplify(self.b) * sqrt2(simplify(self.r))


#import sage.all
#import sage.rings.qqbar as qqbar
#AA = qqbar.AA

def sqrt(n):
    """Bugfix wrapper for .constructible.sqrt(n).
    
    Fixes ZeroDivisionErrors which occasionally arise when computing
    constructible.sqrt(..) on Constructible(..) numbers whose internal
    representations are not sufficiently simplified.
    """
    if isinstance(n, float): return math.sqrt(n)
#    return AA(n).sqrt()
    return sqrt2(simplify(n))

def acos(cos):
    """Allows the string representation of ConstrAngle to be evaluable."""
    return ConstrAngle(cos=cos)

class ConstrAngle(Real):
    """Constructible angles.
    
    The cosines of constructible angles are constructible numbers.
    This implementation supports treating it as a 2-tuple (sgn, cos).
    
    TODO:
        __rmul__, __radd__ for encoding ConstrAngleLengths as (quot*pi + rem)
    
    Attributes:
        sgn (int): 1 if 0 <= theta <= pi, -1 otherwise.
        cos (Constructible): cos(theta).
    """
    
    def __init__(self, sgn=1, cos=sqrt(1), exact=True):
        # Regard -acos(1) as equal to acos(1).
        self.sgn = 1 if cos == 1 else sgn
        # Keep the internal representation as simple as possible (expected to
        # be faster than comparing different representations of the same
        # constructible number with complex structures)
        self.cos = ((-1 if cos < 0 else 1) * sqrt(cos*cos)
                    if exact else float(cos))
        self.exact = exact
    
    @property
    def sin(self):
        """Sine of a constructible angle."""
        cos = self.cos
        return self.sgn * sqrt(1-cos*cos)
    
    @property
    def cos_sin(self):
        """Returns (cos, sin); analogous to AngleVector[] in Wolfram Lang."""
        return self.cos, self.sin
    
    def bisect(self, other):
        """Returns the closest midpoint of two angles.
        
        The algorithm is set up so that the intermediate quantities computed
        are always between -pi and (inclusively) pi.
        """
        min_theta, max_theta = sorted((self, other))
        return max_theta + (min_theta-max_theta)/2
    
    def __getitem__(self, items):
        if items == 0: return self.sgn
        elif items == 1: return self.cos
        raise IndexError()
    
    def __iter__(self):
        return iter((self.sgn, self.cos))

    def __repr__(self):
        return ('-' if self.sgn == -1 else '') + 'acos(%s)' % self.cos
    
    def __add__(self, other):
        """Adds two constructible angles, where pi < self + other <= pi."""
        cos1, sin1 = self.cos_sin
        cos2, sin2 = other.cos_sin
        return ConstrAngle(-1 if sin1*cos2 + cos1*sin2 < 0 else 1,
                           cos1*cos2 - sin1*sin2, self.exact)
    
    def __truediv__(self, other):
        """Divides by 2, where -pi <= self <= pi."""
        if other == 2: return ConstrAngle(self.sgn, sqrt((1+self.cos)/2), self.exact)
        return NotImplemented
    
    def __neg__(self):
        return ConstrAngle(-self.sgn, self.cos, self.exact)
    
    def __lt__(self, other):
        """Returns self < other, where -pi <= self, other <= pi."""
        if self.sgn != other.sgn: return other.sgn == 1
        if self.sgn == 1: return other.cos < self.cos  # both in top hemisphere
        return self.cos < other.cos  # both in bottom hemisphere
    
    def __le__(self, other):    
        return self == other or self < other
    
    def __eq__(self, other):
        if not isinstance(other, ConstrAngle): return NotImplemented
        return (self.sgn == other.sgn) and (self.cos == other.cos
                                            if self.exact else
                                            math.isclose(self, other))
    
    def __abs__(self):
        return ConstrAngle(1, self.cos, self.exact)
    
    def __float__(self):
        return self.sgn*math.acos(self.cos)
    
    # Satisfy the numbers.Real abstract base class requirements.
    # This is safe since the executed strings are static.
    for f in """__div__, __floordiv__, __mod__, __mul__, __pos__,
                __pow__, __radd__, __rdiv__, __rfloordiv__, __rmod__, __rmul__,
                __rpow__, __rtruediv__, __trunc__,
                __ceil__, __floor__, __round__""".split():
        exec(f.strip(',') + ' = lambda *args: NotImplemented')
    del f

pi = ConstrAngle(1, -sqrt(1))
zero = ConstrAngle()
minuspi = ConstrAngle(-1, -sqrt(1))  # special value for comparison purposes
    

class ConstrAngleLength(Real):
    """Large constructible angle lengths.
    
    One can represent such lengths as sgn * (pi, ..., pi, ConstrAngle(1, ..)).
    Here, they are internally represented as sgn * (quot*pi + rem).
    The winding number is analogous to half this quotient.
    
    Constructors:
        ConstrAngleLength(ConstrAngle(..))
        ConstrAngleLength(sgn, quot, rem)
        TODO: (quot*pi + rem)
    
    Attributes:
        sgn (int): -1 if length is negative, 1 otherwise.
        quot (int): quotient; numbers of pi to add to the absolute length.
        rem (ConstrAngle(1,cos)): remainder; final value to add to the
                                  absolute length.
    """

    def __init__(self, *args):
        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, ConstrAngle):
                self.sgn = arg.sgn
                if abs(arg) == pi:
                    self.quot = 1
                    self.rem = zero
                else:
                    self.quot = 0
                    self.rem = abs(arg)
            else: raise TypeError('ConstrAngleLength must be built from \
                                   ConstrAngle')
        else:
            self.sgn, self.quot, self.rem = args
            if self.quot == 0 and self.rem == zero: self.sgn = 1
    
    def __repr__(self):
        return (('-' if self.sgn == -1 else '')
                + '(%s*pi + %s)' % (self.quot, self.rem))
    
    def __add__(self, other):
        if self.sgn == other.sgn:
            rem = self.rem + other.rem
            r = int(rem < zero or rem == pi)
            return ConstrAngleLength(self.sgn,
                                     self.quot + other.quot + r,
                                     pi + rem if r else rem)
        # signs differ
        if self.quot == other.quot:
            rem = self.rem - other.rem
            return ConstrAngleLength(self.sgn*rem.sgn, 0, abs(rem))
        # quotients differ too
        a, b = (self, other) if self.quot < other.quot else (other, self)
        rem = b.rem - a.rem
        r = int(rem < zero or rem == pi)
        return ConstrAngleLength(b.sgn,
                                 b.quot - a.quot - r,
                                 pi + rem if r else rem)
    
    def __neg__(self):
        return ConstrAngleLength(-self.sgn, self.quot, self.rem)
    
    def __lt__(self, other):
        self_sgn = self.sgn
        if self_sgn != other.sgn:
            return self_sgn == -1
        if self.quot != other.quot:
            return (self.quot < other.quot if self_sgn == 1 else
                    other.quot < self.quot)
        return (self.rem < other.rem if self_sgn == 1 else
                other.rem < self.rem)
    
    def __le__(self, other):    
        return self == other or self < other
    
    def __eq__(self, other):
        return (self.rem == other.rem and
                self.quot == other.quot and
                self.sgn == other.sgn)
    
    def __abs__(self):
        return ConstrAngleLength(1, self.quot, self.rem)
        
    # Satisfy the numbers.Real abstract base class requirements.
    for f in """__div__, __float__, __floordiv__, __mod__, __mul__, __pos__,
                __pow__, __radd__, __rdiv__, __rfloordiv__, __rmod__, __rmul__,
                __rpow__, __rtruediv__, __truediv__, __trunc__,
                __ceil__, __floor__, __round__""".split():
        exec(f.strip(',') + ' = lambda *args: NotImplemented')
    del f

twopi = ConstrAngleLength(1, 2, zero)  # or onepi+onepi
onepi = ConstrAngleLength(1, 1, zero)  # or ConstrAngleLength(pi)
zerolength = ConstrAngleLength(zero)


class ConstrSphPoint(object):
    """Constructible spherical points.
    
    Points on the sphere are represented by pairs (lat, lon) indicating the
    latitude and longitude of the point in spherical coordinates.
    
    Attributes:
        lat (ConstrAngle): latitude, -pi/2 <= lat <= pi/2
        lon (ConstrAngle): longitude, -pi < lon <= pi
    """
    
    def __init__(self, lat=zero, lon=zero):
        self.lat = lat
        self.lon = lon
    
    @property
    def hemi(self):
        """Returns the hemisphere a point lies on.
    
        Returns:
            1 if northern hemisphere,
            0 if equator, or
            -1 if southern hemisphere.
        """
        if self.lat == zero: return 0
        return self.lat.sgn  # 1 if self.lat > zero else -1
    
    @property
    def cart(self):
        """Cartesian coordinates of a point, preserving right-handedness."""
        phi = self.lat
        theta = self.lon
        cos_phi = phi.cos
        return (cos_phi*theta.cos, cos_phi*theta.sin, phi.sin)
    
    def ang(self, lon_e=zero):
        """Returns the angle between self and the equator at a given longitude.
        
        Returns the angle self(zero,lon_e)(zero,lon_e+pi/2). If self is below
        the equator, the angle returned is negative. This algorithm assumes
        that (zero_, lon_e) is not antipodal to a,
        i.e. a != (zero_, lon_e + pi). The algorithm is optimised for the
        special case of the angle being computed about the equator.
        
        See derivation file for details.
        
        Args:
            lon_e: The longitude of the point on the equator about which the
                   angle is computed.
        """
#        if self in (ConstrSphPoint(lon=lon_e), ConstrSphPoint(lon=lon_e+pi)):
#            raise ValueError(self, lon_e)
        assert self not in (ConstrSphPoint(lon=lon_e), ConstrSphPoint(lon=lon_e+pi))
        lat, lon = self
        D_lon = lon - lon_e
        cos_ang = lat.cos * D_lon.sin / ConstrAngle(cos=lat.cos*D_lon.cos).sin
        return ConstrAngle(lat.sgn, cos_ang)
        
    def rotate(self, theta, axis=None, forwards=False):
        """Returns the point rotated about an axis point on the sphere.
        
        Links:
            https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
        
        Args:
            theta: The rotation amount.
            axis: The point about which to rotate self; defaults to zero_zero
            forwards: Whether the point is rotated in an anticlockwise
                      direction (right-hand rule). By default, rotate so as to
                      'snap' the point back onto the equator using the output
                      of ConstrSphPoint.ang(..) as theta.
        """
#        print('rotate')
        if not axis: axis = zero_zero
        cos_theta = theta.cos
        sin_theta = ConstrAngle(theta.sgn if forwards else -theta.sgn,
                                cos_theta).sin
        v = self.cart
        k = axis.cart
        kXv = tuple(k[(i+1)%3]*v[(i+2)%3] - k[(i+2)%3]*v[(i+1)%3]
                    for i in range(3))  # cross product
        # Rodrigues formula RHS, 3rd term scalar coefficient
        v_rot3 = sum(k[i]*v[i] for i in range(3)) * (1-cos_theta)
        v_rot = tuple(v[i]*cos_theta + kXv[i]*sin_theta + k[i]*v_rot3
                      for i in range(3))  # Rodrigues formula
        cos_lat_rot = sqrt(1 - v_rot[2]*v_rot[2])
#        print(v_rot[2])
        lat = ConstrAngle(-1 if v_rot[2]<0 else 1, cos_lat_rot)
        lon = ConstrAngle(-1 if v_rot[1]<0 else 1, # default to lon=0 if pole
                          v_rot[0]/cos_lat_rot if cos_lat_rot else 1)
        return ConstrSphPoint(lat, lon)
                

    def linear(a, b, d_ac, d_bc):
        """Solves the linear problem in spherical geometry.
        
        Given two points of a spherical triangle and distances from these
        points to the third point, this function returns the list of possible
        third points c.
        
        See derivation file for details.
        
        For symmetry `a` is used instead of `self`.
        
        Args:
            a, b: The two points.
            d_ac, d_bc: The distances ac, bc.
        """
        if not b.lat.cos:  # if b is on a pole
            a,b = (b,a)
            d_ac,d_bc = (d_bc,d_ac)
        x_a, y_a, z_a = a.cart
        x_b, y_b, z_b = b.cart
        dp = (d_ac.cos, d_bc.cos)
        detM = x_a*y_b-x_b*y_a
        if not detM:  # edge case: either a.lat.cos==0 or (b.lon-a.lon).sin==0
            assert a.lat.cos==0 or (b.lon-a.lon).sin==0
            m = 1 if a.lon == b.lon else -1
            detN = a.lat.sin*b.lat.cos - m*a.lat.cos*b.lat.sin
            # 26/12/20 1:20am: the m* bit took me 3 days to find and insert;
            #                  a true needle in a haystack
            sin_lat_c = (b.lat.cos*dp[0] - m*a.lat.cos*dp[1])/detN
            if sin_lat_c in (-1,1):  # c is on a pole
                return ConstrSphPoint(lat=ConstrAngle(-1 if sin_lat_c<0 else 1,
                                                      0)
                                      )
            cos_lat_c = sqrt(1-sin_lat_c*sin_lat_c)  # guaranteed != zero
            lat_c = ConstrAngle(-1 if sin_lat_c < 0 else 1, cos_lat_c)
            cos_lon_bc = (a.lat.sin*dp[1] - b.lat.sin*dp[0])/(detN*cos_lat_c)
            lons_bc = (ConstrAngle(sgn, cos_lon_bc) for sgn in 
                       ((1,) if cos_lon_bc in (-1,1) else (1,-1)))
            return tuple(ConstrSphPoint(lat_c, b.lon - lon_bc)
                         for lon_bc in lons_bc)
        # general algorithm
        zp = (z_a, z_b)
        alpha__d = tuple(sum(v[i]*dp[i] for i in range(2))# (alpha_yd, alpha_xd)
                         for v in ((y_b,-y_a), (-x_b,x_a)))
        alpha__z = tuple(sum(v[i]*zp[i] for i in range(2))# (alpha_yz, alpha_xz)
                         for v in ((y_b,-y_a), (-x_b,x_a)))
        alpha2 = sum(al*al for al in alpha__z) + detM*detM
        alpha1 = sum(ald*alz for ald, alz in zip(alpha__d, alpha__z))
        alpha0 = sum(al*al for al in alpha__d) - detM*detM
        Delta = alpha1*alpha1 - alpha2*alpha0
        if Delta < 0: return ()
        zs_c = (((alpha1+sqrt(Delta))/alpha2, (alpha1-sqrt(Delta))/alpha2)
                if Delta else (alpha1 / alpha2,))
        xs_c, ys_c = ([sum(v[i]*(dp[i]-z_c*zp[i]) for i in range(2)) / detM
                       for z_c in zs_c] for v in ((y_b,-y_a), (-x_b,x_a)))
        coss_lat_c = tuple(sqrt(1-z_c*z_c) for z_c in zs_c)
        lats_c = tuple(ConstrAngle(-1 if z < 0 else 1, cos_lat)
                       for z, cos_lat in zip(zs_c, coss_lat_c))
        lons_c = tuple(ConstrAngle(-1 if y < 0 else 1,
                                   x / cos_lat if cos_lat else 1)
                       for x, y, cos_lat in zip(xs_c, ys_c, coss_lat_c))
        return tuple(ConstrSphPoint(lat_c, lon_c)
                     for lat_c, lon_c in zip(lats_c, lons_c))
    
    def dist(self, other, exact=True):
        if not exact:  # dist can be expensive for exact values
            self_lon = ConstrAngle(*self.lon, exact=False)
            other_lon = ConstrAngle(*other.lon, exact=False)
            return ConstrAngle(1,
                float(self.lat.sin)*float(other.lat.sin)
                    + float(self.lat.cos)*float(other.lat.cos)*(self_lon-other_lon).cos,
                exact=False)
        # https://en.wikipedia.org/wiki/Great-circle_distance
        return ConstrAngle(1,
            self.lat.sin*other.lat.sin
                + self.lat.cos*other.lat.cos*(self.lon-other.lon).cos)
    
    def __getitem__(self, items):
        if items==0: return self.lat
        elif items==1: return self.lon
        raise IndexError()
    
    def __iter__(self):
        return iter((self.lat, self.lon))
    
    def __repr__(self):
        return repr(tuple(self))
    
    def __eq__(self, other):
        return self.lat == other.lat and self.lon == other.lon

zero_zero = ConstrSphPoint()
zero_pi = ConstrSphPoint(lon=pi)



if __name__ == '__main__':
    a = ConstrSphPoint()
    
