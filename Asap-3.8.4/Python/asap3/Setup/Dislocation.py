# Dislocation.py  --  Set up a dislocation with mixed character.

# Copyright (C) 2000-2004 CAMP and Jakob Schiotz <schiotz@fysik.dtu.dk>

"""Module for setting up a dislocation with mixed screw and edge components.

The main function is `Dislocation`, it returns the displacement field
(see `Asap.Setup.DisplacementField`) of a dislocation.  The displacement
field is constructed as the sum of the fields of the edge and screw
components, except if one of the components is zero, in which case it
is omitted.  The displacement field can be applied to a lattice.
"""

__docformat__ = "restructuredtext en"


from numpy import *
from asap3.Setup.DisplacementField import DisplacementField

def magn(v):
    """The magnitude of a vector."""
    return sqrt(dot(v, v))

def Dislocation(origin, line, Burgers, PoissonRatio=None, debug=None):
    """Set up a dislocation with mixed screw and edge character.

    In the following description, a 3-vector may be a tuple, list or
    numeric array with three elements.
    
    This function will set up a straight dislocation through a
    specified point with a given Burgers vector.  The first argument
    ('origin', a 3-vector) specifies a point on the dislocation line,
    and the second argument ('line', a 3-vector) gives the dislocation
    line (only the direction of the vector matters). The third
    parameter ('Burgers', a 3-vector) specifies the Burgers vector of
    the dislocation.  Finally, the fourth argument is an optional
    parameter which specifies Poisson's ratio of the material.  The
    default value is 1/3.
    
    Returns a `Setup.DisplacementField.DisplacementField` object.
    Displacement fields may be added together (and multiplied by
    scalars).  The dislplacement field can be applied to any list of
    atoms using field.`apply_to(listofatoms)`.  The list of atoms is
    modified using its GetCartesianPositions() and
    SetCartesianPositions() methods.

    Example:

         from Setup.Lattice.FCC111Ortho import FCC111Ortho
         from Structures.ChemicalElements import *
         from Numeric import *
         from Asap.Setup.Dislocation import Dislocation
         lat = FCC111Ortho((12,12,12), Copper, 6.78)
         center = matrixmultiply(lat.GetCoordinateBasis(), (0.5, 0.5, 0.5)) 
         line = lat.MillerToDirection((0,0.2,1))
         Burgers = lat.MillerToDirection((0.5,0.5,0))
         disloc = Dislocation(center, line, Burgers)
         disloc.apply_to(lat)
         plot = lat.GetPlot()         

    """

    line = array(line) / magn(line)
    Burgers = array(Burgers)
    screw_b = dot(Burgers, line)
    edge_b = Burgers - screw_b * line

    if screw_b == 0.0:
	result = EdgeDislocation(origin, line, edge_b, PoissonRatio,
                                 debug=debug)
    elif magn(edge_b) == 0.0:
	result = ScrewDislocation(origin, line, screw_b, debug=debug)
    else:
	result = (EdgeDislocation(origin, line, edge_b, PoissonRatio,
                                  debug=debug) +
		  ScrewDislocation(origin, line, screw_b, cut=edge_b,
                                   debug=debug))
    return result


# Helper functions, including specialized functions for defining screw
# and edge dislocation.

def ScrewDislocation(origin, line, b, cut=None, debug=None):
    """Set up a straight screw dislocation.

    Set up a straight screw dislocation through a specified point with
    a given burgers vector.  The first argument ('origin', a 3-tuple)
    specifies a point on the dislocation line, the second ('line', a
    3-tuple) gives the dislocation line (only the direction of the
    vector matters).  The third ('b', a scalar) gives the Burgers
    vector (positive means parallel to the dislocation line, negative
    means antiparallel).  The third optional argument (cut, a 3-tuple)
    gives a direction that lies in the cutting plane when starting in
    the origin.  If not given, (1.0, 0.0, 0.0) is used unless parallel
    to the dislocation line, in which case (0.0, 1.0, 0.0) is used.

    Returns a Setup.DisplacementField.DisplacementField object.
    Dislplacement fields may be added together (and multiplied by
    scalars).  The dislplacement field can be applied to any list of
    atoms using 'field.apply_to(listofatoms)'.  The list of atoms is
    modified using its 'GetCartesianPositions()' and
    'SetCartesianPositions()' methods.

    Example:

      >>> from Setup.Lattice.FCC111Ortho import FCC111Ortho
      >>> from Structures.ChemicalElements import *
      >>> from Asap.Setup.Dislocation import ScrewDislocation
      >>> from Numeric import *
      >>> lat = FCC111Ortho((12,12,12), Copper, 3.6)
      >>> center = matrixmultiply(lat.GetCoordinateBasis(), (0.5, 0.5, 0.5))
      >>> disloc = ScrewDislocation(center, lat.MillerToDirection((1,1,1)), 3.6*sqrt(3.0))
      >>> disloc.apply_to(lat)
      >>> plot = lat.GetPlot()
    """

    # Make sure the burgers vector is a float, not a vector.
    if not isinstance(b, float):
        raise TypeError, "The 'Burgers vector' argument must be a number (specifying its length)."
    
    # Normalize line and cut, and assure that they are NumPy arrays
    normline = array(line) / magn(line)
    if cut is not None:
        cut = array(cut) / magn(cut)
    else:
        cut = array((1.0, 0.0, 0.0))
        if parallel(normline, cut):
            cut = array((0.0, 1.0, 0.0))
    assert not parallel(normline, cut)
    origin = array(origin)

    # Turn cut to be orthogonal to normline, without changing the plane
    # defined by cut and normline.
    cut = cut - normline * dot(cut, normline)
    cut = cut / magn(cut)
    return DisplacementField(lambda r, f=_screw, origin=origin, line=normline,
                             b=b, cut=cut, debug=debug:
                             f(r,origin,line,b,cut,debug))

def _screw(r, origin, line, b, cut, debug):
    """The screw dislocation displacement field."""
    if debug:
        print "r: ", r
        print "origin: ", origin
    r = r - origin
    if debug:
        print "r: ", r
    x = vdot(r, cut)
    y = vdot(r, cross(line, cut))
    z = vdot(r, line)         # Distance along line
    if debug:
        print "x: ", x
        print "y: ", y
        print "z: ", z
    x = where(equal(x, 0) * equal(y, 0), ones(x.shape, float), x)
    theta = arctan2(y, x)
    if debug:
        print "theta: ", theta
    return line * (b * theta / (2 * pi))[:,newaxis]

def EdgeDislocation(origin, line, Burgers, PoissonRatio=None, debug=None):
    """Sets up a straight edge dislocation.

    Sets up a straight edge dislocation through a specified point with
    a given Burgers vector.  The first argument ('origin', a 3-tuple)
    specifies a point on the dislocation line, and the second argument
    ('line', a 3-tuple) gives the dislocation line (only the direction of the
    vector matters). The third parameter
    ('Burgers', a 3-tuple) specifies the Burgers vector of the dislocation.
    The Burgers vector must be perpendicular to the line vector. Finally, the
    fourth argument is an optional parameter which specifies Poisson's ratio
    (unlike in the case of a screw dislocation, this number is needed in order
    to calculate the displacement field). The default value of Poisson's ratio
    is 1/3. 
    
    Returns a Setup.DisplacementField.DisplacementField object.
    Dislplacement fields may be added together (and multiplied by
    scalars).  The dislplacement field can be applied to any list of
    atoms using 'field.apply_to(listofatoms)'.  The list of atoms is
    modified using its 'GetCartesianPositions()' and
    'SetCartesianPositions()' methods.

    Example:

         from Setup.Lattice.FCC111Ortho import FCC111Ortho
         from Structures.ChemicalElements import *
         from Numeric import *
         from Asap.Setup.Dislocation import EdgeDislocation
         lat = FCC111Ortho((12,12,12), Copper, 6.78)
         center = matrixmultiply(lat.GetCoordinateBasis(), (0.5, 0.5, 0.5)) 
         line = lat.MillerToDirection((0,0,1))
         Burgers = lat.MillerToDirection((0.5,0.5,0))
         disloc = Edge(center, line, Burgers)
         disloc.apply_to(lat)
         plot = lat.GetPlot()         

    """

    # Normalize line and Burgers, and check their orthogonality
    b = magn(Burgers)
    normBurgers = array(Burgers) / b
    normline = array(line) / magn(line)
    assert perp(normline, normBurgers) 
    origin = array(origin)
    if PoissonRatio is None:
        PoissonRatio=1.0/3.0
    return DisplacementField(lambda r, f=_edge, origin=origin, line=normline,
                             Burgers=normBurgers, b=b,
                             PoissonRatio=PoissonRatio, debug=debug:
                             f(r,origin,line,Burgers,b,PoissonRatio,debug))

def _edge(r, origin, line, Burgers, b, PoissonRatio, debug):
    """The edge dislocation displacement field."""
    if debug:
        print "r: ", r
        print "origin: ", origin
    r = r - origin 
    x = vdot(r, Burgers) 
    y = vdot(r, cross(line, Burgers)) 
    z = vdot(r, line) 
    if debug:
        print "x: ", x 
        print "y: ", y
    theta = arctan2(y, x)
    l = x*x+y*y
    epsilon = 0.1 * b**2
    ux = b / (2 * pi) * (theta + cos(theta) * sin(theta) / (2*(1 - PoissonRatio)))
    uy = b / (2 * pi * 4 * (1 - PoissonRatio)) * ((1-2*PoissonRatio)*log(l+epsilon) + cos(2*theta))
    if debug:
        print "theta: ", theta
        print "ux: ", ux
        print "ux: ", uy
    return  Burgers * ux[:,newaxis]- cross(line, Burgers) * uy[:,newaxis]
    


def vdot(a, b):
    """Dot product of an array of vectors with a vector."""
    # a has shape (N, 3), b has shape (3,) or (N, 3)
    assert a.shape[-1] == 3 and b.shape[-1] == 3
    at = transpose(a)
    bt = transpose(b)
    return at[0] * bt[0] + at[1] * bt[1] + at[2] * bt[2]
    
def magn(v):
    """The magnitude of a vector."""
    return sqrt(dot(v, v))

def cross(a, b):
    """The cross product of two vectors."""
    return array((a[1]*b[2] - b[1]*a[2],
                  a[2]*b[0] - b[2]*a[0],
                  a[0]*b[1] - b[0]*a[1]))

def perp(v1, v2, small=1e-6):
    """Tests if two vectors are perpendicular."""
    if dot(v1,v2) < small:
        return 1
    else:
        return 0


def parallel(v1, v2, small=1e-6):
    """Tests if two vectors are parallel."""
    c = cross(v1, v2)
    if dot(c, c) < small:
        return 1
    else:
        return 0
    
