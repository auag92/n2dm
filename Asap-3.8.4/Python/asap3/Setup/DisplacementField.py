# DisplacementField.py  --  Implements displacement fields

# Copyright (C) 1999 CAMP and Jakob Schiotz <schiotz@fysik.dtu.dk>

"""Defines the DisplacementField class.

Displacement fields are returned e.g. by the `Setup.Dislocation`
modules, they are vector fields (see `VectorFunction`) and can be
applied to a lattice using the `apply_to()` method.  Two displacement
fields can be added or subtracted, and a field can be
multiplied/divided by a scalar.  """

__docformat__ = "restructuredtext en"

import numpy

class VectorFunction:
    """Implements vector function objects.

    Function objects are callable objects (i.e. they behave like
    functions), but they can be added and multiplied with scalars
    and with themselves.
    """
    
    def __init__(self,expr):
        self.__call__=expr

    def __add__(self,other):
        if callable(other):
            return self.__class__(lambda x,other=other,oldfunc=self : oldfunc(x)+other(x))
        else:
            return self.__class__(lambda x,const=other,oldfunc=self : oldfunc(x)+const)

    def __radd__(self, other):
        # Addition is commutative
        return self + other

    def __sub__(self,other):
        if callable(other):
            return self.__class__(lambda x,other=other,oldfunc=self : oldfunc(x)-other(x))
        else:
            return self.__class__(lambda x,const=other,oldfunc=self : oldfunc(x)-const)

    def __rsub__(self, other):
        if callable(other):
            return self.__class__(lambda x,other=other,oldfunc=self : other(x)-oldfunc(x))
        else:
            return self.__class__(lambda x,const=other,oldfunc=self : const-oldfunc(x))
        
    def __neg__(self):
        return self.__class__(lambda x,oldfunc=self : -oldfunc(x))

    def __mul__(self,other):
        if callable(other):
            return self.__class__(lambda x,other=other,oldfunc=self : oldfunc(x)*other(x))
        else:
            return self.__class__(lambda x,const=other,oldfunc=self : oldfunc(x)*const)
    
    def __rmul__(self,other):
        # Multiplication is commutative
        return self * other

    def __div__(self,other):
        if callable(other):
            return self.__class__(lambda x,other=other,oldfunc=self : oldfunc(x)/other(x))
        else:
            return self.__class__(lambda x,const=other,oldfunc=self : oldfunc(x)/const)
    
    def __rdiv__(self,other):
        if callable(other):
            return self.__class__(lambda x,other=other,oldfunc=self : other(x)/oldfunc(x))
        else:
            return self.__class__(lambda x,const=other,oldfunc=self : const/oldfunc(x))

class DisplacementField(VectorFunction):
    """A displacement field.

    Displacement fields are returned e.g. by the `Setup.Dislocation`
    modules, they are vector fields (see `VectorFunction`) and can be
    applied to a lattice using the `apply_to()` method.  The
    DisplacementField instance can also be used as a function taking
    an N x 3 array of positions and returning an N x 3 array of
    displacements.

    To create a displacement field, define a function taking an N x 3
    array of positions and returning an N x 3 array of displacements,
    then create the DisplacementField with the function as the sole
    argument.  A 'lambda' construct will often be useful.  For an
    example, see the source code of `Setup.Dislocation.Screw`.
    
    """

    def apply_to(self, list, fixedbox=0, usepositions=None):
        """Apply the displacement field to a list of atoms.

        The list of atoms object is modified to reflect the displacement
        field.

        Two optional arguments are possible.

        fixedpox (default: 0): Adjust the displacement field to keep
        the corners of the computational box fixed by adding a smooth
        displacement field to the specified field.

        usepositions (default: None): If given, it must be an array
        containing the Cartesian positions of all atoms, to be used
        instead of the positions extracted from the list of atoms.  It
        can be useful if several displacement fields are to be applied
        to the same system, by giving the original positions as this
        argument the fields can be applied sequentially while
        obtaining an effect similar to adding the fields and then
        applying the sum.  This is mainly useful if the displacement
        field contains discontinuities, where a previously applied
        field (or dynamics) may cause some atoms to cross the
        discontinuity.
        """
        if fixedbox:
            box = list.get_cell()
            correction = _fixboxfield(self, box)
            field = self + correction
        else:
            field = self
        oldr = list.get_positions()
        if usepositions is None:
            usepositions = oldr
        else:
            if usepositions.shape != oldr.shape:
                raise ValueError, ("The optional argument usepositions has " +
                                   "the wrong shape: "+str(usepositions.shape))
        u = field(usepositions)
        list.set_positions(oldr + u)

def _fixboxfield(field, box):
    """Create a smooth displacement used to keep the computational box fixed.

    INTERNAL USE ONLY!  This function creates a smooth displacement field
    that can be added to an already existing field to keep the corners of the
    computational box fixed.
    """
    def cfield(r, invbox, k0, k1, k2, k3, k12, k13, k23, k123):
        """The correction field."""
        NewAxis = numpy.newaxis 
        (x1, x2, x3) = numpy.transpose(numpy.dot(r, invbox))
        return (k0 + k1*x1[:,NewAxis] + k2*x2[:,NewAxis] + k3*x3[:,NewAxis]
                + k12*(x1*x2)[:,NewAxis] + k13*(x1*x3)[:,NewAxis] +
                k23*(x2*x3)[:,NewAxis] + k123*(x1*x2*x3)[:,NewAxis])
    
    rs = numpy.zeros((2,2,2,3), float)
    for i in (0,1):
        for j in (0,1):
            for k in (0,1):
                rs[i,j,k,0] = i
                rs[i,j,k,1] = j
                rs[i,j,k,2] = k
    rs.shape = (8,3)                    # Reshape to a list of coordinates
    r = numpy.dot(rs, box) # Real space coordinates
    u = - field(r)                      # The desired corrections
    u.shape = (2,2,2,3)

    k0 = u[0,0,0]
    k1 = -u[0,0,0] + u[1,0,0]
    k2 = -u[0,0,0] + u[0,1,0]
    k3 = -u[0,0,0] + u[0,0,1]
    k12 = u[0,0,0] - u[1,0,0] -u[0,1,0] + u[1,1,0]
    k13 = u[0,0,0] - u[1,0,0] -u[0,0,1] + u[1,0,1]
    k23 = u[0,0,0] - u[0,1,0] -u[0,0,1] + u[0,1,1]
    k123 = (- u[0,0,0] + u[1,0,0] + u[0,1,0] + u[0,0,1]
            - u[1,1,0] - u[1,0,1] - u[0,1,1] + u[1,1,1])
    invbox = LinearAlgebra.inverse(box)

    return DisplacementField(lambda r, k0=k0, k1=k1, k2=k2, k3=k3, k12=k12,
                             k13=k13, k23=k23, k123=k123, invbox=invbox,
                             f=cfield:
                             f(r, invbox, k0, k1, k2, k3, k12, k13, k23, k123))

