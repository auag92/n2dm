# AtomsInArrays.py  --  Implements a list of atoms where the data is stored
#                       in multiarrays.

# Copyright (C) 1999 CAMP and Jakob Schiotz <schiotz@fysik.dtu.dk>

"""Implements a list of atoms where the data is stored in multiarrays.

Two versions are implemented, one where all information is stored in
cartesian coordinates, and one where they are stored in scaled space.

In both versions most access is done through __getattr__ to data
specified as arguments to the constructor.  If the object was
constructed with xxxx=array then they can be accessed through
GetXxxx() and SetXxxx().  The exception is positions, momenta and
forces that are accessed through GetCartesianXxxx() and
SetCartesianXxxx().

In Asap version 2.9 and later, this module is only used to create
empty atoms for use in MakeParallelAtoms.
"""

import string
import types
import Numeric
import LinearAlgebra
#import AtomInList
#import DocumentedAttributes
#from Structures.TypesFromAtomicNumbers import TypesFromAtomicNumbers
#from Structures.FakeUnitCell import FakeUnitCell
# Delayed imports:
# Visualization.Avatars.RasMol is imported in AtomsInCmFile.GetPlot()

class ListOfCartesianAtoms:
    """A list of atoms where the storage is in the form of multiarrays.

    The object is created like this:
        ListOfCartesianAtoms(positions=data, momenta=data, .....)
        
    The data is then accessed through the usual access functions, i.e.
    GetCartesianPositions(), GetCartesianMomenta(),
    GetCartesianForces().  Other information can be given when the
    object is initialized, for example stresses can be given as
    stesses=data.  They can then be accessed with GetStresses() and
    changed with SetStresses().
    """

    _specials = {"positions":1, "momenta":1, "forces":1}
    
    def __init__(self, periodicboundaries=None, **kargs):
        # Check that positions are specified, and that they are a 3xN
        # multiarray
        if not kargs.has_key("positions"):
            raise ArgumentError, "positions missing"
        try:
            shape = kargs["positions"].shape
        except AttributeError:
            raise TypeError, "positions should be a multiarray"
        if shape[1] != 3:
            raise ValueError, "positions should be an Nx3 multiarray"

	if periodicboundaries is not None:
	    self.SetPeriodicBoundaries(periodicboundaries)
	
        len = shape[0]
        self._len = len
        self._data = {}
        for arg in kargs.keys():
            if kargs[arg] is not None:
                shape = kargs[arg].shape
                if shape[0] != len:
                    msg = ("%s has shape %t, expected (%d, ?)." %
                           (arg, shape, len))
                    raise ValueError, msg
                self._data[arg] = kargs[arg]

    # List functions
    def __len__(self):
        return self._len

    def __getitem__(self, n):
        if n < 0 or n >= self._len:
            raise IndexError
        return AtomInList.AtomInList(self, n)

    def __getslice__(self, low, high):
        # Select the data to be used in the slice object
        newdata = {}
        for key in self._data.keys():
            newdata[key] = self._data[key][low:high]
        # Now create a new list of atoms object passing the data to the
        # constructor
        return apply(self.__class__, (), newdata)

    def extend(self, other):
	"""Extend the list of atoms by adding atoms from another list.

	The other list must at least contain the same type of data as this
	one.
	"""
	newdata = {}
	for key in self._data.keys():
	    accessfunc = "GetCartesian"+string.capitalize(key)
	    try:
		newdata[key] = getattr(other, accessfunc)()
	    except AttributeError:
		accessfunc = "Get"+string.capitalize(key)
		newdata[key] = getattr(other, accessfunc)()
	for key in self._data.keys():
	    self._data[key] = Numeric.concatenate((self._data[key], newdata[key]))

    def __add__(self, other):
	"""Add two lists of atoms."""
	newdata = {}
	for key in self._data.keys():
	    accessfunc = "GetCartesian"+string.capitalize(key)
	    try:
		otherdata = getattr(other, accessfunc)()
	    except AttributeError:
		accessfunc = "Get"+string.capitalize(key)
		otherdata = getattr(other, accessfunc)()
	    selfdata = self._data[key]
	    newdata[key] = Numeric.concatenate((selfdata, otherdata))
	return apply(self.__class__, (self.GetPeriodicBoundaries(),), newdata)

    def __delitem__(self, n):
	"""Delete an atom."""
	if n < 0 or n >= self._len:
	    raise IndexError, "list deletion index out of range."
	for key in self._data.keys():
	    d = self._data[key]
	    self._data[key] = Numeric.concatenate((d[:n], d[n+1:]))
	self._len = self._len - 1
	
    # We should also implement access functions that change/delete the
    # elements in the list

    # Some access functions are not handled by the generic ones

    def GetPeriodicBoundaries(self):
        """Get the periodic boundary conditions."""
        return self._periodicboundaries

    def SetPeriodicBoundaries(self, pbc):
        """Set the periodic boundary conditions."""
        if not isinstance(pbc, types.TupleType):
            raise TypeError, "Not a tuple"
        if len(pbc) != 3:
            raise ValueError, "Not a tuple of length 3"
        self._periodicboundaries = pbc

    # Visualization
    def GetPlot(self, interactive=0, debug=0):
	"""Get a RasMol avatar plotting the simulation."""
	from Visualization.Avatars.RasMol import RasMol
	return RasMol(self, interactive=interactive, debug=debug)

    # Generic access functions
    def _getdata(self):
        if self._data.has_key(self._attrname):
            return self._data[self._attrname]
        else:
            # Attempting to access a missing "attribute"
            raise AttributeError, self._requested

    def _setdata(self, coordinates):
        a = self._attrname
        if self._data.has_key(a):
            if coordinates.shape != self._data[a].shape:
                raise ValueError, "wrong shape"
        else:
            # Introducing a new sort of data
            if coordinates.shape[0] != self._len:
                raise ValueError, "wrong shape"
        self._data[a] = coordinates

    def _getsingledata(self, n):
        if self._data.has_key(self._attrname):
            return self._data[self._attrname][n]
        else:
            # Attempting to access a missing "attribute"
            raise AttributeError, self._requested

    def _setsingledata(self, n, coordinates):
        if self._data.has_key(self._attrname):
            self._data[self._attrname][n] = coordinates
        else:
            # Attempting to access a missing "attribute"
            raise AttributeError, self._requested

    def _isspecial(self):
        """Checks if access through [GS]etCartesianXXXX is meaningful."""
        if self._attrname == "momentums":
            self._attrname = "momenta"
        if not self._specials.has_key(self._attrname):
            # The access function is not allowed to exist
            raise AttributeError, self._requested
        
    def __getattr__(self, name):
        """Implement the access functions.

        Catches names of the type GetXXXX or GetCartesianXXXX
        or the similar SetXXXX names.
        """
        self._requested = name
        if name[:18] == "GetSingleCartesian":
            self._attrname = string.lower(name[18:])+"s"
            self._isspecial()
            return self._getsingledata
        elif name[:18] == "SetSingleCartesian":
            self._attrname = string.lower(name[18:])+"s"
            self._isspecial()
            return self._setsingledata
        elif name[:9] == "GetSingle":
            self._attrname = string.lower(name[9:])+"s"
            return self._getsingledata
        elif name[:9] == "SetSingle":
            self._attrname = string.lower(name[9:])+"s"
            return self._setsingledata
        if name[:12] == "GetCartesian":
            self._attrname = string.lower(name[12:])
            self._isspecial()
            return self._getdata
        elif name[:12] == "SetCartesian":
            self._attrname = string.lower(name[12:])
            self._isspecial()
            return self._setdata
        elif name[:3] == "Get":
            self._attrname = string.lower(name[3:])
            return self._getdata
        elif name[:3] == "Set":
            self._attrname = string.lower(name[3:])
            return self._setdata
        else:
            raise AttributeError, name

def _lowerplural(str):
    """Returns the plural of an attribute name converted to lower case."""
    try:
        return string.lower(_pluralmap[str])
    except KeyError:
        return string.lower(str+"s")

Attributes = (
    ("Position", "The position."),
    ("Momentum", "The momentum of the atom.", "Momenta"),
    ("Force", "The forces on the atoms."),
    ("Class", "An integer 'classifying' the atoms.", "Classes"),
    ("Type",
     "An object describing the type of the atoms, e.g. the chemical element."),
    ("AtomicNumber", "The atomic number of the atom."),
    ("Energy", "The energy of the atom", "Energies"),
    ("Mass", "The mass of the atom", "Masses"),
    ("Velocity", "The derivative of the position.", "Velocities"))

def _pl(t):
    l = len(t)
    assert l == 2 or l == 3
    if l == 2:
        return (t[0], t[0]+"s")
    elif l == 3:
        return (t[0], t[2])
    else:
        raise AssertionError, "Internal error in _pl, length is not 2 or 3."

Plurals = map(_pl, Attributes)

del _pl

_pluralmap = {}
for _i in Plurals:
    _pluralmap[_i[0]] = _i[1]
del _i


