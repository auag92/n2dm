# Copyright (C) 2005  CAMP
# Please see the accompanying LICENSE file for further information.

"""Subset filter for Asap.

A subset filter will pick out some of the atoms in a ListOfAtoms, and
behave as a smaller ListOfAtoms.

This version has been optimized for large lists of atoms.

It will work with parallel Asap simulations, provided that the
dynamics (or any other object using this filter instead of the
original list of atoms) do not attempt to store local data about the
atoms using the RegisterArray method, as that method is not passed
through this filter.  The RegisterArray method can be extended to work
with this filter, but that will only be done as need arises.

"""

__docformat__ = 'reStructuredText'


import Numeric as num

from ASE.ListOfAtoms import GeneralizedPositionsAndForces
from ASE.Filters.ListOfAtomsFilter import ListOfAtomsFilter

class Subset(GeneralizedPositionsAndForces, ListOfAtomsFilter):
    """Subset filter. for Asap.

    A subset filter will pick out some of the atoms in a ListOfAtoms, and
    behave as a smaller ListOfAtoms.

    This version has been optimized for large lists of atoms.

    It will work with parallel Asap simulations, provided that the
    dynamics (or any other object using this filter instead of the
    original list of atoms) do not attempt to store local data about
    the atoms using the RegisterArray method, as that method is not
    passed through this filter.  The RegisterArray method can be
    extended to work with this filter, but that will only be done as
    need arises.  For this reason, multiple Subset filters cannot be
    chained with parallel Asap simulations.

    """
    def __init__(self, atoms, mask=None, indices=None):
        """Construct a subset filter.

        The subset of the atoms to include can be specified with a
        mask (a list of booleans - one for each atom) or with a list
        of indices of the atoms to include.

        For parallel simulations, the indices are interpreted as local
        indices on this CPU at the time the filter is created.  As
        atoms move, the indices may change or atoms may migrate to
        other processors, but the 'status' of the atom, i.e. if it is
        included or excluded by this filter, follows the atom.
        
        """

        ListOfAtomsFilter.__init__(self, atoms)
        try:
            self._isparallel = self.parallel = atoms.parallel
            # self._isparallel is used internally
            # self.parallel will inform other filters/dynamics that this is a
            #   parallel simulation.
        except AttributeError:
            self._isparallel = False

        self.mask = num.zeros(len(atoms))
        if self._isparallel:
            self.mask = atoms.RegisterArray(self.mask)
            
        if mask is None and indices is not None:
            num.put(self.mask, indices, 1)
        elif mask is not None and indices is None:
            self.mask[:] = num.not_equal(mask, 0) # Force to 1 or 0.
        else:
            raise ValueError, "Either mask or indices must be specified, but not both."

        self.makeidx(force=True)

    def __del__(self):
        if self._isparallel:
            self.atoms.DeregisterArray(self.mask)
            
    def makeidx(self, force=False):
        """Make the index arrays, if necessary."""
        if force or self._isparallel:
            natoms = len(self.atoms)
            self.indices = num.compress(self.mask, num.arange(natoms))
            self.n = len(self.indices)
            self.reverse = num.clip(num.add.accumulate(self.mask)-1, 0, natoms)

    def __len__(self):
        self.makeidx()
        return self.n

    def __getitem__(self, n):
        self.makeidx()
        return self.atoms[self.indices[n]]

    def GetCartesianPositions(self):
        return self.get(self.atoms.GetCartesianPositions)

    def SetCartesianPositions(self, pos):
        self.set(self.atoms.GetCartesianPositions,
                 self.atoms.SetCartesianPositions, pos, 3)

    def GetCartesianMomenta(self):
        return self.get(self.atoms.GetCartesianMomenta)

    def SetCartesianMomenta(self, mom):
        self.set(self.atoms.GetCartesianMomenta,
                 self.atoms.SetCartesianMomenta, mom, 3)

    def GetCartesianVelocities(self):
        return self.get(self.atoms.GetCartesianVelocities)

    def SetCartesianVelocities(self, v):
        self.set(self.atoms.GetCartesianVelocities,
                 self.atoms.SetCartesianVelocities, v, 3)

    def GetCartesianForces(self):
        return self.get(self.atoms.GetCartesianForces)

    def GetMasses(self):
        return self.get(self.atoms.GetMasses)

    def GetAtomicNumbers(self):
        return self.get(self.atoms.GetAtomicNumbers)

    def SetAtomicNumbers(self, z):
        return self.set(self.atoms.GetAtomicNumbers,
                        self.atoms.SetAtomicNumbers, z)

    def GetStresses(self):
        return self.get(self.atoms.GetStresses)

    def GetPotentialEnergies(self):
        return self.get(self.atoms.GetPotentialEnergies)

    def GetTags(self):
        return self.get(self.atoms.GetTags)

    def SetTags(self, t):
        self.set(self.atoms.GetTags, self.atoms.SetTags, t)

    def GetKineticEnergies(self):
        return self.get(self.atoms.GetKineticEnergies)

    def GetKineticEnergy(self):
        return num.sum(self.GetKineticEnergies())

    def SetUnitCell(self, cell, fix=False):
        "Change the unit cell by changing the unit cell of the full list of atoms."
        self.atoms.SetUnitCell(cell, fix)

    def GetUnitCell(self):
        "Get the unit cell, which equals the unit cell of the full list of atoms."
        return self.atoms.GetUnitCell()

    # The hard work is done in the get and set methods
    def get(self, getmethod):
        # Since calling some methods (e.g. GetCartesianForces) may trigger
        # a migration, it must be called before makeidx is called!
        tmp = getmethod()
        self.makeidx()
        return num.take(tmp, self.indices)

    def set(self, getmethod, setmethod, values, dim=0):        
        self.makeidx()
        if dim == 0:
            m = self.mask
        else:
            m = self.mask[:,num.NewAxis] + num.zeros(dim)[num.NewAxis,:]
        newval = num.where(m, num.take(values, self.reverse), getmethod())
        setmethod(newval)
