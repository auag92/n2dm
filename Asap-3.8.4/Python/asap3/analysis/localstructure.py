"""Defines functions and observers analysing the local structure.

The following FUNCTIONS are defined.  They operate on the atoms, and
return an array with the resulting data.  They do not modify the
atoms, but may trigger neighbor list updates in the attached
potential.

CNA: Simple CNA analysis, finds atoms with local FCC and HCP
coordination. 

CoordinationNumbers: Finds the number of neighbors of each atom.

The following CLASSES are defined.  They are can be used as
observers (or base classes for observers) of the atoms by attaching 
their update method to the dynamics.  They can also have their own
observers, so fx a trajectory can be told to write *after* the analysis
is done.

RestrictedCNA: Simple CNA analysis, finds atoms with local FCC and HCP
coordination.  The tags of the atoms are set accordingly.

FullCNA: Full CNA analysis.  Several methods are available for returning
different kind of CNA information.  To use as an observer, the update
method MUST be overridden and augmented to call one of the get_xxx_cna
methods, and use the returned information somehow.

"""

import numpy as np
from asap3 import _asap
import ase.data
from asap3.Internal.ListOfElements import ListOfElements
from asap3.Internal.Subject import Subject

def GuessLatticeConstant(atoms):
    """Gets an estimate of the lattice constant for a system.

    If only a single element is present, the lattice constant for that
    element is returned.

    If more than one element is present, a weighted average of their
    lattice constants is returned.

    In any case, the lattice constant is the 'a' constant if the
    crystal structure is not cubic, except in the case of HCP structure
    where sqrt(2)*a is returned.

    Note that this function is only intended for guessing a cutoff
    radius for e.g. CNA analysis, do not rely on this lattice constant
    being very meaningful.

    If any element does not have a lattice constant, a KeyError is raised.
    """
    refstate = ase.data.reference_states
    z = atoms.get_atomic_numbers()
    if z.min() == z.max():
        # Cast to int is a numpy bug workaround: np.int32 objects are not ints
        if refstate[int(z[0])]['symmetry'] == 'hcp':
            return refstate[int(z[0])]['a'] * np.sqrt(2)
        else:
            return refstate[int(z[0])]['a']
    a = 0.0
    n = 0
    for e in ListOfElements(atoms):
        n1 = np.equal(z, e).sum() # Number of atoms with this element
        if refstate[int(e)]['symmetry'] == 'hcp':
            a += n1 * refstate[int(e)]['a'] * np.sqrt(2)
        else:
            a += n1 * refstate[int(e)]['a']
        n += n1
    return a / n    

def CNA(atoms, rCut=None):
    """Restricted Common Neighbor Analysis on a list of atoms.

    Returns an array of integers (int8) reflecting the local crystal structure.
    FCC is class 0, HCP is class 1 and everything else is class 2.

    It requires a value of rCut between the first and second shell.
    It tries to guess a value for rCut from the atomic numbers if none
    is specified.

    Consider using `Asap.Filters.RestrictedCNA` instead.
    """
    if not rCut:
        rCut = GuessLatticeConstant(atoms) * 0.825
    return _asap.RestrictedCNA(atoms, rCut)


def CoordinationNumbers(atoms, rCut = None):
    """Coordination numbers of the atoms.

    Returns the coordination numbers of the atoms.  The cutoff can be
    given, otherwise it tries to guess a value between the first and
    second coordination shell, based on the atomic numbers.
    """
    if not rCut:
        rCut = GuessLatticeConstant(atoms) * 0.825
    return _asap.CoordinationNumbers(atoms, rCut)

class FullCNA:
    """Full Common Neighbor Analysis on a list of atoms
    
    This CNA calculator can be used as a stand-alone analyzer or
    as an observer.  To calculate the CNA of some atoms once, 
    create the object and call the appropriate get_cna function.
        
    The same object can be used to calculate the CNA multiple times
    when an atoms object is altered (although little if anything is
    saved by reusing the object).  Each time a get_cna_xxx() method
    is called, the CNA are recalculated.  If you use a diffent
    atoms object or want to change the cutoff, call the update
    method.
    
    To use as an observer, derive a subclass and give it a method
    which gets the relevant CNA data and stores it (or whatever 
    the observer should do with the data).
    """
    def __init__(self, atoms=None, rCut=None):
        self.atoms = atoms
        self.rawcna = None
        self.rCut = rCut
        if self.rCut is None and atoms is not None:
            self.rCut = GuessLatticeConstant(self.atoms) * 0.825
        self.tuples = {}
    
    def update(self, atoms=None, rCut=None):
        """Change the atoms or the cutoff."""
        if atoms is not None:
            self.atoms = atoms
        if self.atoms is None:
            raise ValueError("Atoms must be specified when creating the object or when calling update.")
        if rCut is not None:
            self.rCut = rCut
        if not self.rCut:
            self.rCut = GuessLatticeConstant(self.atoms) * 0.825
        
    def get_number_of_pairs(self):
        return len(self.get_raw_cna())

    def get_raw_cna(self):
        """Return raw CNA data as a numpy array."""
        return _asap.FullCNA(self.atoms, self.rCut).get_raw_cna()
    
    def get_cna_tuple(self, cnaint):
        """Convert a CNA integer to a CNA tuple."""
        try:
            return self.tuples[cnaint]
        except KeyError:
            cna = ((cnaint & 0xFF0000) >> 16, (cnaint & 0xFF00) >> 8, cnaint & 0xFF)
            self.tuples[cnaint] = cna
            return cna
    
    def get_normal_cna(self):
        """CNA as a dictionary per atom.
        
        The keys of the dictionaries are the CNA tuples, the values are 
        the number of times they occur.
        """
        return _asap.FullCNA(self.atoms, self.rCut).get_per_atom_cna()
        # The code below implements the same in Python.  It is too slow but is
        # kept for inspiration
#        rawcna = self.get_raw_cna()
#        n = len(self.atoms)
#        result = [{} for i in range(n)]
#        for bond in rawcna:
#            for atom in bond[:2]:
#                assert(atom < n)
#                cna = self.get_cna_tuple(bond[2])
#                try:
#                    result[atom][cna] += 1
#                except KeyError:
#                    result[atom][cna] = 1
#        return result
    
    def get_total_cna(self):
        """Get the total CNA for the system.
        
        A dictionary is returned with the CNA tuples as keys, and the total
        number of bonds with that CNA index as values.
        """
        return _asap.FullCNA(self.atoms, self.rCut).get_total_cna()
        # The code below implements the same in Python.  It is slower but is
        # kept for inspiration
#        rawcna = self.get_raw_cna()
#        temp = {}
#        for x in rawcna[:,2]:
#            try:
#                temp[x] += 1
#            except KeyError:
#                temp[x] = 1
#        result = {}
#        for k, v in temp.iteritems():
#            result[self.get_cna_tuple(k)] = v
#        return result
    
    def get_total_perz_cna(self):
        """Get the total CNA for the system per element.
        
        A dictionary is returned where the keys are the pairs of
        elements encountered in the system, and the keys are
        dictionaries summarizing the CNA (as for get_total_cna())
        for bonds between these elements.  The pairs of elements
        are always given with the smallest atomic number first.
        """

        rawcna = self.get_raw_cna()
        z = self.atoms.get_atomic_numbers()
        temp = {}
        rawcna[:,:2] = z[rawcna[:,:2]]
        # Summarize the data with minimal processing
        for line in rawcna:
            try:
                temp[tuple(line)] += 1
            except KeyError:
                temp[tuple(line)] = 1
        # Now do the processing.
        # First, sum data with z1 > z2 into entries with z1 < z2
        for z1, z2, cna in temp.keys():
            if z2 < z1:
                temp[(z2, z1, cna)] = temp.get((z2, z1, cna), 0) + temp[(z1,z2,cna)]
        result = {}
        for k in temp.iterkeys():
            z1, z2, cna = k
            if (z1 <= z2):
                # Valid entry
                try:
                    d = result[(z1, z2)]
                except KeyError:
                    d = result[(z1, z2)] = {}
                d[self.get_cna_tuple(cna)] = temp[k]
        return result
    
    def get_normal_and_total_cna(self):
        """Return the per-atom and total cna.
        
        Returns a tuple with two elements, the first is
        the same as get_normal_cna() would return, the second
        is what get_total_cna() would return.  However, unlike
        calling the two methods separately, the actual CNA 
        calculation is only done once.
        """
        # This is intended as an example of how to return
        # the CNA in multiple ways without repeating the calculations
        workhorse = _asap.FullCNA(self.atoms, self.rCut)
        normal = workhorse.get_per_atom_cna()
        total = workhorse.get_total_cna()
        return (normal, total)
    
class RestrictedCNA(Subject):
    """Restricted Common Neighbor analysis on a list of atoms.

    Sets the tags of the atoms according to their local crystal
    structure.  FCC is class 0, HCP is class 1 and everything else is
    class 2.
    
    It requires a value of rCut between the first and second shell.
    It tries to guess a value for rCut from the atomic numbers if none
    is specified.
    
    This class is intended to use as an observer, so its analyze()
    function is called automatically by e.g. the dynamics.  It can
    itself act as a subject, so a Plotter or a Trajectory can be
    called just after the calculations.

    Optional parameters:
        
        rCut = None:  The cutoff.  If not given, one is guessed.

        analyze_first = True: Should an analysis be made as soon as
        this object is created?  Leave as true if you plan to attach a
        Trajectory as an observer to this object, and if the
        Trajectory will save the initial state.
        
    """
    def __init__(self, atoms, rCut=None, analyze_first=True):
        Subject.__init__(self)
        self.atoms = atoms
        if rCut:
            self.rCut = rCut
        else:
            self.rCut = GuessLatticeConstant(atoms) * 0.825
        if analyze_first:
            self.analyze()  # There will not be any observers yet.

    def analyze(self):
        "Runs the CNA analysis."
        cna = CNA(self.atoms, self.rCut)
        self.atoms.set_tags(cna)
        self.call_observers()
        
    update = analyze
    
