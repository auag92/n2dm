from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.testtools import *
from numpy import *
import ase.data

def TestLists(nblist, nblist2, name, count=None):
    "Run the tests on a half and a full neighbor list."
    print "Running tests on", name
    print "Testing %s: Lengths of lists." % (name,)
    if count:
        sum = 0
        for lst in nblist:
            sum += len(lst)
        ReportTest("Absolute length of list", sum, count*len(atoms), 0)

    sum1 = sum2 = 0
    for nb in enumerate(nblist):
        sum1 += len(nb)
    for nb in enumerate(nblist2):
        sum2 += len(nb)
        
    ReportTest("Equal total length of lists", sum2, sum1, 0)
    
    print "Testing %s: List 1 atoms on list 2." % (name,)
    for i, nb in enumerate(nblist):
        for jj in nb:
            j = int(jj)
            ReportTest.BoolTest("Atom %d on list %d (forward)" % (j, i),
                                (j in nblist2[i]) != (i in nblist2[j]),
                                silent=True)
        if ReportTest.GetNumberOfErrors() > 10:
            print "*** Too many errors - giving up! ***"
            break
    
    print "Testing %s: List 2 atoms on list 1." % (name,)
    for i, nb in enumerate(nblist2):
        for jj in nb:
            j = int(jj)
            ReportTest.BoolTest("Atom %d on list %d (reverse)" % (j, i),
                                (j in nblist[i]) != (i in nblist[j]),
                                silent=True)
        if ReportTest.GetNumberOfErrors() > 10:
            print "*** Too many errors - giving up! ***"
            break
    

print_version(1)

element = "Cu"
latconst = ase.data.reference_states[ase.data.atomic_numbers[element]]['a']

atoms = FaceCenteredCubic(directions=[[1,0,0],[0,1,1],[0,0,1]], size=(9,7,5),
                          symbol=element, debug=0, pbc=(0,0,0))
print "Symmetry:", atoms.get_pbc()
atoms.set_calculator(EMT())
epot = atoms.get_potential_energy()

nblist1 = NeighborList(latconst * 0.5 * (1/sqrt(2) + 1), atoms, 0.0)
nblist2 = NeighborCellLocator(latconst * 0.5 * (1/sqrt(2) + 1), Atoms(atoms), 0.0)
TestLists(nblist1, nblist2, "nearest-neigbor lists (free)")

atoms = Atoms(atoms, pbc=(1,1,1))
print "Symmetry:", atoms.get_pbc()

nblist1 = NeighborList(latconst * 0.5 * (1/sqrt(2) + 1), atoms, 0.0)
nblist2 = NeighborCellLocator(latconst * 0.5 * (1/sqrt(2) + 1), Atoms(atoms), 0.0)
TestLists(nblist1, nblist2, "nearest-neigbor lists (periodic)", 6)

nblist1 = NeighborCellLocator(4.98409, atoms, 0.0)
nblist2 = NeighborList(4.98409, Atoms(atoms), 0.0)
TestLists(nblist1, nblist2, "long neigbor lists (periodic)", 21)

atoms = Atoms(atoms, pbc=(0,0,0))
print "Symmetry:", atoms.get_pbc()

nblist = NeighborCellLocator(4.98409, atoms, 0.0)
nblist2 = NeighborList(4.98409, Atoms(atoms), 0.0)
TestLists(nblist, nblist2, "long neigbor lists (free)")

atoms1 = Atoms(atoms, pbc=(1,1,1))
atoms2 = Atoms(atoms, pbc=(1,1,1))
nblist = NeighborCellLocator(4.98409, atoms1)
nblist2 = NeighborList(4.98409, atoms2)

x0 = atoms1.get_positions()[2,0]
for i in range(30):
    r = atoms1.get_positions()
    r[2,0] += 0.05
    atoms1.set_positions(r)
    print r[2,0],
    r = atoms2.get_positions()
    r[2,0] += 0.05
    print r[2,0]
    atoms2.set_positions(r)
    nblist.check_and_update(atoms1)
    nblist2.check_and_update(atoms2)
    TestLists(nblist, nblist2, "Translation step "+str(i))

ReportTest.Summary()
