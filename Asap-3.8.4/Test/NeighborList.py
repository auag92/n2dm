from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.testtools import *
from numpy import *
import ase.data

def TestLists(nblist, fnb, name, count=None):
    "Run the tests on a half and a full neighbor list."
    print ""
    if count:
        print "Testing %s: Length of lists" % (name,)
        sum = 0
        for lst in nblist:
            sum += len(lst)
        ReportTest("   Half list", sum, count*len(atoms), 0)

        lfnb = map(len, fnb)
        assert len(lfnb) == len(atoms)
        ReportTest("   Shortest full list", min(lfnb), 2*count, 0)
        ReportTest("   Longest full list", max(lfnb), 2*count, 0)

    print ("Testing %s: Symmetry of full list; full list atoms on half-lists."
           % (name,))
    for i, nb in enumerate(fnb):
        for jj in nb:
            j = int(jj)
            ReportTest.BoolTest("Atom %d on list %d" % (j, i), i in fnb[j],
                                silent=True)
            ReportTest.BoolTest("Exactly one of atoms %d and %d on half-lists"
                                % (j, i),
                                (i in nblist[j]) != (j in nblist[i]),
                                silent=True)
        if ReportTest.GetNumberOfErrors() > 10:
            print "*** Too many errors - giving up! ***"
            break

    print "Testing %s: Half-list atoms on full list." % (name,)
    for i, nb in enumerate(nblist):
        for jj in nb:
            j = int(jj)
            ReportTest.BoolTest("Atom %d on list %d (forward)" % (j, i),
                                j in fnb[i], silent=True)
            ReportTest.BoolTest("Atom %d on list %d (reverse)" % (i, j),
                                i in fnb[j], silent=True)
        if ReportTest.GetNumberOfErrors() > 10:
            print "*** Too many errors - giving up! ***"
            break
    

print_version(1)

element = "Cu"
latconst = ase.data.reference_states[ase.data.atomic_numbers[element]]['a']

atoms = FaceCenteredCubic(directions=[[1,0,0],[0,1,1],[0,0,1]], size=(9,7,5),
                          symbol=element, debug=0)
atoms.set_calculator(EMT(minimum_image=True))
epot = atoms.get_potential_energy()

nblist = atoms.get_calculator().get_neighborlist()
count = {}
for lst in nblist:
    n = len(lst)
    try:
        count[n] += 1
    except KeyError:
        count[n] = 1
# print "Histogram:"
numbers = count.keys()
numbers.sort()
sum = 0
for i in numbers:
    #print i, count[i]
    sum += i*count[i]

ReportTest("Number of neighbors (EMT's NB list)", sum, 21*len(atoms), 0)

nblist = NeighborList(latconst * 0.5 * (1/sqrt(2) + 1), atoms, 0.0)
#nblist = NeighborCellLocator(latconst * 0.5 * (1/sqrt(2) + 1), atoms, 0.0)
fnb = FullNeighborList(latconst * 0.5 * (1/sqrt(2) + 1), Atoms(atoms))
TestLists(nblist, fnb, "nearest-neigbor lists (periodic)", 6)

ReportTest("Energy unperturbed 1", atoms.get_potential_energy(), epot, 1e-11)
atoms.set_positions(atoms.get_positions())
ReportTest("Energy unperturbed 2", atoms.get_potential_energy(), epot, 1e-11)

nblist = NeighborList(4.98409, atoms, 0.0)
fnb = FullNeighborList(4.98409, Atoms(atoms))
TestLists(nblist, fnb, "long neigbor lists (periodic)", 21)

ReportTest("Energy unperturbed 3", atoms.get_potential_energy(), epot, 1e-11)
atoms.set_positions(atoms.get_positions())
ReportTest("Energy unperturbed 4", atoms.get_potential_energy(), epot, 1e-11)

atoms = Atoms(atoms, pbc=(0,0,0))

nblist = NeighborList(latconst * 0.5 * (1/sqrt(2) + 1), atoms, 0.0)
fnb = FullNeighborList(latconst * 0.5 * (1/sqrt(2) + 1), Atoms(atoms))
TestLists(nblist, fnb, "nearest-neigbor lists (non-periodic)")

atoms = Atoms(atoms, pbc=(0,1,0))

nblist = NeighborList(4.98409, atoms, 0.0)
fnb = FullNeighborList(4.98409, Atoms(atoms))
TestLists(nblist, fnb, "long neigbor lists (semi-periodic)")

ReportTest.Summary()
