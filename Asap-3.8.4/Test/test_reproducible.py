"""Usage: python test_reproducible.py [init]

Used to test if the EMT potential still produces the same energies, forces
and stresses as it used to.

Call it the first time with the "init" argument to generate a data file.

Subsequent calls calculate the same quantities, but compare them with the
database.
"""

from numpy import *
from asap3 import *
import cPickle, os, sys, time

# PrintVersion(1)


def main():
    datafile = "test_reproducible2.dat"
    is_data = os.access(datafile, os.R_OK)
    initialize = 0
    if len(sys.argv) == 2 and sys.argv[1] == "init":
        initialize = 1
        if os.access(datafile, os.R_OK):
            print "ERROR: Data file already exists.  Delete it if you want to make a new one."
            print ""
            print __doc__
            sys.exit(1)
    elif len(sys.argv) != 1:
        print "ERROR: Wrong command line."
        print ""
        print __doc__
        sys.exit(1)
    elif not os.access(datafile, os.R_OK):
        print "ERROR: No data file."
        print ""
        print __doc__
        sys.exit(1)

    indata = cPickle.load(file("test_reproducible.in"))
    lattice = Atoms(positions=indata['positions'], cell=indata['cell'], symbols=["Cu"]*len(indata['positions']), pbc=True)
    #_v = Verbose(1)
    # Make a lattice
    #lattice = FCCOrtho(((1,-1,0),(1,1,-2),(1,1,1)), (15,8,14), Copper,
    #                   symmetry=(1,1,1), half=1)
    print "Number of atoms", len(lattice)
    print lattice.get_cell()
    # Perturb the positions
    r = lattice.get_positions()
    dr =  0.05 * sin(arange(3*len(lattice)))
    dr.shape = (-1,3)
    r += dr
    lattice.set_positions(r)

    atoms1 = Atoms(lattice)
    atoms1.set_calculator(EMT())

    print "Total energy:", atoms1.get_potential_energy()
    print "Stress:", atoms1.get_stress()
    #CNA(atoms1)
    #plot = atoms1.GetPlot()
    #time.sleep(10)

    z = lattice.get_atomic_numbers()
    z[0] = z[1] = 47
    lattice.set_atomic_numbers(z)
    atoms2 = Atoms(lattice)
    atoms2.set_calculator(EMT())
    
    if initialize:
        Cu = {}
        n = len(atoms1)
        Cu["energies"] = atoms1.get_potential_energies()
        Cu["forces"] = atoms1.get_forces()
        Cu["stresses"] = atoms1.GetStresses()
        assert Cu["energies"].shape == (n,)
        assert Cu["forces"].shape == (n, 3)
        assert Cu["stresses"].shape == (n, 6)
        AgCu = {}
        AgCu["energies"] = atoms2.GetPotentialEnergies()
        AgCu["forces"] = atoms2.GetCartesianForces()
        AgCu["stresses"] = atoms2.GetStresses()
        data = {"Cu": Cu, "AgCu": AgCu}
        cPickle.dump(data, open(datafile, "w"))
    else:
        data = cPickle.load(open(datafile))
        Cu = data["Cu"]
        AgCu = data["AgCu"]
        print "*** Checking pure copper ***"
        e = 0
        e = e + evaluate("energies", Cu["energies"], atoms1.get_potential_energies())
        e = e + evaluate("forces", Cu["forces"], atoms1.get_forces())
        e = e + evaluate("stresses", Cu["stresses"], atoms1.get_stresses())
        print "*** Checking silver in copper ***"
        e = e + evaluate("energies", AgCu["energies"], atoms2.get_potential_energies())
        e = e + evaluate("forces", AgCu["forces"], atoms2.get_forces())
        e = e + evaluate("stresses", AgCu["stresses"], atoms2.get_stresses())
        if e == 6:
            print "*** All tests passed ***"
        else:
            print "*** THERE WERE ERRORS IN SOME TESTS! ***"
            sys.exit(2)

def evaluate(text, expected, actual):
    diff = max(abs(expected.flat[:] - actual.flat[:]))
    passed = diff < 1e-10
    if passed:
        print "Checking %s: max error = %g   OK" % (text, diff)
    else:
        print "Checking %s: max error = %g   FAILED!     <<<<<<<<<<" % (text, diff)
    #print "   average over all atoms:", sum(actual) / len(actual)
    return passed

    
main()
