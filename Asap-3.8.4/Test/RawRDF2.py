"Test the RawRadialDistributionFunction function."

from numpy import *
from asap3 import *
from asap3 import _asap
from ase.lattice.compounds import NaCl
from asap3.testtools import *
from asap3.Internal.ListOfElements import ListOfElements

print_version(1)

# testtypes: latticeconst, maxRDF, bins, useEMT
testtypes = ((5.1, 6.001, 100, False),
             (5.2, 6.001, 100, True),
             (5.1, 15.001, 100, False),
             (5.15, 15.001, 100, True))

for latconst, maxrdf, nbins, withemt in testtypes:

    atoms = NaCl(directions=[[1,0,0],[0,1,0],[0,0,1]], symbol=("Cu","Au"),
                 size=(10,10,10), latticeconstant=latconst, debug=0)
    natoms = len(atoms)
    ReportTest("Number of atoms", natoms, 8000, 0)

    if withemt:
        atoms.set_calculator(EMT())
        print atoms.get_potential_energy()

    result = _asap.RawRDF(atoms, maxrdf, nbins, zeros(len(atoms), int32), 1,
                          ListOfElements(atoms))
    z = atoms.get_atomic_numbers()[0]

    globalrdf, rdfdict, countdict = result

    print globalrdf
    localrdf = zeros(globalrdf.shape)
    for i in (29,79):
        for j in (29,79):
            x = rdfdict[0][(i,j)]
            print "LOCAL", i, j
            print x
            localrdf += x


    ReportTest("Local and global RDF are identical",
               min( globalrdf == localrdf), 1, 0)
    ReportTest("Atoms are counted correctly",
               countdict[0][29] + countdict[0][79], natoms, 0)
    
    shellpop = [6, 12, 8, 6, 24, -1]
    shell = [sqrt(i+1.0)/2.0 for i in range(6)]
    print shell
    print shellpop
    n = 0

    dr = maxrdf/nbins

    for i in range(nbins):
        if (i+1)*dr >= shell[n] * latconst:
            if shellpop[n] == -1:
                print "Reached the end of the test data"
                break
            ReportTest(("Shell %d (%d)" % (n+1, i)), globalrdf[i],
                       natoms*shellpop[n], 0)
            n += 1
        else:
            ReportTest(("Between shells (%d)" % (i,)), globalrdf[i], 0, 0)

    
ReportTest.Summary()

