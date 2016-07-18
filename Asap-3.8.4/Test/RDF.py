"Test the RawRadialDistributionFunction function."

from numpy import *
from asap3 import *
from asap3.analysis.rdf import RadialDistributionFunction
from ase.lattice.cubic import FaceCenteredCubic
from asap3.testtools import *

print_version(1)

# testtypes: latticeconst, maxRDF, bins, useEMT
testtypes = ((3.6, 6.001, 100, False),
             (3.7, 6.001, 100, True),
             (3.6, 15.001, 100, False),
             (3.65, 15.001, 100, True))

for latconst, maxrdf, nbins, withemt in testtypes:

    atoms = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]], symbol="Cu",
                              size=(10,10,10), latticeconstant=latconst)
    natoms = len(atoms)
    ReportTest("Number of atoms", natoms, 4000, 0)

    if withemt:
        atoms.set_calculator(EMT())
        print atoms.get_potential_energy()

    rdf = RadialDistributionFunction(atoms, maxrdf, nbins)
    z = atoms.get_atomic_numbers()[0]

    globalrdf = rdf.get_rdf()
    localrdf = rdf.get_rdf(elements=(z,z))

    print globalrdf

    ReportTest("Local and global RDF are identical",
               max(abs(globalrdf - localrdf)), 0.0, 1e-6)
        
    shellpop = [12, 6, 24, 12, 24, -1]
    shell = [sqrt(i+1.0)/sqrt(2.0) for i in range(6)]
    print shell
    print shellpop
    n = 0

    dr = maxrdf/nbins

    for i in range(nbins):
        if (i+1)*dr >= shell[n] * latconst:
            if shellpop[n] == -1:
                print "Reached the end of the test data"
                break
            rho = len(atoms) / atoms.get_volume()
            expected = shellpop[n] / (4 * pi * ((i+0.5) * dr)**2 * dr * rho)
            maxerr = shellpop[n] / (8 * pi * ((i+0.5) * dr)**3 * rho)
            print "Shell", n+1, globalrdf[i]
            ReportTest(("Shell %d (%d)" % (n+1, i)), globalrdf[i],
                       expected, maxerr)
            n += 1
        else:
            ReportTest(("Between shells (%d)" % (i,)), globalrdf[i], 0, 0)


ReportTest.Summary()

