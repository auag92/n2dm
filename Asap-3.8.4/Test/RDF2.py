"Test the RawRadialDistributionFunction function."

from numpy import *
from asap3 import *
from asap3.analysis.rdf import RadialDistributionFunction
from ase.lattice.compounds import NaCl
from asap3.testtools import *
import sys

print_version(1)

if len(sys.argv) == 2 and sys.argv[1] == "-exit":
    ReportTest.ExitAfterError()
    
# testtypes: latticeconst, maxRDF, bins, useEMT, save
testtypes = ((5.1, 6.001, 100, False, True),
             (5.2, 6.001, 100, True, False),
             (5.1, 15.001, 100, False, False),
             (5.15, 15.001, 100, True, False)
             )

for latconst, maxrdf, nbins, withemt, save in testtypes:

    atoms = NaCl(directions=[[1,0,0],[0,1,0],[0,0,1]],
                 symbol=("Cu", "Au"), size=(10,10,10),
                 latticeconstant=latconst, debug=0)
    natoms = len(atoms)
    ReportTest("Number of atoms", natoms, 8000, 0)

    if withemt:
        atoms.set_calculator(EMT())
        print atoms.get_potential_energy()

    rdf = RadialDistributionFunction(atoms, maxrdf, nbins)
    if save:
        rdf.output_file("RDF2-rdf")
    z = atoms.get_atomic_numbers()[0]

    globalrdf = rdf.get_rdf()
    localrdf = zeros(globalrdf.shape, float)
    for i in (29,79):
        for j in (29,79):
            localrdf += rdf.get_rdf(elements=(i,j))

    ReportTest("Local and global RDF are identical",
               min( globalrdf == localrdf), 1, 0)
        
    shellpop = [6, 12, 8, 6, 24, -1]
    shell = [sqrt(i+1.0)/2.0 for i in range(6)]
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
            print "expected, error", expected, maxerr, globalrdf[i]-expected
            ReportTest(("Shell %d (%d)" % (n+1, i)), globalrdf[i],
                       expected, maxerr)
            n += 1
        else:
            ReportTest(("Between shells (%d)" % (i,)), globalrdf[i], 0, 0)

    if save:
        rdf = RadialDistributionFunction.load("RDF2-rdf0000.rdf")

        newglobalrdf = rdf.get_rdf()
        ReportTest("Saved and original RDF have same length",
                   len(newglobalrdf), len(globalrdf), 0)
        for i in range(len(globalrdf)):
            ReportTest("Saved and global RDF at position %d" % (i,),
                       newglobalrdf[i], globalrdf[i], 0)
            
        localrdf = zeros(globalrdf.shape, float)
        for i in (29,79):
            for j in (29,79):
                x = rdf.get_rdf(elements=(i,j))
                localrdf += x

        ReportTest("Saved local and global RDF are identical",
                   min( globalrdf == localrdf), 1, 0)

        os.unlink("RDF2-rdf0000.rdf")
    
ReportTest.Summary()

