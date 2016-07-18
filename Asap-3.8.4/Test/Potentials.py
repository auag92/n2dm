# Testing various potentials.

from asap3 import *
from asap3.md.verlet import VelocityVerlet
from asap3.EMT2013Parameters import PtY_parameters
from ase.lattice.cubic import *
from ase.lattice.compounds import *
from numpy import *
from asap3.testtools import ReportTest
try:
    import potResults
except ImportError:
    resultfail = True
else:
    resultfail = False
import os

timeunit = 1.018047e-14             # Seconds
femtosecond = 1e-15 / timeunit      # Femtosecond in atomic units

print_version(1)


def dotest(atoms, nsteps, ampl, name):
    print "Potential energy", atoms.get_potential_energy() / len(atoms)
    r = atoms.get_positions()
    r.flat[:] += ampl * sin(arange(3*len(atoms)))
    atoms.set_positions(r)
    print "Potential energy", atoms.get_potential_energy() / len(atoms)

    print "Running Verlet dynamics (%s)" % (name,)
    dyn = VelocityVerlet(atoms, 2*femtosecond)
    etot1 = (atoms.get_potential_energy() + atoms.get_kinetic_energy())
    dyn.run(nsteps)
    etot2 = (atoms.get_potential_energy() + atoms.get_kinetic_energy())
    ReportTest(("Energy conservation (%s)" % (name,)), etot1, etot2, 1.0)
    print etot1, etot2

    epot = atoms.get_potential_energies()
    stress = atoms.get_stresses()
    if firsttime:
        print "Reporting energies and stresses"
        e = []
        s = []
        j = 0
        for i in range(0, len(atoms), 100):
            e.append(epot[i])
            s.append(stress[i,j])
            j = (j + 1) % 6
        print >> out, "e"+name+" =", repr(e)
        print >> out, "s"+name+" =", repr(s)
    else:
        print "Testing energies and stresses"
        j = 0
        eres=getattr(potResults, "e"+name)
        sres=getattr(potResults, "s"+name)
        for i in range(len(atoms)/100):
            ReportTest(("%s energy %d" % (name, i*100)),
                       epot[i*100], eres[i], 1e-8, silent=True)
            ReportTest(("%s stress %d" % (name, i*100)),
                       stress[i*100, j], sres[i], 1e-8, silent=True)
            j = (j + 1) % 6

firsttime = (len(sys.argv) >= 2 and sys.argv[1] == '--first')
if firsttime and os.path.exists("potResults.py"):
    print "This will overwrite the result file 'potResults.py'."
    print "If you really want to do this, erase it and run this again."
    sys.exit(1)
if resultfail and not firsttime: 
    print "Importing 'potResults.py' failed!"
    print "Maybe you need to create it with the --first option."
    sys.exit(1)
    
if firsttime:
    print "Creating the file 'potResults.py'"
    out = open('potResults.py', "w")    

atoms = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]], size=(15,15,15),
                          symbol="Cu", pbc=(1,0,1), debug=0)
ReportTest("Number of Cu atoms", len(atoms), 13500, 0)

atoms.set_calculator(EMT())
dotest(atoms, 50, 0.1, "EMT_Cu")

atoms = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]], size=(15,15,15),
                          symbol="Cu", pbc=(1,0,1), debug=0)
ReportTest("Number of Cu atoms", len(atoms), 13500, 0)

atoms.set_calculator(EMT(EMTRasmussenParameters()))
dotest(atoms, 50, 0.1, "EMT_Cu_Rasm")

#atoms = BodyCenteredCubic([[1,0,0],[0,1,0],[0,0,1]], size=(15,15,30),
#                          element="Mo", periodic=(1,0,1), debug=0)
#ReportTest("Number of Mo atoms", len(atoms), 13500, 0)
#
#atoms.SetCalculator(MoPotential())
#dotest(atoms, 50, 0.06, "Mo")

atoms = L1_2(directions=[[1,0,0],[0,1,0],[0,0,1]], size=(15,15,15),
             symbol=("Cu", "Au"), latticeconstant=3.95, pbc=(1,0,1), 
             debug=0)
ReportTest("Number of alloy atoms", len(atoms), 13500, 0)
nCu = sum(equal(atoms.get_atomic_numbers(), 29))
nAu = sum(equal(atoms.get_atomic_numbers(), 79))
ReportTest("Number of Cu atoms in alloy", nCu, 13500/4, 0)
ReportTest("Number of Au atoms in alloy", nAu, 3*13500/4, 0)

atoms.set_calculator(EMT())
dotest(atoms, 50, 0.06, "EMT_CuAu3")

atoms = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]], size=(15,15,15),
                          symbol="Pt", pbc=(1,0,1), debug=0)
ReportTest("Number of Pt atoms", len(atoms), 13500, 0)

atoms.set_calculator(EMT2013(PtY_parameters))
dotest(atoms, 50, 0.1, "EMT2013_Pt")

atoms = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]], size=(15,15,15),
                          symbol="Y", latticeconstant=4.97,
                          pbc=(1,0,1), debug=0)
ReportTest("Number of Y atoms", len(atoms), 13500, 0)

atoms.set_calculator(EMT2013(PtY_parameters))
dotest(atoms, 50, 0.1, "EMT2013_Y")

atoms = L1_2(directions=[[1,0,0],[0,1,0],[0,0,1]], size=(15,15,15),
             symbol=("Y", "Pt"), latticeconstant=4.06, pbc=(1,0,1), 
             debug=0)
ReportTest("Number of alloy atoms", len(atoms), 13500, 0)
nY = sum(equal(atoms.get_atomic_numbers(), 39))
nPt = sum(equal(atoms.get_atomic_numbers(), 78))
ReportTest("Number of Cu atoms in alloy", nY, 13500/4, 0)
ReportTest("Number of Au atoms in alloy", nPt, 3*13500/4, 0)

atoms.set_calculator(EMT2013(PtY_parameters))
dotest(atoms, 50, 0.06, "EMT2013_Pt3Y")

if firsttime:
    # Create the "main" routine in the results module.
    print >> out, 'if __name__ == "__main__":'
    print >> out, '    print "This is not a test, but a module containing test results"'
    out.close()
    
ReportTest.Summary()

