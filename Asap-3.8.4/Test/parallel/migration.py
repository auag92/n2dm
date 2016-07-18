from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.testtools import ReportTest
from asap3.EMT2013Parameters import PtY_parameters
from asap3.mpi import world
from asap3.Internal.ParallelListOfAtoms import ParallelAtoms
import numpy as np

#DebugOutput("migration%d.log", nomaster=True)

def pot():
    #return EMT2013(PtY_parameters)
    return EMT()

#set_verbose(1)
master = world.rank == 0

if master:
    atoms0 = FaceCenteredCubic(symbol='Pt', size=(15,15,30))
else:
    atoms0 = None
    
atoms0 = MakeParallelAtoms(atoms0, (1,1,2))
atoms0.set_calculator(pot())

print >>sys.stderr, "*********** FIRST FORCE CALCULATION ************"
print >>sys.stderr, "len(atoms) =", len(atoms0), "   no. atoms =", atoms0.get_number_of_atoms()
f0 = atoms0.get_forces()
perturbation = 0.01 * np.random.standard_normal(atoms0.get_positions().shape)
r = atoms0.get_positions() + perturbation
atoms0.set_positions(r)
print >>sys.stderr, "*********** SECOND FORCE CALCULATION ************"

f1 = atoms0.get_forces()

print >>sys.stderr, "*********** COPYING ATOMS **************"
atoms2 = ParallelAtoms((1,1,2), atoms0.comm, atoms0, distribute=False)
atoms2.set_calculator(pot())
print >>sys.stderr, "*********** THIRD FORCE CALCULATION ************"
f2 = atoms2.get_forces()

#maxdev = abs(f2 - f1).max()
#print maxdev
#ReportTest("Max error 1:", maxdev, 0.0, 1e-6)

#ReportTest.Summary()

print >>sys.stderr, "No crashes - success !!"