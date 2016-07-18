from Asap import *
from Asap.Dynamics.Langevin import Langevin
from Numeric import *
from LinearAlgebra import *
from Asap.testtools import ReportTest
from Asap.Setup.Lattice.Cubic import FaceCenteredCubic

class Check:
    def __init__(self, atoms):
        self.atoms = atoms
        self.invbasis = inverse(atoms.GetUnitCell())
        self.n = 0
    def Update(self):
        self.n += 1
        r = self.atoms.GetWrappedPositions()
        #r = atoms.GetCartesianPositions()
        rs = matrixmultiply(r, self.invbasis).flat
        ReportTest.BoolTest("Underflow (%d)" % (self.n,), min(rs) >= 0.0,
                            silent=True)
        ReportTest.BoolTest("Overflow (%d)" % (self.n,), max(rs) <= 1.0,
                            silent=True)
        print "Ekin =", atoms.GetKineticEnergy()/len(atoms), min(rs), max(rs)
        
PrintVersion(1)

atoms = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]],
                          size = (20,10,4), element='Cu', periodic=(1,1,1))
atoms.SetCalculator(EMT())
dyn = Langevin(atoms, 5*femtosecond, 3000*kB, 0.05)
dyn.Attach(Check(atoms))
dyn.Run(100)

r = matrixmultiply(atoms.GetCartesianPositions(),
                   inverse(atoms.GetUnitCell())).flat
print min(r), max(r)
atoms2 = ListOfAtoms(atoms)
atoms2.SetCartesianPositions(2*atoms2.GetCartesianPositions() - 4.0)
atoms2.SetCalculator(EMT())
r = matrixmultiply(atoms2.GetCartesianPositions(),
                   inverse(atoms2.GetUnitCell())).flat
print min(r), max(r)
Check(atoms2).Update()

ReportTest.Summary()

