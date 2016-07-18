from Asap import *
from Asap.Setup.Lattice.Cubic import *
from Asap.testtools import ReportTest
import ASE.ChemicalElements.mass

print "Testing overflow of neighbor list in Lennard Jones potential (ticket 19)"
atoms = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]], size=(20,20,2),
                          element="Cu")

uc = atoms.GetUnitCell()
uc[2,2] *= 20
atoms.SetUnitCell(uc, fix=True)

elements = [29]
epsilon  = [0.15]
sigma    = [2.7]
masses   = [ASE.ChemicalElements.mass.masses[29]]

for cut in (2.0, 5.0, 10.0):
    print "Cutoff =", cut,
    atoms.SetCalculator(LJPotential(1, elements, epsilon, sigma, masses, cut))
    print "  Energy =", atoms.GetPotentialEnergy()

print "Test passed!"
