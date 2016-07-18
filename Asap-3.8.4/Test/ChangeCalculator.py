from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.testtools import ReportTest

atoms = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]],
                          size=(6,6,6), symbol="Cu")
atoms.set_calculator(EMT())
f1 = atoms.get_forces()
atoms.set_calculator(EMT())
f2 = atoms.get_forces()
maxdev = abs(f2 - f1).max()
print maxdev
ReportTest("Max error 1:", maxdev, 0.0, 1e-6)

atoms2 = Atoms(atoms)
if atoms2.get_calculator() is None:
    # Slightly old ase
    atoms2.set_calculator(atoms.get_calculator())
f2 = atoms2.get_forces()
maxdev = abs(f2 - f1).max()
print maxdev
ReportTest("Max error 2:", maxdev, 0.0, 1e-6)

f2 = atoms.get_forces()
maxdev = abs(f2 - f1).max()
print maxdev
ReportTest("Max error 1:", maxdev, 0.0, 1e-6)

