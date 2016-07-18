from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.testtools import ReportTest
from asap3.md.velocitydistribution import *

#set_verbose(1)

atoms = FaceCenteredCubic(directions=((1,0,0), (0,1,0), (0,0,1)),
                          size=(15,15,15), symbol="Cu", pbc=True)
atoms.set_calculator(EMT())

atoms.get_forces()
atoms.get_forces()
MaxwellBoltzmannDistribution(atoms, 300 * units.kB)
atoms.get_forces()
atoms.get_forces()

atoms = FaceCenteredCubic(directions=((1,0,0), (0,1,0), (0,0,1)),
                          size=(15,15,15), symbol="Cu", pbc=True)
atoms.set_calculator(EMT())
s = atoms.get_stress()
print
print "Stress:", s
s = atoms.get_stress()
print
print "Stress:", s
MaxwellBoltzmannDistribution(atoms, 300 * units.kB)
s = atoms.get_stress()
print
print "Stress:", s
s = atoms.get_stress()
print
print "Stress:", s
MaxwellBoltzmannDistribution(atoms, 300 * units.kB)
s = atoms.get_stress()
print
print "Stress:", s
s = atoms.get_stress()
print
print "Stress:", s

atoms = FaceCenteredCubic(directions=((1,0,0), (0,1,0), (0,0,1)),
                          size=(15,15,15), symbol="Cu", pbc=True)
atoms.set_calculator(EMT())
s = atoms.get_stress()
atoms.get_forces()
atoms.get_forces()
MaxwellBoltzmannDistribution(atoms, 300 * units.kB)
atoms.get_forces()
atoms.get_forces()

atoms = FaceCenteredCubic(directions=((1,0,0), (0,1,0), (0,0,1)),
                          size=(15,15,15), symbol="Cu", pbc=True)
atoms.set_calculator(EMT())
atoms.get_stresses()
atoms.get_stresses()
MaxwellBoltzmannDistribution(atoms, 300 * units.kB)
atoms.get_stresses()
atoms.get_stresses()

print
print
print "No crash: Test passes succesfully!"
