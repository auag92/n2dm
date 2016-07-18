from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.testtools import ReportTest

print "Making atoms"
atoms = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]],
                          size=(25,25,25), symbol="Cu")
memory_usage(atoms)

print "Attaching EMT potential"
atoms.set_calculator(EMT())
memory_usage(atoms)

print "Calculating forces"
atoms.get_forces()
memory_usage(atoms)

print "Test passed!"
