# Checks for Ticket #11

from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.testtools import ReportTest

print "Test for Ticket #11: https://trac.fysik.dtu.dk/projects/Asap/ticket/11"

atoms = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]], size=(6,6,6),
                          symbol="Cu")
atoms.set_calculator(EMT())
r = atoms.get_positions()
print "Orig position", r[-1]

uc = atoms.get_cell()
print uc
r[-1] = 1.51*uc[2]
atoms.set_positions(r)
print atoms.get_potential_energy()

p1 = atoms.get_positions()[-1]
print "p1:", p1

atoms.set_cell(uc, scale_atoms=True)
print atoms.get_potential_energy()
p2  = atoms.get_positions()[-1]
print "p2:", p2

atoms.set_cell(uc, scale_atoms=False)
print atoms.get_potential_energy()
p3 = atoms.get_positions()[-1]
print "p3:", p3

ReportTest("p2 equals p1", p2[2], p1[2], 1e-6)
ReportTest("p3 equals p1", p3[2], p1[2], 1e-6)
ReportTest.Summary()
