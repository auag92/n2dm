from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.testtools import ReportTest

atoms = FaceCenteredCubic(size=(5,5,5), symbol="Cu")
pot = EMT()
atoms.set_calculator(pot)

for i in (1, 2):
    print "*** Pass", i
    ReportTest("Energy required", pot.calculation_required(atoms, ["energy"]), 1, 0)
    ReportTest("Forces required", pot.calculation_required(atoms, ["forces"]), 1, 0)
    ReportTest("Stress required", pot.calculation_required(atoms, ["stress"]), 1, 0)
    ReportTest("Magmom required", pot.calculation_required(atoms, ["magmoms"]), 1, 0)
    e = atoms.get_potential_energy()
    ReportTest("Energy not required", pot.calculation_required(atoms, ["energy"]), 0, 0)
    ReportTest("Forces required (II)", pot.calculation_required(atoms, ["forces"]), 1, 0)
    f = atoms.get_forces()
    ReportTest("Energy not required (II)", pot.calculation_required(atoms, ["energy"]), 0, 0)
    ReportTest("Forces not required", pot.calculation_required(atoms, ["forces"]), 0, 0)
    ReportTest("Energy or forces not required",
               pot.calculation_required(atoms, ["energy", "forces"]), 0, 0)
    ReportTest("Energy or stress required",
               pot.calculation_required(atoms, ["energy", "stress"]), 1, 0)
    s = atoms.get_stress()
    ReportTest("Stress not required", pot.calculation_required(atoms, ["stress"]), 0, 0)

    r = atoms.get_positions()
    r[0,0] += 0.1
    atoms.set_positions(r)

ReportTest.Summary()

