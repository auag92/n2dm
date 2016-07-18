# Checks for Ticket #11

"Test for correct energy calculations with varying unit cell."

from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.testtools import ReportTest

print __doc__

for pbc in (True, False, (1,1,0)):
    for scale in (True, False):
        print "Running test with pbc=%s and scale_atoms=%s" % (pbc, scale)
        
        atoms = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]],
                                  size=(6,6,6), symbol="Cu", pbc=pbc)
        atoms.set_calculator(EMT())
        uc = atoms.get_cell()
        atoms.get_potential_energy()
        for factor in (1.0, 1.01, 1.02, 1.1, 1.5, 1.49, 1.4, 1.4, 1.0, 0.9):
            atoms.set_cell(uc * factor, scale_atoms=scale)
            f = atoms.get_forces()
            e = atoms.get_potential_energy()
            atoms2 = Atoms(atoms)
            atoms2.set_calculator(EMT())
            e2 = atoms2.get_potential_energy()
            f2 = atoms2.get_forces()
            name = "(factor = %.3f  PBC = %s  scale_atoms = %s)" % (factor,
                                                                    pbc, scale)
            ReportTest("Energy "+name, e, e2, 1e-6)
            maxf = max(abs(f.flat[:] - f2.flat[:]))
            ReportTest("Forces "+name, maxf, 0.0, 1e-6)

ReportTest.Summary()
