"""Test the FixAtoms constraint and the Subset filter."""

from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.md.velocitydistribution import MaxwellBoltzmannDistribution
from asap3.md.verlet import VelocityVerlet
from asap3.io.trajectory import PickleTrajectory
from ase.constraints import Filter
from asap3.constraints import FixAtoms
from ase.constraints import FixAtoms as ASE_FixAtoms
from asap3.testtools import ReportTest
import numpy as np

init = FaceCenteredCubic(size=(10,10,10), symbol='Cu', pbc=False)
z = init.get_positions()[:,2]
fixedatoms = np.less(z, 0.501*z.max())
print len(init), sum(fixedatoms)
MaxwellBoltzmannDistribution(init, 2000*units.kB)

print
print "Running simulation with Filter"
atoms1 = Atoms(init)
atoms1.set_calculator(EMT())
atoms1a = Filter(atoms1, mask=np.logical_not(fixedatoms))

dyn = VelocityVerlet(atoms1a, 0.5*units.fs)
dyn.run(50)
r1 = atoms1.get_positions()

print
print "Running simulation with Asap's FixAtoms"
atoms2 = Atoms(init)
atoms2.set_calculator(EMT())
atoms2.set_constraint(FixAtoms(mask=fixedatoms))

dyn = VelocityVerlet(atoms2, 0.5*units.fs)
dyn.run(50)
r2 = atoms2.get_positions()

print
print "Running simulation with ASE's FixAtoms"
atoms3 = Atoms(init)
atoms3.set_calculator(EMT())
# Just to be difficult, convert mask to indices
indx = np.compress(fixedatoms, np.arange(len(atoms3)))
assert len(indx) == fixedatoms.sum()
atoms3.set_constraint(ASE_FixAtoms(indx))

dyn = VelocityVerlet(atoms3, 0.5*units.fs)
dyn.run(50)
r3 = atoms2.get_positions()

err = np.max(np.abs(r1 - r2).flat)
print
print "Filter and Asap's FixAtoms:", err
ReportTest("Identical positions (Filter and Asap's FixAtoms)", err, 0.0, 1e-9)

err = np.max(np.abs(r2 - r3).flat)
print "ASE's and Asap's FixAtoms:", err
ReportTest("Identical positions (ASE's and Asap's FixAtoms)", err, 0.0, 1e-9)

print
print "Running Langevin simulation with Asap's FixAtoms"
atoms4 = Atoms(init)
atoms4.set_calculator(EMT())
atoms4.set_constraint(FixAtoms(mask=fixedatoms))

dyn = Langevin(atoms4, 0.5*units.fs, 1000*units.kB, 0.01)
dyn.run(50)

print
print "Running Langevin simulation with ASE's FixAtoms"
atoms5 = Atoms(init)
atoms5.set_calculator(EMT())
atoms5.set_constraint(ASE_FixAtoms(mask=fixedatoms))

dyn = Langevin(atoms5, 0.5*units.fs, 1000*units.kB, 0.01)
dyn.run(50)

print
sanity = [[atoms1, "Verlet + Filter"],
          [atoms2, "Verlet + Asap's FixAtoms"],
          [atoms3, "Verlet + ASE's FixAtoms"],
          [atoms4, "Langevin + Asap's FixAtoms"],
          [atoms5, "Langevin + ASE's FixAtoms"],
          ]
r_init = init.get_positions()
for a, label in sanity:
    print "Sanity check, %s:", (label,)
    ok = (a.get_positions() == r_init) + np.logical_not(fixedatoms)[:,np.newaxis]
    if ok.all():
        print "  Stationary atoms have not moved: OK"
    else:
        raise RuntimeError("Stationary atoms have moved (%s)" % (label,))
    ok = (a.get_positions() != r_init) + fixedatoms[:,np.newaxis]
    if ok.all():
        print "  Mobile atoms have moved: OK"
    else:
        raise RuntimeError("Mobile atoms have not moved (%s)" % (label,))
    

ReportTest.Summary()




