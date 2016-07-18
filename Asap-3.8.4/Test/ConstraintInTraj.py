from asap3 import *
from asap3.constraints import FixAtoms
from ase.lattice.cubic import FaceCenteredCubic
from ase.io import read, write
from asap3.io.trajectory import PickleTrajectory

fn = "testconstraint.traj"
atoms = FaceCenteredCubic(symbol='Cu', size=(3,3,3))
c = FixAtoms(indices=(2,3))
atoms.set_constraint(c)

write(fn, atoms)

atoms2 = read(fn)
c2 = atoms2.constraints[0]
print "Original class:", c.__class__
print "New class:", c2.__class__

assert c.__class__ is c2.__class__
print "Test passed."
