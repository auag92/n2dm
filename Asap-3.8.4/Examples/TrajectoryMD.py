"""Simple molecular dynamics.

A block of 27 cubic unit cells of Cu is set up, a single atom is given
a significant momentum, and constant energy molecular dynamics is
performed.

"""

from numpy import *
from asap3 import Atoms, EMT, units, PickleTrajectory
from ase.lattice.cubic import FaceCenteredCubic
from asap3.md.verlet import VelocityVerlet

# Create the atoms
atoms = FaceCenteredCubic(size=(3,3,3), symbol="Cu", pbc=False)

# Give the first atom a non-zero momentum
atoms[0].set_momentum(array([0, -11.3, 0]))

# Associate the EMT potential with the atoms
atoms.set_calculator(EMT())

# Now do molecular dynamics, printing kinetic, potential and total
# energy every ten timesteps.
dyn = VelocityVerlet(atoms, 5.0*units.fs)


# Make a trajectory writing output
trajectory = PickleTrajectory("TrajectoryMD-output.traj", "w", atoms)

# Attach it to the dynamics, so it is informed every fifth time a
# timestep is made.
dyn.attach(trajectory, interval=5)

# Now do 1000 timesteps.
dyn.run(1000)
print "The output is in the NetCDF file TrajectoryMD-output.traj"
