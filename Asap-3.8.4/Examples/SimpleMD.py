"""Simple molecular dynamics.

A block of 27 cubic unit cells of Cu is set up, a single atom is given
a significant momentum, and constant energy molecular dynamics is
performed.

"""

from numpy import *
from asap3 import Atoms, EMT, units
from ase.lattice.cubic import FaceCenteredCubic
from asap3.md.verlet import VelocityVerlet

# Create the atoms
atoms = FaceCenteredCubic(size=(3,3,3), symbol="Cu", pbc=False)

# Give the first atom a non-zero momentum
atoms[0].set_momentum(array([0, -11.3, 0]))
print "Kinetic energies of all atoms:"
p = atoms.get_momenta()
kinenergies = 0.5 * (p * p).sum(1) / atoms.get_masses()
print kinenergies

# Associate the EMT potential with the atoms
atoms.set_calculator(EMT())

# Now do molecular dynamics, printing kinetic, potential and total
# energy every ten timesteps.
dyn = VelocityVerlet(atoms, 5.0*units.fs)
print ""
print "Energy per atom:"
print "  %15s %15s %15s" % ("Pot. energy", "Kin. energy", "Total energy")

for i in range(25):
    dyn.run(10)
    epot = atoms.get_potential_energy()/len(atoms)
    ekin = atoms.get_kinetic_energy()/len(atoms)
    print "%15.5f %15.5f %15.5f" % (epot, ekin, epot+ekin)

