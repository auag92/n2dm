import numpy as np

from ase import units
from ase.data import atomic_numbers
from ase.io import PickleTrajectory
from ase.md import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.lattice.cubic import SimpleCubic
from asap3 import Morse, set_verbose
from asap3.testtools import ReportTest
import time

np.random.seed(42)

#set_verbose(2)

# Define constants and calculator
elements = np.array([atomic_numbers['Ru'], atomic_numbers['Ar']])
epsilon = np.array([[5.720, 0.092], [0.092, 0.008]])
alpha = np.array([[1.475, 2.719], [2.719, 1.472]])
rmin = np.array([[2.110, 2.563], [2.563, 4.185]])

calc = Morse(elements, epsilon, alpha, rmin)

# Define atoms
atoms = SimpleCubic('Ar', size=(15,15,15), latticeconstant=50.0)
n = 0
while n < 100:
    i = np.random.randint(len(atoms)-1)
    if atoms[i].number != atomic_numbers['Ru']:
        atoms[i].number = atomic_numbers['Ru']
        n += 1
atoms.set_calculator(calc)

# Set initial momentum
MaxwellBoltzmannDistribution(atoms, 1000*units.kB)

# Run dynamics
startcpu, startwall = time.clock(), time.time()

dyn = VelocityVerlet(atoms, 1.0 * units.fs, logfile='test-energy.dat', loginterval=10)
dyn.run(10)
etot = (atoms.get_potential_energy() + atoms.get_kinetic_energy())/len(atoms)
#obs.close()
for i in range(25):
    if i:
        dyn.run(1000)
    epot = atoms.get_potential_energy()/len(atoms)
    ekin = atoms.get_kinetic_energy()/len(atoms)
    print "%9.5f %9.5f   %9.5f" % (epot, ekin, epot+ekin)
    ReportTest("Step %i." % (i,), epot+ekin, etot, 1e-3, silent=True)

cpu, wall = time.clock() - startcpu, time.time() - startwall
fraction = cpu/wall

print ""
print ""
print "TIMING RESULTS:"
print "Morse potential:   CPU time %.2fs  Wall clock time %.2fs (%.0f%%)" % (cpu, wall, fraction*100)
print ""

ReportTest.Summary()
