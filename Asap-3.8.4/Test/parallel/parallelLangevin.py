#! /usr/bin/env python

"""testLangevin.py - test the Langevin dynamics.

Usage: python testLangevin.py

Tests Langevin dynamics using the EMT Copper potential.
"""

import sys, time
from numpy import *
import asap3
from asap3.testtools import ReportTest
from asap3.md.langevin import Langevin
from asap3.md.verlet import VelocityVerlet
from ase.lattice.cubic import FaceCenteredCubic
from asap3 import *
from asap3.mpi import world

nequil = 1000
nequilprint = 25
nsteps = 20000
nprint = 250
tolerance = 0.01
nminor = 25
timestep = 0.5

#set_verbose(1)

# Set up atoms in a regular simple-cubic lattice.
ismaster = world.rank == 0
isparallel = world.size != 1
if world.size == 1:
    cpulayout = None
elif world.size == 2:
    cpulayout = [1,2,1]
elif world.size == 3:
    cpulayout = [1,3,1]
elif world.size == 4:
    cpulayout = [1,2,2]

print_version(1)

if ismaster:
    atoms = FaceCenteredCubic(size=(2,10,10), symbol="Cu", pbc=False,
                              latticeconstant = 3.5)
else:
    atoms = None

if isparallel:
    atoms = MakeParallelAtoms(atoms, cpulayout)
atoms.set_calculator(asap3.EMT())
natoms = atoms.get_number_of_atoms()
    
ReportTest("Number of atoms", natoms, 800, 0)

# Make a small perturbation of the momenta
atoms.set_momenta(1e-6 * random.random([len(atoms), 3]))
print "Initializing ..."
predyn = VelocityVerlet(atoms, 0.5)
try:
    predyn.run(2500)
except:
    print atoms.arrays['positions']
    print atoms.arrays['momenta']
    print atoms.arrays['momenta'].shape
    print atoms.get_masses()
    print atoms.get_masses().shape
    
    raise

initr = atoms.get_positions()
initp = atoms.get_momenta()


def targetfunc(params, x):
    return params[0] * exp(-params[1] * x) + params[2]

output = file("Langevin.dat", "w")

for temp, frict in ((0.01, 0.001),):
    dyn = Langevin(atoms, timestep, temp, frict)
    print ""
    print "Testing Langevin dynamics with T = %f eV and lambda = %f" % (temp, frict)
    ekin = atoms.get_kinetic_energy()/natoms
    print ekin
    output.write("%.8f\n" % ekin)
    temperatures = [(0, 2.0 / 3.0 * ekin)]
    a = 0.1
    b = frict
    c = temp
    print "Equilibrating ..."
    tstart = time.time()
    for i in xrange(1,nequil+1):
        dyn.run(nminor)
        ekin = atoms.get_kinetic_energy() / natoms
        if i % nequilprint == 0:
            print "%.6f  T = %.6f (goal: %f)" % \
                  (ekin, 2./3. * ekin, temp)
        output.write("%.8f\n" % ekin)
    tequil = time.time() - tstart
    print "This took %s minutes." % (tequil / 60)
    output.write("&\n")
    temperatures = []
    print "Taking data - this takes", nsteps/nequil, "times longer!"
    tstart = time.time()
    for i in xrange(1,nsteps+1):
        dyn.run(nminor)
        ekin = atoms.get_kinetic_energy() / natoms
        temperatures.append(2.0/3.0 * ekin)
        if i % nprint == 0:
            tnow = time.time() - tstart
            tleft = (nsteps-i) * tnow / i
            print "%.6f    (time left: %.1f minutes)" % (ekin, tleft/60)
        output.write("%.8f\n" % ekin)
    output.write("&\n")
    temperatures = array(temperatures)
    mean = sum(temperatures) / len(temperatures)
    print "Mean temperature:", mean, "eV"
    print ""
    print "This test is statistical, and may in rare cases fail due to a"
    print "statistical fluctuation."
    print ""
    ReportTest("Mean temperature:", mean, temp, tolerance*temp)
            
output.close()

