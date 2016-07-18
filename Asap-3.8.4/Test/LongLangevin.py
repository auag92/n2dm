#! /usr/bin/env python

"""longLangevin.py - test the Langevin dynamics.

Usage: python longLangevin.py

Tests Langevin dynamics using the EMT Copper potential.  It takes a while
since it is a statistical test of the temperature control.
"""

import sys, time
import numpy as np
from asap3.testtools import ReportTest
from ase import *
from asap3 import *
from asap3.md.langevin import Langevin
from asap3.md.verlet import VelocityVerlet
from ase.lattice.cubic import FaceCenteredCubic
from Scientific.Statistics import standardDeviation, variance
from Scientific.Functions.LeastSquares import leastSquaresFit

nsteps = 1000000
#nsteps = 5000
nprint = 50000
nequil = 100000
nminor = 25
nequilprint = 200
tolerance = 0.01
timestep = 5 # fs
frict = 0.001
temp = 100 # K
repeats = 5

# Set up atoms in a regular simple-cubic lattice.
initial = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]],
                            size = (1,5,5), symbol="Cu", latticeconstant=3.5,
                            pbc=False, debug=0)
initial.set_calculator(EMT())
    
ReportTest("Number of atoms", len(initial), 100, 0)

# Make a small perturbation of the momenta
initial.set_momenta(1e-6 * np.random.normal(size=[len(initial), 3]))
print "Initializing (1) ..."
predyn = VelocityVerlet(initial, 0.5)
predyn.run(2500)

initial2 = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]],
                             size = (4,4,4), symbol="Cu", debug=0)
initial2.set_calculator(EMT())
# Make a small perturbation of the momenta
initial2.set_momenta(1e-6 * np.random.normal(size=[len(initial2), 3]))
print "Initializing (2) ..."
predyn = VelocityVerlet(initial2, 0.5)
predyn.run(2500)


class TemperatureMonitor:
    def __init__(self, report):
        self.report = report
        self.reports = []

    def New(self, name, atoms, temperature, pas):
        self.name = "%s (pass %d)" % (name, pas+1)
        self.atoms = atoms
        self.n = 0
        self.ncall = 0
        self.sum = 0.0
        self.temperature = temperature

    def Update(self):
        self.ncall += 1
        T = self.atoms.get_kinetic_energy() / (1.5 * kB * len(self.atoms))
        self.sum += T
        self.n += 1
        if self.ncall % self.report == 0:
            tavg = self.sum / self.n
            print "%5d  T = %7.2f K  T_avg = %7.2f" % (self.ncall, T, tavg)
                
    def Check(self):
        T = self.sum / self.n
        result = "%s   T = %.2f K   (expected %.2f K)" % (self.name, T,
                                                          self.temperature)
        print result
        self.reports.append(result)
        ReportTest(self.name, T, self.temperature, 0.01*self.temperature)
        
    def Report(self):
        for r in self.reports:
            print r

def targetfunc(params, x):
    return params[0] * exp(-params[1] * x) + params[2]

def equilibrate(atoms, dyn, temp):
    print "Equilibrating ..."
    tstart = time.time()
    temperatures = []
    for i in xrange(1,nequil/nminor+1):
        dyn.run(nminor)
        ekin = atoms.get_kinetic_energy() / len(atoms)
        temperatures.append((i*nminor*timestep, 2.0/3.0 * ekin))
        if i % nequilprint == 0:
            #data = array(temperatures)
            #print data.shape
            try:
                (a, b, c) = leastSquaresFit(targetfunc, (0.1, 2*frict, temp),
                                            temperatures)[0]
            except OverflowError:
                print "leastSquaresFit failed (this is OK)."
                a = b = c = 0.0
            print "%.6f  T_inf = %.6f (goal: %f)  k = %.6f" % \
                  (ekin, c, temp, b)
    tequil = time.time() - tstart
    print "This took %.1f minutes." % (tequil / 60)
    print "Taking data - this takes", repeats*nsteps/nequil, "times longer!"

monitor = TemperatureMonitor(nprint)

atoms = Atoms(initial)
atoms.set_calculator(EMT())
dyn = Langevin(atoms, timestep*femtosecond, temp*kB, 0.001)
equilibrate(atoms, dyn, kB*temp)
dyn.attach(monitor.Update, 5)
for i in range(repeats):
    monitor.New("Free boundaries, fixcm=True", atoms, temp, i)
    dyn.run(nsteps)
    monitor.Check()

atoms = Atoms(initial)
atoms.set_calculator(EMT())
dyn = Langevin(atoms, timestep*femtosecond, temp*kB, 0.001, fixcm=False)
equilibrate(atoms, dyn, kB*temp)
dyn.attach(monitor.Update, 5)
for i in range(repeats):
    monitor.New("Free boundaries, fixcm=False", atoms, temp, i)
    dyn.run(nsteps)
    monitor.Check()

atoms = Atoms(initial2)
atoms.set_calculator(EMT())
dyn = Langevin(atoms, timestep*femtosecond, temp*kB, 0.001)
equilibrate(atoms, dyn, kB*temp)
dyn.attach(monitor.Update, 5)
for i in range(repeats):
    monitor.New("Periodic boundaries, fixcm=True", atoms, temp, i)
    dyn.run(nsteps)
    monitor.Check()

atoms = Atoms(initial2)
atoms.set_calculator(EMT())
dyn = Langevin(atoms, timestep*femtosecond, temp*kB, 0.001, fixcm=False)
equilibrate(atoms, dyn, kB*temp)
dyn.attach(monitor.Update, 5)
for i in range(repeats):
    monitor.New("Periodic boundaries, fixcm=False", atoms, temp, i)
    dyn.run(nsteps)
    monitor.Check()

monitor.Report()
ReportTest.Summary()


