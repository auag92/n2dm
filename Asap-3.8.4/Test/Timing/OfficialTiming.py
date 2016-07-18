#!/usr/bin/env python

"""OfficialTiming.py - Asap timing script.

Reports timing of various types of dynamics in usec/atom/timestep.
"""

from numpy import *
from asap3 import *
from asap3.md.verlet import VelocityVerlet
from asap3.md.langevin import Langevin
from asap3.md.npt import NPT
from asap3.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.lattice.cubic import FaceCenteredCubic, BodyCenteredCubic
from ase.lattice.compounds import L1_2
#from Asap.Timing import *
import sys, cPickle, time, commands, os, re
import numpy as np
from asap3.testtools import ReportTest

import asap3.Internal.Threads

from StringIO import StringIO
import sys
import commands
import time

nthreads = 1
if len(sys.argv) > 1 and sys.argv[1] == "-t":
    nthreads = AsapThreads()

if len(sys.argv) > 1 and sys.argv[1] == "-T":
    nthreads = AsapThreads(2, force=True)

host = commands.getoutput("hostname")
logfilename = "officialtiming.log"
when = time.strftime("%a %d %b %Y %H:%M", time.localtime(time.time()))
asapversion = get_version()

random.seed([42, 12345])

class Logger(StringIO):
    def __init__(self, filename):
        StringIO.__init__(self)
        self.logfile = filename
    def write(self, x):
        sys.stdout.write(x)
        StringIO.write(self, x)
    def close(self):
        outfile = open(self.logfile, "a")
        outfile.write(self.getvalue())
        outfile.close()
        StringIO.close(self)

logger = Logger(logfilename)

def Timing(name, func, timesteps, natoms):
    print "Preparing timing"
    func(2)
    print "Running timing:", name
    startcpu, startwall = time.clock(), time.time()
    func(timesteps)
    cpu, wall = time.clock() - startcpu, time.time() - startwall
    fraction = 100.0 * cpu/wall
    print "Time:", cpu, "CPU seconds  ", wall, "Wall seconds" 
    cpu *= 1e6 / (timesteps * natoms)
    wall *= 1e6 / (timesteps * natoms)
    logger.write("%s: %.2f usec/atom/timestep CPU-time  %.2f Wall time (%.1f%%)\n"
                 % (name, cpu, wall, fraction))

def temperature(atoms):
    return 2.0/3.0 * atoms.get_kinetic_energy() / len(atoms) / units.kB

def MakeCu(T=300, size=(29,29,30)):
    print "Preparing", T, "K Copper system."
    atoms = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]],
                              symbol='Cu', size=size)
    atoms.set_calculator(EMT())
    MaxwellBoltzmannDistribution(atoms, 2*T * units.kB)
    #dyn = VelocityVerlet(atoms, 5*units.fs)
    dyn = Langevin(atoms, 5*units.fs, T*units.kB, 0.05)
    dyn.run(50)
    print "Done.  Temperature =", temperature(atoms), \
          "K.  Number of atoms: ", len(atoms)
    return atoms

def MakeCu3Ni(T=300):
    print "Preparing", T, "K NiCu3 system."
    atoms = L1_2(directions=[[1,0,0],[0,1,0],[0,0,1]], symbol=('Ni', 'Cu'),
                 latticeconstant=3.61, size=(29,29,30))
    atoms.set_calculator(EMT())
    MaxwellBoltzmannDistribution(atoms, 2*T * units.kB)
    #dyn = VelocityVerlet(atoms, 5*units.fs)
    dyn = Langevin(atoms, 5*units.fs, T*units.kB, 0.05)
    dyn.run(50)
    print "Done.  Temperature =", temperature(atoms), \
          "K.  Number of atoms: ", len(atoms)
    return atoms

# Make sure that CPU speed is revved up.
dummy = L1_2(directions=[[1,0,0],[0,1,0],[0,0,1]], symbol=('Au', 'Cu'),
             latticeconstant=4.08, size=(10,10,10), debug=0)
dummy.set_calculator(EMT())
dyn = Langevin(dummy, 5*units.fs, 300*units.kB, 0.05)
dyn.run(10)
del dummy

logger.write("\n\nRunning on %s %s\n" % (host, when))
modelname = "Unknown CPU model name"
cpumhz = "Unknown CPU speed"
try:
    lines = open("/proc/cpuinfo").readlines()
    for line in lines:
        if line[:10] == "model name":
            modelname = line
            break
    for line in lines:
        if line[:7] == "cpu MHz":
            cpumhz = line
            break
except:
    print "Cannot get CPU info from /proc/cpuinfo"

logger.write(cpumhz)
logger.write(modelname)
logger.write(asapversion+"\n")
logger.write("Number of threads: " + str(nthreads)+"\n")

atoms300 = MakeCu()

atoms = Atoms(atoms300)
atoms.set_calculator(EMT())
dyn = VelocityVerlet(atoms, 5*units.fs)
Timing("Ver300", dyn.run, 50, len(atoms))

atoms = MakeCu(1000)
dyn = VelocityVerlet(atoms, 5*units.fs)
Timing("Ver1000", dyn.run, 50, len(atoms))

atoms = Atoms(atoms300)
atoms.set_calculator(EMT())
dyn = Langevin(atoms, 5*units.fs, 300*units.kB, 0.001)
Timing("Langevin", dyn.run, 50, len(atoms))

atoms = Atoms(atoms300)
atoms.set_calculator(EMT())
dyn = NPT(atoms, 5*units.fs, 300*units.kB, 0, 25*units.fs,
          (75*units.fs)**2 * 140*units.GPa)
Timing("NPT", dyn.run, 50, len(atoms))

atoms = Atoms(atoms300, pbc=(0,0,0))
atoms.set_calculator(EMT())
dyn = VelocityVerlet(atoms, 5*units.fs)
Timing("FreeBC", dyn.run, 50, len(atoms))

atoms = MakeCu3Ni()
atoms.set_calculator(EMT())
dyn = VelocityVerlet(atoms, 5*units.fs)
Timing("Cu3Ni", dyn.run, 50, len(atoms))

atoms = MakeCu(size=(6,6,7))
atoms.set_calculator(EMT())
dyn = VelocityVerlet(atoms, 5*units.fs)
Timing("Tiny", dyn.run, 5000, len(atoms))

elements = [29]
epsilon  = [0.15]
sigma    = [2.34]
rcut = 0.5*(sqrt(3)+2)*1.09*sigma[0]

atoms = Atoms(atoms300)
atoms.set_calculator(LennardJones(elements, epsilon, sigma, rcut, True))
dyn = VelocityVerlet(atoms, 5*units.fs)
Timing("L-J", dyn.run, 50, len(atoms))

logger.close()
