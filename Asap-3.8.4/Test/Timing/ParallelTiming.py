#!/usr/bin/env python

"""OfficialTiming.py - Asap timing script.

Reports timing of various types of dynamics in usec/atom/timestep.
"""

from numpy import *
from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.testtools import ReportTest
from asap3.md.velocitydistribution import MaxwellBoltzmannDistribution
from asap3 import mpi
from asap3.md.npt import NPT
from ase import units
from StringIO import StringIO
import sys
import commands
import time

possibleLayouts = [(1,1,2), (1,2,2), (2,2,2), (2,2,3), (2,2,4), (2,3,3),
                   (3,3,3), (3,3,4), (3,4,4), (4,4,4), (4,4,5), (4,5,5),
                   (4,6,6), (6,6,6)]
systemSizes = [((6,6,7), 1500), ((13,14,14), 150), ((29,29,30), 15), ((63,63,63), 8)]
##systemSizes = [((6,6,7), 100), ((13,14,14), 10)]

if len(sys.argv) > 1:
    postfix = "-" + sys.argv[1]
else:
    postfix = ""
logfilename = "paralleltiming%s.log" % (postfix,)
bigtablename = "paralleltiming%s-big.table" % (postfix,)
fasttablename = "paralleltiming%s-fast.table" % (postfix,)

ismaster = mpi.world.rank == 0

class Logger(StringIO):
    def __init__(self, filename, silent=False):
        StringIO.__init__(self)
        self.logfile = filename
        self.silent = silent
    def write(self, x):
        if not self.silent:
            sys.stdout.write(x)
        StringIO.write(self, x)
    def close(self):
        if ismaster:
            outfile = open(self.logfile, "a")
            outfile.write(self.getvalue())
            outfile.close()
        StringIO.close(self)

logger = Logger(logfilename, silent=(not ismaster))
bigtable = Logger(bigtablename, silent=True)
fasttable = Logger(fasttablename, silent=True)

def SerialTiming(name, func, timesteps, natoms):
    if ismaster:
        print "Preparing timing"
    func(2)
    if ismaster:
        print "Running serial timing:", name
    mpi.world.barrier()
    startcpu, startwall = time.clock(), time.time()
    func(timesteps)
    cpu, wall = time.clock() - startcpu, time.time() - startwall
    fraction = 100.0 * cpu/wall
    if ismaster:
        print "Time:", cpu, "CPU seconds  ", wall, "Wall seconds" 
    wall = mpi.world.sum(wall)
    wall *= 1e6 / (timesteps * natoms * mpi.world.size)
    fraction = mpi.world.min(fraction)
    if ismaster:
        logger.write("Serial %s: %.2f usec/atom/timestep Wall-time  (%.1f%%)\n"
                     % (name, wall, fraction))
    return wall

def ParallelTiming(name, func, timesteps, natoms, serialtime, summary):
    if ismaster:
        print "Preparing timing"
    func(2)
    if ismaster:
        print "Running parallel timing:", name
    mpi.world.barrier()
    startcpu, startwall = time.clock(), time.time()
    func(timesteps)
    mpi.world.barrier()
    cpu, wall = time.clock() - startcpu, time.time() - startwall
    if ismaster:
        print "Time:", cpu, "CPU seconds  ", wall, "Wall seconds"
    wall *= 1e6 / (timesteps * natoms) * mpi.world.size
    efficiency = 100.0 * serialtime / wall
    if ismaster:
        logger.write("Parallel %s: %.2f usec/atom/timestep Wall-time. Efficiency: %.1f%%\n"
                     % (name, wall, efficiency))
        summary.write("%5.2f %3d%% " % (wall, efficiency))

def temperature(atoms):
    return 2.0/3.0 * atoms.get_kinetic_energy() / len(atoms) / units.kB

def MakeCu(size, T=300):
    if ismaster:
        print "Preparing", T, "K Copper system."
    atoms = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]],
                              symbol='Cu', size=size, debug=0)
    atoms.set_calculator(EMT())
    MaxwellBoltzmannDistribution(atoms, 2*T * units.kB)
    #dyn = VelocityVerlet(atoms, 5*units.fs)
    dyn = Langevin(atoms, 5*units.fs, T*units.kB, 0.05)
    dyn.run(50)
    if ismaster:
        print "Done.  Temperature =", temperature(atoms), \
              "K.  Number of atoms: ", len(atoms)
    return atoms


host = commands.getoutput("hostname")
when = time.strftime("%a %d %b %Y %H:%M", time.localtime(time.time()))
asapversion = get_version()

random.seed([42+mpi.world.rank, 12345*mpi.world.rank])

layouts = {}
for l in possibleLayouts:
    layouts[ l[0]*l[1]*l[2] ] = l

if ismaster:
    print "Running on %d processors with %s as master" % (mpi.world.size, host)

cpuLayout = layouts[mpi.world.size]


# Make sure that CPU speed is revved up.
dummy = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]], symbol='Cu',
                          size=(10,10,10), debug=0)
dummy.set_calculator(EMT())
dyn = Langevin(dummy, 5*units.fs, 300*units.kB, 0.05)
dyn.run(10)
del dummy

logger.write("\n\nRunning on %d processors with %s as master (%s)\n"
             % (mpi.world.size, host, when))
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

if ismaster:
    bigtable.write("%-5d " % mpi.world.size)
    fasttable.write("%-5d " % mpi.world.size)

for size, steps in systemSizes:
    atoms = MakeCu(size)
    atomspercpu = len(atoms)
    positions = atoms.get_positions()
    basis = atoms.get_cell()
    atomicnumber = int(atoms.get_atomic_numbers()[0])

    # Run the serial timings
    MaxwellBoltzmannDistribution(atoms, 300 * units.kB)
    name = "Verlet%d" % (atomspercpu,)
    atoms.set_calculator(EMT())
    dyn = VelocityVerlet(atoms, 5*units.fs)
    verlet = SerialTiming(name, dyn.run, steps, atomspercpu)
    name = "Langevin%d" % (atomspercpu,)
    dyn = Langevin(atoms, 5*units.fs, 300*units.kB, 0.001)
    langevin = SerialTiming(name, dyn.run, steps, atomspercpu)
    #name = "NPT%d" % (atomspercpu,)
    #dyn = NPT(atoms, 5*units.fs, 300*units.kB, 0, 25*units.fs,
    #          (75*units.fs)**2 * 140*units.GPa)
    #npt = SerialTiming(name, dyn.run, steps, atomspercpu)

    # A distributed atoms.Repeat
    n = mpi.world.rank
    gridpos = (n % cpuLayout[0], (n / cpuLayout[0]) % cpuLayout[1],
               n / (cpuLayout[0] * cpuLayout[1]))
    newpos = positions + dot(array(gridpos), basis)[newaxis,:]
    atoms = Atoms(positions=newpos,
                  cell = dot(array(cpuLayout), basis),
                  pbc = (1,1,1))
    atoms.set_atomic_numbers(atomicnumber*ones(len(atoms),int))
    atoms = MakeParallelAtoms(atoms, cpuLayout)

    # Run the parallel timings (big systems)
    MaxwellBoltzmannDistribution(atoms, 300 * units.kB)
    name = "Verlet%d-big" % (atomspercpu,)
    atoms.set_calculator(EMT())
    dyn = VelocityVerlet(atoms, 5*units.fs)
    ParallelTiming(name, dyn.run, 4*steps, atoms.get_number_of_atoms(),
                   verlet, bigtable)
    name = "Langevin%d-big" % (atomspercpu,)
    dyn = Langevin(atoms, 5*units.fs, 300*units.kB, 0.001)
    ParallelTiming(name, dyn.run, 4*steps, atoms.get_number_of_atoms(),
                   langevin, bigtable)
    #name = "NPT%d-big" % (atomspercpu,)
    #dyn = NPT(atoms, 5*units.fs, 300*units.kB, 0, 25*units.fs,
    #          (75*units.fs)**2 * 140*units.GPa)
    #ParallelTiming(name, dyn.run, 4*steps, atoms.get_number_of_atoms(),
    #               npt, bigtable)

    # Make parallel systems without repeating
    if ismaster:
        atoms = Atoms(positions=positions,
                      cell = basis,
                      pbc = (1,1,1))
        atoms.set_atomic_numbers(atomicnumber*ones(len(atoms),int))
    else:
        atoms = None
    try:
        atoms = MakeParallelAtoms(atoms, cpuLayout)
        atoms.set_calculator(EMT())
        atoms.get_potential_energy()
    except Exception:
        print "*** SKIPPING TESTS - TOO FEW ATOMS/CPU ***"
        fasttable.write(("%-5s      "*2) % (("`-`",)*2))
        continue

    # Run the parallel timings (fast systems)
    MaxwellBoltzmannDistribution(atoms, 300 * units.kB)
    name = "Verlet%d-fast" % (atomspercpu,)
    atoms.set_calculator(EMT())
    dyn = VelocityVerlet(atoms, 5*units.fs)
    ParallelTiming(name, dyn.run, 4*steps, atoms.get_number_of_atoms(),
                   verlet, fasttable)
    name = "Langevin%d-fast" % (atomspercpu,)
    dyn = Langevin(atoms, 5*units.fs, 300*units.kB, 0.001)
    ParallelTiming(name, dyn.run, 4*steps, atoms.get_number_of_atoms(),
                   langevin, fasttable)
    #name = "NPT%d-fast" % (atomspercpu,)
    #dyn = HooverNPT(atoms, 5*units.fs, 300*units.kB, 0, 25*units.fs,
    #                (75*units.fs)**2 * 140*units.GPa)
    #ParallelTiming(name, dyn.run, 4*steps, atoms.get_number_of_atoms(),
    #               npt, fasttable)

if ismaster:
    bigtable.write("\n")
    fasttable.write("\n")
    
logger.close()
bigtable.close()
fasttable.close()
