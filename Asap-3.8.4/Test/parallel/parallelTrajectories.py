from asap3 import *
from ase.md.verlet import VelocityVerlet
from ase.lattice.cubic import FaceCenteredCubic
from asap3.io.trajectory import *
from numpy import *
import sys, os, time
from asap3.testtools import ReportTest
from asap3.mpi import world

debug = 0
if debug == 1:
    DebugOutput("makeverlet%d.log", nomaster=True)
elif debug == 2:
    time.sleep(world.rank)
    print "PID:", os.getpid()
    time.sleep(20)

print_version(1)
#set_verbose(1)

ismaster = world.rank == 0
isparallel = world.size != 1
if world.size == 1:
    cpulayout = None
elif world.size == 2:
    cpulayout = [2,1,1]
elif world.size == 3:
    cpulayout = [1,3,1]
elif world.size == 4:
    cpulayout = [2,1,2]

delete = True
precision = 1e-8

def maketraj(atoms, t, nstep):
    e = [atoms.get_potential_energy()]
    print "Shape of force:", atoms.get_forces().shape
    dyn = VelocityVerlet(atoms, 5*units.fs)
    for i in range(nstep):
        dyn.run(10)
        energy = atoms.get_potential_energy()
        e.append(energy)
        if ismaster:
            print "Energy: ", energy
        if t is not None:
            t.write()
    return e

def checktraj(t, e, cpus=None):
    i = 0
    for energy in e:
        atoms = t.get_atoms(i, cpus)
        atoms.set_calculator(EMT())
        ReportTest("Checking frame %d / cpus=%s" % (i, str(cpus)),
                   atoms.get_potential_energy(), energy, precision)
        i += 1

if ismaster:
    initial = FaceCenteredCubic(size=(10,10,10), symbol="Cu", pbc=(1,0,0))
else:
    initial = None
if isparallel:
    atoms = MakeParallelAtoms(initial, cpulayout)
else:
    atoms = initial.copy()
    
atoms.set_calculator(EMT())
print "Writing trajectory"
traj = PickleTrajectory("traj1.nc", "w", atoms)
traj.write()
energies = maketraj(atoms, traj, 10)
traj.close()

if ismaster:
    print "Reading trajectory (serial)"
    traj = PickleTrajectory("traj1.nc")
    checktraj(traj, energies)

if isparallel:
    print "Reading trajectory (parallel)"
    traj = PickleTrajectory("traj1.nc")
    checktraj(traj, energies, cpulayout)

print "Repeating simulation"
atoms = traj.get_atoms(5, cpulayout)
atoms.set_calculator(EMT())
energies2 = maketraj(atoms, None, 5)
if ismaster:
    for i in range(5):
        ReportTest("Rerun[%d]" % (i,), energies2[i], energies[i+5], precision)
traj.close()

print "Appending to trajectory"
traj = PickleTrajectory("traj1.nc", "a")
atoms = traj.get_atoms(-1, cpulayout)
atoms.set_calculator(EMT())
traj.set_atoms(atoms)
energies2 = maketraj(atoms, traj, 5)
traj.close()

if ismaster:
    print "Reading longer trajectory"
    traj = PickleTrajectory("traj1.nc")
    checktraj(traj, energies + energies2[1:])

if ismaster and delete:
    print "Deleting trajectory"
    os.unlink("traj1.nc")
    
ReportTest.Summary()
