#import production

from asap3 import *
from asap3.md.verlet import VelocityVerlet
from cPickle import *
from numpy import *
#from Asap.Trajectories.NetCDFTrajectory import *
import sys, os, time
from asap3.testtools import *
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

def report():
    for i in range(world.size):
        if i == world.rank:
            print "Data on processor", i
            for key in atoms.arrays.keys():
                print "  ", key, atoms.arrays[key].shape
            r = atoms.get_positions()
            if len(r):
                print "Limits to the positions:"
                print ("[%.4f, %.4f]  [%.4f, %.4f]  [%.4f, %.4f]" %
                       (min(r[:,0]), max(r[:,0]), min(r[:,1]), max(r[:,1]),
                        min(r[:,2]), max(r[:,2])))
            if world.size > 1:
                print "Ghost data on processor", i
                for key in atoms.ghosts.keys():
                    print "  ", key, atoms.ghosts[key].shape
                r = atoms.ghosts['positions']
                if len(r):
                    print "Limits to the ghost positions:"
                    print ("[%.4f, %.4f]  [%.4f, %.4f]  [%.4f, %.4f]" %
                           (min(r[:,0]), max(r[:,0]), min(r[:,1]), max(r[:,1]),
                            min(r[:,2]), max(r[:,2])))
        world.barrier()


# Ensure same time step as in Asap-2
timeunit = 1.018047e-14             # Seconds
femtosecond = 1e-15 / timeunit      # Marginally different from units.fs

if isparallel:
    print "RUNNING PARALLEL VERSION OF TEST SCRIPT"
else:
    print "RUNNING SERIAL VERSION OF TEST SCRIPT"

if ismaster:
    # First, load a small system
    if os.path.exists("testVerlet.pickle"):
        data = load(file("testVerlet.pickle"))
    else:
        data = load(file("../testVerlet.pickle"))
    init_pos = array(data["initial"])
    init_pos.shape = (-1,3)
    init_box = array(data["box"])
    init_box.shape = (3,3)
    a = Atoms(positions=init_pos, cell=init_box, pbc=(1,1,1))
    a.set_atomic_numbers(29*ones(len(a)))
    #atoms = a.repeat((2,2,3))
    # Ensure same order of atoms as in ASE-2
    atoms = a.repeat((1,1,3)).repeat((1,2,1)).repeat((2,1,1))
    dx = 0.1 * sin(arange(3*len(atoms))/10.0)
    dx.shape = (-1,3)
    atoms.set_positions(atoms.get_positions() + dx)
    del dx
    out = sys.stderr
else:
    atoms = None
    out = file("/dev/null", "w")

if isparallel:
    atoms = MakeParallelAtoms(atoms, cpulayout)
    nTotalAtoms = atoms.get_number_of_atoms()
else:
    nTotalAtoms = len(atoms)

#report()

print "Setting potential"
atoms.set_calculator(EMT())

dyn = VelocityVerlet(atoms, 2 * femtosecond)
# if isparallel:
#     traj = ParallelNetCDFTrajectory("ptraj.nc", atoms, interval=20)
# else:
#     traj = NetCDFTrajectory("ptraj.nc", atoms, interval=20)
# dyn.Attach(traj)
# traj.Update()

print "Number of atoms:", nTotalAtoms

epot = atoms.get_potential_energy() / nTotalAtoms
ekin = atoms.get_kinetic_energy() / nTotalAtoms
etotallist = [epot+ekin]
ekinlist = [ekin]

#report()

if ismaster:
    print "\nE_pot = %-12.5f  E_kin = %-12.5f  E_tot = %-12.5f" % (epot, ekin,
                                                                 epot+ekin)
ReportTest("Initial potential energy", epot, 0.04312, 1e-4)
ReportTest("Initial kinetic energy", ekin, 0.0, 1e-9)

dyn.attach(MDLogger(dyn, atoms, "-", peratom=True), interval=10)

for i in range(40):
    dyn.run(10)
    epot = atoms.get_potential_energy() / nTotalAtoms
    ekin = atoms.get_kinetic_energy() / nTotalAtoms
    etotallist.append(epot+ekin)
    ekinlist.append(ekin)

if ismaster:
    print "Average total energy:", sum(etotallist)/len(etotallist)
    print "Average kinetic energy:", sum(ekinlist)/len(ekinlist)

ReportTest("Agv. total energy", sum(etotallist)/len(etotallist), 0.04309,
           0.0001)
ReportTest("Agv. kinetic energy", sum(ekinlist)/len(ekinlist), 0.02165,
           0.002)
ReportTest.Summary(0)

