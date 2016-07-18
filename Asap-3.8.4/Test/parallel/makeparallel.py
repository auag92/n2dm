from asap3 import *
from asap3.mpi import world
from ase.lattice.cubic import *
from asap3.testtools import ReportTest
import time
import os

debug = 0
if debug == 1:
    DebugOutput("makeparallel%d.log", nomaster=True)
elif debug == 2:
    time.sleep(world.rank)
    print "PID:", os.getpid()
    time.sleep(20)

if world.size == 1:
    layout = None
elif world.size == 2:
    layout = [2,1,1]
elif world.size == 3:
    layout = [1,3,1]
elif world.size == 4:
    layout = [2,1,2]

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


if world.size == 1 or world.rank == 1:
    atoms = FaceCenteredCubic(symbol="Cu", size=(20,20,30))
else:
    atoms = None

if world.size > 1:
    atoms = MakeParallelAtoms(atoms, layout, distribute=False)
    natoms = atoms.get_number_of_atoms()
else:
    natoms = len(atoms)
    
report()

print
print

if world.size > 1:
    atoms.distribute()

report()

print
print
print "Calculating total energy."
#if world.rank == 0:
#    set_verbose(1)

atoms.set_calculator(EMT())
e = atoms.get_potential_energy()

report()


ReportTest("Potential energy", e, -2.40461833565e-3 * natoms/4, 1e-5)

print "Exiting."
world.barrier()
#ReportTest.Summary()
