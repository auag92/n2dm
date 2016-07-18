from asap3 import *
import cPickle as pickle
from numpy import *
import sys, os, time
from asap3.testtools import *
from asap3.mpi import world
from ase.lattice.cubic import *

debug = 0
if debug == 1:
    DebugOutput("moveatomsLJ%d.log", nomaster=True)
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
    
picklefile = "moveatomsLJ.pickle"

if isparallel:
    print "RUNNING PARALLEL VERSION OF TEST SCRIPT"
else:
    print "RUNNING SERIAL VERSION OF TEST SCRIPT"

if ismaster:
    atoms = FaceCenteredCubic(size=(10,10,10), symbol="Cu",
                              latticeconstant=1.09*2.7*1.41,
                              pbc=True)
else:
    atoms = None

#set_verbose(1)

if isparallel:
    atoms = MakeParallelAtoms(atoms, cpulayout, distribute=False)
    atoms.ghosts["ID"] = zeros(0,int)
    print "Distributing"
    atoms.distribute()
    nTotalAtoms = atoms.get_number_of_atoms()
else:
    nTotalAtoms = len(atoms)
    atoms.arrays["ID"] = arange(nTotalAtoms)

print "Setting potential"
atoms.set_calculator(LennardJones([29], [0.15], [2.7], 2.7*1.9, True))

try:
    expected = pickle.load(open(picklefile))
except IOError:
    expected = None
    newexpected = []


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


for i in range(10):
    print "Step number", i
    if i:
        dz = 0.1 * sin(atoms.arrays["ID"] / 10.0)
        r = atoms.get_positions()
        r[:,2] += dz
        atoms.set_positions(r)
    e = atoms.get_potential_energy() / nTotalAtoms
    print "\n%8f" % (e,)
    if expected:
        ReportTest("Step "+str(i), e, expected[i], 1e-9)
    else:
        newexpected.append(e)
    #report()
    #print atoms.get_potential_energies()
    
if expected is None and ismaster:
    print "Saving energies"
    out = open(picklefile, "w")
    pickle.dump(newexpected, out)
    out.close()
    
ReportTest.Summary()

