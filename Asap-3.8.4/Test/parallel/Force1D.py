from asap3 import *
from asap3.mpi import world
from asap3.testtools import ReportTest
from ase.visualize import view
import numpy as np
import time
import os
import sys

N = 20
d = 3.0
w = 10.0

cpulayout = (2,1,1)
isparallel = world.size > 1
ismaster = world.rank == 0

debug = 0
if debug == 1:
    DebugOutput("force1D_%d.log", nomaster=True)
elif debug == 2:
    time.sleep(world.rank)
    print "PID:", os.getpid()
    time.sleep(20)
    
sys.stderr.write("\n\n\n\n")

def main():
    positions = np.zeros((N, 3))
    positions[:,0] = d * np.arange(N)
    if ismaster:
        atoms = Atoms("Cu%i" % (N,), positions=positions, cell=(N*d, w, w),
                      pbc=False)
        atoms.center()
    else:
        atoms = None
    if isparallel:
        atoms = MakeParallelAtoms(atoms, cpulayout)
        id = atoms.get_ids()
    else:
        id = np.arange(N)
    #view(atoms)
    
    atoms.set_calculator(EMT())
    f = atoms.get_forces()
    x = atoms.get_positions()[:,0]
    # time.sleep(world.rank)
    for i in range(len(atoms)):
        print world.rank, x[i], f[i]
        
    # Now running actual tests - these are different on the two processors, but that should
    # not be a problem
    for i in range(len(f)):
        if id[i] > 1 and id[i] < N - 2:
            ftot = (f[i] * f[i]).sum()
            ReportTest("Force on atom %i" % (id[i],), ftot, 0.0, 1e-13)
    world.barrier()        
    ReportTest.Summary()
    
if world.size <= 2:
    main()
else:
    print "Skipping test when running on more than 2 processors."
    