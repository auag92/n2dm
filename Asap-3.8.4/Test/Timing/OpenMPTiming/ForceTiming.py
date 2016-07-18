"""Test force calculations alone with the purpose of optimizing threading.

Usage:

python ForceTiming.py CMD file.pickle

where CMD is one of A (serial Alone), S (serial), O (OpenMP) or M (MPI).

The result is placed in the pickle file.
"""

from asap3 import *
from asap3.Timing import report_timing
from ase.lattice.cubic import FaceCenteredCubic
import time
import sys
import pickle
import numpy as np

print_version(1)

nsteps = 10

if len(sys.argv) != 3:
    print >> sys.stderr, __doc__
    print >> sys.stderr, "ERROR: Requires two arguments."
    sys.exit(1)
cmd = sys.argv[1]
if cmd not in "ASOM":
    print >> sys.stderr, __doc__
    print >> sys.stderr, "ERROR: First argument is incorrect."
    sys.exit(1)

if cmd in "SM":
    from asap3.mpi import world

try:
    omp_num_threads = int(os.environ["OMP_NUM_THREADS"])
except KeyError:
    omp_num_threads = None

print "OMP_NUM_THREADS =", omp_num_threads

size = 50 * np.ones(3, int)
if cmd == "A" or cmd == "S":
    layout = np.ones(3)
    ncpu = 1
elif cmd == "O":
    ncpu = int(os.environ["OMP_NUM_THREADS"])
elif cmd == "M":
    ncpu = world.size

if ncpu == 1:
    layout = (1,1,1)
elif ncpu == 2:
    layout = (2,1,1)
elif ncpu == 4:
    layout = (2,2,1)
elif ncpu == 8:
    layout = (2,2,2)
else:
    raise ValueError("Cannot run on %i CPUs." % ncpu)
   
if cmd == "M":
    if world.rank == 0:
        atoms = FaceCenteredCubic(size=size*layout, symbol='Cu')
    else:
        atoms = None
    atoms = MakeParallelAtoms(atoms, layout)
else:
    atoms = FaceCenteredCubic(size=size*layout, symbol='Cu')

atoms.set_calculator(EMT())
natoms = atoms.get_number_of_atoms()
print "Number of atoms:", natoms
assert natoms == ncpu * 500000
print "Potential energy:", atoms.get_potential_energy()
start = time.time()

d = 0.1

for i in range(nsteps):
    atoms.arrays['positions'][50][0] += d
    d = -d
    f = atoms.get_forces()
    
wall = time.time() - start
if cmd in "SM":
    master = world.rank == 0
    wall = world.sum(wall)
    wall /= world.size
else:
    master = True
    
if master:
    report_timing()
    print "Wall time elapsed:", wall
    f = open(sys.argv[2], "a")
    pickle.dump((cmd, ncpu, wall), f)
    f.close()
