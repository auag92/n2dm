"""Test force calculations alone with the purpose of optimizing threading."""

from asap3 import *
from asap3.Timing import report_timing
from ase.lattice.cubic import FaceCenteredCubic
import time

nsteps = 10

start = time.time()
if len(sys.argv) >= 2 and sys.argv[1] == '-t':
    AsapThreads()   
if len(sys.argv) >= 2 and sys.argv[1] == '-T':
    AsapThreads(6)   

atoms = FaceCenteredCubic(size=(100,100,50), symbol='Cu')
atoms.set_calculator(EMT())
print "Number of atoms:", len(atoms)

d = 0.1
pos = atoms.arrays['positions']  # Nasty!

for i in range(nsteps):
    pos[50][0] += d
    d = -d
    f = atoms.get_forces()
    
report_timing()

print "Wall time elapsed:", time.time() - start


