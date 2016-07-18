"""Test for memory leaks in the CNA code."""

import asap3
from ase.lattice.compounds import L1_2
import numpy as np

steps = 30
msteps = 5

def makeatoms():
    atoms = L1_2(size=(10,10,10), symbol=('Au', 'Cu'), latticeconstant=4.0)
    r = atoms.get_positions()
    r += np.random.normal(0.0, 0.0001, r.shape)
    atoms.set_positions(r)
    return atoms

m0 = asap3.heap_mallinfo()
if m0 < 0:
    print "Memory monitoring not supported, test skipped."
else:
    print "Memory at startup", m0, 'kB'
    
    for type in ('CNA', 'COORD'):
        for haspot in (False, True):
            for replaceatoms in (False, True):
                print ("Running test: Type=%s HasPot=%s ReplaceAtoms=%s"
                       % (type, haspot, replaceatoms))
                leak = 0
                if not replaceatoms:
                    atoms = makeatoms()
                    if haspot:
                        atoms.set_calculator(asap3.EMT())
                        atoms.get_potential_energy()
                for i in range(steps):
                    if replaceatoms:
                        atoms = makeatoms()
                        if haspot:
                            atoms.set_calculator(asap3.EMT())
                            atoms.get_potential_energy()
                    if type == 'CNA':
                        asap3.CNA(atoms)
                    else:
                        asap3.CoordinationNumbers(atoms)
                   
                    m = int(asap3.heap_mallinfo())
                    if i == msteps:
                        m0 = m
                    if i % msteps == 0 and i > msteps:
                        print "    Memory usage:", m, "kB"
                        if m > m0:
                            leak = max(leak, m - m0)
                if leak > 1:
                    print "      *** MEMORY LEAK DETECTED ***"

                            
                            
                            