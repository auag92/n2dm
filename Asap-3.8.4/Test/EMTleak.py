"""Test for memory leaks in the EMT code."""

import asap3
from ase.lattice.compounds import L1_2
import numpy as np

steps = 30
msteps = 5

def makeatoms(causetrouble):
    atoms = L1_2(size=(10,10,10), symbol=('Au', 'Cu'), latticeconstant=4.0)
    r = atoms.get_positions()
    r += np.random.normal(0.0, 0.0001, r.shape)
    if causetrouble:
        r[100] = r[110]
    atoms.set_positions(r)
    return atoms

m0 = asap3.heap_mallinfo()
if m0 < 0:
    print "Memory monitoring not supported, test skipped."
else:
    print "Memory at startup", m0, 'kB'
    
    for throwexception in (False, True):
        for replacepotential in (False, True):
            for replaceatoms in (False, True):
                print ("Running test: ReplacePotential=%s ReplaceAtoms=%s Exception=%s"
                       % (replacepotential, replaceatoms, throwexception))
                leak = 0
                if not replaceatoms:
                    atoms = makeatoms(throwexception)
                if not replacepotential:
                    pot = asap3.EMT()
                for i in range(steps):
                    if replaceatoms:
                        atoms = makeatoms(throwexception)
                    if replacepotential:
                        pot = asap3.EMT()
                    if i == 0 or replaceatoms or replacepotential:
                        atoms.set_calculator(pot)
                    try:
                        e = atoms.get_potential_energy()
                        #f = atoms.get_forces()
                    except asap3.AsapError:
                        pass
                    m = int(asap3.heap_mallinfo())
                    if i == msteps:
                        m0 = m
                    if i % msteps == 0 and i > msteps:
                        print "    Memory usage:", m, "kB"
                        if m > m0:
                            leak = max(leak, m - m0)
                if leak > 1:
                    print "      *** MEMORY LEAK DETECTED ***"

                            
                            
                            