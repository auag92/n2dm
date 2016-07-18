from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from ase.optimize import LBFGS
import numpy as np

calc = EMT()
symbol = 'Cu'

syms = ((True, True, True),
        (False, False, False),
        (True, True, False),
        (True, False, False))

for symmetry in syms:
    print "Periodicity:", symmetry
    for n in (4,6):
        atoms = FaceCenteredCubic(symbol=symbol, size=(4, 8, n), pbc=symmetry)
        atoms.set_calculator(calc)

        dyn = LBFGS(atoms, trajectory=None, logfile=None)
        dyn.run(fmax=0.05)

        e = atoms.get_potential_energy()
        print "Layers, energy: %2i %8.3f" % (n, e)

print "TEST PASSED."
