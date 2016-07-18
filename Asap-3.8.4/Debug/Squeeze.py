from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
import numpy as np

set_verbose(1)

#AsapThreads()

n = 100

atoms = FaceCenteredCubic(symbol='Cu', size=(10,10,10))
atoms.set_calculator(EMT())
uc = atoms.get_cell()
for i in range(n+1):
    factor = float(n - i) / n
    print "Squeeze factor", factor
    uc2 = np.array(uc)
    uc2[2,2] *= factor
    atoms.set_cell(uc2, scale_atoms=True)
    #nbl = NeighborList(10.0, atoms)
    print atoms.get_potential_energy()
    
    