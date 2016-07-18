from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic

pbc = False
atoms = FaceCenteredCubic(symbol='Cu', size=(10,10,10), pbc=pbc)
atoms.set_calculator(EMT())

paramfile = open("configuration.txt", "w")
print >>paramfile, int(pbc)
print >>paramfile, len(atoms)
z = atoms.get_atomic_numbers()
r = atoms.get_positions()
for i in range(len(atoms)):
    print >>paramfile, z[i], r[i,0], r[i,1], r[i,2]
    
paramfile.close()
print "Energy:", atoms.get_potential_energy()
