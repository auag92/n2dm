from asapserial3 import *
emt = EMT()
from ase import *
from ase.lattice.cubic import *
atoms = FaceCenteredCubic(directions=[[1,-1,0], [1,1,-2], [1,1,1]],
                          size=(5,5,5), symbol='Cu', pbc=(1,1,1))

atoms.set_calculator(emt)
e = atoms.get_potential_energy()
print "Energy per atom:", e/len(atoms), "eV"
