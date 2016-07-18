from asap3 import *
from ase.lattice.cubic import *

atoms = FaceCenteredCubic(size=(10, 10, 10), symbol="Cu")
z = atoms.get_atomic_numbers()
z[7]=8
atoms.set_atomic_numbers(z)

try:
    atoms.set_calculator(EMT())
except AsapError:
    print "Test 1 passed."
else:
    print "Test 1 failed."
    
atoms = FaceCenteredCubic(size=(1, 10, 10), symbol="Cu")
atoms.set_calculator(EMT())
try:
    e = atoms.get_potential_energy()
except AsapError:
    print "Test 2 passed."
else:
    print "Test 2 failed."
    

atoms = FaceCenteredCubic(size=(10, 10, 10), symbol="Cu")
r = atoms.get_positions()
r[10] = r[11]
atoms.set_positions(r)
atoms.set_calculator(EMT())
try:
    e = atoms.get_potential_energy()
except AsapError:
    print "Test 3 passed."
else:
    print "Test 3 failed."
    