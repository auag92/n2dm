import ase
from ase import Atoms
from ase.atom import Atom
import sys
from ase.visualize import view
import pickle
f = open(sys.argv[1],'r') #The .amc file
p = pickle.load(f)

positions = p['atomspositions']
atms = Atoms()

for p0 in positions:
	a = Atom('Au',position=p0)
	atms.append(a)

atms.center(vacuum=2)
view(atms)
