from ase.io import read
from asap3 import CoordinationNumbers
from ase.visualize import view
import sys

trajfile  = sys.argv[1]

atoms = read(trajfile)

atoms.set_tags(CoordinationNumbers(atoms))
print "Color Scheme in HUE:"
print "0 : (0,0,0)"
print "1 : (0,0,0)"
print "2 : (0,0,0)"
print "3 : (0,0,0)"
print "4 : Violet ()"
print "5 : (1,0,0)"
print "6 : Orange(1,0.75,0)"
print "7 : (0,1,0)"
print "8 : (0,0.5,0)"
print "9 : (1,0.4,0)"
print "10: (0,0,1)"
print "11: (0,0,0)"
print "12: (1,1,1)"
print "13: (0,0,0)"

view(atoms)
