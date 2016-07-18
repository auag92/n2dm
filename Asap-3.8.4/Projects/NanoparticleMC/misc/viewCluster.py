"""generates atoms object from .smc file and views in ase gui
usage: python viewCluster.py smcfile
where smcfile is the full path to .smc file
"""
import sys
from montecarlo import *
from ase.cluster.cubic import FaceCenteredCubic
import ase
from data import fcc
from atommontecarlodata import AtomMonteCarloData
from ase.visualize import view

if len(sys.argv) < 2:
	print >>sys.stderr, __doc__
    	sys.exit(1)
#Store path to file
pp = str(sys.argv[1])
nconf= int(sys.argv[2])
#Instantiate d as SurfaceM.. object
d = SurfaceMonteCarloData()
#Read data from file to d
d.read(pp)
surfaces = fcc.surface_names
#Construct atoms
atoms = FaceCenteredCubic('Au', surfaces, d[nconf][1], latticeconstant=4.055)
view(atoms) #View atoms
#sys.exit("Script Terminated with code 0")
