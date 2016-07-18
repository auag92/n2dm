### Mail to user
#PBS -m ae

### Queue name (small, medium, long, verylong)
#PBS -q verylong

### Number of nodes
#PBS -l nodes=1:opteron


"""Runs a very long EMT Verlet dynamics on a small system.

Usage:
   python LongEMTsmall.py EMT
or
   python LongEMTsmall.py EMT2013
"""

from asap3 import *
from asap3.EMT2013Parameters import EMT2013_parameters
from asap3.md.langevin import Langevin
from asap3.md.verlet import VelocityVerlet
from asap3.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.lattice.cubic import FaceCenteredCubic
import sys

if len(sys.argv) == 2:
    arg = sys.argv[1]
else:
    arg = 'Wrong number of arguments.'
if arg == 'EMT':
    filename = 'small_EMT.dat'
elif arg == 'EMT2013':
    filename = 'small_EMT2013.dat'
else:
    print __doc__
    sys.exit(1)

    
T = 450
atoms = FaceCenteredCubic(size=(10,10,10), symbol='Cu', pbc=False)
if arg == 'EMT':
    atoms.set_calculator(EMT())
else:
    atoms.set_calculator(EMT2013(EMT2013_parameters))
    
print "Setting temperature:", T, "K"
MaxwellBoltzmannDistribution(atoms, 2*T*units.kB)
Stationary(atoms)

dyn1 = Langevin(atoms, 2*units.fs, temperature=T*units.kB, friction=0.05)
dyn1.run(200)

dyn = VelocityVerlet(atoms, 3*units.fs, logfile=filename, loginterval=25)
dyn.run(10000000) # 10^7 timesteps is 3 ns.

