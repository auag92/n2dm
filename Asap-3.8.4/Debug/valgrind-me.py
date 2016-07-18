# For valgrind, print something before module initialization

print "*****"
print "*****  The script is starting now by importing cPickle and Numeric."
print "*****"
from cPickle import *
from Numeric import *
print "*****"
print "*****  Now, ASAP is being imported."
print "*****"
from Asap import *
from Asap.Dynamics.VelocityVerlet import VelocityVerlet
from Asap.testtools import ReportTest

print "*****"
print "*****  The script is now running.  Valgrind should not complain any more."
print "*****"


PrintVersion(1)

data = load(file("testVerlet.pickle"))
init_pos = array(data["initial"])
init_pos.shape = (-1,3)
init_box = array(data["box"])
init_box.shape = (3,3)
atoms = ListOfAtoms(positions=init_pos, cell=init_box)
atoms.SetAtomicNumbers(47)
atoms.SetCalculator(EMT())

dyn = VelocityVerlet(atoms, 2 * femtosecond)


for i in range(3):
    dyn.Run(10)
    epot = atoms.GetPotentialEnergy() / len(atoms)
    ekin = atoms.GetKineticEnergy() / len(atoms)
    print "E_pot = %-12.5f  E_kin = %-12.5f  E_tot = %-12.5f" % (epot, ekin,
                                                                 epot+ekin)
print "*****"
print "*****  The script is now finished, and Python finalizes."
print "*****"

