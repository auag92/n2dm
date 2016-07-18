print 
print "This test runs Verlet dynamics with the MonteCarloEMT potential instead"
print "of the usual EMT potential.  The result must be the same, but the"
print "performance will be slightly worse."
print 

from Asap import *
from Asap.Dynamics.VelocityVerlet import VelocityVerlet
from cPickle import *
from Numeric import *
from Asap.testtools import ReportTest

PrintVersion(1)

data = load(file("testVerlet.pickle"))
init_pos = array(data["initial"])
init_pos.shape = (-1,3)
init_box = array(data["box"])
init_box.shape = (3,3)
atoms = ListOfAtoms(positions=init_pos, cell=init_box)
atoms.SetAtomicNumbers(47)
atoms.SetCalculator(MonteCarloEMT())

dyn = VelocityVerlet(atoms, 2 * femtosecond)


for i in range(20):
    dyn.Run(10)
    epot = atoms.GetPotentialEnergy() / len(atoms)
    ekin = atoms.GetKineticEnergy() / len(atoms)
    print "E_pot = %-12.5f  E_kin = %-12.5f  E_tot = %-12.5f" % (epot, ekin,
                                                                 epot+ekin)
final_pos = array(data["final"])
diff = max(abs(atoms.GetCartesianPositions().flat - final_pos))
print "Maximal deviation of positions:", diff
ReportTest("Maximal deviation of positions", diff, 0, 1e-9)

diff = max(abs(atoms.GetStresses().flat - array(data["stress"])))

print "Maximal deviation of stresses:", diff
ReportTest("Maximal deviation of stresses", diff, 0, 1e-9)
ReportTest.Summary()

