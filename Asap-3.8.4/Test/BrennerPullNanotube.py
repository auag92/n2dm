from ase.structure import molecule
from ase.structure import nanotube
from ase.optimize import BFGS, MDMin
from asap3 import BrennerPotential
from asap3.testtools import ReportTest
import numpy as np

atoms = nanotube(5, 5, length=10)
atoms.set_calculator(BrennerPotential())
uc = atoms.get_cell()
l0 = uc[2,2]

e = -1e100
for eps in np.linspace(0, 0.5, 201):
    uc[2,2] = l0 * (1 + eps)
    atoms.set_cell(uc, scale_atoms=True)
    dyn = MDMin(atoms, logfile=None, dt=0.05)
    dyn.run(fmax=0.01, steps=100)
    epot = atoms.get_potential_energy()
    if epot > e:
        e = epot
    print "%.3f  %f" % (eps, epot)

print "Maximal energy:", e
ReportTest("Maximal energy", e, -969.0, 5.0)
if abs(e - -620.98) < 1.0:
    print "Bug 42 is back!"
