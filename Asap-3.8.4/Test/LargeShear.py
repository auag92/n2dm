import numpy as np
from ase.all import *
from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.testtools import ReportTest

#set_verbose(1)
atoms = FaceCenteredCubic(symbol='Cu', size=(4,4,4))
atoms.set_calculator(EMT())
uc = atoms.get_cell()
#vals = (0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3)
#vals = linspace(-0.1444, 0.1444, 20)
vals = np.linspace(0, 0.2, 20)

res = []
for n,i in enumerate(vals):
    uc[0,2] = i
    atoms.set_cell(uc, scale_atoms=True)
    epot = atoms.get_potential_energy()
    print n, i, epot
    res.append(epot)

poly = np.polyfit(vals, res, 2)
fits = np.polyval(poly, vals)

rms = np.sqrt((fits - res) * (fits - res))
print rms
maxrms = rms.max()
ReportTest("Worse fit", maxrms, 0.0, 1e-5)
ReportTest.Summary()



    

