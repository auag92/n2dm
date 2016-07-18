from ase.constraints import StrainFilter
from ase.optimize import MDMin
from asap3 import *
from ase.lattice.cubic import *
from asap3.testtools import ReportTest

print_version(1)
size = 5

atoms = FaceCenteredCubic(size=(size,size,size), symbol="Cu", pbc=True)
defsize = atoms.get_volume()
atoms.set_cell(atoms.get_cell() * 1.1, scale_atoms=True)
atoms.set_calculator(EMT())

def printvol(a):
    print "Volume:", a.get_volume(), " energy:",\
          atoms.get_potential_energy() + atoms.get_kinetic_energy(),\
          " stress:", atoms.get_stress()[:3]

f = StrainFilter(atoms, [1, 1, 1, 0, 0, 0])
opt = MDMin(f, logfile="/dev/null", dt=0.01/(atoms.get_cell()[0,0]))
printvol(atoms)
opt.attach(printvol, 10, atoms)
opt.run(0.01)
printvol(atoms)
print "Original vol:", defsize
print
stress = atoms.get_stress()
for i in range(6):
    ReportTest("Stress component %d (T=0)" % (i,), stress[0], 0.0, 1e-5)

atoms = FaceCenteredCubic(size=(size,size,size), symbol="Cu", pbc=True)
atoms.set_calculator(EMT())
dyn = Langevin(atoms, 10*units.fs, 500*units.kB, 0.02)
dyn.attach(printvol, 100, atoms)
dyn.run(500)

print

f = StrainFilter(atoms, [1, 1, 1, 0, 0, 0])
opt = MDMin(f, logfile="/dev/null", dt=0.01/(atoms.get_cell()[0,0]))
printvol(atoms)
opt.attach(printvol, 10, atoms)
opt.run(0.01)
printvol(atoms)
print "Original vol:", defsize
print

stress = atoms.get_stress()
for i in range(6):
    ReportTest("Stress component %d (T=500K)" % (i,), stress[0], 0.0, 1e-5)

ReportTest.Summary()
