"""Test the FixAtoms constraint and the Subset filter."""

from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.md.velocitydistribution import MaxwellBoltzmannDistribution
from asap3.md.verlet import VelocityVerlet
from asap3.io.trajectory import PickleTrajectory
from asap3.constraints import Filter, FixAtoms
from ase.constraints import FixAtoms as ASE_FixAtoms
from asap3.testtools import ReportTest
from asap3.mpi import world
import numpy as np

debug = 0
if debug == 1:
    DebugOutput("parallelconstraints%d.log", nomaster=True)
elif debug == 2:
    time.sleep(world.rank)
    print "PID:", os.getpid()
    time.sleep(20)

print_version(1)
#set_verbose(1)

ismaster = world.rank == 0
isparallel = world.size != 1
if world.size == 1:
    cpulayout = None
elif world.size == 2:
    cpulayout = [2,1,1]
elif world.size == 3:
    cpulayout = [1,3,1]
elif world.size == 4:
    cpulayout = [2,1,2]

if ismaster:
    init = FaceCenteredCubic(size=(10,10,10), symbol='Cu', pbc=False)
    z = init.get_positions()[:,2]
    fixedatoms = np.less(z, 0.501*z.max())
    print len(init), sum(fixedatoms)
    MaxwellBoltzmannDistribution(init, 6000*units.kB)
    init.set_tags(fixedatoms)
else:
    init = None

print
print "Running simulation with Filter"
atoms1 = MakeParallelAtoms(init, cpulayout)
atoms1.arrays['r_init'] = atoms1.get_positions()
atoms1.set_calculator(EMT())
atoms1a = Filter(atoms1, mask=np.logical_not(atoms1.get_tags()))

dyn = VelocityVerlet(atoms1a, 3*units.fs)
dyn.run(100)

print
print "Running simulation with Asap's FixAtoms"
atoms2 = MakeParallelAtoms(init, cpulayout)
atoms2.arrays['r_init'] = atoms2.get_positions()
atoms2.set_calculator(EMT())
atoms2.set_constraint(FixAtoms(mask=atoms2.get_tags()))

dyn = VelocityVerlet(atoms2, 3*units.fs)
dyn.run(100)

print
print "Running Langevin simulation with Asap's FixAtoms"
atoms4 = MakeParallelAtoms(init, cpulayout)
atoms4.arrays['r_init'] = atoms4.get_positions()
atoms4.set_calculator(EMT())
atoms4.set_constraint(FixAtoms(mask=atoms4.get_tags()))

dyn = Langevin(atoms4, 3*units.fs, 3000*units.kB, 0.01)
dyn.run(100)


print
sanity = [[atoms1, "Verlet + Filter"],
          [atoms2, "Verlet + Asap's FixAtoms"],
          [atoms4, "Langevin + Asap's FixAtoms"],
          ]
for a, label in sanity:
    print world.rank, "Sanity check, %s:", (label,)
    r_init = a.arrays['r_init']
    ok = (a.get_positions() == r_init) + np.logical_not(a.get_tags())[:,np.newaxis]
    if ok.all():
        print world.rank, "  Stationary atoms have not moved: OK"
    else:
        raise RuntimeError("Stationary atoms have moved (%s)" % (label,))
    ok = (a.get_positions() != r_init) + a.get_tags()[:,np.newaxis]
    if ok.all():
        print world.rank, "  Mobile atoms have moved: OK"
    else:
        raise RuntimeError("Mobile atoms have not moved (%s)" % (label,))
    

print "ALL TESTS SUCCEEDED!"



