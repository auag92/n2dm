from numpy import *
#from Scientific.Statistics import *
from asap3 import *
from asap3.md.npt import NPT
from asap3.md.langevin import Langevin
from asap3.io.trajectory import PickleTrajectory
from ase.lattice.cubic import FaceCenteredCubic
from asap3.testtools import ReportTest
from asap3.mpi import world
from ase import units
import sys

T_goal = 300 # K
p_goal = 1 # GPa
bulk = 137 # Gpa
cp = 24.43 # J / (K * mol)
step1 = 50
step2 = 25

if len(sys.argv) == 1 or sys.argv[0].endswith("TestAll.py"):
    out1 = "testNPT-1.out"
    out2 = "testNPT-2.out"
    ttime = 25
    ptime = 75
elif len(sys.argv) == 3:
    ttime = float(sys.argv[1])
    ptime = float(sys.argv[2])
    out1 = "testNPT-%d-%d-1.out" % (ttime, ptime)
    out2 = "testNPT-%d-%d-2.out" % (ttime, ptime)
else:
    raise RuntimeError, "Expected zero or two parameter to the script!" + str(sys.argv)

ismaster = world.rank == 0
if world.size == 2:
    cpulayout = [2,1,1]
elif world.size == 4:
    cpulayout = [2,2,1]
elif world.size == 8:
    cpulayout = [2,2,2]
else:
    raise RuntimeError("Cannot run on %i CPUs" % (world.size,))

if ismaster:
    print_version(1)

if not hasattr(NPT, '_npt_version'):
    print "Skipping test: NP dynamics does not work in parallel with this old version of ASE."
else:
    if ismaster:
        atoms = FaceCenteredCubic(size=(30,15,15), symbol="Cu", pbc=True)
    else:
        atoms = None
    atoms = MakeParallelAtoms(atoms, cpulayout)
    atoms.set_calculator(EMT())
    print "Number of atoms:", atoms.get_number_of_atoms()

    print "Heating to %d K using Langevin" % T_goal
    lgv = Langevin(atoms, 5 * units.fs, temperature=2*T_goal*units.kB, friction=0.05)

    while atoms.get_kinetic_energy() < 1.5 * atoms.get_number_of_atoms() * T_goal * units.kB:
        lgv.run(5)
        T = atoms.get_kinetic_energy() / (1.5 * atoms.get_number_of_atoms() * units.kB)
        print "Temperature is now %.2f K" % (T,)
    print "Desired temperature reached!"

    lgv.set_temperature(T_goal*units.kB)

    for i in range(2):
        lgv.run(20)
        s = atoms.get_stress()
        p = -(s[0] + s[1] + s[2])/3.0 / units.GPa
        T = atoms.get_kinetic_energy() / (1.5 * atoms.get_number_of_atoms() * units.kB)
        print "Pressure is %f GPa, desired pressure is %f GPa (T = %.2f K)" % (p, p_goal, T)
        dv = (p - p_goal) / bulk
        print "Adjusting volume by", dv
        cell = atoms.get_cell()
        atoms.set_cell(cell * (1.0 + dv/3.0))

    T = atoms.get_kinetic_energy() / (1.5 * atoms.get_number_of_atoms() * units.kB)
    print "Temperature is now %.2f K" % (T,)

    stressstate = array([-2, -1, 0, 0, 0, 0])*p_goal*units.GPa
    dyn = NPT(atoms, 5 * units.fs, T_goal*units.kB, stressstate,
              ttime*units.fs, (ptime*units.fs)**2 * bulk * units.GPa)
    traj = PickleTrajectory("NPT-atoms.traj", "w", atoms)
    #dyntraj = ParallelHooverNPTTrajectory("NPT-dyn-traj.nc", dyn, interval = 50)
    dyn.attach(traj, interval=50)
    #dyn.Attach(dyntraj)

    out = open(out1, "w")

    temp = []
    pres = []
    vol = []
    for i in xrange(step1):
        dyn.run(5)
        T = atoms.get_kinetic_energy() / (1.5 * atoms.get_number_of_atoms() * units.kB)
        s = atoms.get_stress() / units.GPa
        p = -(s[0] + s[1] + s[2])/3.0
        out.write("%5.2f %5.2f %7.5f      %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f\n" %
                  ((dyn.get_time(), T, p)+tuple(s)))
        out.flush()
        cell = atoms.get_cell()
        v = cell[0,0] * cell[1,1] * cell[2,2]
        temp.append(T)
        pres.append(p)
        vol.append(v)
        if i % 10 == 0:
            print i,"/",step1

    del atoms, dyn, traj
    #del dyntraj

