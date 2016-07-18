#!/usr/bin/env python

#PBS -N testNPT
#PBS -m ae
#PBS -q verylong
#PBS -l nodes=1

from Numeric import *
from Scientific.Statistics import *
from Asap import *
from Asap.Dynamics.NPTdynamics import *
from Asap.Dynamics.Langevin import *
from Asap.Trajectories.NetCDFTrajectory import *
from Asap.Setup.Lattice.FCCOrtho import FCCOrtho
from ASE.ChemicalElements import Element
from ASE.Utilities.Tests import *
import sys

T_goal = 300 # K
p_goal = 1 # GPa
bulk = 137 # Gpa
cp = 24.43 # J / (K * mol)
step1 = 100
step2 = 50

zapfiles = []

if len(sys.argv) == 1 or sys.argv[0].endswith("TestAll.py"):
    out1 = "testNPT-1.out"
    out2 = "testNPT-2.out"
    zapfiles.extend([out1,out2])
    ttime = 25
    ptime = 75
elif len(sys.argv) == 3:
    ttime = float(sys.argv[1])
    ptime = float(sys.argv[2])
    out1 = "testNPT-%d-%d-1.out" % (ttime, ptime)
    out2 = "testNPT-%d-%d-2.out" % (ttime, ptime)
else:
    raise RuntimeError, "Expected zero or two parameter to the script!"

if True:
    atoms = FCCOrtho(((1,0,0), (0,1,0), (0,0,1)),
                     (15,15,15), Element("Cu"), symmetry=(1,1,1))
    atoms.SetCalculator(EMT())
    print "Number of atoms:", len(atoms)

    print "Heating to %d K using Langevin" % T_goal
    lgv = Langevin(atoms, 5 * femtosecond, temperature=2*T_goal*kB, friction=0.01)

    while atoms.GetKineticEnergy() < 1.5 * len(atoms) * T_goal * kB:
        lgv.Run(10)
        T = atoms.GetKineticEnergy() / (1.5 * len(atoms) * kB)
        print "Temperature is now %.2f K" % (T,)
    print "Desired temperature reached!"

    lgv.SetTemperature(T_goal*kB)

    for i in range(2):
        lgv.Run(50)
        s = atoms.GetStress()
        p = -(s[0] + s[1] + s[2])/3.0 / gigapascal
        T = atoms.GetKineticEnergy() / (1.5 * len(atoms) * kB)
        print "Pressure is %f GPa, desired pressure is %f GPa (T = %.2f K)" % (p, p_goal, T)
        dv = (p - p_goal) / bulk
        print "Adjusting volume by", dv
        cell = atoms.GetUnitCell()
        atoms.SetUnitCell(cell * (1.0 + dv/3.0))

    T = atoms.GetKineticEnergy() / (1.5 * len(atoms) * kB)
    print "Temperature is now %.2f K" % (T,)

    stressstate = array([-2, -1, 0, 0, 0, 0])*p_goal*gigapascal
    dyn = HooverNPT(atoms, 5 * femtosecond, T_goal*kB, stressstate,
                    ttime*femtosecond, (ptime*femtosecond)**2 * bulk * gigapascal)
    traj = NetCDFTrajectory("NPT-atoms-traj.nc", atoms, interval = 50,
                            singleprecision = False)
    dyntraj = HooverNPTTrajectory("NPT-dyn-traj.nc", dyn, interval = 50)
    zapfiles.extend(["NPT-atoms-traj.nc", "NPT-dyn-traj.nc"])
    dyn.Attach(traj)
    dyn.Attach(dyntraj)

    out = open(out1, "w")

    temp = []
    pres = []
    vol = []
    for i in xrange(step1):
        dyn.Run(5)
        T = atoms.GetKineticEnergy() / (1.5 * len(atoms) * kB)
        s = atoms.GetStress() / gigapascal
        p = -(s[0] + s[1] + s[2])/3.0
        out.write("%5.2f %5.2f %7.5f      %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f\n" %
                  ((dyn.GetTime(), T, p)+tuple(s)))
        out.flush()
        cell = atoms.GetUnitCell()
        v = cell[0,0] * cell[1,1] * cell[2,2]
        temp.append(T)
        pres.append(p)
        vol.append(v)
        if i % 10 == 0:
            print i,"/",step1

    del atoms, dyn, traj, dyntraj

atoms = NetCDFTrajectory("NPT-atoms-traj.nc").GetListOfAtoms(frame=5)
atoms.SetCalculator(EMT())
dyn = HooverNPTTrajectory("NPT-dyn-traj.nc").GetDynamics(frame=5)
dyn.AttachAtoms(atoms)

out = open(out2, "w")

stress = []
for i in xrange(step2):
    dyn.Run(5)
    T = atoms.GetKineticEnergy() / (1.5 * len(atoms) * kB)
    s = atoms.GetStress() / gigapascal
    p = -(s[0] + s[1] + s[2])/3.0
    out.write("%5.2f %5.2f %7.5f      %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f\n" %
              ((dyn.GetTime(), T, p)+tuple(s)))
    out.flush()
    stress.append(s)
    if i % 10 == 0:
        print i,"/",step2,"(part 2)"

for f in zapfiles:
    print "Deleting", f
    os.unlink(f)
    
