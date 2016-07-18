#!/usr/bin/env python

#PBS -N HooverNPT
#PBS -m ae
#PBS -q long
#PBS -l nodes=1:opteron

from numpy import *
from Scientific.Statistics import *
from asap3 import *
from asap3.md.npt import *
from asap3.md.langevin import *
from ase.md.velocitydistribution import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.testtools import ReportTest
import sys

T_goal = 300 # K
p_goal = 1 # GPa
bulk = 137 # Gpa
cp = 24.43 # J / (K * mol)
beta = 3 * 16.5e-6 # 1/K
#step1 = 10000
step1 = 10000
step2 = 10000

trajfile = "npt.bundle"

if len(sys.argv) == 1:
    out1 = "testNPT-1.out"
    out2 = "testNPT-2.out"
    ttime = 50
    ptime = 75
elif len(sys.argv) == 3:
    ttime = float(sys.argv[1])
    ptime = float(sys.argv[2])
    out1 = "testNPT-%d-%d-1.out" % (ttime, ptime)
    out2 = "testNPT-%d-%d-2.out" % (ttime, ptime)
else:
    raise RuntimeError, "Expected zero or two parameter to the script!"

#ReportTest.Verbose()

atoms = FaceCenteredCubic(directions=((1,0,0), (0,1,0), (0,0,1)),
                          size=(15,15,15), symbol="Cu", pbc=True)
atoms.set_calculator(EMT())
print "Number of atoms:", len(atoms)

print "Heating to %d K using Langevin" % T_goal
lgv = Langevin(atoms, 5 * units.fs, temperature=2*T_goal*units.kB,
               friction=0.01)
MaxwellBoltzmannDistribution(atoms, 0.9*T_goal*units.kB)

while atoms.get_kinetic_energy() < 1.5 * len(atoms) * T_goal * units.kB:
    lgv.run(10)
    T = atoms.get_kinetic_energy() / (1.5 * len(atoms) * units.kB)
    print "Temperature is now %.2f K" % (T,)
print "Desired temperature reached!"

lgv.set_temperature(T_goal*units.kB)

for i in range(4):
    lgv.run(100)
    s = atoms.get_stress()
    p = -(s[0] + s[1] + s[2])/3.0 / units.GPa
    T = atoms.get_kinetic_energy() / (1.5 * len(atoms) * units.kB)
    print "Pressure is %f GPa, desired pressure is %f GPa (T = %.2f K)" % (p, p_goal, T)
    dv = (p - p_goal) / bulk
    print "Adjusting volume by", dv
    cell = atoms.get_cell()
    atoms.set_cell(cell * (1.0 + dv/3.0), scale_atoms=True)
    
T = atoms.get_kinetic_energy() / (1.5 * len(atoms) * units.kB)
print "Temperature is now %.2f K" % (T,)

dyn = NPT(atoms, 5 * units.fs, T_goal * units.kB, p_goal * units.GPa,
          ttime * units.fs, (ptime*units.fs)**2 * bulk * units.GPa)
traj = BundleTrajectory(trajfile, "w", atoms)
dyn.attach(traj, interval=50)

out = open(out1, "w")

temp = []
pres = []
vol = []
for i in xrange(step1):
    dyn.run(5)
    T = atoms.get_kinetic_energy() / (1.5 * len(atoms) * units.kB)
    s = atoms.get_stress() / units.GPa
    p = -(s[0] + s[1] + s[2])/3.0
    out.write("%5.2f  %7.5f      %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f\n" %
              ((T, p)+tuple(s)))
    out.flush()
    cell = atoms.get_cell()
    v = cell[0,0] * cell[1,1] * cell[2,2]
    temp.append(T)
    pres.append(p)
    vol.append(v)
    if i % 10 == 0:
        print i,"/",step1

temp = array(temp[len(temp)/2:])
pres = array(pres[len(pres)/2:])
vol = array(vol[len(vol)/2:])

print "Average temperature: %.5f K  (standard deviation: %.5g K)" % (mean(temp), standardDeviation(temp))
ReportTest("Average temperature", mean(temp), T_goal,
           4*standardDeviation(temp))
CV = cp * len(atoms) / 6.02205e23 / 1.60219e-19 # eV/K
expected_stddev_energy = sqrt(units.kB * T * T * CV)  # eV
expected_stddev_temp = sqrt(units.kB * T * T / CV) # K
#ReportTest("Standard deviation of temperature", standardDeviation(temp),
#           expected_stddev_energy / (3 * len(atoms) * units.kB), 0.05, True)
ReportTest("Standard deviation of temperature",
           standardDeviation(temp[len(temp)/2:]),
           expected_stddev_temp,
           0.5*expected_stddev_temp)

ReportTest("Average pressure", mean(pres), p_goal, 0.01, True)
expected_stddev_pressure = sqrt(units.kB * T * bulk * units.GPa / mean(vol))
expected_stddev_volume = sqrt(units.kB * T * mean(vol) / (bulk * units.GPa))
ReportTest("Standard deviation of pressure",
           standardDeviation(pres[len(pres)/2:]),
           expected_stddev_pressure/units.GPa,
           0.2*expected_stddev_pressure/units.GPa)
ReportTest("Standard deviation of volume",
           standardDeviation(vol),
           expected_stddev_volume,
           0.2*expected_stddev_volume)

dyn = NPT(atoms, 5 * units.fs, T_goal*units.kB,
          array([-2, -1, 0, 0, 0, 0])*p_goal*units.GPa,
          25*units.fs, (75*units.fs)**2 * bulk)

out = open(out2, "w")

stress = []
for i in xrange(step2):
    dyn.run(5)
    T = atoms.get_kinetic_energy() / (1.5 * len(atoms) * units.kB)
    s = atoms.get_stress() / units.GPa
    out.write("%5.2f      %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f\n" %
              ((T, )+tuple(s)))
    out.flush()
    stress.append(s)
    if i % 10 == 0:
        print i,"/",step2,"(part 2)"
    
stress = array(stress[len(stress)/2:])
s0 = mean(stress[:,0])
s1 = mean(stress[:,1])
s2 = mean(stress[:,2])
ReportTest("Stress component xx", s0, -2*p_goal, 0.05*p_goal)
ReportTest("Stress component yy", s1, -p_goal, 0.05*p_goal)
ReportTest("Stress component zz", s2, 0.0, 0.05*p_goal)
ReportTest.Summary()
