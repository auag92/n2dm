#! /usr/bin/env python

from Asap import *
from Asap.Dynamics.VelocityVerlet import VelocityVerlet
from Asap.Dynamics.Langevin import Langevin
from Asap.Setup.Lattice.Cubic import FaceCenteredCubic
import sys, cPickle, time, commands, os, re
import RandomArray
from Numeric import *
from Asap.testtools import ReportTest
import ASE.ChemicalElements.mass

# cpu time:  time.clock().   Wall clock time: time.time()

AsapThreads()

host = commands.getoutput("hostname")
timesteps = 5
dbfilename = "bigtiming.dat"
selfcheckfilename = "bigtiming-selfcheck.dat"
logfilename = "bigtiming.log"
asapversion = GetVersion()
when = time.strftime("%a %d %b %Y %H:%M", time.localtime(time.time()))

RandomArray.seed(42, 12345)

PrintVersion(1)
print "Running ASAP timing on "+host+"."
if re.match("^n\d\d\d.dcsc.fysik.dtu.dk$", host):
    print "    This is a d512 node on Niflheim."
    fullhost = "niflheim-d512/%s" % (host.split(".")[0])
    host = "niflheim-d512"
elif re.match("^[stu]\d\d\d.dcsc.fysik.dtu.dk$", host):
    print "    This is an s50 node on Niflheim."
    fullhost = "niflheim-s50/%s" % (host.split(".")[0])
    host = "niflheim-s50"
else:
    fullhost = host
print "Current time is "+when
print ""

print "Preparing system"
initial = FaceCenteredCubic([[1,0,0],[0,1,0],[0,0,1]], size=(50, 50, 50),
                            element="Cu")
#initial = FaceCenteredCubic([[1,0,0],[0,1,0],[0,0,1]], size=(25, 25, 25),
#                            element="Cu")
print "Number of atoms:", len(initial)
ReportTest("Number of atoms", len(initial), 500000, 0)
r = initial.GetCartesianPositions()
r.flat[:] += 0.14 * sin(arange(3*len(initial)))
initial.SetCartesianPositions(r)



print "Preparing to run Verlet dynamics (periodic boundaries)."
atoms = ListOfAtoms(initial, periodic=(1,1,1))
elements = [29]
epsilon  = [0.15]
sigma    = [2.7]
masses   = [ASE.ChemicalElements.mass.masses[29]]

atoms.SetCalculator(LJPotential(1, elements, epsilon, sigma, masses, -1.0,
                                True, False))
dynamics = VelocityVerlet(atoms, 5*femtosecond)

print "Running Verlet dynamics."
startcpu, startwall = time.clock(), time.time()
dynamics.Run(timesteps)

vpcpu, vpwall = time.clock() - startcpu, time.time() - startwall
vpfraction = vpcpu/vpwall
sys.stderr.write("\n")
print "Verlet dynamics done."
del dynamics, atoms


print ""
print ""
print "TIMING RESULTS:"
print "Verlet (pbc):   CPU time %.2fs  Wall clock time %.2fs (%.0f%%)" % (vpcpu, vpwall, vpfraction * 100)
print ""

