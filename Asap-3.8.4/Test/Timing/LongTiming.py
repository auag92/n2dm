#! /usr/bin/env python

from Asap import *
from Asap.Dynamics.VelocityVerlet import VelocityVerlet
from Asap.Dynamics.Langevin import Langevin
from Asap.Setup.Lattice.Cubic import FaceCenteredCubic
import sys, cPickle, time, commands, os, re
import RandomArray
from Numeric import *
from Asap.testtools import ReportTest

# cpu time:  time.clock().   Wall clock time: time.time()

host = commands.getoutput("hostname")
timesteps = 100
dbfilename = "longtiming.dat"
selfcheckfilename = "timing-selfcheck.dat"
logfilename = "longtiming.log"
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
initial = FaceCenteredCubic([[1,0,0],[0,1,0],[0,0,1]], size=(30, 30, 30),
                            element="Cu")
ReportTest("Number of atoms", len(initial), 108000, 0)
r = initial.GetCartesianPositions()
r.flat[:] += 0.14 * sin(arange(3*len(initial)))
initial.SetCartesianPositions(r)

print "Running self-test."
atoms = ListOfAtoms(initial)
atoms.SetCalculator(EMT())
e = atoms.GetPotentialEnergies()
f = atoms.GetCartesianForces()
if os.access(selfcheckfilename, os.F_OK):
    olde, oldf = cPickle.load(open(selfcheckfilename))
    de = max(fabs(e - olde))
    df = max(fabs(f.flat - oldf.flat))
    print "Maximal deviation:  Energy", de, "  Force", df
    ReportTest("Max force error", df, 0.0, 1e-11)
    ReportTest("Max energy error", de, 0.0, 1e-11)
    del olde, oldf
else:
    print "WARNING: No self-check database found, creating it."
    cPickle.dump((e, f), open(selfcheckfilename, "w"))
del e,f,atoms

ReportTest.Summary(exit=1)

print "Preparing to run Verlet dynamics (free boundaries)."
atoms = ListOfAtoms(initial, periodic=(0,0,0))
atoms.SetCalculator(EMT())
dynamics = VelocityVerlet(atoms, 5*femtosecond)

print "Running Verlet dynamics."
startcpu, startwall = time.clock(), time.time()
dynamics.Run(timesteps)

vfcpu, vfwall = time.clock() - startcpu, time.time() - startwall
vffraction = vfcpu/vfwall
sys.stderr.write("\n")
print "Verlet dynamics done."
del dynamics, atoms

print "Preparing to run Verlet dynamics (periodic boundaries)."
atoms = ListOfAtoms(initial, periodic=(1,1,1))
atoms.SetCalculator(EMT())
dynamics = VelocityVerlet(atoms, 5*femtosecond)

print "Running Verlet dynamics."
startcpu, startwall = time.clock(), time.time()
dynamics.Run(timesteps)

vpcpu, vpwall = time.clock() - startcpu, time.time() - startwall
vpfraction = vpcpu/vpwall
sys.stderr.write("\n")
print "Verlet dynamics done."
del dynamics, atoms

print "Preparing to run Langevin dynamics (periodic boundaries)."
atoms = ListOfAtoms(initial, periodic=(1,1,1))
atoms.SetCalculator(EMT())
dynamics = Langevin(atoms, 5*femtosecond, 400*kB, 0.001)

print "Running Langevin dynamics."
startcpu, startwall = time.clock(), time.time()
dynamics.Run(timesteps)

lcpu, lwall = time.clock() - startcpu, time.time() - startwall
lfraction = lcpu/lwall
sys.stderr.write("\n")
print "Langevin dynamics done."
del dynamics, atoms

print ""
print ""
print "TIMING RESULTS:"
print "Verlet (free):  CPU time %.2fs  Wall clock time %.2fs (%.0f%%)" % (vfcpu, vfwall, vffraction * 100)
print "Verlet (pbc):   CPU time %.2fs  Wall clock time %.2fs (%.0f%%)" % (vpcpu, vpwall, vpfraction * 100)
print "Langevin (pbc): CPU time %.2fs  Wall clock time %.2fs (%.0f%%)" % (lcpu, lwall, lfraction * 100)
print ""

key = (host, timesteps)
if os.access(dbfilename, os.F_OK):
    database = cPickle.load(open(dbfilename))
    founddb = 1
else:
    print "WARNING: No database with previous results ("+dbfilename+") - creating it."
    database = {}
    founddb = 0
    
try:
    best = database[key]
except:
    isbetter = (1,1,1)
    best = {}
    print "No old data to compare with."
else:
    print "Best timing so far on "+host+":"
    print "Verlet (free):  CPU time %.2fs  Wall clock time %.2fs (%.0f%%)" % (best["vfcpu"], best["vfwall"], best["vffraction"] * 100)
    print "Obtained "+best["vfwhen"]
    print best["vfversion"]
    print ""
    print "Verlet (pbc):   CPU time %.2fs  Wall clock time %.2fs (%.0f%%)" % (best["vpcpu"], best["vpwall"], best["vpfraction"] * 100)
    print "Obtained "+best["vpwhen"]
    print best["vpversion"]
    print ""
    print "Langevin (pbc): CPU time %.2fs  Wall clock time %.2fs (%.0f%%)" % (best["lcpu"], best["lwall"], best["lfraction"] * 100)
    print "Obtained "+best["lwhen"]
    print best["lversion"]
    print ""
    vfimp = (best["vfcpu"] - vfcpu) * 100.0 / best["vfcpu"]
    vpimp = (best["vpcpu"] - vpcpu) * 100.0 / best["vpcpu"]
    limp = (best["lcpu"] - lcpu) * 100.0 / best["lcpu"]
    print "Improvement (CPU times):"
    print "  Verlet: %.1f%%, %.1f%%   Langevin: %.1f%%" % (vfimp, vpimp, limp)
    print ""
    isbetter = (vfimp > 0.0, vpimp > 0.0, limp > 0.0)

if (vffraction < 0.95 or vpfraction < 0.95 or lfraction < 0.95):
    print "These data are not quite reliable - got less than 95% of the CPU time."
    print "Data will not be saved, even if they are best."
    isbetter = (0,0,0)

if isbetter != (0,0,0):
    print "Saving improved result"
    value = best
    if isbetter[0]:
        best["vfwhen"] = when
        best["vfcpu"] = vfcpu
        best["vfwall"] = vfwall
        best["vffraction"] = vffraction
        best["vfversion"] = asapversion
    if isbetter[1]:
        best["vpwhen"] = when
        best["vpcpu"] = vpcpu
        best["vpwall"] = vpwall
        best["vpfraction"] = vpfraction
        best["vpversion"] = asapversion
    if isbetter[2]:
        best["lwhen"] = when
        best["lcpu"] = lcpu
        best["lwall"] = lwall
        best["lfraction"] = lfraction
        best["lversion"] = asapversion
    database[key] = best
    newfilename = dbfilename + ".new"
    bakfilename = dbfilename + ".bak"
    cPickle.dump(database, open(newfilename, "w"), 1)
    if founddb:
        os.rename(dbfilename, bakfilename)
    os.rename(newfilename, dbfilename)

now=time.strftime("%Y/%m/%d %H:%M")
version = asapversion.split()[2]
compiler = asapversion.split("'")[1]

percent = 100.0 * (vfcpu+vpcpu+lcpu) / (vfwall+vpwall+lwall) 

logline = ( """%s %s
%3.0f%% cpu=%6.2f wall=%6.2f    %6.2f %6.2f %6.2f
ASAP ver %s  %s
"""  % (now, fullhost,
        percent, vfcpu+vpcpu+lcpu, vfwall+vpwall+lwall, vfcpu, vpcpu, lcpu,
        version, compiler))

print "\nLog line:"
print logline

logfile = open(logfilename, "a")
logfile.write(logline+"\n")
logfile.close()
