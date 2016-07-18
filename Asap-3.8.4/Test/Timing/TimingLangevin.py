#! /usr/bin/env python

#PBS -N Timing
#PBS -l nodes=1:ppn=4:opteron285
#PBS -q small

from numpy import *
from asap3 import *
from asap3.md.verlet import VelocityVerlet
from asap3.md.langevin import Langevin
from ase.lattice.cubic import FaceCenteredCubic
from asap3.Timing import report_timing
import sys, cPickle, time, commands, os, re
import numpy as np
from asap3.testtools import ReportTest

# cpu time:  time.clock().   Wall clock time: time.time()

#set_verbose(1)

usethread = (len(sys.argv) > 1 and
             (sys.argv[1] == "-t" or sys.argv[1] == "-T"))
if usethread:
    if sys.argv[1] == "-t":
        AsapThreads()
    else:
        AsapThreads(4)

host = commands.getoutput("hostname")
timesteps = 100
if usethread:
    dbfilename = "timing-thread.dat"
    logfilename = "timing-thread.log"
else:
    dbfilename = "timing.dat"
    logfilename = "timing.log"
selfcheckfilename = "timing-selfcheck.dat"
asapversion = get_version()
when = time.strftime("%a %d %b %Y %H:%M", time.localtime(time.time()))

randomstate = "randomstate.pickle"
if os.path.isfile(randomstate):
    np.random.set_state(cPickle.load(open(randomstate)))
else:
    print "Saving random state for next call."
    rndfile = open(randomstate, "w")
    cPickle.dump(np.random.get_state(), rndfile)
    rndfile.close()
    

#PrintVersion(1)
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
initial = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]],
                            size=(30, 30, 30),
                            symbol="Cu")
ReportTest("Number of atoms", len(initial), 108000, 0)
r = initial.get_positions()
r.flat[:] += 0.14 * sin(arange(3*len(initial)))
initial.set_positions(r)

print "Running self-test."
atoms = Atoms(initial)
atoms.set_calculator(EMT())
e = atoms.get_potential_energies()
f = atoms.get_forces()
if os.access(selfcheckfilename, os.F_OK):
    olde, oldf = cPickle.load(open(selfcheckfilename))
    de = max(fabs(e - olde))
    df = max(fabs(f.flat[:] - oldf.flat[:]))
    print "Maximal deviation:  Energy", de, "  Force", df
    ReportTest("Max force error", df, 0.0, 1e-11)
    ReportTest("Max energy error", de, 0.0, 1e-11)
    del olde, oldf
else:
    print "WARNING: No self-check database found, creating it."
    cPickle.dump((e, f), open(selfcheckfilename, "w"))
del e,f,atoms

ReportTest.Summary(exit=1)

print "Preparing to run Langevin dynamics."
atoms = Atoms(initial)
atoms.set_calculator(EMT())
dynamics = Langevin(atoms, 5*units.fs, 400*units.kB, 0.001)

print "Running Langevin dynamics."
startcpu, startwall = time.clock(), time.time()
dynamics.run(timesteps)

lcpu, lwall = time.clock() - startcpu, time.time() - startwall
lfraction = lcpu/lwall
sys.stderr.write("\n")
print "Langevin dynamics done."
print "Temperature:", atoms.get_temperature()
#del dynamics, atoms

print ""
print ""
print "TIMING RESULTS:"
print "Langevin: CPU time %.2fs  Wall clock time %.2fs (%.0f%%)" % (lcpu, lwall, lfraction * 100)
print ""


report_timing()
