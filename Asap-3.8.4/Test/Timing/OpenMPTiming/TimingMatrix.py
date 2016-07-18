#!/usr/bin/env python

#PBS -N Timing
#PBS -m ae
#PBS -q long
#PBS -l nodes=1:opteron:ppn=4

from asap3 import *
from asap3.md.verlet import VelocityVerlet
from asap3.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.lattice.cubic import FaceCenteredCubic
import numpy as np
import time
import pickle

print_version(1)

T = 300

threads = (1, 2, 4, 8)
s1 = (100, 1000, 10000, 100000, 1000000)
s2 = (1, 2, 5)
sizes = []
for a in s1:
    for b in s2:
        sizes.append(a*b)
print sizes

targettime = 30.0
laststeps = 5000
lastsize = sizes[0]
lasttime = targettime
results = {}
for nthreads in threads:
    try:
        AsapThreads(nthreads)
    except ValueError:
        break
    maxthread = nthreads
    for natoms in sizes:
        print "Timing with %i atoms (%i threads)" % (natoms, nthreads)
        blocksize = int(np.ceil((natoms/4)**(1./3.)))
        atoms = FaceCenteredCubic(symbol='Cu', size=(blocksize,blocksize,blocksize), pbc=False)
        print "Creating block with %i atoms, cutting to %i atoms" % (len(atoms), natoms)
        atoms = atoms[:natoms]
        assert len(atoms) == natoms
        atoms.set_calculator(EMT())
        MaxwellBoltzmannDistribution(atoms, 2 * T * units.kB)
        dyn = VelocityVerlet(atoms, 5*units.fs)
        ptsteps = int(laststeps * (0.1 * targettime / lasttime) * lastsize / natoms)
        if ptsteps < 100:
            ptsteps = 100
        print "Running pre-timing (%i steps)..." % (ptsteps,)
        t1 = time.time()
        dyn.run(ptsteps - 50)
        MaxwellBoltzmannDistribution(atoms, (2 * T - atoms.get_temperature()) * units.kB)
        dyn.run(50)
        t1 = time.time() - t1
        steps = int(ptsteps * targettime / t1)
        if steps < 200:
            steps = 200 
        print "Temperature is %.1f K" % (atoms.get_temperature(),)
        print "Running main timing (%i steps)" % (steps,)
        MaxwellBoltzmannDistribution(atoms, T * units.kB)
        t1 = time.time()
        dyn.run(steps)
        t1 = time.time() - t1
        lasttime = t1
        print "... done in %.1f s  (T = %.1f K)." % (t1, atoms.get_temperature())
        t1 *= 1e6 / (natoms * steps)
        print "RESULT: %.3f us/atom/step  (%i atoms, %i threads)" % (t1, natoms, nthreads)
        results[(natoms, nthreads)] = t1
        laststeps = steps
        lastsize = natoms
        
out = open("timing.pickle", "w")
pickle.dump(results, out)
pickle.dump(maxthread, out)
pickle.dump(threads, out)
pickle.dump(sizes, out)
out.close()

        
