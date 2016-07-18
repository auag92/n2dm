#!/usr/bin/env python

from Asap import *
from Asap.Setup.Lattice.Cubic import *
from Asap.Timing import ReportTiming
import sys, time

size = (29,29,30)
timesteps = 50
#size = (63,63,63)
#timesteps = 10

print "Timing script, only for calculations."

def runit():
    atoms = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]],
                              size=size, element="Cu")
    natoms = len(atoms)
    print "Number of atoms:", natoms, size
    print "Number of iterations:", timesteps

    atoms.SetCalculator(EMT())
    atoms.GetCartesianForces()
    r = atoms.GetCartesianPositions()

    print "Starting timing."
    startcpu, startwall = time.clock(), time.time()
    for i in range(timesteps):
        atoms.SetCartesianPositions(r)
        atoms.GetCartesianForces()
    cpu, wall = time.clock() - startcpu, time.time() - startwall
    print "Real time:", cpu, "(cpu),   ", wall, "(wall)"
    cpu *= 1e6 / (timesteps * natoms)
    wall *= 1e6 / (timesteps * natoms)

    print
    print
    print "CPU time:", cpu, "usec/atom/timestep"
    print "Wall time:", wall, "usec/atom/timestep"
    print "Ratio:", 100.0 *cpu/wall, "%"
    print
    return wall

wallone = runit()
nthreads = AsapThreads()
wallmany = runit()

print "Threading speedup:", 100.0 * wallone / wallmany, "%"
print "Efficiency:", 100.0 * wallone / wallmany / nthreads, "%"

ReportTiming()
