#!/usr/bin/env python

from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.Setup.Dislocation import Dislocation
from Asap.Filters.TimeAveragedPositions import TimeAveragedPositions
from Asap.Filters.RestrictedCNA import RestrictedCNA
from Asap.Dynamics.Langevin import Langevin
from Asap.Trajectories import NetCDFTrajectory
from ASE.Visualization.PrimiPlotter import *
from ASE.Visualization.FieldPlotter import *
from Numeric import *

PrintVersion(1)

cleanup = True  # Delete files at end.

#size = (50, 88, 35)
size = (30, 25, 7)
n_avg=25
n_output=50
n_total = 2*n_output

atoms = FaceCenteredCubic(directions=((1,1,-2), (-1,1,0), (1,1,1)),
                         size=size, element="Au", debug=0, periodic=(0,0,1))
basis = atoms.GetUnitCell()

center = 0.5 * array([basis[0,0], basis[1,1], basis[2,2]]) + array([0.1, 0.1, 0.1])

disl = Dislocation(center, atoms.MillerToDirection((1,1,0)),
                   atoms.MillerToDirection((1,1,0))/2.0)

#atoms = ListOfAtoms(slab)
disl.ApplyTo(atoms)

atoms.SetCalculator(EMT(EMTRasmussenParameters()))

# Constant temperatur dynamics
dyn = Langevin(atoms, 5 * femtosecond, 300 * kB, 0.002)

# Run Common Neighbor Analysis on each n_output timestep, using the
# positions averaged over the last n_avg timesteps, and writing the
# result to the output file.
avg = TimeAveragedPositions(atoms, n_avg, n_output)
dyn.Attach(avg)
cna = RestrictedCNA(avg, verbose=1)
avg.Attach(cna)

# We output and plot the instantaneous positions.  Replace atoms with
# avg to output and plot the averaged positions instead.
plotlog = open("plot.log", "w")
plotter = PrimiPlotter(atoms, verbose=1)
plotter.SetOutput(PostScriptFile("cna_plot"))
plotter.SetOutput(GifFile("cna_plot"))
plotter.SetLog(plotlog)

# Only plot atoms in the dislocation.

def invis(a):
    c = a.GetTags()
    return equal(c, 0)

plotter.SetInvisibilityFunction(invis)
plotter.SetColors({0:"red", 1:"yellow", 2:"blue"})
cna.Attach(plotter)

# We output and plot the instantaneous positions.  Replace atoms with
# avg to output and plot the averaged positions instead.
plotter = FieldPlotter(atoms, atoms.GetTags, verbose=1)
plotter.SetDataRange('plot')
plotter.SetOutput(PostScriptFile("fieldplot"))
plotter.SetOutput(GifFile("fieldplot"))
plotter.SetLog(plotlog)

cna.Attach(plotter)

# We output and plot the instantaneous positions.  Replace atoms with
# avg to output and plot the averaged positions instead.
plotter = FieldPlotter(atoms, atoms.GetTags, verbose=1)
plotter.SetBlackWhiteColors(reverse=True)
plotter.SetDataRange('plot')
plotter.SetOutput(PostScriptFile("bwplot"))
plotter.SetOutput(GifFile("bwplot"))
plotter.SetLog(plotlog)

cna.Attach(plotter)

# Run the simulation
dyn.Run(n_total)

plotlog.close()

if cleanup:
    print "Removing output files"
    print "Deleting plot.log"
    os.unlink("plot.log")
    for name in ("cna_plot", "fieldplot", "bwplot"):
        for num in ("0000", "0001"):
            for ending in (".ps", ".gif"):
                fn = name+num+ending
                print "Deleting", fn
                os.unlink(fn)

print
print "Test passed!"
