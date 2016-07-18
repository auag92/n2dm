from Asap import *
from Asap.Dynamics.Langevin import *
from Asap.Filters.TimeAveragedPositions import TimeAveragedPositions
from Asap.Filters.RestrictedCNA import RestrictedCNA
from Asap.Trajectories.NetCDFTrajectory import *
from ASE.Visualization.PrimiPlotter import *
from Numeric import *
import Scientific.MPI as MPI
import commands
import sys, os, time

infile = "initial.nc"
infofile = "relax-01.dat"
avgfile = "relax-avg-01.nc"
outfile = "relax-01.nc"
debugfile = "/scratch/schiotz/relax-n%d.err"
n_avg = 100
n_output = 5000
n_summary = 100
n_total = 10 * n_output
timestep = 5.0   # fs
lgvfrict = 0.001
temperature = 300 # Kelvin
kB = 8.61734e-5

def now():
    return time.ctime(time.time())

nodes = (5, 5, 2)

DebugOutput(debugfile, nomaster=True)
PrintVersion(1)
ismaster = MPI.world.rank == 0

if ismaster:
    info = open(infofile, "w")
else:
    info = sys.stderr  # Which has been redirected to a debugging file.

info.write("Simulation started: %s\n" % now())
info.flush()
atoms = ParallelNetCDFTrajectory(infile).GetListOfAtoms(nodes)
nAtoms = atoms.GetTotalNumberOfAtoms()
atoms.SetCalculator(EMT(EMTRasmussenParameters()))
atoms.SetCartesianMomenta(zeros((len(atoms),3), Float))

dyn = Langevin(atoms, timestep * femtosecond,
               temperature * 15 * kB, 3*lgvfrict)

i = 0
while atoms.GetKineticEnergy() < 1.5 * temperature * kB * nAtoms:
    i = i + 1
    dyn.Run(10)
    info.write("Energy per atom (step %d): Epot = %.3f eV   Ekin = %.3f eV\n" %
               (i, atoms.GetPotentialEnergy()/nAtoms, atoms.GetKineticEnergy()/nAtoms))
    info.flush()

info.write("Temperature has reached %d K: %s\n" % (temperature, now()))
info.flush()
dyn.SetTemperature(temperature * kB)
dyn.SetFriction(lgvfrict)

avg = TimeAveragedPositions(atoms, n_avg, n_output, modifyAtoms=1, verbose=1)
dyn.Attach(avg)

cna = RestrictedCNA(atoms, verbose=1)  # Runs CNA on the modified atoms
avg.Attach(cna)

writer = ParallelNetCDFTrajectory(outfile, atoms, interval=n_output, mode="w")
dyn.Attach(writer)  

avgwriter = ParallelNetCDFTrajectory(avgfile, atoms, interval=1, mode="w")
# By attaching to cna, the atoms will hold their averaged positions when called.
cna.Attach(avgwriter)

print atoms.nCells
print avg.nCells
print cna.nCells

def invis(a):
    c = a.GetTags()
    z = a.GetCartesianPositions()[:,2]
    basis = a.GetUnitCell()
    return not_equal(c, 1)

p = ParallelPrimiPlotter(atoms)
p.SetInvisibilityFunction(invis)
p.SetColors({0:"red", 1:"yellow", 2:"blue"})
#p.SetOutput(X11Window())
p.SetOutput(GifFile("plot/test"))
cna.Attach(p)

# Cause a first output (not averaged)
avg.ForceNotify()
# Also of the atoms
writer.Update()

for i in xrange(n_total / n_summary):
    dyn.Run(n_summary)
    info.write("Energy per atom (step %d): Epot = %.3f eV   Ekin = %.3f eV\n" %
               (i, atoms.GetPotentialEnergy()/nAtoms,
                atoms.GetKineticEnergy()/nAtoms))
    info.write("    Stress: %g %g %g\n            %g %g %g\n" %
               tuple(atoms.GetStress()))
    info.flush()
    
print "avg.counter is", avg.counter
