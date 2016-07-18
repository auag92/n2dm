#!/usr/bin/env python

from Asap.Setup.Lattice.FCCOrtho import *
from Asap.Setup.Dislocation import Dislocation
from ASE.ChemicalElements import Element
from Asap import *
from Asap.Trajectories import NetCDFTrajectory
from ASE.Visualization.PrimiPlotter import *
from Numeric import *

splitting = 5
size = (25, 35, 10)

# Create a slab of Gold
Gold = Element("Au")
atoms = FCCOrtho(((1,1,-2), (-1,1,0), (1,1,1)), size, Gold)
basis = atoms.GetUnitCell()
print "Number of atoms:", len(atoms)

# Create the dislocation
center = 0.5 * array([basis[0,0], basis[1,1], basis[2,2]]) + array([0.1, 0.1, 0.1])
offset = 0.5 * splitting * atoms.MillerToDirection((-1,0,1))

d1 = Dislocation(center - offset, atoms.MillerToDirection((-1,-1,0)),
                 atoms.MillerToDirection((-2,-1,1))/6.0)
d2 = Dislocation(center + offset, atoms.MillerToDirection((1,1,0)),
                 atoms.MillerToDirection((1,2,1))/6.0)

# Insert the dislocation in the slab
(d1+d2).ApplyTo(atoms)

# Write the configuration to a file
traj = NetCDFTrajectory("initial.nc", atoms)
traj.Update()
traj.Close()

# Analyse the local crystal structure prior to plotting.
print "Now running Common Neighbor Analysis"
atoms.SetCalculator(EMT())
CNA(atoms)

# Now make two plots of the system.  The surfaces at the back and
# front are removed, and all atoms in perfect FCC order are removed
# from the plot.
c = atoms.GetTags()
y = array(atoms.GetCartesianPositions()[:,1])
z = array(atoms.GetCartesianPositions()[:,2])
invis1 = equal(c, 0) + less(y, 5) + greater(y, basis[1,1] - 5)
invis2 = equal(c, 0) + less(z, 5) + greater(z, basis[2,2] - 5)

p1 = PrimiPlotter(atoms)
p1.SetRotation([-88,0,0])
p1.SetInvisible(invis1)
p1.SetOutput(GifFile("sideview"))

p2 = PrimiPlotter(atoms)
p2.SetInvisible(invis2)
p2.SetOutput(GifFile("topview"))

for p in (p1,p2):
    p.SetColors({0:"red", 1:"yellow", 2:"blue"})
    p.SetOutput(X11Window())
    p.Plot()



