from Asap import *

from Setup.Dislocation.Screw import *
from Setup.Lattice.FCC111Ortho import *
from Structures.ChemicalElements import Copper
import Structures.ChemicalElements.AtomicWeight
from Structures.IonDynamics import VelocityVerlet

from Numeric import *
from LinearAlgebra import *

# Make a chunk with (111) in the z direction.
latticeconstant = 6.58
lattice = FCC111Ortho((36,36,18), Copper, latticeconstant, symmetry=(1,1,0))

# Place the dislocations
dimensions = lattice.GetCoordinateBasis()
offset = lattice.MillerToDirection((-1,1,0)) * 1.2
center1 = matrixmultiply(dimensions, array((0.5001, 0.5001, 0.5001))) + offset
center2 = matrixmultiply(dimensions, array((0.5001, 0.5001, 0.5001))) - offset
line = lattice.MillerToDirection((1,1,0))
b = latticeconstant/sqrt(2)
#line = lattice.MillerToDirection((1,1,1))
#b = latticeconstant*sqrt(3)
print "Dislocation 1:  position =", center1, " direction =", line
print "Dislocation 2:  position =", center2, " direction =", line

disloc1 = Screw(center1, line, b)
disloc2 = Screw(center2, line, -b)

# Apply dislocation dislplacement field to the lattice, keeping the corners
# fixed to minimize trouble at the boundaries.
(disloc1 + disloc2).ApplyTo(lattice, fixedbox=1)

# Make a simulation object
a = Atoms(lattice, potential=EMTPotential())
#CNA(a)
Color(a, a.GetStresses()[:,4])
plot = a.GetPlot()
plot.SetRotation([-70, 20, -171])

# Set momenta to zero (this should be the default).
a.SetCartesianMomenta(zeros(a.GetCartesianPositions().shape, Float))
mover = VelocityVerlet(a, lambda:None, 0.5)
p = a.GetCartesianMomenta()
print "    Ekin =", 0.5 * sum((p*p).flat) / Copper.GetProperty("AtomicWeight"), "eV"

# raise RuntimeError

for i in range(10):
    for i in range(5):
        mover(5)
        p = a.GetCartesianMomenta()
        print "    Ekin =", 0.5 * sum((p*p).flat) / Copper.GetProperty("AtomicWeight"), "eV"        
    #CNA(a)
    Color(a, a.GetStresses()[:,4])
    print "Updating plot ..."
    plot.Update(a)
    print " done"

