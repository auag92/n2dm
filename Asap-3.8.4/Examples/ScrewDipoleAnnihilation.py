from Asap import *
from Numeric import *

a = EMTPotential().GetLatticeConstant() # lattice constant
b = a / sqrt(2)                         # length of Burgers vector
basis = array([( b / 2, a / 2, b / 2),  # fcc basis vectors
               (-b / 2, a / 2, b / 2),
               (0, 0, b)])
box = array((17 * b, 13 * a, 5 * b))    # box to be filled woth atoms
center = 0.5 * box

from Upset import *
region = Rectangle(box, center)
positions = Fill(region, basis)         # a NumPy array with positions
superCell = SuperCell(box, [1, 1, 1])   # a supercell with periodic boundary
                                        # conditions in all directions
atoms = Atoms(positions, superCell, EMTPotential())
print "Perfect FCC - the energy should be zero:", atoms.GetEnergy()

# Set up a screw dislocation dipole (or a 2D array of dipoles):
from Setup.Dislocation.Screw import Screw
line = (0, 0, 1)                        # this is a [110] direction
screw1 = Screw(center + ( 3 * b, a / 4, 0), line, b) # first screw
screw2 = Screw(center + (-3 * b, a / 4, 0), line, b) # second screw
(screw1 - screw2).ApplyTo(atoms)

Writer("ScrewDipole.nc", atoms).Write()
CNA(atoms)
plot = atoms.GetPlot()

print "Dipole introduced - 20 x 20 relaxation steps:"

from Structures.IonDynamics import *

dynamics=QuickminAllCoordinates(atoms, lambda:None, timestep=0.2)
majors, minors = 20, 20

for j in range(majors):
    dynamics(iterations=minors)
    CNA(atoms)
    plot.Update()
    plot.ColorByClasses(atoms.GetColors())
    print "%2d %9.4f" % (j, atoms.GetEnergy())

forces = atoms.GetCartesianForces()
print "Average force:", sqrt(dot(forces.flat, forces.flat) / len(forces))

stresses = atoms.GetStresses()
p = -(stresses[:,0] + stresses[:,3] + stresses[:,5])
Color(atoms, p)
plot2 = atoms.GetPlot()
