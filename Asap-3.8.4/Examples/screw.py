from Asap import *
from Numeric import *
from Upset import *
from Visualization.Avatars.RasMol import *
from Setup.Dislocation.Screw import Screw

pot = EMTPotential()
a = pot.GetLatticeConstant()
b = a / sqrt(2)
h = b * sqrt(3) / 2
h3 = a / sqrt(3)
D = array((0.0, 0.0, 0.0))
A = array((h / 3, h3, b / 2))
B = array((h, 0.0, b / 2))
C = array((0.0, 0.0, b))
alpha = (B + C + D) / 3.0
beta = (A + C + D) / 3.0
gamma = (A + B + D) / 3.0
delta = (A + B + C) / 3.0

basis = array([B, A, C])

Lz = 5 * b
L2x = 48 * h
L2y = 48 * h3
rect = Rectangle([L2x, L2y, Lz])
pos = Fill(rect, basis)
superCell = SuperCell([L2x, L2y, Lz], [0, 0, 1])
atoms = Atoms(pos, superCell)
screw = Screw(delta, C - D, b)
screw.ApplyTo(atoms)
atoms.SetPotential(pot)

atoms.CNA()
p = RasMol(atoms)
p.ColorByClasses({0 : "red", 1 : "white", 2 : "blue"})

from Structures.IonDynamics import *

dynamics=QuickminAllCoordinates(atoms, lambda:None, timestep=0.2)
majors, minors = 100, 20

for j in range(majors):
    dynamics(iterations=minors)
    print atoms.GetEnergy()

