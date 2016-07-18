from QC import *
from Numeric import *
from IPAF import IPAF
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

Lx = 8.1 * h
Ly = 8.1 * h3
Lz = 16 * b
rect1 = Rectangle([Lx, Ly, Lx])
pos = Fill(rect1, basis)
L2x = 48 * h
L2y = 48 * h3
rect2 = Rectangle([L2x, L2y, Lz])
qcpos = Fill(rect2 - rect1, 4 * basis)
superCell = SuperCell([L2x, L2y, Lz], [0, 0, 0])
atoms = QCAtoms(Atoms(pos, superCell), qcpos)
screw = Screw(delta, C - D, b)
#screw.ApplyTo(atoms)

tetrahedra = IPAF(atoms, 0.49 * Lz)

lattice = SuperCell(basis)
atoms.SetTetrahedra(tetrahedra, lattice)

atoms.SetPotential(pot)
stop
from Structures.IonDynamics import *

dynamics=QuickminAllCoordinates(atoms, lambda:None, timestep=0.2)
majors, minors = 10, 4

for j in range(majors):
    dynamics(iterations=minors)
    print atoms.GetEnergy()

