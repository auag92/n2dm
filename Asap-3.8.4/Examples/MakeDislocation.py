#!/usr/bin/env python

from numpy import *
from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.Setup.Dislocation import Dislocation
from ase.visualize.primiplotter import *

print_version(1)

splitting = 5
size = (50, 88, 35)
#size = (30, 25, 7)

Gold = "Au"
slab = FaceCenteredCubic(directions=((1,1,-2), (-1,1,0), (1,1,1)),
                         size=size, symbol=Gold)
basis = slab.get_cell()
print basis
print "Number of atoms:", len(slab)

center = 0.5 * array([basis[0,0], basis[1,1], basis[2,2]]) + array([0.1, 0.1, 0.1])
offset = 0.5 * splitting * slab.miller_to_direction((-1,0,1))
print center

d1 = Dislocation(center - offset, slab.miller_to_direction((-1,-1,0)),
                 slab.miller_to_direction((-2,-1,1))/6.0)
d2 = Dislocation(center + offset, slab.miller_to_direction((1,1,0)),
                 slab.miller_to_direction((1,2,1))/6.0)

atoms = Atoms(slab)
(d1+d2).apply_to(atoms)
del slab

# Write the file in case CNA fails
traj = PickleTrajectory("initial.traj", "w", atoms)
traj.close()

print "Now attempting CNA"

atoms.set_calculator(EMT(EMTRasmussenParameters()))
c = CNA(atoms)
atoms.set_tags(c)


for i in range(3):
    print i, sum(equal(c, i))
    
y = array(atoms.get_positions()[:,1])
z = array(atoms.get_positions()[:,2])
invis1 = equal(c, 0) + less(y, 5) + greater(y, basis[1,1] - 5)
invis2 = equal(c, 0) + less(z, 5) + greater(z, basis[2,2] - 5)

# Overwrite the file, this time including CNA data
traj = PickleTrajectory("initial.traj", "w", atoms)
traj.close()

p1 = PrimiPlotter(atoms)
p1.set_rotation([-88,0,0])
p1.set_invisible(invis1)

p2 = PrimiPlotter(atoms)
p2.set_invisible(invis2)

for p in (p1,p2):
    p.set_colors({0:"red", 1:"yellow", 2:"blue"})
    p.set_output(X11Window())
    p.set_output(GifFile("init"))
    p.plot()



