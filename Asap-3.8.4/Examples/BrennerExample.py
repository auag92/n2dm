# Copyright (C) 2007  CAMP
# Please see the accompanying LICENSE file for further information.

from ASE import Atom
from Asap import *
from Asap.Dynamics.VelocityVerlet import VelocityVerlet
from Asap.Utilities import BrennerImport

# read the coordinates and the types of atoms from a
# brenner dot d file:
brenner = BrennerImport.BrennerImport()
brenner.initCoord("coord.d")
a=brenner.GetListOfAtoms()

# initialize the Brenner Potential with automatically adapting cube size
potential = BrennerPotential(True)
# set accuracy
potential.SetRll(2.5)
#define the cube size (You can only set the cube size, if you created the potential using BrennerPotential(False))
#potential.SetCubeSize(brenner.GetCubeSize())

a.SetCalculator(potential)
#step size 0.2 femtoseconds
dyn = VelocityVerlet(a, 0.2*femtosecond)

step = 0;
etot = a.GetPotentialEnergy() + a.GetKineticEnergy()

for i in range(6):
    dyn.Run(100)
    step+=100
    print "Total energy after %4.d steps: %.8f"%(step,a.GetPotentialEnergy()+a.GetKineticEnergy())
    #write xmol.d log file
    potential.WriteXmol(step, 0.2*femtosecond);
    

