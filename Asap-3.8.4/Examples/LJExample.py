# Copyright (C) 2007  CAMP
# Please see the accompanying LICENSE file for further information.

# Illustrates how to use the Lennard Jones Potential

import ASE.ChemicalElements.mass
from ASE import Atom
from ASE.ChemicalElements import Element

from Asap import *
from Asap import ListOfAtoms
from Asap.Setup.Lattice.Cubic import *
from Asap.Dynamics.VelocityVerlet import VelocityVerlet
from Asap.Setup.Lattice.Cubic import *

from Numeric import *

PrintVersion(1)

nsteps = 100 # the number of steps to be performed
stepsize=1*femtosecond # the size of one step
size=10      # the cube size: here 10*10*10
rCut = -1    # cut-off radius (optional parameter, default: -1)
v0 = True    # add v0 during calculation (optional parameter, default: True)

elements = [A, B, C]                #<-- add here the atomic numbers of your elements 
epsilon  = [eps_AA,    0,      0,   #<-- add here the epsilons 
            eps_AB, eps_BB,    0,   #<-- add here the epsilons 
            eps_AC, eps_BC, eps_CC] #<-- add here the epsilons 
sigma    = [sig_AA,    0,      0,   #<-- add the sigmas  
            sig_AB, sig_BB,    0,   #<-- add the sigmas  
            sig_AC, sig_BC, sig_CC] #<-- add the sigmas  

masses   = [ASE.ChemicalElements.mass.masses[A],  #<-- add here the atomic numbers of your elements again
            ASE.ChemicalElements.mass.masses[B],
            ASE.ChemicalElements.mass.masses[C]]

atoms = L1_2([[1,0,0],[0,1,0],[0,0,1]], size=(size,size,size),
                          element=("A", "B", "C"),  #<-- add here the names of the elements
                          periodic=(1,0,1), 
                          debug=0,
                          latticeconstant=3.95)     #<-- adapt the lattice constant

atoms.SetCalculator(LJPotential(3, elements, epsilon, sigma, masses, rCut, v0))

epot = atoms.GetPotentialEnergy()

print "Potential energy: ", epot
