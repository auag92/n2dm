# Script which tests the EMTv102.EMT() calculators get_potential_energy() and get_forces() function.
# The test is done by creating an atoms object consisting of two Cu atoms placed at a distance r from each other.
# The get_potential_energy() and get_forces() functions are called and the results are stored, then r is incremented 
# and the process is repeated. The Energy and Force of the system is then plottet as a function of r up to just 
# beyond the cutoff distance of approx 4.8. The numerical deriviative of the energy is then calculated and compared to the # force.

# Written by: Rasmus E. Christiansen

# Initialization
import numpy
from ase import EMT
import asap3
import matplotlib.pyplot as ppl
from ase import Atoms 

# The cutoff distance for cobber is defined, r_cut_Cu = 4.7695, rounded to 0 decimals
r_cut_Cu = 5
# The distance r is defined
rstart = r = 1.5
# The incrementation is defined
inc = 0.01
# The number of steps is defined
num_steps = int(round(r_cut_Cu/inc))

# The atoms object is created
atomsNew = Atoms(positions=[[0,0,0], [0,0,r]], symbols='CuCu', cell=(10,10,10))
atomsOld = Atoms(positions=[[0,0,0], [0,0,r]], symbols='CuCu', cell=(10,10,10))


# The calculator is initialized
BOHR = 0.5292

emt = asap3.EMT2013({29: {'eta2': 1.652 * 1/BOHR,'lambda': 1.906 * 1/BOHR,'kappa': 2.740 * 1/BOHR,'e0': -3.51,'V0': 2.476,'seq': 2.67*BOHR,'neq': 0.0091 * 1/(BOHR * BOHR * BOHR),'mass': 63.546}})

# The Energy is calculated using the Effective Medium Theory by attaching the EMT calculator to the atoms object.
atomsNew.set_calculator(emt)
atomsOld.set_calculator(EMT())


# Calculations

# The calculations for the forces and energy are performed
Enew = numpy.zeros([num_steps])
Eold = numpy.zeros([num_steps])
Fnew = numpy.zeros([num_steps,2,3])
Fold = numpy.zeros([num_steps,2,3])
for i in range(num_steps):
    # Energy and force are calculated and stored
    Enew[i] = atomsNew.get_potential_energy()
    Fnew[i,:] = atomsNew.get_forces()
    Eold[i] = atomsOld.get_potential_energy()
    Fold[i,:] = atomsOld.get_forces()
    # The distance r is incremented
    r += inc
    # The atoms object is updated
    atomsNew.set_positions([[0,0,0],[0,0,r]])
    atomsOld.set_positions([[0,0,0],[0,0,r]])
    


# Sets up the x-axis for the plots
x = rstart+inc*numpy.array(range(num_steps))

# Calculates the minimum for the energy
EminNew = min(Enew)
EminOld = min(Eold)
# Calculates the distance for the minimum energy
rminNew = rstart+inc*Enew.argmin()
rminOld = rstart+inc*Eold.argmin()


# Calculates the numerical deriviative of the energy in the z-direction for atom 1:
FnumNew = numpy.zeros([num_steps,2])
FnumOld = numpy.zeros([num_steps,2])
for i in range(num_steps-1):
    FnumNew[i] = [(Enew[i+1] - Enew[i]) / inc , -(Enew[i+1] - Enew[i]) / inc]
    FnumOld[i] = [(Eold[i+1] - Eold[i]) / inc , -(Eold[i+1] - Eold[i]) / inc]


# Plots
ppl.figure(1)
ppl.plot(x,Enew)
ppl.plot(x,Eold)
ppl.figure(2)
ppl.plot(x,Fnew[:,0,2])
ppl.plot(x,FnumNew[:,0])
ppl.plot(x,Fold[:,0,2])
ppl.plot(x,FnumOld[:,0])
ppl.figure(3)
ppl.plot(x,Fnew[:,1,2])
ppl.plot(x,FnumNew[:,1])
ppl.plot(x,Fold[:,1,2])
ppl.plot(x,FnumOld[:,1])
ppl.show()







