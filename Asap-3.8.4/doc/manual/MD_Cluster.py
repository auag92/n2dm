# Import the necessary modules

from Numeric import *
from Asap import ListOfAtoms, EMT, kB, femtosecond
from Asap.Dynamics.Langevin import Langevin
from Asap.Trajectories.NetCDFTrajectory import NetCDFTrajectory

# Make the positions
Z = 29            # Copper
latconst = 3.61   # Lattice constant of copper
N = 5
pos = []

for i in range(N):
    for j in range(N):
        for k in range(N):
            pos.append((i,j,k))
            pos.append((i,j+0.5,k+0.5))
            pos.append((i+0.5,j,k+0.5))
            pos.append((i+0.5,j+0.5,k))
pos = array(pos) * latconst
supercell = [N*latconst, N*latconst, N*latconst]

# Create the atoms object
atoms = ListOfAtoms(positions=pos, cell=supercell)
atoms.SetAtomicNumbers(Z)

# The dynamics
atoms.SetCalculator(EMT())  # Use the EMT potential
dyn = Langevin(atoms, 5*femtosecond, 800*kB, 0.002)

# Output to NetCDF file
trajectory = NetCDFTrajectory("output.nc", atoms, interval=2000)
dyn.Attach(trajectory)

# Do the dynamics
for i in range(100):
    if i == 50:
        dyn.SetTemperature(2000*kB)
        print "Setting temperature to 2000K"
    dyn.Run(100)
    print "Ekin = %.5f   Epot = %.5f  Etot = %.5f" % \
          (atoms.GetKineticEnergy()/len(atoms),
           atoms.GetPotentialEnergy()/len(atoms),
           (atoms.GetKineticEnergy()+atoms.GetPotentialEnergy())/len(atoms))

           
