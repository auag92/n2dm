# First a documentation string to describe the example
"""An example script demonstrating an NVT MD simulation with Asap.

This example runs a molecular dynamics simulation on a small copper
cluster in the NVT ensemble, i.e. with constant particle number,
constant volume and constant temperature.  It uses Langevin dynamics
to get the constant temperature.
"""

# Import ASAP and relatives
from Asap import *
from Asap.Dynamics.Langevin import Langevin
from Asap.Trajectories import NetCDFTrajectory

# Import other essential modules
from Numeric import *

# Set up an fcc crystal of copper.
element = 29
lattice_constant = 3.59 # Angstrom
size = (10, 10, 10)
fcc_basis = array([[0.0, 0.5, 0.5],
                   [0.5, 0.0, 0.5],
                   [0.5, 0.5, 0.0]])
# Now make an array which looks like
# [ [0, 0, 0]
#   [0, 0, 1]
#      ...
#   [0, 0, 9]
#   [0, 1, 0]
#      ...
#   [9, 9, 9] ]
#
# When it is multiplied with fcc_basis it gives a small crystal with 1000
# atoms.
#
# The array is made with a simple loop, although fancy Numeric array
# operations could do it faster (but it would be close to
# unreadable).
pos = []
for i in range(size[0]):
    for j in range(size[1]):
        for k in range(size[2]):
            pos.append([i,j,k])
positions = matrixmultiply(pos, fcc_basis) * lattice_constant

# Now the atoms are created from the positions and a sufficiently
# large unit cell.  The unit cell is uncritical as there are no
# periodic boundaries, and the atoms are able to move outside the unit
# cell, so in principle it should work with the default
# 1-Angstrom-cube unit cell.

atoms = ListOfAtoms(positions=positions, cell=[max(positions[:,0]),
                                               max(positions[:,1]),
                                               max(positions[:,2])])
atoms.SetAtomicNumbers(element)

# Now the standard EMT Calculator is attached
atoms.SetCalculator(EMT())

# Make an object doing Langevin dynamics at a temperature of 800 K
dyn = Langevin(atoms, dt=5*femtosecond, temperature=800*kB, friction=0.005)

# Make a trajectory
traj = NetCDFTrajectory('MD_Cluster.nc', atoms, interval=500)
dyn.Attach(traj)

# Write the initial configuration
traj.Update()

# Print the energies
def printenergy(a, step=[0,]):
    n = len(a)
    ekin = a.GetKineticEnergy() / n
    epot = a.GetPotentialEnergy() / n
    print ("%4d: E_kin = %-9.5f  E_pot = %-9.5f  E_tot = %-9.5f  T = %.1f K" %
           (step[0], ekin, epot, ekin+epot, 2.0/3.0*ekin/kB))
    step[0] += 1
          
printenergy(atoms)

# Now do the dynamics, doing 5000 timesteps, writing energies every 50 steps
for major in range(100):
    dyn.Run(50)
    printenergy(atoms)

# Now increase the temperature to 2000 K and continue
dyn.SetTemperature(2000 * kB)
for major in range(100):
    dyn.Run(50)
    printenergy(atoms)
