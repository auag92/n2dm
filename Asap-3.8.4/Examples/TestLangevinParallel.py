from Numeric import *
from Simulations.Asap import *
import Scientific.MPI
from Setup.Lattice.FCC111Ortho import *
from Structures.ChemicalElements import Copper
import Structures.ChemicalElements.AtomicWeight
from Simulations.Asap.IonDynamics import Langevin
import sys, os

PrintVersion(1)

sys.stdout = sys.stderr
#newerror = os.open("stderr."+str(Scientific.MPI.world.rank), os.O_WRONLY|os.O_SYNC|os.O_CREAT|os.O_TRUNC, 0660)
#os.dup2(newerror, sys.stderr.fileno())
#v=Verbose(3)

v=Verbose(0)

initial = FCC111Ortho((21,20,24), Copper, 3.59164)
totalatoms = len(initial)
basis = initial.GetCoordinateBasis()
supercell = SuperCell(basis, [1,1,1])

if Scientific.MPI.world.rank == 0:
    print "Master:", len(initial), "atoms."
    atoms = MakeParallelAtoms(initial, [2,1,1], supercell)
else:
    atoms = MakeParallelAtoms(None, [2,1,1], supercell)
print Scientific.MPI.world.rank, ":", len(atoms), "atoms."
atoms.SetPotential(EMT())
print Scientific.MPI.world.rank, ":", len(atoms), "atoms (2)."
atoms.SetCartesianMomenta(zeros((len(atoms),3), Float))
print Scientific.MPI.world.rank, ": Momenta are now zero."
atoms.CalculateEnergiesAndForces()
dynamics = Langevin(atoms, atoms.CalculateEnergiesAndForces, 0.3, 0.1, 0.01)
print Scientific.MPI.world.rank, ": Langevin dynamics created."

if Scientific.MPI.world.rank == 0:
    print 0, atoms.GetKineticEnergy()/totalatoms
else:
    dummy = atoms.GetKineticEnergy()/totalatoms

#for i in range(500):
for i in range(200):
    if i == 100:
        print "RESETTING TEMPERATURE"
        dynamics.SetTemperature(0.1)
    dynamics(10)
    if Scientific.MPI.world.rank == 0:
        print i+1, atoms.GetKineticEnergy()/totalatoms
        sys.stdout.flush()
    else:
        dummy = atoms.GetKineticEnergy()/totalatoms
print "Done.!"
del dynamics
print "Dynamics deleted"
print "References to atoms:", sys.getrefcount(atoms) - 1
print atoms.__del__
del atoms
print "Atoms deleted"
