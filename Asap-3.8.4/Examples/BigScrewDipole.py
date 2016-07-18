"""BigScrewDipole.py   -   MD simulation of two screw dipoles annihilating.

This example shows the annihilation of two screw dipoles.  It can use the
old or the new EMT potential, it can plot or not, and it can run in serial or
in parallel.

Usage: [mpi]python BigScrewDipole.py {newpot|oldpot} {plot|noplot}
                                     {serial|parallel}

If the third parameter is parallel, two or four nodes must be reserved.
"""

from Simulations.Asap import *
from Numeric import *
from Structures.IonDynamics import *
from Simulations.Asap.Upset import *
from Setup.Dislocation.Screw import Screw
from Structures.ChemicalElements import Copper
import Structures.ChemicalElements.AtomicWeight
from Structures.AtomsInArrays import ListOfCartesianAtoms
import sys
from sys import argv
import time, os

sys.stdout = sys.stderr


#newerror = os.open("stderr."+str(Scientific.MPI.world.rank), os.O_WRONLY|os.O_SYNC|os.O_CREAT|os.O_TRUNC, 0660)
#os.dup2(newerror, sys.stderr.fileno())

if len(argv) != 4:
    print __doc__
    raise RuntimeError, "There must be three arguments."

if argv[1] == "newpot":
    newpot = 1
elif argv[1] == "oldpot":
    newpot = 0
else:
    print __doc__
    raise RuntimeError, "First argument must be newpot or oldpot."

if argv[2] == "plot":
    doplot = 1
elif argv[2] == "noplot":
    doplot = 0
else:
    print __doc__
    raise RuntimeError, "Second argument must be plot or noplot."

if argv[3] == "serial":
    parallel = 0
    ismaster = 1
    filename = "BigScrewDipole.serial.energy"
elif argv[3] == "parallel":
    from Scientific.MPI import world
    parallel = world.size
    ismaster = (world.rank == 0)
    mpirank = world.rank
    filename = "BigScrewDipole.parallel%d.energy" % (parallel,)
    if parallel != 2 and parallel != 4:
        raise RuntimeError, "Parallel runs must use 2 or 4 nodes, not "+str(world.size)
else:
    print __doc__
    raise RuntimeError, "Third argument must be serial or parallel."

if ismaster:
    print __doc__
    print ""
    print "New potential: ", newpot
    print "Plotting:      ", doplot
    print "Parallel nodes:", parallel
    print "Output file:   ", filename

if doplot and parallel:
    raise RuntimeError, "Cannot plot a parallel simulation."

def ekin(a):
    p = a.GetCartesianMomenta()
    return 0.5 * sum((p*p).flat) / Copper.GetProperty("AtomicWeight")

if ismaster:
    def printenergies(a, outfile = None, offset=0, n = [0,0.0]):
        try:
            l = a.GetTotalNumberOfAtoms()
        except AttributeError:
            l = len(a)
        ep = a.GetEnergy()/l
        if offset:
            n[1] = ep
        ep = ep - n[1]
        ek = a.GetKineticEnergy()/l
        print "%4d Epot =%10.6f    Ekin =%10.6f    Etot =%10.6f" % (n[0], ep,
                                                                    ek, ep+ek)
        if outfile is not None:
            outfile.write("%4d %10.6f %10.6f %10.6f\n" % (n[0], ep, ek, ep+ek))
        n[0] = n[0] + 1
else:
    def printenergies(a, outfile = None, offset=0, n = [0,0.0]):
        a.GetEnergy()  # Make sure the nodes stay in sync.
        a.GetKineticEnergy()

#verb = Verbose(3)
verb = Verbose(ismaster)
#verb = Verbose(0)


if newpot:
    a = 6.789483
else:
    a = EMTPotential(29).GetLatticeConstant() # lattice constant
b = a / sqrt(2)                         # length of Burgers vector
basis = array([( b / 2, a / 2, b / 2),  # fcc basis vectors
               (-b / 2, a / 2, b / 2),
               (0, 0, b)])
box = array((46 * b, 35 * a, 5 * b))    # box to be filled woth atoms
center = 0.5 * box

region = Rectangle(box, center)
if ismaster:
    positions = Fill(region, basis)         # a NumPy array with positions
else:
    positions = zeros((0,3), Float)

initialatoms = ListOfCartesianAtoms(positions=positions,
                                    momenta=zeros(positions.shape, Float))
initialatoms.SetAtomicNumbers(Copper.GetProperty("AtomicNumber") *
                              ones((len(initialatoms),)))

superCell = SuperCell(box, [1, 1, 1])   # a supercell with periodic boundary
                                        # conditions in all directions
if newpot:
    potential = EMT(EMTAtomicUnitParameters())
else:
    potential = EMTPotential(29)

if parallel:
    if parallel == 2:
        cells = [2,1,1]
    elif parallel == 4:
        cells = [2,2,1]
    print "Node %d: Preparing to create parallel list of atoms (%d atoms)." % (mpirank, len(initialatoms))
    atoms = PAtoms(initialatoms, superCell, nCells=cells)
    print "Node %d: Distribution done, now %d atoms."  % (mpirank, len(atoms))
else:
    atoms = Atoms(initialatoms, superCell)

print "Setting potential."
atoms.SetPotential(potential)
print "Calculating energy."
e = atoms.GetEnergy()
print "The system contains",len(atoms),"atoms."
print "Perfect FCC - the energy should be zero:", e
printenergies(atoms, offset=1)

# Set up a screw dislocation dipole (or a 2D array of dipoles):
line = (0, 0, 1)                        # this is a [110] direction
screw1 = Screw(center + ( 3 * b, a / 4, 0), line, b) # first screw
screw2 = Screw(center + (-3 * b, a / 4, 0), line, b) # second screw
(screw1 - screw2).ApplyTo(atoms)

# Writer("ScrewDipole.nc", atoms).Write()
CNA(atoms)
if doplot:
    plot = atoms.GetPlot()
    plot("zoom 120")

majors, minors = 30, 50
print "Dipole introduced - %d x %d timesteps:" % (majors, minors)

dynamics=VelocityVerlet(atoms, lambda:None, timestep=0.2)
if ismaster:
    outfile = open(filename, "w")
else:
    outfile = None
    
printenergies(atoms, outfile)

for j in range(1,majors+1):
    dynamics(iterations=minors)
    #print "%2d %9.4f" % (j, atoms.GetEnergy())
    printenergies(atoms, outfile)
#    print "Node %d: #atoms = %d, max(Z) = %d, min(Z) = %d, max(m) = %d, min(m) = %d" % (Scientific.MPI.world.rank, len(atoms), max(atoms.GetAtomicNumbers()), min(atoms.GetAtomicNumbers()), max(atoms.GetMasses()), min(atoms.GetMasses()))

    if doplot and (j % 5 == 0):
        print "Plotting step", j
        CNA(atoms)
        plot.Update()
        plot.ColorByClasses(atoms.GetColors())

forces = atoms.GetCartesianForces()
print "Average force:", sqrt(dot(forces.flat, forces.flat) / len(forces))

if doplot:
    print "Exiting in 10 seconds."
    time.sleep(10)

#stresses = atoms.GetStresses()
#p = -(stresses[:,0] + stresses[:,3] + stresses[:,5])
#Color(atoms, p)
#plot2 = atoms.GetPlot()
