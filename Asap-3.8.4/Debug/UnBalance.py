#PBS -N UnBalance
#PBS -m ae
#PBS -q long
#PBS -l nodes=1:opteron:ppn=2

"""Test handling of extreme load-unbalancing."""

from asap3 import *
from asap3.md import MDLogger
from ase.lattice.cubic import FaceCenteredCubic
import numpy as np
from asap3.mpi import world

#DebugOutput("UnBalance.%d.out")
#set_verbose(1)
print_version(1)

fast = False

#AsapThreads()

cpulayout = (1,1,2)

element = 'Pt'
size = (20,20,100)

master = world.rank == 0

if master:
    atoms = FaceCenteredCubic(symbol=element, size=size, pbc=(True, True, False))
    atoms.center(vacuum=10.0, axis=2)
    atoms.set_momenta(np.zeros((len(atoms),3)))
    # Select an atom to get a kick
    r = atoms.get_positions()
    uc = atoms.get_cell()
    x = r[:,0] - 0.5 * uc[0,0]
    y = r[:,1] - 0.5 * uc[1,1]
    z = r[:,2]
    zprime = z - 0.01 * (x * x + y * y)
    n = np.argmax(zprime)
    #a = atoms[n]
    #dp = np.sqrt(2 * a.mass * 1000.0)
    #a.momentum = np.array([0, 0, dp])
    t = np.zeros(len(atoms), int)
    t[n] = 1
    atoms.set_tags(t)
else:
    atoms = None
atoms = MakeParallelAtoms(atoms, cpulayout)
print len(atoms), atoms.get_number_of_atoms()
atoms.set_calculator(EMT())

traj = PickleTrajectory("UnBalance.traj", "w", atoms)
if fast:
    atoms.get_forces()
    traj.write()
    for i in range(50):
        print "\n\n\n\n*** STEP %i ***\n\n\n\n\n" % (i,)
        r = atoms.get_positions()
        r += atoms.get_tags().reshape((-1,1)) * np.array([[0, 0, 20.0],])
        atoms.set_positions(r)
        atoms.get_forces()
        traj.write()
else:
    dyn = VelocityVerlet(atoms, 5*units.fs)
    logger = MDLogger(dyn, atoms, 'UnBalance.log', stress=True, peratom=True)
    dyn.attach(logger, interval=10)
    dyn.attach(traj, interval=100)
    dyn.run(10000)


    