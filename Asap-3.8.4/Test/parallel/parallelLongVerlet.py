#import production

### Job name
#PBS -N LongVerlet

### Mail to user
#PBS -m ae

### Queue name (small, medium, long, verylong)
#PBS -q verylong

### Number of nodes
#PBS -l nodes=1:xeon:ppn=8


from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from cPickle import *
from numpy import *
import sys, os, time
from asap3.testtools import *
from asap3.mpi import world

debug = 1
if debug == 1:
    DebugOutput("longverlet%d.log", nomaster=True)
elif debug == 2:
    time.sleep(world.rank)
    print "PID:", os.getpid()
    time.sleep(20)

writetraj = True
#set_verbose(1)


ismaster = world.rank == 0
isparallel = world.size != 1
if ismaster:
    print_version(1)
    
if world.size == 1:
    cpulayouts = [None]
elif world.size == 2:
    cpulayouts = [[2,1,1], [1,2,1]]
elif world.size == 3:
    cpulayouts = [[3,1,1], [1,1,3]]
elif world.size == 4:
    cpulayouts = [[4,1,1], [2,1,2]]
elif world.size == 8:
    cpulayouts = [[8,1,1], [4,2,1], [2,2,2]]
elif world.size == 16:
    cpulayouts = [[8,1,2], [4,2,2]]
else:
    raise RuntimeError("No cpu layout for %d CPUs" % (world.size,))

e_start_dict = {}
for cpulayout in cpulayouts:
    for boundary in ((1,1,1), (0,0,0), (0,1,1), (1,0,0)):
        print ("CPU Layout: %s.  Periodic boundary conditions: %s."
               % (str(cpulayout), str(boundary)))
        if ismaster:
            atoms = FaceCenteredCubic(size=(160,20,20), symbol="Cu",
                                      pbc=boundary, latticeconstant=3.61*1.04)
        else:
            atoms = None
        if isparallel:
            atoms = MakeParallelAtoms(atoms, cpulayout)
        natoms = atoms.get_number_of_atoms()
        atoms.set_calculator(EMT())
        MaxwellBoltzmannDistribution(atoms, 3000*units.kB)
        if ismaster:
            print "Initializing"
        e_st_pot = atoms.get_potential_energy()/natoms
        try:
            e_start_pot = e_start_dict[boundary]
        except KeyError:
            e_start_pot = e_start_dict[boundary] = e_st_pot
        else:
             ReportTest("Initial energy ", e_st_pot, e_start_pot, 1e-4,
                        silent=True)
        dyn = VelocityVerlet(atoms, logfile="-", dt=3*units.fs, loginterval=1)
        dyn.run(50)
        e_start = (atoms.get_potential_energy() 
                   + atoms.get_kinetic_energy())/natoms
        if ismaster:
            print "Running"
        dyn = VelocityVerlet(atoms, dt=5*units.fs)
        logger = MDLogger(dyn, atoms, '-', peratom=True)
        logger()
        dyn.attach(logger, interval=25)
        if writetraj:
            if cpulayout is None:
                tname = "longVerlet-serial-%d-%d-%d.traj" % tuple(boundary)
            else:
                tname = ("longVerlet-%d-%d-%d--%d-%d-%d.traj" %
                         (tuple(cpulayout) + tuple(boundary)))
            traj = PickleTrajectory(tname, "w", atoms)
            traj.write()
            dyn.attach(traj, interval=1000)
        
        for i in range(100):
            dyn.run(100)
            if i % 5 == 4:
                print i+1, "%"
            e = (atoms.get_potential_energy() + atoms.get_kinetic_energy())/natoms
            ReportTest("Step "+str(i), e, e_start, 1e-4, silent=True)

ReportTest.Summary()

