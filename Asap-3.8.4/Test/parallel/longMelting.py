#PBS -N longMelting
#PBS -q long
#PBS -l nodes=1:ppn=4:opteron4
#PBS -m ae

from asap3 import *
from asap3.md.langevin import Langevin
from asap3.md.nptberendsen import Inhomogeneous_NPTBerendsen
from asap3.md import MDLogger
from asap3.testtools import ReportTest
from ase.lattice.cubic import FaceCenteredCubic
import numpy as np
from asap3.mpi import world
import random
import time
import shutil
import os

quickie = False  # Runs quickly, but fails.
cleanup = True  # Remove output directory when done.

if world.size == 2:
    cpulayout1 = (1,2,1)
    cpulayout2 = (1,1,2)
elif world.size == 4:
    cpulayout1 = (1,2,2)
    cpulayout2 = (1,1,4)
elif world.size == 8:
    cpulayout1 = (2,1,4)
    cpulayout2 = (1,2,4)
elif world.size == 16:
    cpulayout1 = (2,2,4)
    cpulayout2 = (1,2,8)
else:
    raise RuntimeError("Not prepared to run on %s CPUs.  Supported: 2, 4, 8, 16." 
                       % (world.size,))

size = (10,10,50)
if quickie:
    steps1 = 1000
    steps2 = 1000
    takedata = 500
else:
    steps1 = 30000  # Phase 1
    steps2 = 1000000 # Phase 2
    takedata = 250000

friction = 0.0005

master = world.rank == 0


class Thermometer:
    def __init__(self, atoms, filename, blocks=50):
        self.atoms = atoms
        self.output = open(filename, "w")
        self.blocks = blocks

    def update(self):
        zmax = self.atoms.get_cell()[2,2]
        z = self.atoms.get_positions()[:,2]
        group = np.round(z/zmax * self.blocks).astype(int)
        group = np.clip(group, 0, self.blocks-1)
        p = self.atoms.get_momenta()
        ekin = (p*p).sum(axis=1) / (2 * self.atoms.get_masses())
        ekinsum = np.zeros(self.blocks)
        blocksum = np.zeros(self.blocks, int)
        for i in range(self.blocks):
            thisblock = np.equal(group, i)
            ekinsum[i] = ekin[thisblock].sum()
            blocksum[i] = thisblock.sum()
        world.sum(ekinsum)
        world.sum(blocksum)
        if master:
            np.clip(blocksum, 1, blocksum.max(), out=blocksum)
            T = ekinsum / ( 1.5 * units.kB * blocksum)
            for i in range(self.blocks):
                self.output.write("%i %.1f\n" % (i, T[i],))
            self.output.write("&\n")
            self.output.flush()
            
class AveragingThermometer:
    def __init__(self, atoms):
        self.atoms = atoms
        self.Tsum = 0
        self.n = 0
        
    def update(self):
        self.Tsum += self.atoms.get_temperature()
        self.n += 1
        
    def get_temp(self):
        return self.Tsum / self.n

def main(element, T_up, T_low, T_expect, bulkmodulus, makepot):
    def set_temp_fric(dyn, atoms):
        z = atoms.get_positions()[:,2]
        zmax = atoms.get_cell()[2,2]
        zmin = 0
        T = T_low + (T_up - T_low) * (z - zmin) / (zmax - zmin)
        dyn.set_temperature(T * units.kB)

    rnd = world.sum(random.randint(0,2**32)) # Insures the same number on all nodes.
    prefix = "melt-%s-%s" % (element, hex(rnd)[2:])
    if master:
        print "Using output directory", prefix
        os.mkdir(prefix)
    else:
        while not os.path.exists(prefix):
            time.sleep(1)
    if master:
        atoms = FaceCenteredCubic(symbol=element, size=size, pbc=(True, True, False))
        atoms.center(vacuum=10.0, axis=2)
    else:
        atoms = None
    atoms = MakeParallelAtoms(atoms, cpulayout1)
    print world.rank, '-', len(atoms), atoms.get_number_of_atoms()
    atoms.set_calculator(makepot())
    #view(atoms)
    
    dyn = Langevin(atoms, 5*units.fs, 0.0, friction)
    logger = MDLogger(dyn, atoms, prefix+'/Melt-phase1.log', stress=True, peratom=True)
    dyn.attach(logger, interval=100)
    set_temp_fric(dyn, atoms)
    unstress = Inhomogeneous_NPTBerendsen(atoms, 50*units.fs, 0, taup=500*units.fs,
                                          compressibility=1/bulkmodulus)
    dyn.attach(unstress.scale_positions_and_cell, interval=10)
    #traj = PickleTrajectory("melting-phase1.traj", "w", atoms)
    fn1 = prefix + "/Melt-phase1.bundle" 
    if master:
        print "Opening file", fn1
    traj = BundleTrajectory(fn1, "w", atoms, split=True)
    dyn.attach(traj, interval=500)
    
    therm = Thermometer(atoms, prefix+"/TempProfile-phase1.dat")
    dyn.attach(therm.update, interval=500)
    
    for i in range(steps1 / 100):
        dyn.run(100)
        set_temp_fric(dyn, atoms)
    if master:
        print "Closing file", fn1
    traj.close()
    del traj, atoms
    
    # Reread the bundle
    atoms = BundleTrajectory(fn1).get_atoms(-1, cpulayout2)
    atoms.set_calculator(makepot())
    dyn = VelocityVerlet(atoms, 5*units.fs)
    logger = MDLogger(dyn, atoms, prefix+'/PtMelt-phase2.log', stress=True, peratom=True)
    dyn.attach(logger, interval=1000)
    unstress = Inhomogeneous_NPTBerendsen(atoms, 50*units.fs, 0,
                                          taup=5000*units.fs,
                                          compressibility=1/bulkmodulus)
    dyn.attach(unstress.scale_positions_and_cell, interval=10)
    fn2 = prefix + "/Melt-phase2.bundle" 
    if master:
        print "Opening file", fn2
    traj = BundleTrajectory(fn2, "w", atoms)
    dyn.attach(traj, interval=steps2/100)
    
    therm = Thermometer(atoms, prefix+"/TempProfile-phase2.dat")
    dyn.attach(therm.update, interval=steps2/10)
    
    dyn.run(steps2 - takedata)
    therm2 = AveragingThermometer(atoms)
    dyn.attach(therm2.update, interval=100)
    dyn.run(takedata)  # Run the last part
    T = therm2.get_temp()
    if master:
        print "Melting temperature:", T
    ReportTest("Melting temperature", T, T_expect, 3)
    if cleanup:
        world.barrier()
        # Cleanup may fail due to the directory not being updated from
        # writes on the other nodes.  Wait for NFS file system.
        time.sleep(60)
        if master:
            shutil.rmtree(prefix)
    world.barrier()
    ReportTest.Summary()    

if __name__ == '__main__':
    main(element = 'Pt',
         T_up = 2000.0,  #K
         T_low = 1000.0, #K
         T_expect = 1357, #K
         bulkmodulus = 230.0e9 / 1.0e5,  # Bulk modulus of Pt in bar
         makepot = EMT
         )
