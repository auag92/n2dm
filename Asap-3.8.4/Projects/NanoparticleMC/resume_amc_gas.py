#PBS -l nodes=20:ppn=4:opteron4
#PBS -q verylong
#PBS -N amc_n100_conv1
#PBS -m ae
import os
from montecarlo import SurfaceMonteCarloData
from ase.cluster.cubic import FaceCenteredCubic
from ase.cluster import data
from asap3.MonteCarlo.Metropolis import Metropolis
from asap3.MonteCarlo.Moves import SurfaceMove
from asap3 import EMT
import numpy as np
from ase.visualize import view
from resizecluster import resizecluster
from ase.io.trajectory import PickleTrajectory
import sys
from atommontecarlodata import AtomMonteCarloData
from ase.parallel import world
from AdsCalc import adscalc
from time import time,sleep
import pickle
#Change note: Added gas option, check for indentation tab vs. spaces error.
#Added resume option.
#Arguments:
filename = sys.argv[1]
temperature = float(sys.argv[2])
nsteps = int(sys.argv[3])
outdir= sys.argv[4]

tgas = float(sys.argv[5])
pgas = float(sys.argv[6])

species = "AuCO"

def read_and_do_montecarlo(filename,use_gas):
    d = SurfaceMonteCarloData()
    d.read(filename)
    print "Starting "+str(len(d))+" sims."
    surfaces = data.fcc.surface_names
            
    #for n in range(0,len(d)):
    for n in range(world.rank,len(d),world.size):
        file = outdir+"/a%05i.amc.gz" % n
        if not os.path.exists(file):
            layers = d[n][1]  # Really d[n]["layers"]
            atoms = FaceCenteredCubic(d.atomic_number,
            surfaces, layers,latticeconstant=d.lattice_constant)
            resizecluster(atoms, d.fitsize)
            print "Resized number of atoms:", len(atoms)
            do_monte_carlo(atoms,n,outdir,use_gas)
    world.barrier()#Let the cpu's wait until all in same state.


def do_monte_carlo(atoms,iteration,outdir,use_gas):
    
    tempcalc = EMT()
    if use_gas==True:
        atoms.set_calculator(adscalc(tempcalc,temperature=tgas,pressure=pgas,species=species))
    else:
        atoms.set_calculator(tempcalc)
    Esmc = atoms.get_potential_energy()
    mc = Metropolis(atoms=atoms,log=None)
    surfmove =SurfaceMove()
    mc.attach_move(surfmove)
    
    outfilename = "a%05i.amc" % iteration
    
    amcd = AtomMonteCarloData(atoms=atoms,surfmove=surfmove,temp=temperature,filename=outfilename,Esmc=Esmc)
    mc.attach_observer(amcd.accept_move) #Because default event is at acceptmove
    mc.attach_observer(amcd.reject_move, attime='reject')
    amcd.accept_move() #We need to write the first configuration, 
    
    mc.run(nsteps, temp=temperature)
    
    amcd.write(os.path.join(outdir,outfilename))

        
read_and_do_montecarlo(filename,True)
