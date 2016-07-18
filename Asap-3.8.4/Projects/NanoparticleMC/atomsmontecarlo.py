#PBS -l nodes=14:ppn=4:opteron4
#PBS -q verylong
#PBS -N amc_n300_conv3
#PBS -m ae

"""Starts amc simulations.
Usage:
       asap-qsub atomsmontecarlo.py smc_file sim_temp steps (tgas pgas)/--gas=False
Where: smc_file is the path to .smc file to start these sims. from, 
       sim_temp is the fictious temperature in Kelvin, 
       steps is number of simulation steps for one smc configuration,
       tgas is the real gas temperature in Kelvin.
       pgas is the gas pressure in Pa.
       --gas=False is used to calculate in vacuum, it must replace tgas and pgas to work.
"""
import os
from asap3.nanoparticle_mc.montecarlo import SurfaceMonteCarloData
from ase.cluster.cubic import FaceCenteredCubic
from ase.cluster import data
from asap3.MonteCarlo.Metropolis import Metropolis
from asap3.MonteCarlo.Moves import SurfaceMove
from asap3 import EMT, MonteCarloEMT, MonteCarloAtoms
import numpy as np
from ase.visualize import view
from asap3.nanoparticle_mc.resizecluster import resizecluster
from ase.io.trajectory import PickleTrajectory
import sys
from asap3.nanoparticle_mc.atommontecarlodata import AtomMonteCarloData
from ase.parallel import world
from asap3.nanoparticle_mc.AdsCalc import adscalc
from ase.cluster.cluster import Cluster #as ase_Cluster
from time import time,sleep
from asap3.nanoparticle_mc.Logger import *

#Check if user asked for help with -h or --help or -doc
for a in sys.argv:
    if a=="-h" or a == "--help" or a == "-doc" or a == "--documentation":
        print >> sys.stderr, __doc__
        sys.exit(0)

#Arguments:
filename = sys.argv[1]
temperature = float(sys.argv[2])
nsteps = int(sys.argv[3])
species = "AuCO"

#standard constants:
pgas=1E2 #1mbar
tgas = 300 #600 K

# Make sure that we use the Asap Cluster class instead of the ASE one.
assert hasattr(FaceCenteredCubic, "Cluster")
class MonteCarloCluster(MonteCarloAtoms,Cluster):
    pass
FaceCenteredCubic.Cluster = MonteCarloCluster

#Logging initialization:
log = Log(name ="amc")
log.dumpstring("Atoms Monte Carlo Log:")
params = {"Smc Filename":filename,"Simulation Temperature": str(temperature),"MC Steps":str(nsteps)}

if len(sys.argv)>=5 and (sys.argv[4]=="--gas=false" or sys.argv[4]=="--gas=False"
 or sys.argv[4]=="--gas=0"):
    #print "Not using gas"
    use_gas = False
    if world.rank==0:
        log.dumpstring("Vacuum Simulation")
else:
    #print "Using gas"
    if world.rank==0:
        log.dumpstring("Gas Simulation on "+str(species))
    
    use_gas = True
    if len(sys.argv) >= 5:
        tgas = float(sys.argv[4])
    if len(sys.argv) >= 6:
        pgas = float(sys.argv[5])
    
    params['Gas Species'] = species
    params['Gas Temperature'] = tgas
    params['Gas Pressure'] = pgas


def read_and_do_montecarlo(filename,use_gas):
    d = SurfaceMonteCarloData()
    d.read(filename)
    print "Starting "+str(len(d))+" sims."
    surfaces = data.fcc.surface_names
    #Only one worker should create the filename
    outdir = determine_create_dirname(filename)
    
        
    #for n in range(0,len(d)):
    for n in range(world.rank,len(d),world.size):
        layers = d[n][1]  # Really d[n]["layers"]
        multiplicity = d.get_multiplicity(n)
        atoms = FaceCenteredCubic(d.atomic_number,
                                  surfaces, layers,latticeconstant=d.lattice_constant)
        print "Number of atoms:", len(atoms)
        resizecluster(atoms, d.fitsize)
        print "Resized number of atoms:", len(atoms)
        do_monte_carlo(atoms,n,outdir,use_gas,multiplicity)
    world.barrier()#Let the cpu's wait until all in same state.



def do_monte_carlo(atoms,iteration,outdir,use_gas,multiplicity):
    tempcalc = MonteCarloEMT()
    if use_gas==True:
        atoms.set_calculator(adscalc(tempcalc,temperature=tgas,pressure=pgas,species=species))
    else:
        atoms.set_calculator(tempcalc)
    Esmc = atoms.get_potential_energy()
        
    mc = Metropolis(atoms=atoms,log=None)
    surfmove =SurfaceMove()
    mc.attach_move(surfmove)
    
    outfilename = "a%05i.amc" % iteration
    
    amcd = AtomMonteCarloData(atoms=atoms, surfmove=surfmove,
                  temp=temperature, filename=outfilename,
                  Esmc=Esmc, multiplicity=multiplicity,total_steps = nsteps)
    mc.attach_observer(amcd.accept_move) #Because default event is at acceptmove
    mc.attach_observer(amcd.reject_move, attime='reject')
    amcd.accept_move() #We need to write the first configuration, 
    
    mc.run(nsteps, temp=temperature)
    
    amcd.write(os.path.join(outdir,outfilename))
    

def determine_create_dirname(filename):
    outdir, ext = os.path.splitext(filename)
    assert(ext == ".smc")
    suffix = 1
    if use_gas:
        outdirtmp = outdir+"_amc_gas"
    else:
        outdirtmp = outdir+"_amc"
    while os.path.exists(outdirtmp):
        if use_gas:
            outdirtmp = outdir+"_amc_gas_"+str(suffix)
        else:
            outdirtmp = outdir+"_amc_"+str(suffix)
        suffix+=1
    
    #Create this new directory:
    world.barrier()
    if world.rank == 0:
        os.mkdir(outdirtmp)
    else:
        while not os.path.exists(outdirtmp):
            sleep(1)
    #Now we have a unique directory, return it.
    return outdirtmp
#Dump params to log:
if world.rank==0:
    log.dumpparamdict(params)
read_and_do_montecarlo(filename,use_gas)
