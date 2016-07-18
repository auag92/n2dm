#!/usr/bin/env python

#PBS -N Cu1000_se
#PBS -e se.err
#PBS -o se.log
#PBS -m ae
#PBS -q long
#PBS -l nodes=1:ppn=1:opteron4
"""Prepares a finished surface MC simulation for atoms MC by filtering and symmetry elimination.
Usage: cd simparentfolder 
python surface_end.py direc Tmc Nleft smc_log_file
where Tmc is the temperature used in the surface monte carlo simulation.
where direc is the path to finished surfaceMC simulations.
where Nleft is the number of cfg. to Boltzmann-choose.
where smc_log_file is the filename of the smc_Start log.
"""
from ase import units
import cPickle as pickle
import os
import numpy as np
from ase import *
from asap3.nanoparticle_mc.montecarlo import *
import sys
from asap3.nanoparticle_mc.Logger import *
#Check wheter user asked for doc:
for a in sys.argv:
    if a=="-h" or a=="--help" or a=="-doc" or a=="--documentation":
        print  __doc__
        sys.exit(0)
#Check for inputs:
print(str(len(sys.argv)))

if len(sys.argv)<3:
    print >> sys.stderr, __doc__
    sys.exit(1)

def Boltzmann_Filtering(d,Tmc):
    #Allowed deviation in Nleft:
    Ntolerance = 1000
    #Find minimum energy
    T = Tmc
    s = d.get_sizes()
    energies = (d.get_energies())*s.mean() #To have energies not per atom.
    min_energy = min(energies)
    
    #Now Adjust T until the sum of P is 10000 e.g.
    found = False #used in while loop
    while (not found):
        probabilities = np.exp(-(energies-min_energy)/units.kB/T)
        N = sum(probabilities)
        print N
        if N< Nleft+Ntolerance and N>Nleft:
            log.dumpstring("Boltzmann Filtering: Ended with Sum over probabilities "+str(N))
            log.dumpstring("The temperature was here: "+str(T))
            found = True    
        if N<0:
            raise RuntimeError("Oops, the fictious temperature changed too fast")
        if not found: 
            if N>=Nleft+Ntolerance:
                T-=0.5
            else:
                T+=2.0
    
    print "Last temperature in loop"+str(T)    
    #Now given the probabilities do the monte carlo filtering:
    random_floats = np.random.random(len(energies))
    chosen_idx = []
    for cf in range(len(energies)):
        if probabilities[cf] >= random_floats[cf]:
            chosen_idx.append(cf)
    #Now given the chosen indices make the new array to return
    print "Got to replacement of arrays"
    for name, a in d.arrays.items():
            d.arrays[name] = a[chosen_idx]
    return d



#Default values
Tmc = 4000.0
Nleft = 10000

dir = str(sys.argv[1])
Tmc = float(sys.argv[2])
Nleft = int(sys.argv[3]) # The number of configurations to end up with
#Logging
smc_log = sys.argv[4]



#Load the smc log to merge them:
log = Log(name="smc_end")
log.load2append(smc_log) #Write the smc_start.py log into the new log now.
params = {"Smc_Directory":dir,"T_mc_filter":Tmc,"Count_Filter_Number":Nleft}
log.dumpparamdict(params)


f = open(dir + '_data.pickle')
filenames = pickle.load(f)
f.close()

d = SurfaceMonteCarloData()
d.read(dir + '/' + filenames[0])
d1 = SurfaceMonteCarloData()
for file in filenames[1:]:
    print "reading"+str(file)
    d1.read(dir + '/' + file)
    d.extend(d1)

d.write(dir + '.smc')
log.dumpstring("I have written to dir+.smc")
d.read(dir + '.smc')
log.dumpstring("Now I've read dir+.smc")
#Now filter Boltzmann probabilities
d = Boltzmann_Filtering(d,Tmc)
#
log.dumpstring("I have filtered by Boltzmann prob. current n.o. configurations: " +str(len(d.get_energies())))
d.filter_surfaces()
log.dumpstring("I have filtered by symmetry, begin writing the output file+_surf_end.smc")

d.write(dir + '_N'+ str(Nleft)+ '_surf_end.smc')
