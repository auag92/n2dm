#PBS -l nodes=1:ppn=8:xeon8
#PBS -q verylong
#PBS -N amc_end
#PBS -m ae
"""Ends atom simulation and gives result in file.
Usage:
       asap-qsub endatomsimulation.py amcfolder e_filter T_real (result_log_file)
Where: e_filter is the maximal accepted energy w.r.t. the global min such that a cfg. counts.
       T_real is the real gas temperature in K
       result_log_file is the desired name of result file. sys.stdout if not specified. 
"""
import sys
import os
from asap3.nanoparticle_mc.atommontecarlodata import AtomMonteCarloData
import numpy as np
from ase import units #used to get kB in eV/K
from mc_result import *
#Check if user asked for help with -h or --help or -doc
for a in sys.argv:
    if a=="-h" or a == "--help" or a == "-doc" or a == "--documentation":
        print >> sys.stderr, __doc__
        sys.exit(0)

fdir = sys.argv[1] #Directory of .amc files to join
e_filter = float(sys.argv[2])
T = float(sys.argv[3]) #The real temperature e.g. 300 K 

#First redirect all stdout to log file:
if len(sys.argv)>=5:
    logfile = str(sys.argv[4])
    old_stdout = sys.stdout
    if os.path.exists(logfile):
        os.system("mv "+logfile+" "+logfile+".bck")
    log_file = open(logfile,"w")
    sys.stdout=log_file
    print "Atom Monte Carlo Simulation Summary \n -------------------------------------------"



files = []
os.chdir(fdir)

for f in os.listdir("./"):
    if f.endswith(".amc.gz"):
        files.append(f)
#Now we know the file names. Read merge the informations.
d=AtomMonteCarloData()
dtemp = AtomMonteCarloData()
#Init arrays
relaxed_energies = np.array([])
unrelaxed_energies = np.array([])
mc_energies = np.array([])
real_energies = np.array([])
counts = np.array([])
types = np.array([])
coordinations = None
vac_coordinations=None
vac_types = np.array([])
multiplicities = np.array([], int)
eglobalmin=None
eglobalmin_id = None

def delmultiple(listtohandle,idstorem):
    offset = 0
    for i in idstorem:
        del listtohandle[i-offset]
        offset+=1
    return listtohandle

for f in files: #Now cut half the MC steps off, to only consider Thermal equilibrium.
    #At the same time find the global min. energy
    dtemp.read(f) #loads pickle to dtemp!
    id_r = dtemp.id_relax
    relaxed_energies = np.append(relaxed_energies,dtemp.lists['relaxed_energies'][id_r:])
    
    unrelaxed_energies = np.append(unrelaxed_energies,dtemp.lists['unrelaxed_energies'][id_r:])
    
    mc_energies = np.append(mc_energies,dtemp.lists['energy'][id_r:])
    
    #Calculate the real energies:
    r2append = np.array(dtemp.lists['energy'][id_r:]) + np.array(dtemp.lists['relaxed_energies'][id_r:])-np.array(dtemp.lists['unrelaxed_energies'][id_r:])
    
    real_energies = np.append(real_energies,r2append)

    
    
    counts = np.append(counts,dtemp.lists['counts'][id_r:])

    nn = len(dtemp.lists['relaxed_energies'])
    multiplicities = np.append(multiplicities, dtemp.multiplicity * np.ones(nn, int)[id_r:])
    
    types =  np.append(types,dtemp.lists['types'][id_r:])
    #print dtemp.lists['coordinations'][int((len(dtemp.lists['coordinations'])-1)/2):]
    #print "ID Relax " + str(id_r)
    coord = dtemp.lists['coordinations'][id_r:]
    
    if coordinations is None:
        coordinations = np.zeros((0, len(coord[0])), int) 
    coordinations = np.vstack([coordinations, coord])
    del coord
    
    
    #If lower energy is found
    emin = real_energies.min()
    if eglobalmin is None or emin<eglobalmin:
        eglobalmin = emin #The new global min is stored.
        eglobalmin_id = real_energies.argmin() #indexof equivalent
        eglobalmin_file = f
    #Store temperature:
    Tsim = dtemp.temp
    
    #Store latticeconstant:
    lc = dtemp.lattice_constant
    if lc is None:
        lc = 4.055 #The Au fcc-EMT lattice constant

print "Found Global Min " +str(eglobalmin)+ "eV"
print "Global min lies in file: "+str(eglobalmin_file)+" with index "+str(eglobalmin_id)
print "\n ---------------------------------\n"

print "Tsim = "+str(Tsim)
print "Tcorrect = "+str(T)
print "E_filter = "+str(e_filter)       
#Now we know which file contains the min.en config and which energy aswell as the id in
# the combined array(energies)
# Use the big arrays to filter out.
badidx = [] #Store the non-relevant configurations ids.
for i in range(len(real_energies)):
    if real_energies[i] > eglobalmin+e_filter:
        badidx.append(i)

print "counts min,max,shape:", counts.min(), counts.max(), counts.shape
    
#Now remove these indices from the arrays
relaxed_energies = np.delete(relaxed_energies,badidx)
unrelaxed_energies = np.delete(unrelaxed_energies,badidx)
mc_energies = np.delete(mc_energies,badidx)
real_energies = np.delete(real_energies,badidx)
counts = np.delete(counts,badidx)
multiplicities = np.delete(multiplicities,badidx)
types = np.delete(types,badidx)
coordinations = np.delete(coordinations,badidx,axis=0)
#vac_coordinations = delmultiple(vac_coordinations,badidx)
#vac_types = np.delete(vac_types,badidx)
print "Counting configurations" + str(len(counts))
print "Found "+str(len(badidx))+ " irrelevant configurations \n ---------------------------------\n"

print "counts min,max,shape:", counts.min(), counts.max(), counts.shape



#Now compute Boltzmann averages with the modified BZ-factor
p = counts*multiplicities*np.exp(-(real_energies-eglobalmin)/(units.kB*T)+
                   (mc_energies-mc_energies.min())/(units.kB*Tsim))
#Compute the normalization factor:
factnorm = float(p.sum())

#Normalize the weights.
if not factnorm==0.:
    p = p/factnorm
print "Sum of weights: "+str(p.sum())+"\n --------------------------------- \n"

#Compute Bz averages:
emean = (real_energies*p).sum()
#Make np.array hold the bz average of coordination 0 on [0] and 1 on [1]...
#c = [c_0,c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9,c_10,c_11,c_12]
#cormean = [(cm*p).sum() for cm in c]
print coordinations.shape
cormean = [(coordinations[:,i]*p).sum() for i in range(13)]
ix=0
for cc in cormean:
    print "Found average "+str(ix)+"-coordinated atoms: "+str(cc)
    ix+=1
print "\n --------------------------------- \n"
print "Sum of average coordinated atoms counted: "+str(sum(cormean))

print "number of atoms: "+str(coordinations[0].sum())

print "Found mean energy "+str(emean)+" eV"
natoms = float(coordinations[0].sum())
nunitcells = natoms/4.0
diam = lc*(6*nunitcells/np.pi)**(1.0/3)

print "Particle Diameter Under Spherical assumption: "+str(diam) + " Ang"

#Dump result to pickle here an example of vacuum, if gas you could put in Tgas,Pgas:
#filename=None makes the result determine filename upon pickling.
res = MC_Result(av_coordi = cormean,av_energy=emean,tgas=T,tamc=Tsim,natoms=natoms,diameter=diam)
res.PickleIt()
