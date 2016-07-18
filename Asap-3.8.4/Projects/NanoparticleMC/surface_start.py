#!/usr/bin/env python
"""Starts a number of surface-MC simulations of Au.

Usage: cd /DesiredParentWorkingDir/
       python /ScriptDirectory/surface_start.py natoms temp steps nsimul queue [tgas] [pgas]

Where: 
       natoms is the number of atoms, 
       temp is the fake simulation temperature in Kelvin, 
       steps is number of simulation steps,
       nsimul is the number of simulations to start.
       queue is the niflheim queue parameter i.e. verylong long medium or small.
       tgas is the temperature used for calculating gas adsorption.
       pgas is the gas pressure.

The file AdsorptionParameters.py must be present in the running directory.
Copy it from Asap/Projects/NanoparticleMC.

"""
import cPickle as pickle
import subprocess
import os
import sys
from asap3.nanoparticle_mc.Logger import *
#Initialize the log file
log = Log(name="smc_log")

#Check if user asked for help with -h or --help or -doc
for a in sys.argv:
    if a=="-h" or a == "--help" or a == "-doc" or a == "--documentation":
        print >> sys.stderr, __doc__
        sys.exit(0)

if len(sys.argv)<6 and len(sys.argv)>=2:
    print >> sys.stderr, __doc__ #Throw error if too few parameters
    sys.exit(1)

#Some constants which can be changed:
tgas = 600 #600 K
pgas = 1E2 # 1 mbar
LC = 4.055 #lattice constant
species="Au"
gas_species = "CO"



if len(sys.argv)==1: #Makes arguments optional
    size = 500
    temp = 6000
    steps = 200000
    number = 6
    queue = "verylong" #Default queue
else:
    size = int(sys.argv[1])
    temp = int(sys.argv[2])
    steps = int(sys.argv[3])
    number = int(sys.argv[4])
    queue = str(sys.argv[5])
    if len(sys.argv)>6 and sys.argv[6] is not "--gas=False" and sys.argv[6] is not "--gas=0" and sys.argv[6]=="--gas=false":
        tgas = float(sys.argv[6])
        if len(sys.argv)>7:
            pgas = float(sys.argv[7])
    
#Add gas based on input parameters, default is to add gas:
calc_gas = True #by defualt true
if len(sys.argv)==7 and (sys.argv[6]=="--gas=False" or sys.argv[6]=="--gas=0" or sys.argv[6]=="--gas=false"):
    calc_gas = False
    
    

#Check queue param for validity
valid_queue_str = ["verylong","long","medium","small"]
if (queue in valid_queue_str)==False :
    sys.exit("Error: \n queue must be one of the following: verylong, long, medium or small")

dir = 't%is%ia%i_smc' % (temp, steps,size) #MJ: Changed to include Natoms as name of directory.
n = 1
while os.path.exists(dir):
    if calc_gas:
        dir = 't%is%ia%i_smc_gas_%i' % (temp, steps,size, n)
    else:
        dir = 't%is%ia%i_smc_%i' % (temp, steps,size, n)
    n += 1
print "Working in directory", dir
os.mkdir(dir)
os.chdir(dir)
os.symlink('../AdsorptionParameters.py', 'AdsorptionParameters.py')

template = """#!/usr/bin/env python

#PBS -N Au##SIZE##_smc_##NUMBER##
#PBS -e cluster_##NUMBER##.err
#PBS -o cluster_##NUMBER##.log
#PBS -m n
#PBS -q ##QUEUE##
#PBS -l nodes=1:ppn=1:opteron4


from asap3.nanoparticle_mc.montecarlo import *
from asap3 import EMT
from asap3.nanoparticle_mc.AdsCalc import adscalc


a = ##LC##
size = ##SIZE##
temp = ##TEMP##
steps = ##STEPS##

layers = get_layers_from_size(size)
file = 'cluster_##NUMBER##'

smc = SurfaceMonteCarlo("Au", layers, size, a, debug=0)
tempcalc = EMT()
##CALC##

smc.run(steps, temp, file + '.txt', file + '.smc')
smc.data.write()
smc.data.filter_surfaces()
smc.data.write(file + '_sym.smc')

"""

oops = open("oops.sh", "w")
oops.write("#!/bin/bash\n\n")
filenames = []

for i in range(number):
    file = "cluster_%02i.py" % (i,)
    filenames.append('cluster_%02i_sym.smc' % (i,))
    script = open(file, "w")
    s = template.replace("##NUMBER##", "%02i" % (i,))
    s = s.replace("##SIZE##", str(size))
    s = s.replace("##TEMP##", str(temp))
    s = s.replace("##STEPS##", str(steps))
    s = s.replace("##QUEUE##", queue)
    s = s.replace("##LC##",str(LC))
    #s = s.replace("##SPECIES##",species)
    if calc_gas == True:
        s = s.replace("##CALC##","smc.set_calculator(adscalc(tempcalc,temperature="+str(tgas)+",pressure="+str(pgas)+"))")
    else:
        s = s.replace("##CALC##","smc.set_calculator(tempcalc)")
    script.write(s)
    script.close()

    qsub = ["qsub", file]
    qsubproc = subprocess.Popen(qsub, stdout=subprocess.PIPE, close_fds=True)
    (out, err) = qsubproc.communicate()
    errcode = qsubproc.wait()
    if errcode:
        print >>sys.stderr, "qsub failed with error code", str(errcode)
        print >>sys.stderr, "Command line:", qsub
        sys.exit("qsub failed")
    print "JOB ID:", out
    oops.write("qdel %s\n" % (out,))

oops.close()
os.system("chmod +x oops.sh")
os.chdir('..')
f = open(dir + '_data.pickle', 'w')
pickle.dump(filenames, f)
f.close()
#Now log the parameters
params = {"Name":"SurfaceMonteCarlo", "latticeconstant":LC,
"MonteCarlo_Temperature": temp,"N Atoms":size,"Nsteps":steps,"Fcc_species":species}

if calc_gas==True:
    params['Gas_Temperature']=tgas
    params['Gas_Pressure'] = pgas
    params['Gas_species'] = gas_species
log.dumpparamdict(params)
