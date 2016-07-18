#!/usr/bin/env python

import os
import sys
import time
import subprocess
import numpy as np

smc_file = sys.argv[1]
id_digits = 6
v0 = 1
n = 2
name = 'AuO_Amc'

outdir = smc_file.split('.smc')[0]+"_amc_gas" #string to split is e.g. t4000s500000a300_smc_e1.0_surf_end.smc and folder is named t4000s500000a300_smc_e1.0_surf_end_amc_gas

script = '~/Bachelor/developed-scripts/mcscripts/atomsmontecarlo.py '+smc_file+' 1200 1500000 300 100'
resumescript = '~/Bachelor/developed-scripts/mcscripts/resume_amc_gas.py '+smc_file+' 1200 1500000 '+outdir+' 300 100'

def submit(command):
    qsubproc = subprocess.Popen(command, shell=True,
                                stdout=subprocess.PIPE)
    (out, err) = qsubproc.communicate()
    errcode = qsubproc.wait()
    if errcode:
        print >>sys.stderr, "qsub failed with error code", str(errcode)
        print >>sys.stderr, "Command line:", command
        sys.exit("qsub failed")
    print "JOB ID:", out
    return out

oopsname = time.strftime('oops-%d-%b-%Y--%H-%M-%S.sh')
oopsfile = None

for i in range(v0,n+1):
    if i == v0:
        jobid = submit('asap-qsub -N %s_v%i %s ' % (name, i,script))
    else:
        jobid = submit('asap-qsub -N %s_v%i -W depend=afternotok:%s %s ' % (name, i, old_jobid,resumescript))
    if oopsfile is None:
        oopsfile = open(oopsname, 'a', 1)
    
    #print jobid.split("JOB ID: ")[1]
        
    old_jobid = jobid.split("JOB ID: ")[1].split(".")[0]
    oopsfile.write('qdel %s  # %s\n' % (jobid, dir))

if oopsfile is not None:
    oopsfile.close()
    print "An oops file was written:", oopsname
