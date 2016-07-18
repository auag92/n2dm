"""Check that good parameters have been found by more than one process.

Usage:

python check_result.py dir [dir2 ...]

"""

import production
from asap3.Tools.ParameterOptimization.GetParameters import *
import sys
import glob
import os
import numpy as np

dirs = sys.argv[1:]

if not dirs:
    print __doc__
    sys.exit(1)

files = []
for d in dirs:
    files.extend(glob.glob(d+"/fit-*.dat"))
plur = ["y", "ies"]
print ("Found %i files in %i director%s."  
       % (len(files), len(dirs), plur[len(dirs) != 1]))
 
element = None
parameters = []
errfuncs = []
okfiles = []
for f in files:
    #elem, params, errfunc = get_parameters(f)
    try:
        param_dict, errfunc = get_parameters(f)
    except ValueError:
        continue
    elem = param_dict.keys()
    if element is None:
        element = elem
    elif elem != element:
        print ("Inconsistent element specification: %s != %s"
               % (elem, element))
        sys.exit(2)
    params = []
    for e in elem:
        params.extend(param_dict[e])
    parameters.append(np.array(params))
    errfuncs.append(errfunc)
    okfiles.append(f)

order = np.argsort(errfuncs)

def diff(p1, p2):
    d = (p1 - p2) / p2
    return np.sqrt((d*d).sum()/len(d))

i0 = order[0]
n = 0
for i in order:
    n += 1
    print ("%3i %40s: %8.5f %10.2f %6.2f"
           % (n, okfiles[i], diff(parameters[i], parameters[i0]),
              errfuncs[i], errfuncs[i] - errfuncs[i0]))
