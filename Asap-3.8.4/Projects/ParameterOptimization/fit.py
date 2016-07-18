#!/usr/bin/env python

#PBS -N EMT-${symbol}
#PBS -m n
#PBS -q long
#PBS -l nodes=1:ppn=8:xeon8

import numpy as np

from asap3.Tools.ParameterOptimization.SearchParallel import ParameterSearch
from asap3.Tools.ParameterOptimization.EMT import EMT2011Fit, EMTStdParameters

symbol = '${symbol}'

# Parameter names: ['eta2', 'lambda', 'kappa/delta', 'E0', 'V0', 'S0', 'n0']
initparam = {(symbol,symbol): EMTStdParameters(symbol, 'delta'),
             }

varparam = {(symbol,symbol): [True, True, True, True, True, True, False]}

quantities = [('lattice_constant_a', 'fcc', symbol, ${a}),
              ('bulk_modulus', 'fcc', symbol, ${B}),
              ('elastic_anisotropy', 'fcc', symbol, ${A}),
              #('elastic_constant_C11', 'fcc', symbol, ${C11}),
              #('elastic_constant_C12', 'fcc', symbol, ${C12}),
              ('elastic_constant_C44', 'fcc', symbol, ${C44}),
              ('cohesive_energy', 'fcc', symbol, ${Ecoh}),
              ('surface_energy', 'fcc111', symbol, ${E111}),
              #('surface_energy', 'fcc100', symbol, ${E100}),
              ('surface_ratio', 'fcc111-fcc100', symbol, ${E111_100}),
              ('force_match', '../'+symbol+'_force_PBEsol.traj', symbol, 0, ${fmatch}), 
              ('stacking_fault', 'fcc', symbol, ${Esf}),
              ]

latticeconstants = [['fcc', symbol, ${lc}],
                    ]

opt = ParameterSearch([symbol], EMT2011Fit, initparam, varparam, quantities,
                      latticeconstants, ('delta',), False)
opt.run(0.1, 0.144, '12.5h', xtol=1e-3, ftol=5e-2, delta=1e-2, dat='fit.dat',
        log='fit.log', err='fit.err')

# First parameter: Gauss width for MC part
# Second parameter: MC Temperature
# Third parameter: run time, in h or m.  Or the number of steps if integer.
# xtol and ftol are simplex convergence criterial.  
# delta is the initial simplex size, relatively.

