#!/usr/bin/env python

#PBS -N EMT-alloys
#PBS -m n
#PBS -q long
#PBS -l nodes=2:ppn=4:opteron4

import numpy as np

from asap3.Tools.ParameterOptimization import ParameterPerformance
from asap3.Tools.ParameterOptimization.Optimization import ParameterOptimization
from asap3.Tools.ParameterOptimization.SearchParallel import ParameterSearch
from asap3.Tools.ParameterOptimization.EMT import EMT2011Fit, EMTStdParameters

# Parameter names: ['eta2', 'lambda', 'kappa/delta', 'E0', 'V0', 'S0', 'n0']
initparam = {${initparam}}

varparam = {${varparam}}

quantities = [${quantities}]

latticeconstants = [${latticeconstants}]

calculator = EMT2011Fit(${symbols}, initparam, 'delta')
ParameterPerformance(calculator, quantities, latticeconstants, debug=False)

#opt = ParameterOptimization(${symbols}, EMT2011Fit, initparam, varparam,
#                            quantities, latticeconstants, ('delta',), False)
#(error, optpar, optval) = opt.fit(xtol=1e-3, ftol=5e-2, delta=5e-4,
#                                  log='fit.log', err='fit.err')
#print 'Optimization ended with error function %.6e' % (error,)
#opt.write_result(initparam, optpar, optval)

opt = ParameterSearch(${symbols}, EMT2011Fit, initparam, varparam, quantities,
                      latticeconstants, ('delta',), False)
opt.run(0.1, 0.144, '8.0h', xtol=1e-3, ftol=5e-2, delta=1e-3, dat='fit.dat',
        log='fit.log', err='fit.err')

