import numpy as np

from asap3.Tools.ParameterOptimization.Optimization import ParameterOptimization
from asap3.Tools.ParameterOptimization.Search import ParameterSearch

class ParameterSampling(ParameterSearch):
    def __init__(self, elements, calc, initparam, varparam, stdparam, points,
                 quantities, latticeconstants, calc_args=None, debug=False):
        self.points = points
        ParameterSearch.__init__(self, elements, calc, initparam, varparam,
                                 stdparam, quantities, latticeconstants,
                                 calc_args, debug=0)

    def run(self, xtol=5e-4, ftol=2e0, delta=5e-2, dat=None, txt=None):
        self.nsteps = 0
        self.numberlist = []
        self.errorlist = []
        self.fitparlist = []
        self.optparlist = []
        self.optvallist = []

        fitparam_list = []
        self.sampling(fitparam_list)
        #print len(fitparam_list)
        #for pars in fitparam_list[163:175]:
        #    for (k, p) in pars.items():
        #        print k, p
        #    print "-----"
        #return

        for fitpar in fitparam_list:
            self.nsteps += 1

            # Run optimization
            opt = ParameterOptimization(self.elements, self.calc, fitpar,
                                        self.varparam, self.values,
                                        self.latticeconstants, self.calc_args,
                                        debug=self.debug)
            (error, optpar, optval) = opt.fit(xtol, ftol, delta, txt)

            self.save(error, fitpar, optpar, optval)

        self.write(dat, 'sampling')

    def sampling(self, output, tempparam=None, n=1):
        # Which parameter should be altered
        varcount = 0
        for (key, var) in self.varparam.items():
            for i, v in enumerate(var):
                if v:
                    varcount += 1
                    if varcount == n:
                        k = key
                        m = i
        #print "Varying %s at position %i (%i out of %i)" % (k, m, n, varcount)

        # What is the "standard diviation"
        if isinstance(self.stdparam[k], list):
            std = self.stdparam[k][m]
        else:
            std = self.stdparam[k]

        # Vary the parameter
        ps = self.initparam[k][m] * np.linspace(1.0 - std, 1.0 + std, self.points)
        for p in ps:
            newparam = {}
            for key in self.initparam.keys():
                if tempparam is None:
                    newparam[key] = self.initparam[key][:]
                else:
                    newparam[key] = tempparam[key][:]
            newparam[k][m] = p
            if n < varcount:
                self.sampling(output, newparam, n+1)
            else:
                output.append(newparam)


