import sys
import time
import numpy as np

from asap3.Tools.ParameterOptimization.Optimization import ParameterOptimization

class ParameterSearch:
    def __init__(self, elements, calc, initparam, varparam, stdparam,
                 quantities, latticeconstants, calc_args=None, debug=False):
        self.elements = elements
        self.calc = calc
        self.initparam = initparam
        self.varparam = varparam
        self.stdparam = stdparam
        self.values = quantities
        self.latticeconstants = latticeconstants
        self.calc_args = calc_args
        self.debug = debug

    def get_new_parameters(self, parameters):
        new_parameters = {}
        for key in parameters.keys():
            std = self.stdparam[key]
            mask = np.array(self.varparam[key], float)
            change = np.random.uniform(1.0 - std, 1.0 + std, len(mask))
            change = mask * change - (mask - 1.0)
            new_parameters[key] = (np.array(parameters[key]) * change).tolist()
        return new_parameters

    def run(self, timer, type='follow', xtol=1e-1, ftol=5e-2, delta=5e-2,
            dat=None, txt=None):
        self.nsteps = 0
        self.numberlist = []
        self.errorlist = []
        self.fitparlist = []
        self.optparlist = []
        self.optvallist = []

        while True:
            self.nsteps += 1
            # Make a random set of parameters close to the
            if type == 'stay':
                if self.nsteps == 1:
                    fitpar = self.initparam
                else:
                    fitpar = self.get_new_parameters(self.initparam)
            elif type == 'follow':
                if self.nsteps == 1:
                    fitpar = self.initparam
                else:
                    fitpar = self.get_new_parameters(optpar)
            else:
                raise ValueError("Cannot determine the search type %s" % (type,))

            # Run optimization
            opt = ParameterOptimization(self.elements, self.calc, fitpar,
                                        self.varparam, self.values,
                                        self.latticeconstants, self.calc_args,
                                        debug=self.debug)
            (error, optpar, optval) = opt.fit(xtol, ftol, delta, txt)

            self.save(error, fitpar, optpar, optval)

            if isinstance(timer, int):
                if self.nsteps >= timer:
                    break
            else:
                if timer.done():
                    break

        self.write(dat)

    def save(self, error, fitpar, optpar, optval):
        # Sort and save the 10 best results
        length = len(self.errorlist)
        if length == 0:
            self.numberlist.append(self.nsteps)
            self.errorlist.append(error)
            self.fitparlist.append(fitpar)
            self.optparlist.append(optpar)
            self.optvallist.append(optval)
        else:
            new_min = -1
            new_rank = -1
            for i, err in enumerate(self.errorlist):
                op1 = np.empty(0)
                op2 = np.empty(0)
                for key in self.initparam.keys():
                    mask = np.array(self.varparam[key])
                    op1 = np.append(op1, np.array(optpar[key])[mask])
                    op2 = np.append(op2, np.array(self.optparlist[i][key])[mask])

                ov1 = np.empty(0)
                ov2 = np.empty(0)
                for v1, v2 in zip(optval, self.optvallist[i]):
                    ov1 = np.append(ov1, v1[3])
                    ov2 = np.append(ov2, v2[3])

                ov_dev_max = (2 * np.abs(ov1 - ov2) / np.abs(ov1 + ov2)).max()
                op_dev_max = (2 * np.abs(op1 - op2) / np.abs(op1 + op2)).max()

                if ov_dev_max < 2e-2 and op_dev_max < 5e-2 and new_min < 0 and i > 0:
                    new_min = i

                if error < err and new_rank < 0 and i > 0:
                    new_rank = i

            #print new_min, new_rank, length
            if new_min < 0:
                # New minimum found
                if new_rank > 0:
                    # Better error function
                    self.numberlist.insert(new_rank, self.nsteps)
                    self.errorlist.insert(new_rank, error)
                    self.fitparlist.insert(new_rank, fitpar)
                    self.optparlist.insert(new_rank, optpar)
                    self.optvallist.insert(new_rank, optval)
                    if length >= 21:
                        del self.numberlist[-1]
                        del self.errorlist[-1]
                        del self.fitparlist[-1]
                        del self.optparlist[-1]
                        del self.optvallist[-1]
                else:
                    # Worst error function
                    if length < 21:
                        self.numberlist.append(self.nsteps)
                        self.errorlist.append(error)
                        self.fitparlist.append(fitpar)
                        self.optparlist.append(optpar)
                        self.optvallist.append(optval)
            else:
                # Old minimum found
                if  new_rank > 0 and new_rank <= new_min:
                    # Better error function
                    del self.numberlist[new_min]
                    del self.errorlist[new_min]
                    del self.fitparlist[new_min]
                    del self.optparlist[new_min]
                    del self.optvallist[new_min]
                    self.numberlist.insert(new_rank, self.nsteps)
                    self.errorlist.insert(new_rank, error)
                    self.fitparlist.insert(new_rank, fitpar)
                    self.optparlist.insert(new_rank, optpar)
                    self.optvallist.insert(new_rank, optval)

    def write(self, txt, name="search"):
        if txt is None:
            f = sys.stdout
        else:
            f = open(txt, 'a')
        T = time.localtime()
        f.write("Parameter %s %02i-%02i-%i " % (name, T[2], T[1], T[0]) +
                "at %02i:%02i:%02i\n" % (T[3], T[4], T[5]) +
                "-" * 90 + "\n" +
                "%i optimizations has been made. " % (self.nsteps,) +
                "Below the initial and the 20 best are listed:\n\n")

        for num, error, fitpar, optpar, optval in zip(self.numberlist, self.errorlist,
                                                      self.fitparlist, self.optparlist,
                                                      self.optvallist):
            f.write("Optimization no. %i with error function %.3f:\n" % (num, error))
            for key in fitpar.keys():
                f.write("  %s%s parameters (Initial, Optimal, Variable):\n" % (key[0], key[1]))
                for pi, po, v in zip(fitpar[key], optpar[key], self.varparam[key]):
                    f.write("  %12.9f  %12.9f  %s\n" % (pi, po, v))
            f.write("\n  Fitting values (Quantity, Structure: " +
                    "Target, Optimal, Weight, Deviation):\n")
            devsum = 0
            for ival, oval in zip(self.values, optval):
                if len(oval[2]) > 1:
                    struct = "%s%s-%s" % (oval[2] + (oval[1],))
                else:
                    struct = "%s-%s" % (oval[2] + (oval[1],))
                name = "%s, %s:" % (oval[0], struct)
                dev = np.abs(oval[3] / ival[3] - 1) * 100
                devsum += dev
                f.write("  %-45s  %7.3f  " % (name, ival[3]) +
                        "%10.6f  %5.3f  %7.2f%%\n" % (oval[3], ival[4], dev))
            f.write("%84.2f%%\n\n" % (devsum / len(self.values),))
        f.close()


