import sys
import time
import numpy as np

from ase.parallel import world

from asap3.Tools.ParameterOptimization.Optimization import ParameterOptimization

def now():
    "Time as a string."
    return time.strftime("%H:%M:%S")

class ParameterSearch(ParameterOptimization):
    """Parallel parameter optimization."""
    def __init__(self, elements, calc_class, initparam, varparam, quantities,
                 latticeconstants, calc_args=None, debug=False):
        ParameterOptimization.__init__(self, elements, calc_class, initparam,
                                       varparam, quantities, latticeconstants,
                                       calc_args, debug)

    def random_parameters(self, parameters, std):
        change = np.random.normal(1.0, std, len(parameters))
        return (np.array(parameters) * change).tolist()

    def run(self, std, temp, maxstep=100, xtol=1e-2, ftol=5e-2, delta=5e-3,
            ptol=1e-2, dat=None, log=None, err=None):
        """Run the optimization.
        
        First parameter: Gauss width for MC part
        Second parameter: MC Temperature
        Third parameter: run time, in h or m.  Or the number of steps if integer.
        xtol and ftol are simplex convergence criterial.  
        delta is the initial simplex size, relatively.
        """
        self.ptol = ptol
        self.nsteps = 0
        self.numberlist = []
        self.errorlist = []
        self.countlist = []
        self.fitparlist = []
        self.optparlist = []
        self.optvallist = []

        self.maxstep = None
        self.endtime = None
        self.done(maxstep)

        if log != None:
            s = log.rfind('.')
            log = open(log[:s] + '-' + str(world.rank) + log[s:], 'a')

        if err != None:
            s = err.rfind('.')
            err = open(err[:s] + '-' + str(world.rank) + err[s:], 'w')

        if world.rank == 0:
            done = np.zeros(world.size)
            request_s = [0] * world.size
            request_r = [0] * world.size
            result = [0] * world.size
        else:
            done = False
            request_s = None
            request_r = None
            result = np.empty(1)

        accepts = 0
        old_error = 1e6
        old_optpar = None
        while not self.done():
            # Make a random set of parameters close to the
            if old_optpar == None:
                initparam = self.parameters_dict2list()
                if world.rank == 0:
                    fitpar = initparam
                else:
                    fitpar = self.random_parameters(initparam, std * 2.0)
            else:
                fitpar = self.random_parameters(old_optpar, std)

            # Run optimization
            (error, optpar, optval) = self.fit(fitpar, xtol, ftol, delta,
                                               log, err, True)

            self.nsteps += 1
            if error < old_error or temp < 1e-6:
                prop = 1.0
            else:
                derr = (error - old_error) / old_error
                prop = np.exp(-derr/temp)

            if prop > np.random.uniform():
                accepts += 1
                old_error = error
                old_optpar = optpar[:]
                self.save_opt(error, fitpar, optpar, optval)
         
        # Gather results and Write output
        sys.stderr.write("Task %i [%s]: Finished optimizing, writing results.\n"
                         % (world.rank, now()))
        self.write_local(dat, std, temp, accepts)
        if dat != None:
            sys.stderr.write("Task %i [%s]: Communicating.\n" % (world.rank, now()))
            a = np.zeros(21)
            alen = min(len(self.errorlist), 21)
            a[:alen] = np.array(self.errorlist)[:alen]
            if world.rank == 0:
                errors = np.zeros(world.size * len(a))
                tasks = np.arange(world.size * len(a)) / 21
                if world.size > 1:
                    world.gather(a, 0, errors)
                else:
                    errors = a
                tasks = tasks[errors > 1e-20]
                errors = errors[errors > 1e-20]
                tasks = tasks[errors.argsort()]
                errors = errors[errors.argsort()]
                f = open(dat, 'a')
                f.write('\nParameter search at ' + time.asctime() + '\n')
                f.write('--------------------------------------------\n')
                f.write('Rank  Task  Error function\n')
                for i, (n, e) in enumerate(zip(tasks, errors)):
                    f.write('%-4i  %-4i  %.5e\n' % (i, n, e))
                f.close()
            else:
                world.gather(a, 0)
        sys.stderr.write("Task %i [%s]: Done.\n" % (world.rank, now()))


    def done(self, maxstep=None, hard=False):
        if maxstep != None:
            if isinstance(maxstep, int):
                self.maxstep = maxstep
            elif isinstance(maxstep, str):
                multiplier = {'h': 3600, 'm': 60}
                elapse = float(maxstep[:-1]) * multiplier[maxstep[-1]]
                self.endtime = time.time() + elapse
            else:
                raise ValueError('Cannot understand the given maxstep.')
        else:
            if self.maxstep == None:
                return time.time() >= self.endtime
            else:
                if hard:
                    return self.nsteps > self.maxstep * 1.2
                else:
                    return self.nsteps >= self.maxstep

    def write_local(self, dat, std, temp, accepts):
        s = dat.rfind('.')
        dat = open(dat[:s] + '-' + str(world.rank) + dat[s:], 'a')

        T = time.localtime()
        dat.write("Parameter search %02i-%02i-%i " % (T[2], T[1], T[0]) +
                  "at %02i:%02i:%02i\n" % (T[3], T[4], T[5]) +
                  "-" * 90 + "\n" +
                  "%i optimizations has been made and " % (self.nsteps,) +
                  "%i have been accepted.\n" % (accepts,) +
                  "Below the initial and the 20 best are listed.\n\n" +
                  "Gauss width: %.3f\nTemperature: %.3f\n\n" % (std, temp))

        for num, err, count, fitpar, optpar, optval in zip(self.numberlist,
                                                           self.errorlist,
                                                           self.countlist,
                                                           self.fitparlist,
                                                           self.optparlist,
                                                           self.optvallist):
            dat.write("Optimization no. %i with " % (num,) +
                      "error function %.5e (%i):\n" % (err, count))
            self.write_result(fitpar, optpar, optval, dat)
        dat.close()

    def save_opt(self, error, fitpar, optpar, optval):
        # Sort and save the 20 best results
        length = len(self.errorlist)
        if length == 0:
            self.numberlist.append(self.nsteps)
            self.errorlist.append(error)
            self.countlist.append(1)
            self.fitparlist.append(fitpar)
            self.optparlist.append(optpar)
            self.optvallist.append(optval)
        else:
            new_min = -1
            new_rank = -1
            for i, err in enumerate(self.errorlist):
                op1 = np.array(optpar)
                op2 = np.array(self.optparlist[i])

                #ov1 = np.array(optval)
                #ov2 = np.array(self.optvallist[i])

                #ov_dev_max = np.max(2 * np.abs(ov1 - ov2) / np.abs(ov1 + ov2))
                op_dev_max = np.max(2 * np.abs(op1 - op2) / np.abs(op1 + op2))

                if op_dev_max < self.ptol and new_min < 0 and i > 0: #ov_dev_max < 2e-2 and
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
                    self.countlist.insert(new_rank, 1)
                    self.fitparlist.insert(new_rank, fitpar)
                    self.optparlist.insert(new_rank, optpar)
                    self.optvallist.insert(new_rank, optval)
                    if length >= 21:
                        del self.numberlist[-1]
                        del self.errorlist[-1]
                        del self.countlist[-1]
                        del self.fitparlist[-1]
                        del self.optparlist[-1]
                        del self.optvallist[-1]
                else:
                    # Worst error function
                    if length < 21:
                        self.numberlist.append(self.nsteps)
                        self.errorlist.append(error)
                        self.countlist.append(1)
                        self.fitparlist.append(fitpar)
                        self.optparlist.append(optpar)
                        self.optvallist.append(optval)
            else:
                # Old minimum found
                if  new_rank > 0 and new_rank <= new_min:
                    # Better error function
                    del self.numberlist[new_min]
                    del self.errorlist[new_min]
                    count = self.countlist.pop(new_min)
                    del self.fitparlist[new_min]
                    del self.optparlist[new_min]
                    del self.optvallist[new_min]
                    self.numberlist.insert(new_rank, self.nsteps)
                    self.errorlist.insert(new_rank, error)
                    self.countlist.insert(new_rank, count + 1)
                    self.fitparlist.insert(new_rank, fitpar)
                    self.optparlist.insert(new_rank, optpar)
                    self.optvallist.insert(new_rank, optval)
                else:
                    self.countlist[new_min] += 1


if __name__ == '__main__':
    from asap3.Tools.ParameterOptimization.EMT import EMT2011Fit

    # Parameter names: ['eta2', 'lambda', 'kappa', 'E0', 'V0', 'S0', 'n0']
    initparam = {('Pt','Pt'): [3.4242, 4.1423, 5.9432, -5.85, 4.067, 1.5346, 0.05412]}
    varparam = {('Pt','Pt'): [True, True, True, True, True, True, False]}

    quantities = [('lattice_constant_a', 'fcc', 'Pt', 3.92,      0.001),
                  #('lattice_constant_a', 'hcp', 'Pt', 2.77,      0.05),
                  #('lattice_ratio_ca', 'hcp', 'Pt', 4.78 / 2.77, 0.05),
                  ('bulk_modulus', 'fcc', 'Pt', 278.3,           0.01),
                  #('elastic_constant_C11', 'fcc', 'Pt', 346.7,   0.1),
                  #('elastic_constant_C12', 'fcc', 'Pt', 250.7,  0.1),
                  #('elastic_constant_C44', 'fcc', 'Pt', 76.5,    0.02),
                  ('cohesive_energy', 'fcc', 'Pt', 5.84,         0.01),
                  #('phase_energy', 'fcc-hcp', 'Pt', -0.05,       0.02),
                  #('surface_energy', 'fcc111', 'Pt', 0.631,      0.02),
                  #('surface_ratio', 'fcc111-fcc100', 'Pt', 0.631 / 0.892, 0.01),
                  ]

    latticeconstants = [['fcc', 'Pt', 3.9],
                        ['hcp', 'Pt', (2.8, 4.8)]]

    opt = ParameterSearch(['Pt'], EMT2011Fit, initparam, varparam, quantities,
                          latticeconstants, ('kappa',), True)
    opt.run(0.01, 0.0, '0.1m', dat='SearchParallel.dat', log='SearchParallel.log',
            err='SearchParallel.err')


