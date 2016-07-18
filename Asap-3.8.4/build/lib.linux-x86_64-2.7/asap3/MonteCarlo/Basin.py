import numpy as np

from asap3.MonteCarlo.Base import GlobalOptimizer
from ase.units import kB
from ase.optimize import LBFGS

class BasinHopping(GlobalOptimizer):
    def __init__(self, atoms, log='-', traj=None, localopt=LBFGS, fmax=1e-2):
        GlobalOptimizer.__init__(self, atoms, log, traj)
        self.localoptimizer = localopt
        self.fmax = fmax
        self.Eo = atoms.get_potential_energy()
        self.En = 1.0e32
        self.P = 0.0

        self.log_header = ('%-15s  %-12s  %-6s  %-7s\n' %
                           ('Energy', 'Probability', 'Accept', 'OptStep'))

    def initialize(self, temp=1000.0):
        GlobalOptimizer.initialize(self)
        self.temp = temp
        self.log_info = '\nTemperature: %i\n' % (temp,)

    def evaluate_move(self):
        dyn = self.localoptimizer(self.atoms, logfile=None)
        dyn.run(fmax=self.fmax, steps=1000)

        self.En = self.atoms.get_potential_energy()
        self.P = min(1.0, np.exp((self.Eo - self.En) / (self.temp * kB)))
        accept = self.P > np.random.uniform()

        self.log_string = ('%-15.6f  %-12.6f  %-6s  %-7i\n' %
                           (self.En, self.P, accept, dyn.nsteps))

        return accept

    def accept_move(self):
        self.Eo = self.En
        GlobalOptimizer.accept_move(self)

