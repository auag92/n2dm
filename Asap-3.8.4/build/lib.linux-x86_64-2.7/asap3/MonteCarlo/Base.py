import sys
import time
import numpy as np

from ase.io import PickleTrajectory

from asap3.Tools.Timing import Timer

class GlobalOptimizer:
    def __init__(self, atoms, log='-', traj=None):
        self.atoms = atoms
        self.move = None
        self.moves = []
        self.moves_weight = []
        self.observers = []
        self.nsteps = {}
        self.naccept = {}
        self.timer = Timer()

        self.log_init = ''
        self.log_header = ''
        self.log_string = ''
        self.log_final = ''

        if isinstance(traj, str):
            traj = PickleTrajectory(traj, 'w', atoms, master=True)
            traj.write()
            self.attach_observer(traj.write)

        if isinstance(log, str):
            if log == '-':
                self.logfile = sys.stdout
            else:
                self.logfile = open(log, 'a')
        else:
            self.logfile = None

    def attach_move(self, move, weight=1.0):
        if not hasattr(move, '__call__'):
            raise ValueError("Attached move is not callable.")

        if hasattr(move, 'set_atoms'):
            move.set_atoms(self.atoms)

        if hasattr(move, 'set_optimizer'):
            move.set_optimizer(self)

        self.moves.append(move)
        self.moves_weight.append(weight)

    def random_move(self):
        return self.moves[(np.random.uniform() < self.moves_weight).argmax()]

    def attach_observer(self, observer, interval=1, attime='accept', *args, **kwargs):
        if not hasattr(observer, '__call__'):
            raise ValueError("Attached observer is not callable.")

        if hasattr(observer, 'set_atoms'):
            observer.set_atoms(self.atoms)

        if hasattr(observer, 'set_optimizer'):
            observer.set_optimizer(self)

        self.observers.append((observer, interval, attime, args, kwargs))

    def call_observers(self, thistime):
        for observer, interval, attime, args, kwargs in self.observers:
            if self.get_nsteps() % interval == 0 and attime == thistime:
                observer(*args, **kwargs)

    def initialize(self):
        if not len(self.moves) > 0:
            raise ValueError("You must attach at least one move.")

        w = np.array(self.moves_weight, dtype=float)
        self.moves_weight = np.add.accumulate(w/w.sum())

        for m in self.moves:
            self.nsteps[m.get_name()] = 0
            self.naccept[m.get_name()] = 0

    def finalize(self):
        pass

    def run(self, steps, *args, **kwargs):
        self.timer.start('Initialize')
        self.initialize(*args, **kwargs)
        self.timer.stop('Initialize')

        if self.logfile is not None:
            name = self.__class__.__name__
            T = time.localtime()
            self.logfile.write("\n\n%s global optimization, " % (name,) +
                               "%02i-%02i-%i\n" % (T[2], T[1], T[0]) +
                               "*" * 75 + "\n" +
                               "%-15s%-7s\n" % ("Move", "Weight"))

            for m, w in zip(self.moves, self.moves_weight):
                self.logfile.write("%-15s%-7.5f\n" % (m.get_name(), w))

            self.logfile.write(self.log_init)
            self.logfile.write('\n%-8s  %-7s  %-11s  ' %
                               ('Time', 'Steps', 'Move') + self.log_header)
            self.logfile.flush()

        for i in range(steps):
            self.timer.start('Move')
            self.move = self.random_move()
            self.move(self.atoms)
            self.timer.stop('Move')

            self.call_observers('move')

            self.timer.start('Evaluate')
            accept = self.evaluate_move()
            self.timer.stop('Evaluate')

            self.call_observers('evaluate')

            if accept:
                self.accept_move()
            else:
                self.reject_move()

            self.call_observers('after')

            if self.logfile is not None:
                name = self.move.get_name()
                T = time.localtime()
                self.logfile.write('%02d:%02d:%02d  %-7i  %-11s  ' %
                                   (T[3], T[4], T[5], self.get_nsteps(),
                                    name[:-4]) + self.log_string)
                self.logfile.flush()

        self.finalize()

        if self.logfile is not None:
            self.logfile.write('\n%-13s  %-15s  %-15s\n' %
                               ('Move', 'Steps', 'Accepts'))
            self.logfile.write('=' * 60 + '\n')
            for m in self.moves:
                name = m.get_name()
                ns = self.get_nsteps(name)
                fs = 1.0 * ns / self.get_nsteps()
                na = self.get_naccept(name)
                fa = 1.0 * na / ns
                self.logfile.write('%-13s  %-7i (%5.3f)  %-7i (%5.3f)\n' %
                                   (name, ns, fs, na, fa))
            self.logfile.write('=' * 60 + '\n')
            ns = self.get_nsteps()
            na = self.get_naccept()
            self.logfile.write('%-13s  %-7i (%5.3f)  %-7i (%5.3f)\n' %
                               ('Total', ns, 1.0, na, 1.0 * na / ns))
            self.logfile.write('=' * 60 + '\n')
            self.logfile.write(self.log_final)
            self.timer.write(self.logfile)
            self.logfile.flush()

    def evaluate_move(self):
        return True

    def accept_move(self):
        self.move.accept()
        self.call_observers('accept')
        self.naccept[self.move.get_name()] += 1
        self.nsteps[self.move.get_name()] += 1

    def reject_move(self):
        self.move.reject()
        self.call_observers('reject')
        self.nsteps[self.move.get_name()] += 1

    def get_nsteps(self, name='total'):
        if name == 'total':
            return sum(self.nsteps.values())
        else:
            return self.nsteps[name]

    def get_naccept(self, name='total'):
        if name == 'total':
            return sum(self.naccept.values())
        else:
            return self.naccept[name]

