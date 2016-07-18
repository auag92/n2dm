# Copied from GPAW (https://trac.fysik.dtu.dk/projects/gpaw/)

import time
import math
import sys

from ase.parallel import world

def function_timer(func, *args, **kwargs):
    out = kwargs.pop('timeout', sys.stdout)
    t1 = time.time()
    r = func(*args, **kwargs)
    t2 = time.time()
    print >> out, t2 - t1
    return r

class Timer:
    def __init__(self, print_levels=1000):
        self.timers = {}
        self.t0 = time.time()
        self.running = []
        self.print_levels = print_levels
        
    def start(self, name):
        names = tuple(self.running + [name])
        self.timers[names] = self.timers.get(names, 0.0) - time.time()
        self.running.append(name)
        
    def stop(self, name=None):
        if name is None: name = self.running[-1]
        names = tuple(self.running)
        running = self.running.pop()
        if name != running:
            raise RuntimeError('Must stop timers by stack order.  '
                               'Requested stopping of %s but topmost is %s'
                               % (name, running))
        self.timers[names] += time.time()
            
    def get_time(self, *names):
#        print self.timers, names
        return self.timers[names]
                
    def write(self, out=sys.stdout):
        while self.running:
            self.stop()
        if len(self.timers) == 0:
            return

        t0 = time.time()
        tot = t0 - self.t0

        n = max([len(names[-1]) + len(names) for names in self.timers]) + 1
        out.write('\n%-*s    incl.     excl.\n' % (n, 'Timing:'))
        out.write('%s\n' % ('=' * 60))
        tother = tot
        
        inclusive = self.timers.copy()

        exclusive = self.timers
        keys = exclusive.keys()
        keys.sort()
        for names in keys:
            t = exclusive[names]
            if len(names) > 1:
                if len(names) < self.print_levels + 1:
                    exclusive[names[:-1]] -= t
            else:
                tother -= t
        exclusive[('Other',)] = tother
        inclusive[('Other',)] = tother
        keys.append(('Other',))
        for names in keys:
            t = exclusive[names]
            tinclusive = inclusive[names]
            r = t / tot
            p = 100 * r
            i = int(40 * r + 0.5)
            if i == 0:
                bar = '|'
            else:
                bar = '|%s|' % ('-' * (i - 1))
            level = len(names)
            if level > self.print_levels:
                continue
            name = (level - 1) * ' ' + names[-1] + ':'
            out.write('%-*s%9.3f %9.3f %5.1f%% %s\n' %
                      (n, name, tinclusive, t, p, bar))
        out.write('%s\n' % ('=' * 60))
        out.write('%-*s%9.3f %5.1f%%\n' % (n + 10, 'Total:', tot, 100.0))
        out.write('%s\n' % ('=' * 60))
        out.write('date: %s\n' % time.asctime())

    def add(self, timer):
        for name, t in timer.timers.items():
            self.timers[name] = self.timers.get(name, 0.0) + t


class DebugTimer(Timer):
    def __init__(self, print_levels=1000, comm=world, txt=sys.stdout):
        Timer.__init__(self, print_levels)
        ndigits = 1 + int(math.log10(comm.size))
        self.srank = '%0*d' % (ndigits, comm.rank)
        self.txt = txt

    def start(self, name):
        Timer.start(self, name)
        t = self.timers[tuple(self.running)] + time.time()
        self.txt.write('T%s >> %s (%7.5fs) started\n' % (self.srank, name, t))

    def stop(self, name=None):
        if name is None: name = self.running[-1]
        t = self.timers[tuple(self.running)] + time.time()
        self.txt.write('T%s << %s (%7.5fs) stopped\n' % (self.srank, name, t))
        Timer.stop(self, name)


class StepTimer(Timer):
    """Step timer to print out timing used in computation steps.
    
    Use it like this::

      from gpaw.utilities.timing import StepTimer
      st = StepTimer()
      ...
      st.write_now('step 1')
      ...
      st.write_now('step 2')

    The parameter write_as_master_only can be used to force the timer to
    print from processess that are not the mpi master process.
    """
    
    def __init__(self, out=sys.stdout, name=None, write_as_master_only=True):
        Timer.__init__(self)
        if name is None:
            name = '<%s>' % sys._getframe(1).f_code.co_name
        self.name = name
        self.out = out
        self.alwaysprint = not write_as_master_only
        self.now = 'temporary now'
        self.start(self.now)


    def write_now(self, mark=''):
        self.stop(self.now)
        if self.alwaysprint or world.rank == MASTER:
            print >> self.out, self.name, mark, self.gettime(self.now)
        self.out.flush()
        del self.timers[self.now]
        self.start(self.now)


