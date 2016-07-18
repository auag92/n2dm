"""An observer keeping track of running average positions of atoms.

"""

from asap3.Internal.Subject import Subject

class RunningAverage:
    """Calculate running average positions of atoms.
    
    Attach this observer to the atoms to make average positions available
    e.g. for CNA analysis.  The filter is called often (maybe every 5 time
    steps).  After it has been called interval times (or a multiple hereof)
    the average positions are calculated, and the attached analysis 
    functions are called.  numavg specifies over how many calls are used
    to calculate the average positions.  If numavg==interval, all calls
    are used to collect data, and the average is over all times since the
    last analysis.  If numavg<interval then only the last numavg calls
    before the analysis is used to collect data.  Currently, numavg>interval
    is not supported.
    
    During analysis, positions of the atoms are set to the average positions.
    After analysis, observers can be called e.g. to store the results.  
    Observers may be called before or after the positions are restored to
    their instantaneous values.
    
    Warning: As the atoms are changed before and after analysis, migration
    may be triggered in parallel simulations.

    Parameters:

    atoms: 
        The atoms being observed.  Temporary data is stored in the atoms.
    
    interval:
        How often analysis is done.
    
    numavg:
        How many frames should be used for the averaging.
    
    Important methods:

    update(): 
        Update the running average, and do analysis if needed.  This
        is the function that should be attached to the dynamics.

    __call__(): 
        Same as update()
    
    attach_analysis(callable, *args, **kwargs): 
        Attach an analysis function.  Arguments are passed on.
        
    attach_observer(callable, interval=1, *args, **kwargs):
        Call an external observer for example to store the result
        of the analysis.  The observer is called every interval'th
        time the analysis is perfomed, it is rare to set interval
        different from 1.  The observer is called after the 
        positions have been restored to their instantaneous
        values.
        
    attach_observer_avgpos(callable, interval=1, *args, **kwargs):
        As attach_observer(), except that the observer is called while
        the positions still have their average value.
    """
    
    def __init__(self, atoms, interval, numavg):
        self.atoms = atoms
        self.interval = interval
        self.numavg = numavg
        if numavg > interval:
            raise NotImplementedError("The case numavg>interval has not yet been implemented.")
        self.firstdata = self.interval - self.numavg
        self.analysers = Subject()
        self.observers = Subject()
        self.avgobservers = Subject()
        self.extra_data_buffer = []
        self.n = 0
        self.m = 0
        
    def attach_analysis(self, func, *args, **kwargs):
        "Attach an analysis method"
        self.analysers.attach(func, 1, *args, **kwargs)
        
    def attach_observer(self, func, interval=1, *args, **kwargs):
        "Attach an observer called after analysis"
        if self.extra_data_buffer and hasattr(func, 'set_extra_data'):
            for name, source, once in self.extra_data_buffer:
                func.set_extra_data(name, source, once)
        self.observers.attach(func, interval, *args, **kwargs)
        
    def attach_observer_avgpos(self, func, interval=1, *args, **kwargs):
        "Attach an observer called after analysis but while positions remain averaged."
        self.avgobservers.attach(func, interval, *args, **kwargs)
        
    def reset(self):
        "Reset internal variables after an analysis"
        self.n = 0
        self.m = 0
        if self.atoms.has('avgpositions'):
            self.atoms.set_array('avgpositions', None)
        if self.atoms.has('sumpositions'):
            # This should not normally happen.
            self.atoms.set_array('sumpositions', None)
        
    def update(self):
        self.n += 1
        if self.n <= self.firstdata:
            return  # Do nothing
        # Take data, and include in average
        if self.atoms.has('sumpositions'):
            sp = self.atoms.get_array('sumpositions') + self.atoms.get_positions()
        else:
            sp = self.atoms.get_positions()
        self.atoms.set_array('sumpositions', sp)
        self.m += 1
        if self.n == self.interval:
            assert self.m == self.numavg
            self.atoms.set_array('avgpositions', sp * (1.0/self.m))
            self.atoms.set_array('sumpositions', None)
            # Now save the instantanous positions and set the positions to
            # their average
            self.atoms.set_array('instant_positions', self.atoms.get_positions())
            self.atoms.set_positions(self.atoms.get_array('avgpositions'))
            # Call analysis and then observers (avg pos)
            self.analysers.call_observers()
            self.avgobservers.call_observers()
            # Reset positions to instantaneous value
            self.atoms.set_positions(self.atoms.get_array('instant_positions'))
            self.atoms.set_array('instant_positions', None)
            self.observers.call_observers()
            self.reset()
            
    __call__ = update
        
    def set_extra_data(self, name, source=None, once=False):
        """Forward to observers that extra data should be stored.

        This function is intended to make sure that if a BundleTrajectory
        is attached to an NPT dynamics through a RunningAverage object,
        then information about the dynamics is still saved in the 
        BundleTrajectory.

        It works by forwaring the call to set_extra_data to any observers
        having a set_extra_data method - whether they are attached before
        or after calling this function.
        """
        self.extra_data_buffer.append((name, source, once))
        for obs in self.observers.observers:
            f = obs[0]
            if hasattr(f, "set_extra_data"):
                f.set_extra_data(name, source, once)
