"""Time averages of the positions of a list of atoms.

The `TimeAveragedPositions` object is a Filter giving the averaged
positions of the atoms in a list of atoms.  To do this, it acts as an
Observer of the dynamics object moving the atoms, and it can itself be
observed by e.g. a `RestrictedCNA` filter or a Trajectory.
"""

__docformat__ = "restructuredtext en"

from ASE.Utilities.Subject import Subject
import Numeric

class TimeAveragedPositions(Subject):
    """Time averages of the positions of a list of atoms.

    The TimeAveragedPositions object is a Filter giving the averaged
    positions of the atoms in a list of atoms.  To do this, it acts as
    an Observer of the dynamics object moving the atoms, and it can
    itself be observed by e.g. a CNA filter or a Trajectory.  It must
    be attached to a dynamics to work properly.

    Created with the following arguments:
    
    atoms
        The list of atoms being studied.

    avgtime
        The number of steps over which the average should be calculated.

    avgfrequency
        How often the average should be calculated.  Must be >= avgtime

    modifyAtoms (default False)
        If true, the positions of the original atoms
        are modified to reflect the averaged positions immediately before
        any observers are called, and changed back immediately after.
        Probably only useful for parallel simulations.

    verbose (default False)
        Print debugging info.
    
    """
    def __init__(self, atoms, avgtime, avgfrequency, modifyAtoms=False,
                 verbose=False):
        Subject.__init__(self)
        self.atoms = atoms
        self.avgtime = avgtime
        self.frequency = avgfrequency
        self.modifyAtoms = modifyAtoms
        self.verbose = verbose
        if avgfrequency < avgtime:
            raise NotImplementedError, "Cannot yet handle observations more frequent than the averaging time."
        self.counter = 0
        self.n = 0
        self.avgpos = Numeric.zeros((len(atoms),3), Numeric.Float)
        self.ready = 0
        # Handle parallel simulations
        try:
            self.avgpos = atoms.RegisterArray(self.avgpos)
        except AttributeError:
            self.register = 0
        else:
            self.register = 1

    def __del__(self):
        if self.register:
            self.atoms.DeregisterArray(self.avgpos)

    def Update(self):
        """Called by the dynamics every time the atoms have changed."""
        self.counter += 1
        if self.ready:
            self.n = 0
            self.ready = 0
        # Only do something if it is time to accumulate an average
        if self.counter > self.frequency - self.avgtime:
            newpos = self.atoms.GetScaledPositions()
            self.n += 1
            if self.n == 1:
                if self.verbose:
                    print "Initializing average."
                self.avgpos[:] = newpos
            else:
                if self.verbose:
                    print "Accumulating average."
                self.avgpos += newpos
        if self.counter == self.frequency:
            if self.verbose:
                print "Calculating average (%d steps)." % self.n
            self.ready = 1
            self.avgpos /= self.n
            if self.modifyAtoms:
                if self.register:
                    newpos = self.atoms.RegisterArray(newpos)
                self.atoms.SetScaledPositions(self.avgpos)
                self.atoms.CheckNeighborLists()
            self.Notify()
            if self.modifyAtoms:
                self.atoms.SetScaledPositions(newpos)
                self.atoms.CheckNeighborLists()
                if self.register:
                    self.atoms.DeregisterArray(newpos)
            self.counter = 0

    def ForceNotify(self):
        """Force an average to be calculated, and notify observers.

        Mostly useful to force an output at the beginning of a calculation.
        In that case, no averaging is done, but the rest of the filter chain
        is activated.
        """
        self.counter = self.frequency - 1
        self.Update()

    def GetScaledPositions(self):
        """Get the time-averaged positions, if they are ready."""
        if self.verbose:
            print "Accessing average."
        if self.ready:
            return self.avgpos
        else:
            raise RuntimeError, "Averaged positions are not ready."

    def GetCartesianPositions(self):
        """Get the time-averaged positions, if they are ready."""
        return Numeric.matrixmultiply(self.GetScaledPositions(),
                                      self.atoms.GetUnitCell())

    def GetGeneralizedPositions(self):
        """Get the time-averaged positions as generalized coordinates, if they are ready."""
        return self.GetCartesianPositions().flat

    def SetCartesianPositions(self, pos):
        raise AttributeError, "It makes no sense to set the time-averaged positions."

    def SetGeneralizedPositions(self, pos):
        raise AttributeError, "It makes no sense to set the time-averaged positions."

    # All other calls are passed onto the atoms.
    # Is this a good idea, or should there be a "whitelist" ??    XXXXX

    def __getattr__(self, attr):
        return getattr(self.atoms, attr)

