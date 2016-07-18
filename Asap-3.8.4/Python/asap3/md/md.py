"""Molecular Dynamics."""

import numpy as np

try:
    from ase.optimize.optimize import Dynamics
except ImportError:
    # Fallback to old placement
    from ase.optimize import Dynamics
from ase.data import atomic_masses
from ase.md import MDLogger

_counter = 0  # Prevents identical prefixes


class MolecularDynamics(Dynamics):
    """Base-class for all MD classes."""
    def __init__(self, atoms, timestep, trajectory, logfile=None,
                 loginterval=1):
        Dynamics.__init__(self, atoms, logfile=None, trajectory=trajectory)
        self.dt = timestep
        # Store data locally except on parallel simulations
        self._localdata = not getattr(atoms, "parallel", False)
        if logfile:
            self.attach(MDLogger(dyn=self, atoms=atoms, logfile=logfile),
                        interval=loginterval)

    def run(self, steps=50):
        """Integrate equation of motion."""
        f = self.atoms.get_forces()

        if not self.atoms.has('momenta'):
            self.atoms.set_momenta(np.zeros_like(f))

        for step in xrange(steps):
            f = self.step(f)
            self.nsteps += 1
            self.call_observers()

    def set(self, name, var):
        """Set a local variable.

        If the local variable is a scalar, it is stored locally.  If
        it is an array it is stored on the atoms.  This allows for
        parallel Asap simulations, where such arrays will have to
        migrate among processors along with all other data for the
        atoms.
        """
        if self._localdata or getattr(var, "shape", ()) == ():
            setattr(self, name, var)
        else:
            if hasattr(self, name):
                delattr(self, name)
            self.atoms.set_array(self.prefix+name, var)

    def get(self, name):
        try:
            return getattr(self, name)
        except AttributeError:
            return self.atoms.get_array(self.prefix+name, copy=False)

class ParallelMolDynMixin:
    def __init__(self, prefix, atoms):
        global _counter
        self.prefix = prefix+str(_counter)+"_"
        _counter += 1
        self._uselocaldata = not getattr(atoms, "parallel", False)
        self._localdata = {}
        
    def set(self, name, var):
        """Set a local variable.

        If the local variable is a scalar, it is stored locally.  If
        it is an array it is stored on the atoms.  This allows for
        parallel Asap simulations, where such arrays will have to
        migrate among processors along with all other data for the
        atoms.
        """
        if self._uselocaldata or getattr(var, "shape", ()) == ():
            self._localdata[name] = var
        else:
            if name in self._localdata:
                del self._localdata[name]
            self.atoms.set_array(self.prefix+name, var)

    def get(self, name):
        try:
            return self._localdata[name]
        except KeyError:
            return self.atoms.get_array(self.prefix+name, copy=False)
