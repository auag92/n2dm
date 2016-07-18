from ase import Atom, Atoms
import numpy as np

_mc_err = "Cannot access the position(s) and atomic_number(s) attribute of MC optimized atoms."

class MC_Atom(Atom):
    """A single atom in a Monte Carlo optimized Atoms object."""
    #__slots__ = Atom.__slots__
    
    def get_position(self): 
        p = self.get('position')
        p.flags.writeable = False
        return p

    def set_position(self, position):
        Atom.set(self, 'position', position)
        self._mark()

    def get_atomic_number(self):
        return self.get('number')
    
    def set_atomic_number(self, number):
        Atom.set(self, 'number', number)
        self._mark()
        
    def _mark(self):
        if self.atoms is not None:
            mco = self.atoms.mc_optim
            mco[0] = n = mco[0] + 1
            if n < 100:
                mco[n] = self.index

    def _oops(self):
        raise NotImplementedError(_mc_err)

    number = property(get_atomic_number, set_atomic_number, doc='The atomic number')
    position = property(get_position, set_position, doc='The position (gettting it this way is a read-only array')    

class MonteCarloAtoms(Atoms):
    """A list of atoms class optimized for Monte Carlo simulations.

    This is intended to be used together with the MonteCarloEMT calculator.
    """
    def __init__(self, *args, **kwargs):
        self.mc_optim = np.zeros(101, np.intc)
        Atoms.__init__(self, *args, **kwargs)
        
    def __getitem__(self, i):
        if isinstance(i, int):
            natoms = len(self)
            if i < -natoms or i >= natoms:
                raise IndexError('Index out of range.')

            return MC_Atom(atoms=self, index=i)
        else:
            return Atoms.__getitem__(self, i)
    
    def set_calculator(self, calc=None):
        """Attach calculator object."""
        Atoms.set_calculator(self, calc)
        self.mc_optim[0] = 10000  # MC optim invalid

    def get_potential_energy(self):
        epot = Atoms.get_potential_energy(self)
        self.mc_optim[0] = 0
        return epot

    def get_potential_energies(self):
        epot = Atoms.get_potential_energies(self)
        self.mc_optim[0] = 0
        return epot

    def set_positions(self, newpositions):
        """Set positions."""
        Atoms.set_positions(self, newpositions)
        self.mc_optim[0] = 10000  # MC optim invalid

    def get_test(self):
        return True
        
