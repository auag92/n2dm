import numpy as np
from asap3.nanoparticle_mc import data
from ase.atom import Atom, names as oldnames
from ase.data import atomic_numbers, chemical_symbols, reference_states
from asap3.Internal.MonteCarloAtoms import MC_Atom
class ClusterAtom(MC_Atom):
    """Cluster Atom"""
    names = oldnames.copy()
    names['neighbors'] = ('neighbors', None)
    names['coordination'] = ('coordinations', 0)
    names['type'] = ('types', 0)

    _data = []

    def __init__(self, symbol='X', position=(0.0, 0.0, 0.0), atoms=None, index=None):
        self.atoms = atoms
        self.index = index

        if atoms is None:
            if isinstance(symbol, str):
                self.number = atomic_numbers[symbol]
            else:
                self.number = symbol
 
            self.position = np.array(position, float)

    def __repr__(self):
        output = 'ClusterAtom(%s, %s' % (self.symbol, self.position.tolist())

        for name in self._data:
            if name != 'number' and name != 'position' and self._get(name) is not None:
                output += ', %s=%s' % (name, self._get(name))

        return output + ')'

    def get_raw(self, name):
        """Get attribute, return None if not explicitely set."""
        if name == 'symbol':
            return chemical_symbols[self.get_raw('number')]

        if self.atoms is None:
            return self.data[name]
        
        plural = self.names[name][0]
        if plural in self.atoms.arrays:
            return self.atoms.arrays[plural][self.index]
        else:
            return None

    def get(self, name):
        """Get attribute, return default if not explicitely set."""
        value = self.get_raw(name)
        if value is None:
            if name == 'mass':
                value = atomic_masses[self.number]
            else:
                value = self.names[name][1]
        return value



    #def _get_copy(self, name):
    #    return self._get(name, copy=True)

                
    def set(self, name, value, copy=False):
        """Set attribute."""
        if name == 'symbol':
            name = 'number'
            value = atomic_numbers[value]

        if self.atoms is None or copy:
            assert name in self.names
            self.data[name] = value
        else:
            plural, default = self.names[name]
            if plural in self.atoms.arrays:
                array = self.atoms.arrays[plural]
                if name == 'magmom' and array.ndim == 2:
                    assert len(value) == 3
                array[self.index] = value
            else:
                if name == 'magmom' and np.asarray(value).ndim == 1:
                    array = np.zeros((len(self.atoms), 3))
                elif name == 'mass':
                    array = self.atoms.get_masses()
                else:
                    default = np.asarray(default)
                    array = np.zeros((len(self.atoms),) + default.shape,
                                     default.dtype)
                array[self.index] = value
                self.atoms.new_array(plural, array)


    def has(self, name):
        return name in self._data

    def cut_reference_to_atoms(self):
        for name, a in self.atoms.arrays.items():
            self.set(self.atoms.names[name][0], a[self.index].copy(), True)
        self.atoms = None
        self.index = None

    def get_symbol(self): return chemical_symbols[self.get('number')]
    def get_neighbors(self): return self.get('neighbors')
    def get_type(self): return self.get('type')
    def get_coordination(self): return self.get('coordination')
    #def get_(self): return self._get('')

    def set_symbol(self, value): self.set('number', atomic_numbers[value])
    def set_neighbors(self, value): self.set('neighbors', np.array(value, int))
    def set_type(self, value): self.set('type', value)
    def set_coordination(self, value): self.set('coordination', value)
    #def get_(self, value): return self._set('', value)

    symbol = property(get_symbol, set_symbol, doc='Chemical symbol')
    neighbors = property(get_neighbors, set_neighbors, doc='List of nearest neighbors')
    type = property(get_type, set_type, doc='Atom type')
    coordination = property(get_coordination, set_coordination, doc='Atom coordination')
    # We need to repeat these as MC_Atom has blocked them.
    # Lambda expression needed as _get_position is in a base class :-(
    position = property(lambda self: self._get_position(),
                        lambda self, x : self.set_position(x),
                        doc='XYZ-coordinates')
    number = property(lambda self: self.get_atomic_number(),
                      lambda self, x: self.set_atomic_number(x),
                      doc='Atomic number')
