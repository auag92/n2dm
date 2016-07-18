from ase.constraints import FixAtoms as ASE_FixAtoms
from ase.constraints import Filter as ASE_Filter
import numpy as np

class ConstraintMixin:
    """Mixin class FixAtoms and Filter"""
         
    def prepare_for_asap(self, atoms):
        """Prepare this constraint for optimized Asap dynamics
        
        This function must be called once the atoms are known by
        all dynamics supporting parallel MD.
        """
        # Store the arrays of the atoms, not the atoms, to prevent cyclic references.
        self.atoms_arrays = atoms.arrays
        assert self.indexname not in self.atoms_arrays
        idx = self.index
        self.asap_ready = True
        self.index = idx
        del self._index
        assert self.indexname in self.atoms_arrays
        
    def set_index(self, idx):
        if self.asap_ready:
            natoms = len(self.atoms_arrays['positions'])
            if idx.dtype == bool:
                # Boolean - must be a mask
                assert len(idx) == natoms
            else:
                # Must be a list of indices.  Convert to a mask
                idx2 = np.zeros(natoms, bool)
                idx2[idx] = True
                idx = idx2
            self.atoms_arrays[self.indexname] = idx
        else:
            self._index = idx
            
    def get_index(self):
        if self.asap_ready:
            return self.atoms_arrays[self.indexname]
        else:
            return self._index
        
class FixAtoms(object, ConstraintMixin, ASE_FixAtoms):
    def __init__(self, indices=None, mask=None):
        self.pre_init()
        ASE_FixAtoms.__init__(self, indices, mask)
        
    def pre_init(self):
        self.indexname = "FixAtoms_index"
        self.asap_ready = False
        
    def copy(self):
        if self.index.dtype == bool:
            return FixAtoms(mask=self.index.copy())
        else:
            return FixAtoms(indices=self.index.copy())
        
    def __get_state__(self):
        return {'data': self.index,
                'version': 1}
        
    def __set_state__(self, state):
        try:
            assert(state['version'] == 1)
        except KeyError:
            print state
            raise
        self.pre_init()
        self.index = state['data']
        
    index = property(ConstraintMixin.get_index, ConstraintMixin.set_index)
    
    
class Filter(object, ConstraintMixin, ASE_Filter):
    def __init__(self, atoms, indices=None, mask=None):
        """Filter atoms.

        This filter can be used to hide degrees of freedom in an Atoms
        object.

        Parameters
        ----------
        indices : list of int
           Indices for those atoms that should remain visible.
        mask : list of bool
           One boolean per atom indicating if the atom should remain
           visible or not.
        """

        self.pre_init()
        ASE_Filter.__init__(self, atoms, indices, mask)
        self.prepare_for_asap(atoms)
        
    def pre_init(self):
        self.indexname = "FixAtoms_index"
        self.asap_ready = False
        
    index = property(ConstraintMixin.get_index, ConstraintMixin.set_index)  
    
def check_asap_constraints(atoms, allowed=None):
    """Check that atoms only have allowed constraints.  Return True if so, otherwise False.
    
    An optional second parameter can be a tuple of allowed constraints.
    """
    if allowed is None:
        allowed = (FixAtoms,)
        
    if len(atoms.constraints) == 0:
        return True
    if len(atoms.constraints) > 1:
        return False
    return isinstance(atoms.constraints[0], allowed)

