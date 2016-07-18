"""Asap module Collector.

Defines a filter collecting all information about parallel atoms on the master.
"""

import asap3
import numpy as np

class Collector:
    "Atoms-like filter collecting information on the master node."

    def __init__(self, atoms, master=None):
        self.atoms = atoms
        self.comm = atoms.get_comm()
        if master is None:
            master = (self.comm.rank == 0)
        self.master = master
        self.constraints = self.atoms.constraints

    def __len__(self):
        n = self.atoms.get_number_of_atoms()
        if self.master:
            return n 
        else:
            return 0

    def get_number_of_atoms(self):
        return self.atoms.get_number_of_atoms()

    def get_positions(self):
        return self.collect(self.atoms.get_positions)

    def get_forces(self):
        return self.collect(self.atoms.get_forces)

    def get_momenta(self):
        return self.collect(self.atoms.get_momenta)

    def get_atomic_numbers(self):
        return self.collect(self.atoms.get_atomic_numbers)
    
    def get_tags(self):
        return self.collect(self.atoms.get_tags)

    def get_potential_energy(self):
        return self.atoms.get_potential_energy()

    def get_cell(self):
        return self.atoms.get_cell()

    def get_calculator(self):
        return self.atoms.get_calculator()

    def get_stress(self):
        return self.atoms.get_stress()

    def get_pbc(self):
        return self.atoms.get_pbc()
    
    def get_info(self):
        return self.atoms.info
    
    def get_charges(self):
        raise NotImplementedError
    
    def get_array(self, label):
        return self.collect(lambda a=self.atoms, l=label: a.get_array(l))

    def has(self, name):
        """Check for existance of array.

        name must be one of: 'tags', 'momenta', 'masses', 'magmoms',
        'charges'.
        """
        if name in ['positions', 'tags', 'momenta', 'numbers']:
            return self.atoms.has(name) 
        else:
            return False
    
    def collect(self, method):
        "Collect data from all cpus onto the master."
        ids = self.atoms.get_ids()
        data = method()
        n = self.atoms.get_number_of_atoms()
        if self.master:
            shape = (n,) + data.shape[1:]
            result = np.zeros(shape, data.dtype)
            for cpu in range(self.comm.size):
                if cpu != 0:
                    # Receive from cpu
                    nrecv = np.zeros(1, int)
                    self.comm.receive(nrecv, cpu)
                    nrecv = nrecv[0]
                    ids = np.zeros(nrecv, int)
                    data = np.zeros((nrecv,) + result.shape[1:], result.dtype)
                    self.comm.receive(ids, cpu)
                    self.comm.receive(data, cpu)
                result[ids] = data
            return result
        else:
            assert(len(data) == len(ids))
            nsend = np.array([len(ids)])
            self.comm.send(nsend, 0)
            self.comm.send(ids, 0)
            self.comm.send(data, 0)
            return np.zeros((0,)+data.shape[1:], dtype=data.dtype)

    def _cant_set_pbc(self, pbc):
        "Fake set_pbc method."
        raise NotImplementedError("Cannot change PBC of a Collector instance.")
    
    pbc = property(get_pbc, _cant_set_pbc, "The boundary conditions attribute")
    
    def _cant_set_numbers(self, z):
        "Fake set_atomic_numbers method."
        raise NotImplementedError(
            "Cannot change atomic numbers of a Collector instance.")

    numbers = property(get_atomic_numbers, _cant_set_numbers,
                       "Atomic numbers as a property")
    
    def _cant_set_info(self, info):
        "Cannot set info attribute of Collector instance"
        raise NotImplementedError("Cannot set info attribute of Collector instance")
    
    info = property(get_info, _cant_set_info, "The info dictionary")
    