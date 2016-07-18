'''
Created on Mar 16, 2011

@author: schiotz
'''

from ase.md.md import MolecularDynamics
from ase.md.verlet import VelocityVerlet as VelocityVerlet_ASE
import ase
import asap3
import numpy as np
import sys

class VelocityVerlet_Asap(MolecularDynamics):
    def __init__(self, atoms, dt, trajectory=None, logfile=None,
                 loginterval=1):
        MolecularDynamics.__init__(self, atoms, dt, trajectory, logfile,
                                   loginterval)
        if not atoms.has('momenta'):
            atoms.set_momenta(np.zeros((len(atoms), 3), float))
        self.calculator = atoms.get_calculator()
        self.asap_md = asap3._asap.VelocityVerlet(atoms, self.calculator, dt)
        # Identify FixAtoms constraint
        if atoms.constraints:
            assert len(atoms.constraints) == 1
            constraint = atoms.constraints[0]
            assert isinstance(constraint, asap3.constraints.FixAtoms)
            constraint.prepare_for_asap(atoms)
            mask = constraint.index
            assert mask.shape == (len(atoms),) and mask.dtype == bool
            mult = np.logical_not(mask).astype(float)
            self.atoms.arrays["FixAtoms_mult_double"] = mult
            
    def run(self, steps):
        assert(self.calculator is self.atoms.get_calculator())
        self.asap_md.run(steps, self.observers, self)
        
def VelocityVerlet(atoms, dt, trajectory=None, logfile=None, loginterval=1):
    if isinstance(atoms, ase.Atoms) and asap3.constraints.check_asap_constraints(atoms):
        sys.stderr.write("Using Asap-optimized Verlet algorithm\n")
        return VelocityVerlet_Asap(atoms, dt, trajectory, logfile, loginterval)
    else:
        sys.stderr.write("Using ASE Verlet algorithm\n")
        return VelocityVerlet_ASE(atoms, dt, trajectory, logfile, loginterval)

