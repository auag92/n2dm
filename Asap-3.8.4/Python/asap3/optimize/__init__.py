try:
    from ase.optimize.optimize import Optimizer as ase_Optimizer
except ImportError:
    # Fall back to old placement
    from ase.optimize import Optimizer as ase_Optimizer
from asap3 import parallelpossible as asap_parallel
if asap_parallel:
    from asap3.mpi import world


class Optimizer(ase_Optimizer):
    def converged(self, forces=None):
        """Did the optimization converge?"""
        if forces is None:
            forces = self.atoms.get_forces()
        mxf2 = (forces**2).sum(axis=1).max()
        if asap_parallel:
            mxf2 = world.max(mxf2)
        return mxf2 < self.fmax**2
