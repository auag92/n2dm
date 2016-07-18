import asap3
from ase.io.trajectory import PickleTrajectory as _PickleTrajectory

class _GetAtoms:
    "Mixin class implementing the get_atoms method."
    def get_atoms(self, n, cpus=None, extra=None):
        """Return Atoms object number n from the trajectory.

        traj.get_atoms(n) is very similar to traj[n], but with the following
        differences:

        In serial simulations are traj.get_atoms(n) and traj[n] equivalent.

        In parallel simulations:
        traj[n] returns all atoms on all processors

        traj.get_atoms(n) returns all atoms on the master, and None on
        all other processors.

        traj.get_atoms(n, cpus) returns a ParallelAtoms object with
        atoms distributed among the processors, and processor layout
        given by cpus (three integers).
        
        Optional parameters:
        extra (BundleTrajectory only): 
            A list of names of extra arrays that should be read from the 
            input file.  The data is available with atoms.get_array.
        """
        if self.master:
            atoms = self[n]
            if extra:
                for d in extra:
                    atoms.set_array(d, self.read_extra_data(d, n))
        else:
            atoms = None
        if cpus is None:
            return atoms
        else:
            return asap3.MakeParallelAtoms(atoms, cpus)


class PickleTrajectory(_PickleTrajectory, _GetAtoms):
    write_forces = False
    write_charges = False
    write_magmoms = False

    def __init__(self, filename, mode='r', atoms=None, master=None,
                 backup=True):
        if master is not None:
            if atoms.get_comm().rank != 0:
                raise NotImplementedError("It is required that the cpu with rank 0 is the master")
        _PickleTrajectory.__init__(self, filename, mode, atoms, master,
                                   backup=backup)

    def set_atoms(self, atoms=None):
        if atoms is not None and getattr(atoms, "parallel", False):
            atoms = asap3.Collector(atoms, self.master)
            self.sanitycheck = False
        _PickleTrajectory.set_atoms(self, atoms)
        
    def write(self, atoms=None):
        if atoms is not None and getattr(atoms, "parallel", False):
            atoms = asap3.Collector(atoms, self.master)
            self.sanitycheck = False
        _PickleTrajectory.write(self, atoms)

