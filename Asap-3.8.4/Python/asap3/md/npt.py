from ase.md.npt import NPT as aseNPT
from asap3.md.md import ParallelMolDynMixin
from asap3.mpi import world
import sys

if hasattr(aseNPT, '_npt_version'):
    class NPT(object, aseNPT, ParallelMolDynMixin):        
        def __init__(self, atoms, 
                     timestep, temperature, externalstress, ttime, pfactor,
                     mask=None, trajectory=None, logfile=None, loginterval=1):
            ParallelMolDynMixin.__init__(self, 'NPT', atoms)
            aseNPT.__init__(self, atoms, 
                     timestep, temperature, externalstress, ttime, pfactor,
                     mask=mask, trajectory=trajectory, logfile=logfile,
                     loginterval=loginterval)

        def _getnatoms(self):
            """Get the number of atoms.

            In a parallel simulation, this is the total number of atoms on all
            processors.
            """
            return self.atoms.get_number_of_atoms()    

        def get_center_of_mass_momentum(self):
            "Get the center of mass momentum."
            cm = self.atoms.get_momenta().sum(0)
            world.sum(cm)
            return cm

        @classmethod
        def read_from_trajectory(cls, trajectory,  frame=-1, cpus=None, extra=None):
            """Read dynamics and atoms from trajectory (Class method).

            Simultaneously reads the atoms and the dynamics from a BundleTrajectory,
            including the internal data of the NPT dynamics object (automatically
            saved when attaching a BundleTrajectory to an NPT object).

            Arguments::

            trajectory 
                The filename or an open asap3.io.BundleTrajectory object.

            frame (optional)
                Which frame to read.  Default: the last. 

            cpus (optional)
                If specified, the atoms are read in parallel and distributed
                for a parallel simulation.  This argument should then be the
                cpu layout (a tuple of three integers).
            """
            if isinstance(trajectory, str):
                if trajectory.endswith('/'):
                    trajectory = trajectory[:-1]
                if trajectory.endswith('.bundle'):
                    from asap3.io.bundletrajectory import BundleTrajectory
                    trajectory = BundleTrajectory(trajectory)
                else:
                    raise ValueError("Cannot open '%': unsupported file format" % trajectory)
            atoms = trajectory.get_atoms(frame, cpus, extra=extra)
            # The original read_from_trajectory method does the rest of the
            # work.  As it is a class method, the unbound method is not available
            # in the usual way, but as the im_func attribute.  Change to __func__
            # in Python 3.X
            return aseNPT.read_from_trajectory.im_func(cls, trajectory, 
                                                       frame=frame, atoms=atoms)   

        # Store these arrays on the atoms
        q = property(lambda s: s.get("q"), lambda s, x: s.set("q", x))
        q_past = property(lambda s: s.get("q_past"), lambda s, x: s.set("q_past", x))
        q_future = property(lambda s: s.get("q_future"), lambda s, x: s.set("q_future", x))

else:  
    # hasattr(aseNPT, '_npt_version'):
    sys.stderr.write("Warning: Asap.NPT based on old ASE version.")
    class NPT_OldAse(object, aseNPT):
        def __init__(self, atoms, *args, **kwargs):
            aseNPT.__init__(self, atoms, *args, **kwargs)
            if getattr(atoms, "parallel", False):
                raise NotImplementedError("NPT dynamics does not work in parallel with this version of ASE.  Please use the developer version.")
    NPT = NPT_OldAse

