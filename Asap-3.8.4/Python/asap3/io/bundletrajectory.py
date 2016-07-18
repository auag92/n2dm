import os
import time
import numpy as np
import asap3
from asap3.mpi import world
from asap3.io.trajectory import _GetAtoms
try:
    from ase.io.bundletrajectory import BundleTrajectory as _BundleTrajectory
    from ase.io.bundletrajectory import PickleBundleBackend
except ImportError:
    class _BundleTrajectory:
        "This version of ASE does not support BundleTrajectory."
        def __init__(self, *args, **kwargs):
            raise NotImplementedError(self.__doc__)
    class PickleBundleBackend:
        "This version of ASE does not support BundleTrajectory."
        def __init__(self, *args, **kwargs):
            raise NotImplementedError(self.__doc__)

    
class BundleTrajectory(_BundleTrajectory, _GetAtoms):
    """Reads and writes atoms into a .bundle directory.

    The BundleTrajectory is an alternative way of storing
    trajectories, intended for large-scale molecular dynamics
    simulations, where a single flat file becomes unwieldy.  Instead,
    the data is stored in directory, a 'bundle' (the name bundle is
    inspired from bundles in Mac OS, which are really just directories
    the user is supposed to think of as a single file-like unit).

    Parameters:

    filename:
        The name of the directory.  Preferably ending in .bundle.

    mode (optional):
        The file opening mode.  'r' means open for reading, 'w' for
        writing and 'a' for appending.  Default: 'r'.  If opening in
        write mode, and the filename already exists, the old file is
        renamed to .bak (any old .bak file is deleted), except if the
        existing file is empty.

    atoms (optional):
        The atoms that will be written.  Can only be specified in
        write or append mode.  If not specified, the atoms must be
        given as an argument to the .write() method instead.

    backup (optional):
        Use backup=False to disable renaming of an existing file.

    split (optional):
        If set to True or False, determines whether a split file
        format is used instead of the normal one.  In the split
        format, each processor in a parallel simulation writes its own
        files inside the BundleTrajectory, instead of leaving all I/O
        to the master.  In not specified, a split format is used if
        more than one million atoms.  Ignored in serial simulations.

    iolimit (optional):
        Limits the number of MPI tasks performing I/O simultaneously,
        to prevent overloading the NFS server.  Only enforced if the
        number of tasks is somewhat larger than the limit.
        
    """
    def __init__(self, filename, mode='r', atoms=None,
                 backup=True, split=None, iolimit=10):
        if split is None:
            # Decide if subtype should be split based on number of atoms
            split = (atoms is not None) and atoms.get_number_of_atoms() > 1000000
        # Never use subtype split for serial simulations
        if not getattr(atoms, "parallel", False):
            split = False
        # Use the collector object to join all data on master if subtype is normal
        # and the simulation is parallel.
        if not split and atoms is not None and getattr(atoms, "parallel", False):
            atoms = asap3.Collector(atoms)
        if split:
            self.subtype = 'split'
        else:
            self.subtype = 'normal'
        # self.iolimit may be needed if reading a split bundle.
        if world.size < 1.5 * iolimit:
            self.iolimit = None
        else:
            self.iolimit = iolimit
        
        _BundleTrajectory.__init__(self, filename, mode, atoms,
                                   backup=backup)
        if self.subtype == 'split':
            self.set_extra_data('ID')  # So the atoms can be sorted when read.
        
    def _set_defaults(self):
        subtype = self.subtype   # Preserve it
        _BundleTrajectory._set_defaults(self)
        self.subtype = subtype
        self.datatypes['forces'] = False

    def _set_backend(self, backend=None):
        """Set the backed doing the actual I/O."""
        if backend is not None:
            self.backend_name = backend
        if self.backend_name == 'pickle':
            if self.subtype == 'normal':
                # Use the standard ASE backend
                self.backend = PickleBundleBackend(self.master)
            elif self.subtype == 'split':
                self.backend = PickleSplitBundleBackend(self.master,
                                                        self.iolimit)
        else:
            raise NotImplementedError(
                "This version of ASE cannot use BundleTrajectory with backend '%s'"
                % self.backend_name)

    def write(self, atoms=None):
        if self.subtype == 'normal' and atoms is not None and getattr(atoms, "parallel", False):
            atoms = asap3.Collector(atoms)
        _BundleTrajectory.write(self, atoms)

    def _make_bundledir(self, filename):
        """Make the main bundle directory.

        Since all MPI tasks might write to it, all tasks must wait for
        the directory to appear.

        For performance reasons, the first frame directory is created immediately.
        """
        assert not os.path.isdir(filename)
        world.barrier()
        if self.master:
            self.log("Making directory "+filename)
            os.mkdir(filename)
            framedir = os.path.join(filename, "F0")
            self.log("Making directory "+ framedir)    
            os.mkdir(framedir)
        else:
            i = 0
            while not os.path.isdir(filename):
                time.sleep(1)
                i += 1
            if i > 10:
                self.log("Waiting %d seconds for %s to appear!"
                         % (i, filename))

    def _make_framedir(self, frame):
        """Make subdirectory for the frame.

        For a split bundle, all MPI tasks write to the frame
        directory.  The slaves must therefore wait until it becomes
        available.  To minimize the waiting time, frames are
        pre-created.
        """
        if self.subtype == 'split':
            numdirs = 10
        else:
            numdirs = 1
        if self.master:
            for i in range(frame, frame+numdirs):
                framedir = os.path.join(self.filename, "F"+str(i))
                if not os.path.exists(framedir):
                    self.log("Making directory " + framedir)
                    os.mkdir(framedir)
        framedir = os.path.join(self.filename, "F"+str(frame))
        # Wait for the directory to appear
        i = 0
        while not os.path.isdir(framedir):
            time.sleep(1)
            i += 1
        if i > 10:
            self.log("Waiting %d seconds for %s to appear!"
                     % (i, framedir))
        return framedir

    def close(self):
        """Clean up when closing."""
        if self.state == 'write' and self.master:
            i = self.nframes
            while True:
                fname = os.path.join(self.filename, "F" + str(i))
                if not os.path.exists(fname):
                    break
                self.log("Closing, removing empty directory "+fname)
                os.rmdir(fname)
                i += 1
        _BundleTrajectory.close(self)

    def __del__(self):
        self.close()

class PickleSplitBundleBackend(PickleBundleBackend):
    """A special backend for writing split bundles (ASAP only)."""
    def __init__(self, master, iolimit):
        # Store if this backend will actually write anything
        self.writesmall = master
        self.writelarge = True
        self.writenonarray = master
        self.iolimit = iolimit
        if iolimit:
            self.iostart = np.round(
                np.linspace(0, world.size, iolimit+1)).astype(int)
            self.iotag = 413151
        self.lastwritedir = None
            
    def write_small(self, framedir, smalldata):
        "Write small data to be written jointly."
        smalldata['fragments'] = world.size
        PickleBundleBackend.write_small(self, framedir, smalldata)
        
    def write(self, framedir, name, data):
        "Write data to separate file."
        if hasattr(data, "shape"):
            # We need to store which kind of data was written in this frame 
            # so NFS synchronization is possible when closing file.
            if framedir != self.lastwritedir:
                self.lastwritedir = framedir
                self.writenames = []
            self.writenames.append(name)
            # As expected, we are writing a NumPy array
            self.iosync_start()
            name = "%s_%d" % (name, world.rank)
            PickleBundleBackend.write(self, framedir, name, data)
            self.iosync_end()
        elif self.writenonarray:
            # If the data is not a NumPy array, only the master writes.
            PickleBundleBackend.write(self, framedir, name, data)

    def read(self, framedir, name):
        "Read data from separate file."
        self.iosync_start()
        x = PickleBundleBackend.read(self, framedir, name)
        self.iosync_end()
        return x
    
    def iosync_start(self):
        "Prevents too many simultaneous IO tasks from trashing server."
        if self.iolimit and world.rank not in self.iostart:
            # I must wait.
            token = np.zeros(1, int)
            world.receive(token, world.rank-1, self.iotag)

    def iosync_end(self):
        if self.iolimit and world.rank+1 not in self.iostart:
            # Another task is waiting for me.
            token = np.zeros(1, int)
            world.send(token, world.rank+1, self.iotag)
            
    def close(self, log=None):
        """Make sure that all data is available on disk for all MPI tasks."""
        if self.lastwritedir:
            for name in self.writenames:
                for part in range(world.size):
                    fname = os.path.join(self.lastwritedir, "%s_%d.pickle" % (name, part))
                    if not os.path.exists(fname):
                        if log:
                            log.write("Task %i is waiting for '%s' to appear.\n" %
                                      (world.rank, fname))
                        for i in range(20):
                            time.sleep(5)
                            if os.path.exists(fname):
                                break
                        if not os.path.exists(fname) and log:
                            log.write("WARNING: Task %i gave up waiting for '%s'.\n" %
                                      (world.rank, fname))

