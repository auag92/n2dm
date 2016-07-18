"""Asap module ParallelListOfAtoms.

Defines the parallel list of atoms object (`ParallelAtoms`), and a factory
method for creating them (`MakeParallelAtoms`).

Importing this module also installs a Python exit function causing
MPI_Abort to be called if an uncaught exception occurs.
"""


__docformat__ = "restructuredtext en"

import ase
import ase.units
from asap3 import _asap
import asap3.mpi
import numpy as np
import cPickle, cStringIO
import sys, time
import ase.parallel

class ParallelAtoms(ase.Atoms):
    """Atoms class for parallel Asap simulations.

    It is recommended to create ParallelAtoms objects using
    `MakeParallelAtoms`.
    """
    parallel = 1
    def __init__(self, nCells, comm, atoms, cell=None, pbc=None,
                 distribute=True):
        """Create a ParallelAtoms object.

        WARNING: ParallelAtoms object should normally not be created
        explicitly.  Use MakeParallelAtoms instead.
        """
        self.nCells = np.array(nCells, int)
        self.comm = comm
        if self.nCells.shape != (3,):
            raise ValueError, "ParallelAtoms: nCells must be 3 integers."
        # Extract all data from the atoms object.  Not nice :-)
        self.arrays = {}
        for name in atoms.arrays.keys():
            self.arrays[name] = atoms.arrays[name].copy()
        assert self.arrays["positions"].dtype == np.dtype(float)
        assert self.arrays["positions"].shape == (len(atoms), 3)
        assert self.arrays["numbers"].shape ==  (len(atoms),)

        if cell is not None:
            self.cell = cell
        else:
            self.cell = atoms.get_cell()
        if pbc is not None:
            self.pbc = pbc
        else:
            self.pbc = atoms.get_pbc()
        self.adsorbate_info = {}
        self.info = {}

        self.ghosts = {}
        self.ghosts["positions"] = np.zeros((0,3), float)
        self.ghosts["numbers"] = np.zeros(0, self.arrays["numbers"].dtype)

        # Now make the IDs
        mynatoms = np.array([len(self)])
        natoms_all = np.zeros(self.comm.size, int)
        self.comm.all_gather(mynatoms, natoms_all)
        if not self.arrays.has_key("ID"):
            firstID = sum(natoms_all[:self.comm.rank])
            self.arrays["ID"] = np.arange(firstID, firstID+len(atoms))
        self.total_number_of_atoms = sum(natoms_all)

        # Atoms should have no constraints
        self.set_constraint(None)
        
        if distribute:
            self.distribute()

    def distribute(self):
        _asap.DistributeAtoms(self)

    def get_number_of_atoms(self):
        n = len(self)
        return self.comm.sum(n)
        
    def get_list_of_elements(self):
        """Get a list of elements.

        The list is cached to prevent unnecessary communication.
        """
        try:
            return self.listofelements
        except AttributeError:
            z = self.get_atomic_numbers()
            present = np.zeros(100, int)
            zmax = z.max()
            zmin = z.min()
            present[zmin] = present[zmax] = 1
            for i in range(zmin+1, zmax):
                if np.equal(z, i).any():
                    present[i] = 1
            self.comm.sum(present)
            self.listofelements = []
            for i, p in enumerate(present):
                if p:
                    self.listofelements.append(i)
            return self.listofelements

    def set_atomic_numbers(self, numbers):
        """Set the atomic numbers."""
        try:
            # Discard the cached list of elements
            del self.listofelements
        except AttributeError:
            pass
        ase.Atoms.set_atomic_numbers(self, numbers)

    def get_ids(self):
        """Get the atom IDs in a parallel simulation."""
        return self.arrays["ID"].copy()

    def is_master(self):
        """Return 1 on the master node, 0 on all other nodes."""
        return (self.comm.rank == 0)

    def get_comm(self):
        return self.comm
    
    def wrap_calculator(self, calc):
        "Make an ASAP calculator compatible with parallel simulations."
        try:
            parallelOK = calc.supports_parallel()
        except AttributeError:
            parallelOK = False
        if not parallelOK:
            raise ValueError, "The calculator does not support parallel ASAP calculations."
        return _asap.ParallelPotential(calc)

    def set_calculator(self, calc, wrap=True):
        """Sets the calculator in a way compatible with parallel simulations.
        
        calc: 
            The Calculator to be used.  Normally only Asap calculators will work.
            
        wrap (optional, default=True):
            Indicates if a calculator should be wrapped in a ParallelCalculator object.  
            Wrapping is the default, and should almost always be used, the only exception
            being if the Calculator is implemented as a Python object wrapping an Asap
            calculator, in which case the Asap calculator should first be wrapped in
            a ParallelCalculator object (use atoms.wrap_calculator) and this one should then
            be used by the Python calculator.  The Python calculator is then attached
            without being wrapped again.
        """
        if wrap:
            parcalc = self.wrap_calculator(calc)
        else:
            parcalc = calc
        ase.Atoms.set_calculator(self, parcalc)

    def get_kinetic_energy(self):
        local_ekin = ase.Atoms.get_kinetic_energy(self)
        return self.comm.sum(local_ekin)
    
    def get_temperature(self):
        """Get the temperature. in Kelvin"""
        ekin = self.get_kinetic_energy() / self.get_number_of_atoms()
        return ekin / (1.5 * ase.units.kB)

    def get_ghost_positions(self):
        return self.ghosts['positions'].copy()
    
    def get_ghost_atomic_numbers(self):
        return self.ghosts['numbers'].copy()
    

def MakeParallelAtoms(atoms, nCells, cell=None, pbc=None,
                      distribute=True):
    """Build parallel simulation from serial lists of atoms.

    Call simultaneously on all processors.  Each processor having
    atoms should pass a list of atoms as the first argument, or None
    if this processor does not contribute with any atoms.  If the
    cell and/or pbc arguments are given, they must be given on
    all processors, and be identical.  If it is not given, a supercell
    is attempted to be extracted from the atoms on the processor with
    lowest rank.

    This is the preferred method for creating parallel simulations.
    """
    import cPickle, cStringIO

    mpi = asap3.mpi
    #comm = mpi.world.duplicate()
    comm = mpi.world

    # Sanity check: is the node layout reasonable
    nNodes = nCells[0] * nCells[1] * nCells[2]
    if nNodes != comm.size:
        raise RuntimeError("Wrong number of CPUs: %d != %d*%d*%d" %
                           (comm.size, nCells[0], nCells[1], nCells[2]))
    t1 = np.zeros((3,))
    t2 = np.zeros((3,))
    comm.min(t1)
    comm.max(t2)
    if (t1[0] != t2[0] or t1[1] != t2[1] or t1[2] != t2[2]):
        raise RuntimeError, "CPU layout inconsistent."

    # If pbc and/or cell are given, they may be shorthands in need of
    # expansion.
    if pbc:
        try:
            plen = len(pbc)
        except TypeError:
            # It is a scalar, interpret as a boolean.
            if pbc:
                pbc = (1,1,1)
            else:
                pbc = (0,0,0)
        else:
            if plen != 3:
                raise ValueError, "pbc must be a scalar or a 3-sequence."
    if cell:
        cell = array(cell)  # Make sure it is a numeric array.
        if cell.shape == (3,):
            cell = array([[cell[0], 0, 0],
                          [0, cell[1], 0],
                          [0, 0, cell[2]]])
        elif cell.shape != (3,3):
            raise ValueError, "Unit cell must be a 3x3 matrix or a 3-vector."

    # Find the lowest CPU with atoms, and let that one distribute
    # which data it has.  All other CPUs check for consistency.
    if atoms is None:
        hasdata = None
        mynum = comm.size
    else:
        hasdata = {}
        for name in atoms.arrays.keys():
            datatype = np.sctype2char(atoms.arrays[name])
            shape = atoms.arrays[name].shape[1:]
            hasdata[name] = (datatype, shape)
        mynum = comm.rank
        if pbc is None:
            pbc = atoms.get_pbc()
        if cell is None:
            cell = atoms.get_cell()
    root = comm.min(mynum)   # The first CPU with atoms
    # Now send hasdata, cell and pbc to all other CPUs
    package = cPickle.dumps((hasdata, cell, pbc), 2)
    package = comm.broadcast_string(package, root)
    rootdata, rootcell, rootpbc = cPickle.loads(package)
    if rootdata is None or len(rootdata) == 0:
        raise ValueError, "No data from 'root' atoms.  Empty atoms?!?"
    
    # Check for consistent cell and pbc arguments
    if cell is not None:
        if rootcell is None:
            raise TypeError, "Cell given on another processor than the atoms."
        if (cell.ravel() - rootcell.ravel()).max() > 1e-12:
            raise ValueError, "Inconsistent cell specification."
    else:
        cell = rootcell   # May still be None
    if pbc is not None:
        if rootpbc is None:
            raise TypeError, "PBC given on another processor than the atoms."
        if (pbc != rootpbc).any():
            raise ValueError, "Inconsistent pbc specification."
    else:
        pbc = rootpbc

    # Check for consistent atoms data
    if hasdata is not None:
        if hasdata != rootdata:
            raise ValueError, "Atoms do not contain the sama data on different processors."
    if "positions" not in rootdata:
        raise ValueError, "Atoms do not have positions!"
    
    # Create empty atoms
    if atoms is None:
        atoms = ase.Atoms(cell=cell, pbc=pbc)
        for name in rootdata.keys():
            if atoms.arrays.has_key(name):
                assert np.sctype2char(atoms.arrays[name]) == rootdata[name][0]
                assert len(atoms.arrays[name]) == 0
            else:
                shape = (0,) + rootdata[name][1]
                atoms.arrays[name] = np.zeros(shape, rootdata[name][0])
        
    return ParallelAtoms(nCells, comm, atoms, cell=cell, pbc=pbc, 
                         distribute=distribute)



# A cleanup function should call MPI_Abort if python crashes to
# terminate the processes on the other nodes.
ase.parallel.register_parallel_cleanup_function()

# _oldexitfunc = getattr(sys, "exitfunc", None)
# def _asap_cleanup(lastexit = _oldexitfunc, sys=sys, time=time,
#                   comm = asap3.mpi.world):
#     error = getattr(sys, "last_type", None)
#     if error:
#         sys.stdout.flush()
#         sys.stderr.write("ASAP CLEANUP (node " + str(comm.rank) +
#                          "): " + str(error) +
#                          " occurred.  Calling MPI_Abort!\n")
#         sys.stderr.flush()
#         # Give other nodes a moment to crash by themselves (perhaps
#         # producing helpful error messages).
#         time.sleep(3)
#         comm.abort(42)
#     if lastexit:
#         lastexit()
# sys.exitfunc = _asap_cleanup
        
# END OF PARALLEL STUFF
