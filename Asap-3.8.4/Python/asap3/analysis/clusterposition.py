"""Locate the position of a cluster (its center of mass) in a simulation.

The module calculates the center of mass of a connected cluster of atoms,
even in the case where the cluster straddles the periodic boundaries.  It 
is OK that there are other atoms not connected to the cluster, as long
as the atom with the highest coordination number is part of the cluster,
or an atom being part of the cluster is specified manually.
"""

from asap3 import FullNeighborList
import numpy as np
import sys
from ase.io import PickleTrajectory, BundleTrajectory

class ClusterCenter:
    def __init__(self, cutoff=3.0):
        self.nblist = FullNeighborList(cutoff)
        
    def calculate_center(self, atoms, startatom=None):
        """Calculate the center of mass position of a cluster of atoms.
        
        An atom belonging to the cluster can optionally be specified, if
        not specified one of the highest coordinated atoms is used.
        """
        natoms = len(atoms)
        self.nblist.check_and_update(atoms)
        if startatom is None:
            coordnum = [len(self.nblist[i]) for i in range(natoms)]
            startatom = np.argmax(coordnum)
        self.cluster = []
        self.sumpos = np.zeros(3)
        self.summass = 0
        self.isincluster = np.zeros(natoms, bool)
        self.masses = atoms.get_masses()
        self.add_atom_to_cluster(startatom, atoms[startatom].position)
        com = self.sumpos / self.summass
        return com
    
    def add_atom_to_cluster(self, n, pos):
        "Add an atom and all its neighbors to the cluster."
        self.isincluster[n] = True
        self.sumpos += pos * self.masses[n]
        self.summass += self.masses[n]
        neighbors, reldists, sqdists = self.nblist.get_neighbors(n)
        for i, relpos in zip(neighbors, reldists):
            if not self.isincluster[i]:
                # Add neighboring atom to cluster, using the relative
                # position so periodic boundaries are handled correctly.
                self.add_atom_to_cluster(i, pos + relpos)
            
    def calculate_from_trajectory(self, traj, startatom=None, selector=None):
        """Calculate the center of mass for a cluster in a trajectory file.
        
        traj: The trajectory object, or a file name.
        
        startatom (optional): Specifies an atom guaranteed to be in the cluster.
            If not specified, the atom with the highest coordination number is
            used (if there is only one cluster this should work).
            
        selector (optional): A function defining which atoms should be 
            considered.  The function is called with one argument, the atoms
            object, and should either return an array of booleans, one per
            atom, indicating if they should be included, or return an array
            of integers interpreted as the indices of the atoms to be included.
            This can e.g. be used to select a cluster sitting on a substrate.
            
        This method returns an array of center-of-mass positions, one for each
        frame in the trajectory.
        """
        if isinstance(traj, str):
            if traj.endswith('.traj'):
                traj = PickleTrajectory(traj)
            elif traj.endswith('.bundle'):
                traj = BundleTrajectory(traj)
            else:
                raise ValueError("Cannot handle a file name not ending in .traj or .bundle: " + traj)
        result = []
        for atoms in traj:
            if selector is not None:
                idx = selector(atoms)
                atoms = atoms[idx]
            result.append(self.calculate_center(atoms, startatom))
        return result
    
if __name__ == '__main__':
    # Demo example: Track a cluster of Ni atoms on a surface of something else.
    def sel(a):
        return np.equal(a.get_atomic_numbers(), 28)
    res = ClusterCenter().calculate_from_trajectory(sys.argv[1], selector=sel)
    print np.array(res)
    
