"""nanocrystalline.py - set up a nanocrystalline sample."""

from ase.lattice.cubic import BodyCenteredCubic, FaceCenteredCubic
from ase.visualize import view
from ase.optimize import FIRE
from ase import units
from asap3.md.nptberendsen import Inhomogeneous_NPTBerendsen
from asap3 import NeighborCellLocator
import numpy as np

def random_rot():
    """Returns a random rotation matrix."""
    # The random rotation matrix should most elegantly be produced from three random
    # Euler angles (or similar), with suitable non-uniform distributions to compensate
    # for the "north pole singularity".  However, it is easier (but less elegant) to
    # get a random rotation by chaining a large number of not quite uniformly distributed
    # random rotations.

    g = np.identity(3)  # Initial rotation matrix
    N = 100             # Number of rotations
    for i in range(N):
        phi, theta, psi = np.random.random(3) * 2.0 * np.pi
        # First Euler rotation about z in matrix form
        D = np.array(((np.cos(phi), np.sin(phi), 0.),
                      (-np.sin(phi), np.cos(phi), 0.),
                      (0., 0., 1.)))
        # Second Euler rotation about x:
        C = np.array(((1., 0., 0.),
                      (0., np.cos(theta), np.sin(theta)),
                      (0., -np.sin(theta), np.cos(theta))))
        # Third Euler rotation, 2nd rotation about z:
        B = np.array(((np.cos(psi), np.sin(psi), 0.),
                      (-np.sin(psi), np.cos(psi), 0.),
                      (0., 0., 1.)))
        # Total Euler rotation
        A = np.dot(B, np.dot(C, D))
        # Apply rotation
        g = np.dot(g, A)
    return g
    

def bcc_grains(size, perturbation=0.0):
    "Creates a grain layout based on a perturbed BCC lattice."
    graincenters = BodyCenteredCubic(symbol='H', latticeconstant=1.0, size=size)
    graincenters = graincenters.get_positions() + 0.25
    graincenters /= size
    pert = np.random.standard_normal(graincenters.shape) * perturbation
    graincenters += pert
    rotations = [random_rot() for g in graincenters]
    return graincenters, rotations


class make_nanocrystal:
    """Create a nanocrystalline sample.

    Parameters:
      size:      Size of the system, in Angstrom (either a number or three numbers).
      centers:   Positions of grains in scaled coordinates.
      rotations: Rotation matrices for the grains.
      delta (optional):
                 Grains are cut delta Angstrom before the mathematical boundary
                 between them (default: 0).
      min (optional):
                 If two atoms are closer than this distance, one of them is removed
                 (default: 2.0).

    In addition, either the 'unit' or the 'symbol' parameter must be specified.
      unit (optional):
                 A unit cell for building the crystals (Atoms object).
                 MUST be orthorhombic!
      symbol (optional):
                 If unit is None, then an FCC crystal of this element is used.
      latticeconstant (optional):
                 If symbol is specified, this overrides the default lattice constant.

    Returns:
      A factory that creates a nanocrystalline sample when its run() method is called.
    """
    def __init__(self, size, centers, rotations, delta=0.0, min=2.0,
                 unit=None, symbol=None, latticeconstant=None):
        try:
            assert len(size) == 3
        except TypeError:
            # Not a sequence
            size = (size, size, size)
        self.size = np.array(size)
        assert self.size.shape == (3,)
        self.centers = centers
        self.rotations = rotations
        self.delta = delta
        self.min = min
        assert len(centers) == len(rotations)
        self.centers *= self.size
        if unit is None:
            unit = FaceCenteredCubic(symbol=symbol, size=(1,1,1), pbc=True,
                                     latticeconstant=latticeconstant)
        self.unit = unit

    def run(self):
        "Create the nanocrystalline sample"
        result = None
        grainnumber = 0
        for center, rotation in zip(self.centers, self.rotations):
            assert center.shape == (3,)
            assert rotation.shape == (3,3)
            # Make a chunk of crystal sufficiently large, centered at (0,0,0)
            atoms = self.make_bulk()
            # Rotate it
            pos = atoms.get_positions()
            atoms.set_positions(np.dot(atoms.get_positions(), rotation))
            # Cut a grain out of it
            atoms = self.cut_bulk(atoms, center)
            atoms.set_tags(grainnumber * np.ones(len(atoms), int))
            # Add it to the sample
            if result is None:
                result = atoms
                atoms.set_cell(self.size, scale_atoms=False)
                atoms.set_pbc(True)
            else:
                result.extend(atoms)
        self.remove_close(result)
        return result

    def make_bulk(self):
        "Create a sufficiently large bulk sample"
        minsize = np.sqrt((self.size * self.size).sum())
        uc = self.unit.get_cell().diagonal()
        mult = [int(np.ceil(minsize / x)) * 2 + 1 for x in uc]
        atoms = self.unit.repeat(mult)
        atoms.set_pbc(False)
        center = atoms.get_cell().diagonal() / 2.0
        atoms.set_positions(atoms.get_positions() - center)
        return atoms

    def cut_bulk(self, atoms, center):
        "Center 'atoms' at 'center' and cut of any parts closer to any other center."
        # Replicate grain centers (PBC) into 26 neighboring cells.
        #print "grain at %s:" % (str(center),)
        a = (-1, 0, 1)
        ctr = []
        for i in a:
            for j in a:
                for k in a:
                    ctr.append(self.centers + np.array([i,j,k]) * self.size)
        centers = np.concatenate(ctr)
        # Find directions to other grain centers
        rcent = centers - center
        d2 = (rcent * rcent).sum(axis=1)
        assert d2.shape == (len(centers),)
        order = np.argsort(d2)
        # Vectors and distances to neighboring grains
        rvecs = rcent[order[1:]]
        dists = np.sqrt(d2[order[1:]])
        # Now center these atoms at this grain
        atoms.set_positions(atoms.get_positions() + center)
        # For each other grain we now cut away atoms that are beyond the halfway plane.
        for rvec, dist in zip(rvecs, dists):
            r = atoms.get_positions() - center
            r_projected = np.inner(r, rvec) * (1.0 / dist)
            keep = r_projected < (0.5 * dist - self.delta)
            #print "  Cut toward %s: Keeping %i atoms of %i" % (str(rvec), keep.sum(), len(atoms))
            atoms = atoms[keep]
        print "grain at %s: %i atoms." % (str(center), len(atoms))
        return atoms
    
    def remove_close(self, atoms):
        "Remove atoms that are too close"
        print "Removing atoms that are too close (d_cut = %.3f Angstrom)" % (self.min,)
        remove = {}
        nblist = NeighborCellLocator(self.min, atoms)
        # Loop over all pairs of too-close atoms
        for i in range(len(atoms)):
            for j in nblist[i]:
                # Pair (i,j) is too close
                if not (i in remove or j in remove):
                    if len(remove) % 2:
                        remove[i] = None
                    else:
                        remove[j] = None
        print "Removing %i atoms out of %i" % (len(remove), len(atoms))
        del atoms[remove.keys()]

def minimize_energy(atoms, nstep, pressure_interval=10):
    bulkmodulus = 140e9 / 1e5  # Bulk modulus of typical metal (Cu), in bar.
    dyn = FIRE(atoms)
    unstress = Inhomogeneous_NPTBerendsen(atoms, 50*units.fs, 0,
                                          taup=5000*units.fs,
                                          compressibility=1/bulkmodulus)
    dyn.attach(unstress.scale_positions_and_cell, interval=10)
    dyn.run(steps=nstep)
        
if __name__ == "__main__":
    from asap3 import EMT
    from asap3.analysis.localstructure import RestrictedCNA, CoordinationNumbers

    centers, rotations = bcc_grains((2,2,2), 0.01)
    #centers, rotations = bcc_grains((1,1,1), 0.)
    fact = make_nanocrystal(50.0, centers, rotations, symbol='Cu', delta=0.35, min=2.0)
    atoms = fact.run()
    atoms.set_scaled_positions(atoms.get_scaled_positions())
    print "Running CNA"
    RestrictedCNA(atoms)
    coord = CoordinationNumbers(atoms)
    print "Total number of atoms", len(atoms)
    for i in range(coord.min(), coord.max()+1):
        print "CN=%-2i: %i" % (i, (coord==i).sum())
    view(atoms)
    vol_before = atoms.get_volume()
    print "Energy minimization."
    atoms.set_calculator(EMT())
    minimize_energy(atoms, 10000)
    vol_after = atoms.get_volume()
    print "Volume change: %.1f -> %.1f  (%.1f%%)" % (vol_before, vol_after, 
                                                    100*(vol_after/vol_before - 1))
    print "Running CNA again"
    RestrictedCNA(atoms)
    coord = CoordinationNumbers(atoms)
    print "Total number of atoms", len(atoms)
    for i in range(coord.min(), coord.max()+1):
        print "CN=%-2i: %i" % (i, (coord==i).sum())
    view(atoms)
