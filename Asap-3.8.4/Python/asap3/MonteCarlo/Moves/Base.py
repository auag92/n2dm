import numpy as np

from ase.units import kB, fs
from ase.md import Langevin

from asap3 import CoordinationNumbers

#from asap3.MonteCarlo.Moves.Surface import SurfaceMove

class Move:
    def __init__(self):
        self.optimizer = None
        self.atoms = None

    def __call__(self, atoms=None):
        self.set_atoms(atoms)
        #self.backup()

    def get_name(self):
        return self.__class__.__name__

    def set_optimizer(self, optimizer=None):
        if optimizer is not None:
            self.optimizer = optimizer

    def set_atoms(self, atoms):
        if atoms is not None:
           self.atoms = atoms

    def backup(self, idx=None):
        """Save positions of atoms to be moved.
        
        idx should be None to save all atoms, or either an integer or an 
        array/list of intergers specifying which atoms to save. 
        """
        if idx is None:
            self.backup_positions = self.atoms.get_positions()
        else:
            self.backup_positions = self.atoms.get_positions()[idx]
        self.backup_idx = idx

    def accept(self):
        pass

    def reject(self):
        if self.backup_idx is None:
            self.atoms.set_positions(self.backup_positions)
        elif isinstance(self.backup_idx, int):
            self.atoms[self.backup_idx].position = self.backup_positions
        else:
            # Set atoms one by one to trigger possible MC optimizations
            # in potential.
            for i, r in zip(self.backup_idx, self.backup_positions):
                self.atoms[i].position = r

    def log(self):
        return ""

class SingleMove(Move):
    """Displaces one atom randomly within a sphere around it self."""
    def __init__(self, rmax):
        self.rmax = rmax
        Move.__init__(self)

    def __call__(self, atoms):
        self.set_atoms(atoms)

        natoms = len(self.atoms)
        n = RandChoice(np.arange(natoms))
        d = RandSphere(1, 0.0, self.rmax)

        self.backup(n)
        self.atoms[n].position += d[0]

class ShakeMove(Move):
    """Displaces all atoms randomly within a sphere around them self."""
    def __init__(self, rmax):
        self.rmax = rmax
        Move.__init__(self)

    def __call__(self, atoms):
        self.set_atoms(atoms)
        self.backup()
        natoms = len(atoms)
        positions = atoms.get_positions()
        d = RandSphere(natoms, 0.0, self.rmax)
        positions += d
        atoms.set_positions(positions)

class BallMove(Move):
    """Moves one atom randomly within the cluster modelled by an
    ellipsoid."""
    def __call__(self, atoms):
        self.set_atoms(atoms)

        natoms = len(atoms)
        positions = atoms.get_positions()
        cn = CoordinationNumbers(atoms)

        # Determine surface atoms
        mask = (cn < 11) * (cn > 5)
        if mask.sum() < 12:
            return positions   # XXXX This cannot be right!

        # Fit surface to ellipsoid
        (A, C, Q) = FitEllipsoid(positions[mask])

        # Move atom
        n = RandChoice(np.arange(natoms))
        d = RandEllipsoid(1, (0, A[0]), (0, A[1]), (0, A[2]))
        self.backup(n)
        atoms[n].position = np.dot(Q, d + C)

class ShellMove(Move):
    """Moves n low coordinated atoms randomly around the surface of the
    cluster. The surface is defined as an ellipsiod."""
    def __init__(self, width, steps=1):
        self.width = width
        self.steps = steps
        Move.__init__(self)

    def __call__(self, atoms):
        self.set_atoms(atoms)
        self.backup()
        
        natoms = len(atoms)
        positions = atoms.get_positions()
        cn = CoordinationNumbers(atoms)

        # Determine surface atoms
        mask = (cn < 11) * (cn > 5)
        if mask.sum() < 12:
            return positions

        # Fit surface to ellipsoid
        (A, C, Q) = FitEllipsoid(positions[mask])
        Amin = A - self.width / 2
        Amax = A + self.width / 2
        a = (Amin[0], Amax[0])
        b = (Amin[1], Amax[1])
        c = (Amin[2], Amax[2])

        # Move atoms
        n = 1
        mask = (cn <= n)
        while mask.sum() < self.steps:
            n += 1
            mask = (cn <= n)

        moved = []
        while len(moved) < self.steps:
            n = RandChoice(np.arange(natoms)[mask])
            if n not in moved:
                moved.append(n)
                d = RandEllipsoid(1, a, b, c)
                positions[n] = np.dot(Q, d + C)

        atoms.set_positions(positions)

class CompressMove(Move):
    def __init__(self, cmax=0.3, cmin=0.2, steps=2, rotate=True):
        if cmin > cmax:
            raise ValueError("cmax must be greater than cmin.")
        self.cmax = cmax
        self.cmin = cmin
        self.steps = steps
        self.rotate = rotate
        Move.__init__(self)

    def __call__(self, atoms):
        self.set_atoms(atoms)
        self.backup() 

        natoms = len(atoms)
        positions = atoms.get_positions()
        center = atoms.get_center_of_atoms()

        for i in range(self.steps):
            if self.rotate:
                positions = RotatePositions(positions, RandSphere(),
                                            [0.0, 0.0, 1.0], center)
            c = np.random.uniform(1.0 - self.cmax, 1.0 - self.cmin, natoms)
            positions[:,2] = c * (positions[:,2] - center[2]) + center[2]

        atoms.set_positions(positions)

class BondMove(Move):
    def __call__(self):
        pass

class BrownianMove(Move):
    """Makes a small Langevin simulation."""
    def __init__(self, steps, timestep, temperature, friction):
        self.steps = steps
        self.timestep = timestep * fs
        self.temp = temperature * kB
        self.friction = friction
        Move.__init__(self)

    def __call__(self, atoms):
        self.set_atoms(atoms)
        self.backup()

        # Start med en temperatur!
        dyn = Langevin(atoms, self.timestep, self.temp, self.friction)
        dyn.run(self.steps)

class ExchangeMove(Move):
    """Exchange two atoms with different type."""
    def __call__(self, atoms):
        self.set_atoms(atoms)

        natoms = len(atoms)
        numbers = atoms.get_atomic_numbers()
        positions = atoms.get_positions()

        bin = np.bincount(numbers)
        num = RandChoice(np.arange(len(bin))[bin > 0])

        mask1 = (numbers == num)
        mask2 = np.logical_not(mask1)

        if not (mask2 == True).any():
            raise Warning("Cannot use exchange move with only one element.")

        n1 = RandChoice(np.arange(natoms)[mask1])
        n2 = RandChoice(np.arange(natoms)[mask2])

        p1 = positions[n1].copy()
        p2 = positions[n2].copy()
        self.backup([n1,n2])
        atoms[n1].position = p2
        atoms[n2].position = p1

### Helper functions ###
def RandSphere(length=1, rmin=0.0, rmax=1.0):
    return RandEllipsoid(length, (rmin, rmax), None, None)

def RandEllipsoid(length=1, a=(0.9, 1.1), b=(0.9, 1.1), c=(0.9, 1.1)):
    if isinstance(a, (list, tuple)):
        (amin, amax) = a
        a = np.random.uniform(amin, amax, length)
    else:
        a = a * np.ones(length)

    if b is None:
        b = a
    elif isinstance(b, (list, tuple)):
        (bmin, bmax) = b
        b = np.random.uniform(bmin, bmax, length)
    else:
        b = b * np.ones(length)

    if c is None:
        c = a
    elif isinstance(c, (list, tuple)):
        (cmin, cmax) = c
        c = np.random.uniform(cmin, cmax, length)
    else:
        c = c * np.ones(length)

    t = np.random.uniform(0.0, np.pi, length)
    p = np.random.uniform(0.0, 2*np.pi, length)

    x = a * np.sin(t) * np.cos(p)
    y = b * np.sin(t) * np.sin(p)
    z = c * np.cos(t)

    xyz = np.zeros((length, 3))
    xyz[:,0] = x
    xyz[:,1] = y
    xyz[:,2] = z

    if length > 1:
        return xyz
    else:
        return xyz[0]

def RotatePositions(positions, u, v, center=[0, 0, 0]):
        # Rotate the vector u into v
        center = np.array(center, float)
        pos = positions - center

        u = np.array(u, float)
        u /= np.linalg.norm(u)

        v = np.array(v, float)
        v /= np.linalg.norm(v)

        c = np.dot(u, v)
        w = np.cross(u, v)
        s = np.linalg.norm(w)
        # In case *u* and *v* are parallel, np.cross(u, v) vanish
        # and can't be used as a rotation axis. However, in this
        # case any rotation axis perpendicular to v2 will do.
        eps = 1e-7
        if s < eps:
            w = np.cross((0, 0, 1), v)
            s = np.linalg.norm(w)
        if s < eps:
            w = np.cross((1, 0, 0), v)
            s = np.linalg.norm(w)
        assert s >= eps
        if s > 0:
            w /= s
        return (c * pos - np.cross(pos, s * w) +
                np.outer(np.dot(pos, w), (1.0 - c) * w) + center)

def RandChoice(a):
    n = np.random.random_integers(0, len(a) - 1)
    return a[n]

def FitEllipsoid(xyz):
    """Fitting an ellipsoid based on Qingde Li et al., "Least Squares
    Ellipsoid Specific Fitting", 2004."""
    # Set up matrices
    C = np.array([[-1, 1, 1, 0, 0, 0],
                  [1, -1, 1, 0, 0, 0],
                  [1, 1, -1, 0, 0, 0],
                  [0, 0, 0, -4, 0, 0],
                  [0, 0, 0, 0, -4, 0],
                  [0, 0, 0, 0, 0, -4]])
    Ci = np.linalg.inv(C)

    X = np.array([xyz[:,0]**2, xyz[:,1]**2, xyz[:,2]**2,
                  2*xyz[:,1]*xyz[:,2], 2*xyz[:,0]*xyz[:,2],
                  2*xyz[:,0]*xyz[:,1],
                  2*xyz[:,0], 2*xyz[:,1], 2*xyz[:,2],
                  np.ones(len(xyz))])

    XXT = np.dot(X, X.T)
    S11 = XXT[:6,:6].copy()
    S12 = XXT[:6,6:].copy()
    S12T = XXT[6:,:6].copy()
    S22i = np.linalg.inv(XXT[6:,6:])

    # Solve fitting equation
    eig = np.linalg.eig(np.dot(Ci, S11 - np.dot(S12, np.dot(S22i, S12T))))
    v = eig[1][:, eig[0].argmax()]
    c = np.append(v, -np.dot(np.dot(S22i, S12T), v))

    # Reduce the quadratic equation
    A = np.array([[c[0], c[5], c[4]],
                  [c[5], c[1], c[3]],
                  [c[4], c[3], c[2]]])
    B = np.array([c[6], c[7], c[8]])
    J = c[9]

    eig = np.linalg.eig(A)
    A1 = eig[0]
    Q = eig[1]
    B1 = np.dot(B, Q)
    K = (B1**2 / A1).sum() - J

    # Return semi-axes, center and rotation matrix
    return (np.sqrt(K/A1), -B1/A1, Q)

