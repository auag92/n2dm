import numpy as np

from asap3.MonteCarlo.Moves.Base import RandSphere

from ase.atoms import Atoms, symbols2numbers
from ase.data import covalent_radii

def RandomBoxAtoms(symbols, cell=[5.0, 5.0, 5.0]):
    numbers = symbols2numbers(symbols)
    scaled_positions = np.random.uniform(0.0, 1.0, (len(numbers), 3))
    return Atoms(numbers, scaled_positions=scaled_positions, cell=cell)

def RandomSphereAtoms(symbols, covalent=None, radius=None, vacuum=0.0):
    numbers = symbols2numbers(symbols)
    if covalent is None:
        covalent = 0.8 * np.mean(covalent_radii[numbers])
    if radius is None:
        radius = covalent * float(len(numbers))**(1.0/3.0)
    cell = 3.0 * radius * np.ones(3, float) + vacuum
    positions = radius * np.ones((len(numbers), 3))
    for i in range(len(numbers)):
        pos = radius * RandSphere(1) + 0.5 * cell
        while np.any(np.sum((positions - pos)**2, axis=1) < covalent**2):
            pos = radius * RandSphere(1) + 0.5 * cell
        positions[i] = pos.copy()
    #positions = radius * RandSphere(len(numbers)) + 0.5 * cell
    return Atoms(numbers, positions=positions, cell=cell)

def AlloyAtomsRandom(symbols, atoms):
    numbers_new = symbols2numbers(symbols)
    numbers_old = atoms.get_atomic_numbers()
    natoms = len(numbers_old)

    if len(numbers_old) != len(numbers_new):
        raise ValueError("The number of atoms in symbols (%i) " %(len(numbers_new),) +
                         "must be equal to that in atoms (%i)!" % (len(atoms),))

    numbers = []
    indexes = []
    for n in np.unique(numbers_new + numbers_old.tolist()):
        d = (numbers_new == n).sum() - (numbers_old == n).sum()
        if d > 0: # Add
            while d > 0:
                i = np.random.random_integers(0, len(numbers))
                numbers.insert(i, n)
                d -= 1
        elif d < 0: # Remove
            ind = np.arange(natoms)[numbers_old == n].tolist()
            while d < 0:
                i = np.random.random_integers(0, len(ind) - 1)
                indexes.append(ind.pop(i))
                d += 1

    numbers_old[indexes] = numbers
    atoms.set_atomic_numbers(numbers_old)
    return atoms
