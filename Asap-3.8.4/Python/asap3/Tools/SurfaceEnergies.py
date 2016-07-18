import numpy as np

from ase.optimize import polyfit, LBFGS, LBFGSLineSearch
from ase.lattice.surface import fcc100, fcc110, fcc111, hcp0001

def SurfaceEnergy(images, natoms, calc, fmax=0.01, debug=False):
    """Calculate the surface energy from a list of slab images.

    Parameters
    ----------
    images: List of slab images which the calculation is based on. The x and
    y dimensions of the images unit cell must be conserved.

    natoms: Number of atoms in a atomic layer in the slabs.

    calc: Calculator object that can be attached to the images.
    """
    ucell = images[0].get_cell()

    layers = np.zeros(len(images))
    energies = np.zeros(len(images))

    if debug:
        print "Layers   Energy   LBFGS steps"

    for i, atoms in enumerate(images):
        cell = atoms.get_cell()
        if (ucell[0] != cell[0]).any() or (ucell[1] != cell[1]).any():
            raise ValueError("The x and y dimensions of the unit cell must be conserved.")

        atoms.set_calculator(calc)

        dyn = LBFGS(atoms, logfile=None, trajectory=None)
        dyn.run(fmax=fmax, steps=1000)
        assert dyn.converged(), "LBFGS not converged in 100 steps!"

        layers[i] = len(atoms) / natoms
        energies[i] = atoms.get_potential_energy()

        if debug:
            print "%4i%12.3f%10i" % (layers[i], energies[i], dyn.get_number_of_steps())

    p = polyfit(layers.reshape((-1,1)), energies, 1)
    surface_energy = p.c[0] / (2 * natoms)
    assert surface_energy > 0.0

    return surface_energy

def get_surface_slabs(symbol, surfacetype, surfacesize, layers, latticeconstants,
                      vacuum=10.0):
    images = []
    if surfacetype == 'fcc100' or surfacetype == 'fcc111':
        nlayeratoms = surfacesize[0] * surfacesize[1]
        if isinstance(latticeconstants, (list, tuple)):
            if not len(latticeconstants) == 1:
                raise ValueError("Too many latticeconstants")
            a = latticeconstants[0]
        else:
            a = latticeconstants
    elif surfacetype == 'hcp0001':
        nlayeratoms = surfacesize[0] * surfacesize[1]
        if isinstance(latticeconstants, (list, tuple)):
            if not len(latticeconstants) == 2:
                raise ValueError("Latticeconstants must be on the form [a, c]")
            a = latticeconstants[0]
            c = latticeconstants[1]
        else:
            raise ValueError("Latticeconstants must be on the form [a, c]")
    else:
        raise ValueError("Surface type '%s' is not supported" % (surfacetype,))


    for n in layers:
        size = surfacesize + (n,)
        if surfacetype == 'fcc100':
            images.append(fcc100(symbol=symbol, size=size, a=a, vacuum=vacuum))
        elif surfacetype == 'fcc111':
            images.append(fcc111(symbol=symbol, size=size, a=a, vacuum=vacuum,
                                 orthogonal=True))
        elif surfacetype == 'hcp0001':
            images.append(hcp0001(symbol=symbol, size=size, a=a, c=c,
                                  vacuum=vacuum, orthogonal=True))

    return images, nlayeratoms

if __name__ == '__main__':
    from asap3 import EMT
    from ase.data import reference_states

    debug = False
    layers = np.arange(5, 14, 2)

    print "Calculating (100) surface energies with EMT"
    images, nlayeratoms = get_surface_slabs('Cu', 'fcc100', (6, 6), layers,
                                            reference_states[29]['a'])
    surface_energy = SurfaceEnergy(images, nlayeratoms, EMT(), debug=debug)
    print "Cu:", surface_energy
    images, nlayeratoms = get_surface_slabs('Ag', 'fcc100', (6, 6), layers,
                                            reference_states[47]['a'])
    surface_energy = SurfaceEnergy(images, nlayeratoms, EMT(), debug=debug)
    print "Ag:", surface_energy
    images, nlayeratoms = get_surface_slabs('Au', 'fcc100', (6, 6), layers,
                                            reference_states[79]['a'])
    surface_energy = SurfaceEnergy(images, nlayeratoms, EMT(), debug=debug)
    print "Au:", surface_energy

    print "\nCalculating (111) surface energies with EMT"
    images, nlayeratoms = get_surface_slabs('Cu', 'fcc111', (6, 6), layers,
                                            reference_states[29]['a'])
    surface_energy = SurfaceEnergy(images, nlayeratoms, EMT(), debug=debug)
    print "Cu:", surface_energy
    images, nlayeratoms = get_surface_slabs('Ag', 'fcc111', (6, 6), layers,
                                            reference_states[47]['a'])
    surface_energy = SurfaceEnergy(images, nlayeratoms, EMT(), debug=debug)
    print "Ag:", surface_energy
    images, nlayeratoms = get_surface_slabs('Au', 'fcc111', (6, 6), layers,
                                            reference_states[79]['a'])
    surface_energy = SurfaceEnergy(images, nlayeratoms, EMT(), debug=debug)
    print "Au:", surface_energy


