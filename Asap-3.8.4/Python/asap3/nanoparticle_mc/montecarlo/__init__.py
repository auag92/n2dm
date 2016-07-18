from atom import AtomMonteCarlo, AtomMonteCarloData
from surface import SurfaceMonteCarlo, SurfaceMonteCarloData
#from atom_old import AtomMonteCarlo as AtomMonteCarloOld

def get_layers_from_size(size):
    n1 = int(round(pow(1.8 * size, 1./3.) / 2.))
    n2 = 2 * n1
    n3 = int(3./2. * n1)

    return [n1] * 6 + [n2] * 12 + [n3] * 8
