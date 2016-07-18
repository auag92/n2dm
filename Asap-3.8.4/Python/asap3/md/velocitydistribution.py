# VelocityDistributions.py -- set up a velocity distribution

"""Module for setting up e.g. Maxwell-Boltzmann velocity distributions.

Currently, only one function is defined, MaxwellBoltzmannDistribution,
which sets the momenta of a list of atoms according to a
Maxwell-Boltzmann distribution at a given temperature.
"""

import numpy as np
import asap3.mpi
from ase.md.velocitydistribution import _maxwellboltzmanndistribution


def MaxwellBoltzmannDistribution(atoms, temp, force_temp=False):
    """Sets the momenta to a Maxwell-Boltzmann distribution. temp should be
    fed in energy units; i.e., for 300 K use temp=300.*units.kB. If
    force_temp is set to True, it scales the random momenta such that the
    temperature request is precise.
    """
    class fakecomm:
        "Do-nothing communicator for GPAW-like communication."
        def broadcast(self, data, root):
            pass
    momenta = _maxwellboltzmanndistribution(atoms.get_masses(), temp,
                                            fakecomm())
    atoms.set_momenta(momenta)
    if force_temp:
        temp0 = atoms.get_kinetic_energy() / atoms.get_number_of_atoms() / 1.5
        gamma = temp / temp0
        atoms.set_momenta(atoms.get_momenta() * np.sqrt(gamma))


def Stationary(atoms):
    "Sets the center-of-mass momentum to zero."
    p = atoms.get_momenta()
    p0 = np.sum(p, 0)
    # We should add a constant velocity, not momentum, to the atoms
    m = atoms.get_masses()
    mtot = np.sum(m)
    if getattr(atoms, "parallel", False):
        data = np.zeros(4, float)
        data[:3] = p0
        data[3] = mtot
        asap3.mpi.world.sum(data)
        p0 = data[:3]
        mtot = data[3]
    v0 = p0/mtot
    p -= v0*m[:,np.newaxis]
    atoms.set_momenta(p)
    
