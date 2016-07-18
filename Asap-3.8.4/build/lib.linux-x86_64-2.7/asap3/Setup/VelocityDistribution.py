# VelocityDistributions.py -- set up a velocity distribution

"""Module for setting up e.g. Maxwell-Boltzmann velocity distributions.

Currently, only one function is defined, MaxwellBoltzmannDistribution,
which sets the momenta of a list of atoms according to a
Maxwell-Boltzmann distribution at a given temperature.
"""

import Numeric
import RandomArray

def _maxwellboltzmanndistribution(masses, temp):
    xi = RandomArray.standard_normal(shape=(len(masses),3))
    momenta = xi * Numeric.sqrt(masses * temp)[:,Numeric.NewAxis]
    return momenta

def MaxwellBoltzmannDistribution(atoms, temp):
    """Sets the momenta to a Maxwell-Boltzmann distribution."""
    momenta = _maxwellboltzmanndistribution(atoms.GetMasses(), temp)
    atoms.SetCartesianMomenta(momenta)

    
