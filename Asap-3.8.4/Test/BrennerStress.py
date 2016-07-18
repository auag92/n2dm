#!/usr/bin/env python
# encoding: utf-8

"""
Tests that stresses are calculated correctly by Asap

Name: testStress.py

Description: Part of the Asap test suite.  Tests that stresses are
    calculated correctly by calculating various elastic constants from
    strain and stress and comparing them with the same constants
    calculated using strains and energies.  Can also be imported as a
    module, and used to test advanced calculation methods (QC,
    parallel).

Usage: python testStress.py

Expected result: Some elastic constants for Copper and Silver,
followed by the text 'ALL TESTS SUCCEEDED'.

The elastic constants are calculated by straining the crystal in
various modes, and fitting to the energies or the stresses.
Generally, the energies seems to be most sensitive to numerical noise,
and requires a rather large strain interval (1% or above), whereas the
stresses are much less sensitive to this.  On the other hand,
unlinearities influence the stress fits for large strain intervals.  A
strain interval of [-1%, 1%] is a good compromise, where both methods
work.

C11 and C12 calculated 'directly', i.e. in uniaxial strain, is
sensitive to the strain interval when using the energy to fit C11.
Fitting to the stresses work much better.

C11 and C12 can be calculated in an alternative way using a
volume-conserving deformation and fitting to the energies.

All of the above-mentioned calculations are performed.

"""

from asap3 import *
from ase import data, Atoms
from ase.lattice.cubic import Diamond
from ase.lattice.hexagonal import Graphite
#from asap3.md.verlet import VelocityVerlet
#from asap3.md.langevin import Langevin
from asap3.testtools import ReportTest
import numpy as np

defaultstrains = 0.01 * np.array((-1, -0.75, -0.5, -0.25, 0.0,
                                         0.25, 0.5, 0.75, 1.0))
                                       
book = {}
book['Diamond'] = {"bookbulk": (442.0, 441.4),
                   "bookc11": (1079.0, 1068.9),
                   "bookc12": (124.0, 131.5),
                   "bookc44": (578.0, 735.8),
                   'symbol': 'C'}

# Alexey Bosak and Michael Krisch
# European Synchrotron Radiation Facility, BP 220, F-38043 Grenoble Cedex, France

# Marcel Mohr, Janina Maultzsch, and Christian Thomsen
# Institut für Festkörperphysik, Technische Universität Berlin, Hardenbergstrasse 36, 10623 Berlin, Germany

# Received 22 November 2006; revised 11 January 2007; published 30 April 2007

# The five independent elastic moduli of single-crystalline graphite
# are determined using inelastic x-ray scattering. At room temperature
# the elastic moduli are, in units of GPa, C11=1109, C12=139, C13=0,
# C33=38.7, and C44=4.95. Our experimental results are compared with
# predictions of ab initio calculations and previously reported
# incomplete and contradictory data sets. We obtain an upper limit of
# 1.1 TPa for the on-axis Young’s modulus of homogeneous carbon
# nanotube, thus providing important constraints for further
# theoretical advances and quantitative input to model elasticity in
# graphite nanotubes.

# URL:
# http://link.aps.org.globalproxy.cvt.dk/doi/10.1103/PhysRevB.75.153408
# DOI:
# 10.1103/PhysRevB.75.153408

book['Graphite'] = {"bookbulk": (134.3, 140.77),
                    "bookc11": (1109, 172.75),
                    "bookc12": (139, 124.32),
                    "bookc13": (0, 0),
                    "bookc33": (38.7, 38.7),
                    "bookc44": (75.7, 81.17),
                    'symbol': 'C'}

book['Silicon'] = {"bookbulk": (100.0, 97.8),
                   "bookc11": (124.0, 142.6),
                   "bookc12": (93.7, 75.5),
                   "bookc44": (46.12, 118.9),
                   'symbol': 'Si'}

book['Germanium'] = {"bookbulk": (103.8, 99.27),
                     "bookc11": (124.0, 124.45),
                     "bookc12": (93.7, 86.27),
                     "bookc44": (46.12, 53.44),
                     'symbol': 'Ge'}



defaultstrains = 0.005 * np.array((-1, -0.75, -0.5, -0.25, 0.0,
                                         0.25, 0.5, 0.75, 1.0))

def polynomialLeastSquaresFit(parameters, data, max_iterations=None,
                              stopping_limit = 0.0001):
    """Least-square fit to a polynomial.
    
    Least-squares fit to a polynomial whose order is defined by
    the number of parameter values.

    This is a wrapper function replacing the similar function in
    Scientific.Functions.LeastSquares
    """
    order = len(parameters)
    return (np.polyfit(data[:,0], data[:,1], order)[::-1], None)

def makefits(atoms, strains, indices, shear=0):
    """Do the deformations, and get fits to energies and stresses.

    atoms is a list of atoms.
    strains are the strains as floating point numbers
    indices is a list of indices for the strain components to be used.
        For bulk modulus, it will be ((0,0), (1,1), (2,2)) etc.
    """
    energies = []
    basis = atoms.get_cell()
    vol = np.linalg.det(basis)
    for epsilon in strains:
        if shear:
            adjustment = np.zeros((3,3), np.float)
            for idx in indices:
                adjustment[idx[0]] += idx[2] * epsilon * basis[idx[1]]
            atoms.set_cell(adjustment + basis, scale_atoms=True)
        else:
            scaling = np.ones((3,3), np.float)
            for idx in indices:
                scaling[idx[0]] += idx[1]*epsilon
            atoms.set_cell(scaling * basis, scale_atoms=True)
        energy = atoms.get_potential_energy()
        #print ""
        #print epsilon, energy/len(atoms)
        energies.append((epsilon, energy/vol))
    atoms.set_cell(basis, scale_atoms=True)
    energies = np.array(energies)

    #print "Energies:", energies
    energyfit = polynomialLeastSquaresFit((0.0, 0.0, 0.0), energies)
    #print "EnergyFit:", energyfit
    return energyfit

def findlatticeconst(atoms, latticeconstant):
    """Adjust the volume so the atoms have their lowest energy."""
    basis = atoms.get_cell()
    strains = 0.01 * np.array((-0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1,
                               0.15, 0.2))
    for i in range(5):
        energyfit = makefits(atoms, strains,
                             (((0,0), 1), ((1,1), 1), ((2,2),1)))
        pressurefit = -energyfit[0]
        dilation = -pressurefit[1]/(2*pressurefit[2])
        print "Optimizing lattice constant:", latticeconstant, "->", latticeconstant*(1+dilation/3)
        latticeconstant = latticeconstant*(1+dilation/3)
        basis = (1+dilation/3) * basis
        atoms.set_cell(basis, scale_atoms=True)

def elasticconstants(atoms, name, bookbulk, bookc11, bookc12, bookc44,
                     bookc13=None, bookc33=None,
                     fitfact=1.0, fitfact2=1.0):
    """Check the elastic constants."""
    if bookc13 is None:
        bookc13 = bookc12
    if bookc33 is None:
        bookc33 = bookc11
    energyfit = makefits(atoms, defaultstrains,
                         (((0,0), 1), ((1,1), 1), ((2,2),1)))
    bm = 2.0/9.0 * energyfit[0][2] / units.GPa
    print ""
    print "Calculation for", name
    print "  Bulk modulus from energies:", bm
    print "  Textbook value:", bookbulk[0]
    print ""
    ReportTest("Bulk modulus (%s, energies)" % name, bm, bookbulk[1],
               0.1*fitfact)
    energyfit = makefits(atoms, defaultstrains, (((0,0),1),))
    c11 = 2.0 * energyfit[0][2] / units.GPa
    energyfit = makefits(atoms, defaultstrains,
                         (((0,0),1), ((1,1), -0.5), ((2,2), -0.5)))
    # C11 - C12
    c11mc12 = 4.0/3.0 * energyfit[0][2] / units.GPa
    c12 = c11 - c11mc12
    print ""
    print "Calculation for", name
    print "  C_11 from energies:", c11
    print "  C_12 from energies:", c12
    print "    Textbook values: C_11 = %.1f; C_12 = %.1f" % (bookc11[0],
                                                             bookc12[0])
    
    # B =  (C11 + 2 C12)/3
    altc11 = (3*bm + 2*c11mc12) / 3.0
    altc12 = (3*bm - c11mc12) / 3.0
    print "  C_11 from alternative energies:", altc11
    print "  C_12 from alternative energies:", altc12
    print "    Bulk modulus from C_11 and C_12:", (altc11 + 2 * altc12) / 3.0

    energyfit = makefits(atoms, defaultstrains, (((2,2),1),))
    c33 = 2.0 * energyfit[0][2] / units.GPa
    # C11 - C12
    c11mc12 = 4.0/3.0 * energyfit[0][2] / units.GPa
    c12 = c11 - c11mc12
    print "  C_33 from energies:", c11
    print "    Textbook values: C_33 = %.1f" % (bookc33[0],)
    print ""
    
    ReportTest("C11 from energies", c11, bookc11[1], 0.1*fitfact)
    ReportTest("C11 from alt. energies", altc11, c11, 3.0*fitfact)
    ReportTest("C12 from energies", c12, bookc12[1], 0.5*fitfact)
    ReportTest("C12 from alt. energies", altc12, c12, 3.0*fitfact)
    ReportTest("C33 from energies", c33, bookc33[1], 0.1*fitfact)
    print ""

    
    energyfit = makefits(atoms, defaultstrains,
                         (((2,1), (2,2), 1), ((1,2), (1,1), 1)), shear=1)
    c44 = 0.5 * energyfit[0][2] / units.GPa
    print ""
    print "Calculation for", name
    print "  C_44 from energies:", c44
    print "  Textbook value:", bookc44[0]
    print ""
    ReportTest("C44 from energies", c44, bookc44[1], 0.1*fitfact)
    print ""
    print "Testing other shear modes:"
    energyfit = makefits(atoms, defaultstrains,
                         (((2,0), (2,2), 1), ((0,2), (0,0), 1)), shear=1)
    c44x = 0.5 * energyfit[0][2] / units.GPa
    ReportTest("C44(alt) from energies", c44x, c44, 0.2*fitfact2)
    energyfit = makefits(atoms, defaultstrains,
                         (((1,0), (1,1), 1), ((0,1), (0,0), 1)), shear=1)
    c44x = 0.5 * energyfit[0][2] / units.GPa
    ReportTest("C44(alt) from energies", c44x, c44, 0.1*fitfact)
    print ""

print_version(1)

for element in ('Diamond', 'Graphite', 'Silicon'):
    symbol = book[element]['symbol']
    print "*** Running test on %s ***" % element
    z = data.atomic_numbers[symbol]
    if element == 'Graphite':
        latconst = 2.46
        intial = Graphite(directions=[[2,-1,-1,0], [0,1,-1,0], [0,0,0,1]],
                          symbol=symbol, size=(7,7,5),
                          latticeconstant = {'a':2.46,  'c':6.71})
        #view(initial)
        ReportTest("Number of atoms", len(initial), 7*7*5*8, 0)
    else:
        latconst = data.reference_states[z]['a']
        initial = Diamond(directions=[[1,0,0],[0,1,0],[0,0,1]],
                          symbol=symbol, size=(6,6,6))
        ReportTest("Number of atoms", len(initial), 6*6*6*8, 0)
    atoms = Atoms(initial)
    #view(atoms)
    atoms.set_calculator(BrennerPotential())
    atoms.set_momenta(np.zeros((len(atoms),3), np.float))
    findlatticeconst(atoms, latconst)
    del book[element]['symbol']
    elasticconstants(atoms, element, **book[element])

ReportTest.Summary()
