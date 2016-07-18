#!/usr/bin/env python
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

import numpy as np
import sys
from asap3 import *
from asap3.testtools import ReportTest

# When not imported as a module, the following are imported later:
#from Setup.Lattice.FCC111Ortho import *
#from Structures.ChemicalElements import Copper, Silver
#import Structures.ChemicalElements.AtomicWeight
#import Structures.ChemicalElements.CrystalStructure
#from Structures.IonDynamics import VelocityVerlet, Langevin

defaultstrains = 0.01 * np.array((-1, -0.75, -0.5, -0.25, 0.0,
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
    stresses = []
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
        stress = atoms.get_stress()
        #print ""
        #print epsilon, energy/len(atoms)
        #print stress
        energies.append((epsilon, energy/vol))
        stresses.append((epsilon,) + tuple(stress))
    atoms.set_cell(basis, scale_atoms=True)
    energies = np.array(energies)
    stresses = np.array(stresses)

    #print "Energies:", energies
    energyfit = polynomialLeastSquaresFit((0.0, 0.0, 0.0), energies)
    #print "EnergyFit:", energyfit
    stressfits = []
    for i in range(6):
        stressfits.append(polynomialLeastSquaresFit((0.0, 0.0),
                                                    np.take(stresses,
                                                            (0, i+1), 1)))
    return (energyfit, stressfits)

def findlatticeconst(atoms, latticeconstant):
    """Adjust the volume so the atoms have their lowest energy."""
    basis = atoms.get_cell()
    strains = 0.01 * np.array((-0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1,
                               0.15, 0.2))
    for i in range(5):
        (energyfit, stressfit) = makefits(atoms, strains,
                                          (((0,0), 1), ((1,1), 1), ((2,2),1)))
        pressurefit = -(np.array(stressfit[0][0]) +
                        np.array(stressfit[1][0]) +
                        np.array(stressfit[2][0])) / 3.0
        dilation = -pressurefit[0]/pressurefit[1]
        print "Optimizing lattice constant:", latticeconstant, "->", latticeconstant*(1+dilation/3)
        latticeconstant = latticeconstant*(1+dilation/3)
        basis = (1+dilation/3) * basis
        atoms.set_cell(basis, scale_atoms=True)

def elasticconstants(atoms, name, bookbulk, bookc11, bookc12, bookc44,
                     fitfact=1.0, fitfact2=1.0):
    """Check the elastic constants."""
    (energyfit, stressfits) = makefits(atoms, defaultstrains,
                                       (((0,0), 1), ((1,1), 1), ((2,2),1)))
    bm = 2.0/9.0 * energyfit[0][2] / units.GPa
    bms = []
    for i in range(3):
        bms.append(stressfits[0][0][1] / 3.0 / units.GPa)
    avgbms = np.sum(bms)/3.0
    print ""
    print "Calculation for", name
    print "  Bulk modulus from energies:", bm
    print "  Bulk modulus from stresses:", np.sum(bms)/3.0
    print "  Textbook value:", bookbulk[0]
    print ""
    ReportTest("Bulk modulus (%s, energies)" % name, bm, bookbulk[1],
               0.1*fitfact)
    ReportTest("Bulk modulus (%s, pressure)" % name, avgbms, bm, 0.1*fitfact)
    ReportTest("Bulk modulus (%s, s_xx)" % name, bms[0], avgbms,
               0.00001*avgbms)
    ReportTest("Bulk modulus (%s, s_yy)" % name, bms[1], avgbms,
               0.00001*avgbms)
    ReportTest("Bulk modulus (%s, s_zz)" % name, bms[2], avgbms,
               0.00001*avgbms)
    (energyfit, stressfits) = makefits(atoms, defaultstrains, (((0,0),1),))
    c11 = 2.0 * energyfit[0][2] / units.GPa
    c11s = stressfits[0][0][1] / units.GPa
    c12s = stressfits[1][0][1] / units.GPa
    print ""
    print "Calculation for", name
    print "  C_11 from energies:", c11
    print "  C_11 from stress:", c11s
    print "  C_12 from stress:", c12s
    print "    Bulk modulus from C_11 and C_12:", (c11s + 2 * c12s) / 3.0
    print "    Textbook values: C_11 = %.1f; C_12 = %.1f" % (bookc11[0],
                                                             bookc12[0])
    (energyfit, stressfits) = makefits(atoms, defaultstrains,
                                       (((0,0),1), ((1,1), -0.5),
                                        ((2,2), -0.5)))
    # C11 - C12
    c11mc12 = 4.0/3.0 * energyfit[0][2] / units.GPa
    # B =  (C11 + 2 C12)/3
    altc11 = (3*bm + 2*c11mc12) / 3.0
    altc12 = (3*bm - c11mc12) / 3.0
    print "  C_11 from alternative energies:", altc11
    print "  C_12 from alternative energies:", altc12
    print "    Bulk modulus from C_11 and C_12:", (altc11 + 2 * altc12) / 3.0
    print ""
    
    ReportTest("C11 from stress", c11s, bookc11[1], 0.1*fitfact)
    ReportTest("C11 from energies", c11, c11s, 0.1*fitfact)
    ReportTest("C12 from stress", c12s, bookc12[1], 0.1*fitfact)
    ReportTest("B from C11 and C12", (c11s + 2 * c12s) / 3.0, bookbulk[1],
               0.5*fitfact)
    ReportTest("C11 from alt. energies", altc11, c11s, 0.5*fitfact)
    ReportTest("C12 from alt. energies", altc12, c12s, 0.5*fitfact)
    print ""

    (energyfit, stressfits) = makefits(atoms, defaultstrains,
                                       (((2,1), (2,2), 1),
                                        ((1,2), (1,1), 1)), shear=1)
    c44 = 0.5 * energyfit[0][2] / units.GPa
    c44s = 0.5 * stressfits[3][0][1] / units.GPa
    print ""
    print "Calculation for", name
    print "  C_44 from energies:", c44
    print "  C_44 from stresses:", c44s
    print "  Textbook value:", bookc44[0]
    print "    (please do not expect good agreement for the usual EMT potential)"
    print ""
    ReportTest("C44 from energies", c44, bookc44[1], 0.1*fitfact)
    ReportTest("C44 from stresses", c44s, c44, 0.2*fitfact)
    nulls = (("C14", 0), ("C24", 1), ("C34", 2), ("C45", 4), ("C46", 5))
    print "Testing that the non-existing elastic constants are indeed zero:"
    for name, idx in nulls:
        cnull = stressfits[idx][0][1] / units.GPa
        #print cnull
        ReportTest(name+" from stresses", cnull, 0.0, 0.00001)
    print ""
    print "Testing other shear modes:"
    (energyfit, stressfits) = makefits(atoms, defaultstrains,
                                       (((2,0), (2,2), 1),
                                        ((0,2), (0,0), 1)), shear=1)
    c44x = 0.5 * energyfit[0][2] / units.GPa
    c44sx = 0.5 * stressfits[4][0][1] / units.GPa
    ReportTest("C44(alt) from energies", c44x, c44, 0.2*fitfact2)
    ReportTest("C44(alt) from stresses", c44sx, c44s, 0.2*fitfact2)
    nulls = (("C14", 0), ("C24", 1), ("C34", 2), ("C45", 3), ("C46", 5))
    for name, idx in nulls:
        cnull = stressfits[idx][0][1] / units.GPa
        ReportTest(name+"(alt) from stresses", cnull, 0.0, 0.00001)
    (energyfit, stressfits) = makefits(atoms, defaultstrains,
                                       (((1,0), (1,1), 1),
                                        ((0,1), (0,0), 1)), shear=1)
    c44x = 0.5 * energyfit[0][2] / units.GPa
    c44sx = 0.5 * stressfits[5][0][1] / units.GPa
    ReportTest("C44(alt) from energies", c44x, c44, 0.1*fitfact)
    ReportTest("C44(alt) from stresses", c44sx, c44s, 0.1*fitfact)
    nulls = (("C14", 0), ("C24", 1), ("C34", 2), ("C45", 3), ("C46", 4))
    for name, idx in nulls:
        cnull = stressfits[idx][0][1] / units.GPa
        ReportTest(name+"(alt) from stresses", cnull, 0.0, 0.00001)
    print ""
                                       
book = {}
book['Cu'] = {"bookbulk": (134.3, 133.52),
              "bookc11": (168.3, 171.65),
              "bookc12": (122.1, 114.82),
              "bookc44": (75.7, 88.45)}

book['Copper_Rasmussen'] = {"bookbulk": (134.3, 140.17),
                            "bookc11": (168.3, 172.75),
                            "bookc12": (122.1, 124.32),
                            "bookc44": (75.7, 80.39)}

book['Ag'] = {"bookbulk": (103.8, 99.01),
              "bookc11": (124.0, 124.45),
              "bookc12": (93.7, 86.27),
              "bookc44": (46.12, 53.25)}

if __name__ == "__main__":
    print "This it not a test, but a module imported from a few tests."
    
