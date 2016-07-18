#!/usr/bin/env python
"""
Tests that atoms are counted correctly

Name: AtomCounting.py

Description: Part of the Asap test suite.  Tests indirect addressing
by atomic number (might fail with newer Numeric versions and/or on 64
bit machines)

Usage: python AtomCounting.py

Expected result: The output shoul end with 'ALL TESTS SUCCEEDED'.
"""

import sys
from Numeric import *
from LinearAlgebra import *
from Asap import *
from Asap.testtools import ReportTest
from Numeric import *
from Asap.Setup.Lattice.Compounds import L1_2

atoms = L1_2(directions=((1,0,0),(0,1,0),(0,0,1)), size=(5,5,5),
             element=("Cu", "Au"), latticeconstant=4.0, debug=0)
natoms = len(atoms)
n = zeros(93)
z = atoms.GetAtomicNumbers()
for i in z:
    n[i] += 1

expected = zeros(93)
expected[29] = natoms / 4
expected[79] = 3 * natoms / 4
for i in range(1,93):
    ReportTest(('Number of atoms with Z=%d' % i), n[i], expected[i], 0)

    
    
ReportTest.Summary()
