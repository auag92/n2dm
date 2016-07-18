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

import sys
sys.path.append("..")
from StressModule import *
from ase import data, Atoms
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.verlet import VelocityVerlet
from ase.md.langevin import Langevin
from asap3.mpi import world

debug = 0
if debug == 1:
    DebugOutput("parallelStress-%d.log", nomaster=True)
elif debug == 2:
    time.sleep(world.rank)
    print "PID:", os.getpid()
    time.sleep(20)

print_version(1)

#set_verbose(1)

ismaster = world.rank == 0
isparallel = world.size != 1
if world.size == 1:
    cpulayout = None
elif world.size == 2:
    cpulayout = [2,1,1]
elif world.size == 3:
    cpulayout = [1,3,1]
elif world.size == 4:
    cpulayout = [2,1,2]

defaultstrains = 0.01 * np.array((-1, -0.75, -0.5, -0.25, 0.0,
                                         0.25, 0.5, 0.75, 1.0))
                                       
book = {}
book['Cu'] = {"bookbulk": (134.3, 133.99),
              "bookc11": (168.3, 171.65),
              "bookc12": (122.1, 114.82),
              "bookc44": (75.7, 89.00)}

book['Copper_Rasmussen'] = {"bookbulk": (134.3, 140.77),
                            "bookc11": (168.3, 172.75),
                            "bookc12": (122.1, 124.32),
                            "bookc44": (75.7, 81.17)}

book['Ag'] = {"bookbulk": (103.8, 99.27),
              "bookc11": (124.0, 124.45),
              "bookc12": (93.7, 86.27),
              "bookc44": (46.12, 53.44)}


Copper = "Cu"
Silver = "Ag"
Copper_Rasm = "Copper_Rasmussen"
#for element in (Copper, Copper_Rasm, Silver):
for element in (Copper, Copper_Rasm, Silver):
    symb = element
    if (symb == "Copper_Rasmussen"):
        symbol = "Cu"
    else:
        symbol = symb
    print ""
    print "*** Running test on %s ***" % symb
    z = data.atomic_numbers[symbol]
    latconst = data.reference_states[z]['a']
    initial = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]],
                                symbol=symbol, size=(10,10,10))
    ReportTest("Number of atoms", len(initial), 4000, 0)
    if ismaster:
        atoms = MakeParallelAtoms(initial, cpulayout)
    else:
        atoms = MakeParallelAtoms(None, cpulayout)
    if symb.find("Rasmussen") >= 0:
        print "Using the special EMT parameters by T. Rasmussen"
        atoms.set_calculator(EMT(EMTRasmussenParameters()))
    else:
        atoms.set_calculator(EMT())
    atoms.set_momenta(np.zeros((len(atoms),3), np.float))
    findlatticeconst(atoms, latconst)
    elasticconstants(atoms, symb, **book[symb])

ReportTest.Summary()
