from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.testtools import ReportTest
import numpy as np
from OpenKIM_modelname import openkimmodel

#set_verbose(1)
step = 500

if OpenKIMsupported:
    pbc_candidates = ["NEIGH_RVEC_H", 
                      "NEIGH_PURE_H",
                      "NEIGH_RVEC_F", 
                      "NEIGH_PURE_F",
                      "MI_OPBC_H", 
                      "MI_OPBC_F",
                      "CLUSTER"
                      ]
else:
    pbc_candidates = []
    print "OpenKIM support is not compiled into Asap."

for nbltype in pbc_candidates:
    if nbltype.startswith("MI_OPBC"):
        boundaries = ((1,1,1),)
    elif nbltype == "CLUSTER":
        boundaries = ((0,0,0),)
    else:
        boundaries = ((1,1,1), (0,0,0), (0,0,1))
        
    for pbc in boundaries:
        txt = nbltype + (" PBC=%1i%1i%1i " % pbc)
        # Test that EMT reimported through OpenKIM gives the right results.
        atoms_kim = FaceCenteredCubic(size=(10,10,10), symbol='Cu')
        #atoms_kim = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]],
        #                    size=(30, 30, 30),
        #                    symbol="Cu")
        natoms = len(atoms_kim)
        atoms_kim.set_pbc(pbc)
        r = atoms_kim.get_positions()
        r.flat[:] += 0.1 * np.sin(np.arange(3*natoms))
        atoms_kim.set_positions(r)
        atoms_emt = atoms_kim.copy()
        kim = OpenKIMcalculator(openkimmodel, allowed=nbltype)
        emt = EMT()
        emt.set_subtractE0(False)
        atoms_kim.set_calculator(kim)
        atoms_emt.set_calculator(emt)
        ek = atoms_kim.get_potential_energy()
        ee = atoms_emt.get_potential_energy()
        ReportTest(txt+"Total energy", ek, ee, 1e-8)
        ek = atoms_kim.get_potential_energies()
        ee = atoms_emt.get_potential_energies()
        for i in range(0, natoms, step):
            ReportTest(txt+"Energy of atom %i" % (i,), ek[i], ee[i], 1e-8)
        fk = atoms_kim.get_forces()
        fe = atoms_emt.get_forces()
        n = 0
        for i in range(0, natoms, step):
            n = (n + 1) % 3
            ReportTest(txt+"Force(%i) of atom %i" % (n, i), fk[i, n], fe[i, n], 1e-8)
        sk = atoms_kim.get_stress()
        se = atoms_emt.get_stress()
        for i in range(6):
            ReportTest(txt+"Stress(%i)" % (i,), sk[i], se[i], 1e-8)
        sk = atoms_kim.get_stresses()
        se = atoms_emt.get_stresses()
        for i in range(0, natoms, step):
            n = (n + 1) % 6
            # Volume per atom is not defined the same way: greater tolerance needed
            ReportTest(txt+"Stress(%i) of atom %i" % (n, i), sk[i, n], se[i, n], 1e-3)

    
ReportTest.Summary()
