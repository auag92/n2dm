from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.testtools import ReportTest
from asap3.EMT2013Parameters import PtY_parameters

def MakeLJ():
    "Lennard-Jones potential."
    p = LennardJones([78], [0.15], [2.35], -1.0, False)
    return p

def MakeEMT2013():
    "EMT2013 potential."
    return EMT2013(PtY_parameters)

for pot in (EMT, MakeLJ, MakeEMT2013):
    print "Testing:", pot.__doc__.split("\n")[0]
    
    atoms1 = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]],
                              size=(6,6,6), symbol="Pt")
    emt = pot()
    atoms1.set_calculator(emt)
    f1 = atoms1.get_forces()


    atoms2 = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]],
                              size=(5,5,5), symbol="Pt")
    atoms2.set_calculator(emt)
    f3 = atoms1.get_forces()
    maxdev = abs(f3 - f1).max()
    print maxdev
    ReportTest("Max error 1:", maxdev, 0.0, 1e-6)

    atoms1 = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]],
                              size=(6,6,6), symbol="Pt")
    emt = pot()
    atoms1.set_calculator(emt)
    f1 = atoms1.get_forces()


    atoms2 = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]],
                              size=(5,5,5), symbol="Pt")
    atoms2.set_calculator(emt)
    atoms3 = atoms2.copy()
    atoms3.set_calculator(pot())
    f0 = atoms3.get_forces()
    f2 = atoms2.get_forces()
    f3 = atoms1.get_forces()
    maxdev = abs(f2 - f0).max()
    print maxdev
    ReportTest("Max error 2:", maxdev, 0.0, 1e-6)
    maxdev = abs(f3 - f1).max()
    print maxdev
    ReportTest("Max error 3:", maxdev, 0.0, 1e-6)

ReportTest.Summary()
