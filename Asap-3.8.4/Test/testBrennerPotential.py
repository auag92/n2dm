#!/usr/bin/python
# Copyright (C) 2007  CAMP
# Please see the accompanying LICENSE file for further information.

# Tests the BrennerPotential implementation

from ASE import Atom
from Asap.Trajectories.NetCDFTrajectory import *
from Asap import *
from Asap.Dynamics.VelocityVerlet import VelocityVerlet
from Asap.Dynamics.PredictorCorrector import PredictorCorrector
from Asap.testtools import ReportTest
from Asap.Utilities import BrennerImport
from cPickle import *
from Numeric import *
from time import time, localtime, strftime

PrintVersion(1)


def setParamTest():
  print "\nRunning brenner setParameter test"
  brenner = BrennerImport.BrennerImport()
  brenner.initCoord("coord.d")
  brenner.initInput("input.d")

  atoms=brenner.GetListOfAtoms()
  potential = BrennerPotential()
  atoms.SetCalculator(potential)
  
  potential.SetMinRangeTree(50)
  ReportTest(("SetMinRangeTree: "), 50, potential.GetMinRangeTree(), 0, silent=True)

  potential.SetKFlag(50)
  ReportTest(("SetKFlag: "), 50, potential.GetKFlag(), 0, silent=True)

  potential.SetSqueeze(50)
  ReportTest(("SetSqueeze: "), 50, potential.GetSqueeze(), 0, silent=True)

  potential.SetVolumeScaleDir(50)
  ReportTest(("SetVolumeScaleDir: "), 50, potential.GetVolumeScaleDir(), 0, silent=True)

  potential.SetRll(100)
  ReportTest(("rll: "), 100, potential.GetRll(), 0, silent=True)
   
#Tests the brenner dynamics with the original files provided by the fortran
#implementation
def brennerDynamicsTest():
    print  "\nRunning brenner dynamics test ..."
    brenner = BrennerImport.BrennerImport()
    brenner.initCoord("coord.d")
    brenner.initInput("input.d")

    a=brenner.GetListOfAtoms()
    potential = BrennerPotential(False)
    potential.SetRll(2.5)
    potential.SetCubeSize(brenner.GetCubeSize())
    a.SetCalculator(potential)
 
    print  "cube dimension: ",brenner.GetCubeSize()
    print  "stepsize: ",brenner.GetStepSize()

    dyn = PredictorCorrector(a, brenner.GetStepSize(), brenner.GetMovableParam()) 
    #dyn = VelocityVerlet(a, 0.5*femtosecond) 

    dyn.Run(250)
    step = 250;
    etot = a.GetPotentialEnergy() + a.GetKineticEnergy()
    for i in range(2000/250):
        dyn.Run(250)
        step+=250
        print "Total energy after %4.d steps: %.8f"%(step,a.GetPotentialEnergy()+a.GetKineticEnergy())
        potential.WriteXmol(step, brenner.GetStepSize());
        #print i, a.GetPotentialEnergy()/len(a), atoms.GetKineticEnergy()/len(a), (a.GetPotentialEnergy() + a.GetKineticEnergy())/len(a)
    
    forces = a.GetCartesianForces()
    f = []
    j = 0
    for i in range(0, len(a), 100):
        f.append(forces[i])
        j = (j + 1) % 6
    print "f"+"Cu"+" =", repr(f)
    end = time();
    ReportTest("Approximate energy conservation",
           a.GetPotentialEnergy() + a.GetKineticEnergy(),
           etot, 0.3, silent=True)

#Another carbon nano tube test
def nanoTubeTest():
    print "\nLoading configuration with sp2 bonded carbon (nanotube)"
    a = NetCDFTrajectory("sp2carbon-nanotube.nc").GetListOfAtoms()
    a.SetCartesianMomenta(zeros((len(a),3), Float))
    
    potential = BrennerPotential()
    a.SetCalculator(potential)
    
    print "Calculating forces"
    f = a.GetCartesianForces()
    
    print "Loading correct forces"
    target = load(open("sp2carbon-nanotube-forces.pickle"))
    
    print "Comparing forces"
    for i in range(len(a)):
        for j in (0,1,2):
            ReportTest("Force[%d,%d]" % (i,j), f[i,j], target[i,j], 1e-3, silent=True)

    print "Testing energy conservation"

    #dyn = PredictorCorrector(a, 0.2) #bad result for energy conservation
    dyn = VelocityVerlet(a, 0.2*femtosecond)
    dyn.Run(100)
    step=100
    etot = a.GetPotentialEnergy() + a.GetKineticEnergy()
    print "Total energy after %4.d steps: %.8f"%(step,a.GetPotentialEnergy()+a.GetKineticEnergy())
    for i in range(10):
      dyn.Run(100)
      step+=100
      print "Total energy after %4.d steps: %.8f"%(step,a.GetPotentialEnergy()+a.GetKineticEnergy())
    ReportTest("Approximate energy conservation",
           a.GetPotentialEnergy() + a.GetKineticEnergy(),
           etot, 0.3, silent=True)

start = time();
setParamTest()
brennerDynamicsTest()
nanoTubeTest()

end = time();
print "&&&& The test time for the brenner potential was %d seconds &&&&" %(end-start)
ReportTest.Summary()
