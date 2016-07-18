# Copyright (C) 2007  CAMP
# Please see the accompanying LICENSE file for further information.

"""Imports a brenner *.d coordinate file and/or a input.d file"""
"""This file is used with the BrennerPotential only"""
"""TODO: Add nice error messages"""

from ASE.ChemicalElements.symbol import symbols
from ASE import Atom

from Asap import *
from Asap import ListOfAtoms

import string

class BrennerImport:

  title = "Default parameter set (overwrites default set of brenner implementation)"
  atoms = [] 
  movable = []  #1=movable; 2=non-movable
  step = 0.5    #femtoseconds
  numAtoms = 0  #number of atoms
  kFlag = 1     #brenner k-flag
  rll = 2.5     #some radius that has quite some influence

  #the following parameters cannot be set in adaption of the brenner potential in asap:
  cubedimension = [200, 200, 200] #calculated automatically if not deactivated in BrennerPotential constructor: therefore ignored
  lennardJonesParameter = []      #parameters for lennardJones: ignored
  starttime = 0                   #starttime: ignored
  nxmol = 500                     #every x steps write output file: ignored: 
                                  #but the user can call potential.WriteXmol()

  femtosecond=1e-15

  def __init__(self):
    """Imports coord.d and input.d"""
  
  #Returns a predictor/corrector object
  #TODO: NOT TESTED!
  def getPredictorCorrector(self, atoms):
        pd = PredictorCorrector(atoms, self.GetStepSize(), self.GetMovableParam())
        return pd

  #Sets all parameter that can be set in the potential (a BrennerPotential)
  #TODO: NOT TESTED - especially !
  def setPotentialParameters(self, potential):
        potential.SetKFlag(self.GetKFlag())
        potential.SetRll(self.GetRll())
        potential.SetMinRangeTree(100)


  #To maintain compatibility with the brenner source and the 
  #c implementation of the brenner source we read it as done
  #in the c-source. (See bren2.c/initkt() for details!)
  #Parameter that don't have a self as prefix are IGNORED!!!!
  def initInput(self, filepath):
        #first line:
        f=open(filepath, 'r')
	[kuc, maxkb, self.kFlag, self.nxmol]= self.ReadInt(f, 4)
        self.ReadLine(f) #ignore the rest of the 1st line
        [pseed, self.rll, temp, pbl]= self.ReadFloat(f, 4)
        self.ReadLine(f) #ignore the rest of the 2nd line
        self.ipot = self.ReadInt(f, 1)[0] #=1 REBO (C,H,Si,Ge), =2 tight-binding for C: 
        self.ReadLine(f) #ignore the rest of the 3rd line
        self.lennardJonesParameter = []
        while 1:
          result = self.ReadFloat(f, 4)
          self.ReadLine(f); #ifgnore the rest of the line
          if result==[]:
             break;
          [a,b,c,d] = result;
          if a>0:
             self.lennardJonesParameter.append(a)
             self.lennardJonesParameter.append(b)
             self.lennardJonesParameter.append(c)
             self.lennardJonesParameter.append(d)
        

  def initCoord(self, filepath):
    """Imports the coordinates (the rest is ignored) of a brenner *.d file:
    # are comments that do not belong to the file:
      
      # textual description:
      diamond lattice       
      # atoms:                                     
      63  
      # start time and stepsize:
      0.20000000000E+01   0.50000000000E+00          
      # cube dimension:            
      0.20000000000E+03   0.20000000000E+03   0.20000000000E+03  
      #number,  element number, x,y,z coordinates, 1=movable 2=unmovable atom
      1    6   0.12239212990E+01  -0.37850461006E+01  -0.76280212402E+00  2 #
      2    6  -0.11044621468E+01  -0.38742597103E+01  -0.95745104551E+00  2
      ...
      63   1   0.72413903475E+00   0.10444240570E+02  -0.26598501205E+01  2

    """
    f=open(filepath, 'r')
    self.title = self.ReadLine(f)
    self.numAtoms = self.ReadInt(f, 1)[0]
    self.starttime, self.step = self.ReadFloat(f, 2)
    self.cubedimension = self.ReadFloat(f, 3)
    self.movable = []
    atoms = []
    for i in range(self.numAtoms):
      atomID = self.ReadInt(f, 1)[0]
      symbol = symbols[(self.ReadInt(f, 1)[0])]
      coord  = self.ReadFloat(f, 3)
      movable= self.ReadInt(f, 1)
      atoms.append(Atom(symbol, coord))
      self.movable.append(movable[0])
    self.atoms = ListOfAtoms(atoms, [[self.cubedimension[0], 0, 0],[0, self.cubedimension[1], 0],[0, 0, self.cubedimension[2]]])
    
  def GetListOfAtoms(self):
   return self.atoms;

  def GetMovableParam(self):
   return self.movable;

  def GetCubeSize(self):
   return self.cubedimension;

  #in FEMTOSECONDS
  def GetStepSize(self):
   return self.step;

  #in FEMTOSECONDS
  def GetStartTime(self):
   return self.starttime;
  
  def GetKFlag(self):
   return self.kFlag

  def GetNxmol(self):
   return self.nxmol

  def GetRll(self):
   return self.rll

  #contains tupels of [natom, xmad, epstd, sigtd]
  def GetLennardJonesParameter(self):
   return self.lennardJonesParameter;


     
  #######################################
  #          PRIVATE FUNCTIONS          #
  #######################################
  #reads a line
  def ReadLine(self, f):
    return f.readline();

  #reads 'numbers' integers 
  def ReadInt(self, f, numbers):
    result = []
    for i in range(numbers):
      number = ''
      ch = f.read(1)
      if ch=="": #EOF
        return result
      while ch in string.whitespace:
        ch = f.read(1)
      if ch=="": #EOF
        return result
      while ch not in string.whitespace:
        number += ch
        ch = f.read(1)
      result.append(int(number))
    return result

  #reads 'numbers' floats 
  def ReadFloat(self, f, numbers):
    result = []
    for i in range(numbers):
      number = ''
      ch = f.read(1)
      if ch=="": #EOF
        return result
      while ch in string.whitespace:
        ch = f.read(1)
      while ch not in string.whitespace:
        number += ch
        ch = f.read(1)
      result.append(float(number))
    return result
