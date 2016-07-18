#!/usr/bin/python
# Copyright (C) 2007  CAMP
# Please see the accompanying LICENSE file for further information.
#
#
# This file converts the content of a NetCDFTrajectory file
# to a *.d file that can be used as coord.d input file for the
# brenner potential (c-implementation)

from Asap import *
from Asap.Trajectories.NetCDFTrajectory import *
from cPickle import *
from Asap.testtools import ReportTest


filename = "newcoord.d"
f = open(filename,"w")

print "Loading configuration with sp2 bonded carbon (nanotube)"
nc = NetCDFTrajectory("sp2carbon-nanotube.nc")
atoms = nc.GetListOfAtoms()
cube = nc.Get('UnitCell')


r = atoms.GetCartesianPositions()
m = atoms.GetMasses()
num=atoms.GetAtomicNumbers()

f.write(" %s"%"Some title\n" )
f.write( "   %3.d\n"%len(atoms)   )  #number of atoms
f.write( "   %.11E   %.11E\n"%(0.0, 0.5) )  #start time and timestep
f.write( "   %.11E %.11E %.11E\n"%(cube[0][0], cube[1][1], cube[2][2])) #cube size
for i in range(len(atoms)):
  f.write("  %3.d  %3.d  %.11E  %.11E  %.11E 1\n"%((i+1), num[i], r[i][0], r[i][1], r[i][2])) #1 means movable
#acc, forces, dx3dt3 are just put as zeroes
for i in range(3):
  for j in range(len(atoms)):
    f.write("  %3.d  %.11E  %.11E  %.11E\n"%(j+1, 0.0, 0.0, 0.0))


f.close()
