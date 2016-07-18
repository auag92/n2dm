#! /usr/bin/env python

"""Usage: AsapFileToTrajectory.py oldfile newfile [frame1 [frame2]]

Converts an Asap version 1.x netcdf file to an ASE version 2
trajectory file.  If a single frame number is given, only that frame
is converted.  If two frame numbers are given, all frame from the
first to the last are converted (-1 means the last frame).  If no
frame numbers are given, all frames are converted.

"""

import Numeric as num
from Scientific.IO.NetCDF import NetCDFFile
import sys
import types

old_names = {
    'cartesianPositions': 'CartesianPositions',
    'cartesianMomenta': 'CartesianMomenta',
    'basisVectors': 'UnitCell',
    'classes': 'Tags',
    'atomicNumbers': 'AtomicNumbers',
    'periodic': 'BoundaryConditions'
    }

new_names = {
    #    name                 shape        typecode    once    units
    # -----------------------------------------------------------------
    'CartesianVelocities': (('natoms', 3), num.Float, False, (2, -0.5)),
    'CartesianPositions':  (('natoms', 3), num.Float, False, (1, 0)),
    'CartesianMomenta':    (('natoms', 3), num.Float, False, (2, -0.5)),
    'CartesianForces':     (('natoms', 3), num.Float, False, (-1, 1)),
    'Stress':              ((3, 3),        num.Float, False, (-3, 1)),
    'UnitCell':            ((3, 3),        num.Float, False, (1, 0)),
    'BoundaryConditions':  ((3,),          num.Int,   True , (0, 0)),
    'PotentialEnergy':     ((),            num.Float, False, (0, 1)),
    'AtomicNumbers':       (('natoms',),   num.Int,   True,  (0, 0)),
    'MagneticMoments':     (('natoms',),   num.Float, True,  (0, 0)),
    'Tags':                (('natoms',),   num.Int,   True,  (0, 0))}

def normalize(fr, nfr, default):
    if fr == None:
        fr = default
    if fr < 0:
        fr = nfr + fr
    if fr < 0:
        raise ValueError, "Frame number before beginning of file."
    if fr > nfr:
        raise ValueError, "Frame number after end of file: " + str(fr) + " > " + str(nfr)
    return fr

def AsapFileToTrajectory(oldfile, newfile, firstframe=None, lastframe=None):
    # Check if input file is a filename or a NetCDF file
    if isinstance(oldfile, types.StringTypes):
        oldfile = NetCDFFile(oldfile)

    pos = oldfile.variables['cartesianPositions']  # Must be present
    (nframes, natoms, three) = pos.shape
    print natoms, three, nframes
    firstframe = normalize(firstframe, nframes, 0)
    lastframe = normalize(lastframe, nframes, -1)
    if lastframe < firstframe:
        raise ValueError, "No frames to copy, giving up."

    print "Preparing to copy frames", firstframe, "to", lastframe
    # Now open the output file, and define the variables.
    if isinstance(newfile, types.StringTypes):
        newfile = NetCDFFile(newfile, "w")
    oncevars = []
    manyvars = []
    for v in oldfile.variables.keys():
        try:
            newname = old_names[v]
        except KeyError:
            print "WARNING: Skipping data named", v
            continue
        if new_names[newname][2]:
            shape = new_names[newname][0]
            oncevars.append((v, newname))
        else:
            shape = ("unlim",) + new_names[newname][0]
            manyvars.append((v, newname))
        shape2 = []
        for d in shape:
            if isinstance(d, types.IntType):
                n = d
                d = str(d)
            elif d == 'natoms':
                n = natoms
            elif d == 'unlim':
                n = None
            else:
                raise RuntimeError, "Unknown dimension "+str(d)
            if not newfile.dimensions.has_key(d):
                newfile.createDimension(d, n)
            shape2.append(d)
        print v, "-->", newname, " shape", shape2
        var = newfile.createVariable(newname, oldfile.variables[v].typecode(),
                                     tuple(shape2))
        var.once = new_names[newname][2]
        var.units = new_names[newname][3]
        
    # Now copy the data
    print "Copying global data"
    newfile.history = 'ASE trajectory'
    newfile.version = '0.1'
    newfile.lengthunit = 'Ang'
    newfile.energyunit = 'eV'
    for oldname, newname in oncevars:
        newfile.variables[newname][:] = oldfile.variables[oldname][:]
    
    for n in range(firstframe, lastframe+1):
        print "Copying frame", n
        for oldname, newname in manyvars:
            newfile.variables[newname][n] = oldfile.variables[oldname][n]
    newfile.close()
    
if __name__ == "__main__":
    # sys.tracebacklimit = 0
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        print __doc__
        raise TypeError, "Wrong number of arguments."
    infile = sys.argv[1]
    outfile = sys.argv[2]
    first = last = None
    if len(sys.argv) > 3:
        first = last = int(sys.argv[3])
    if len(sys.argv) > 4:
        last = int(sys.argv[4])
    AsapFileToTrajectory(infile, outfile, first, last)
    
