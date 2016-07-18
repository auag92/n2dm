# Reading old ASAP files (version 1.X NetCDF files).

"""Asap.IO.ReadOldFiles reads NetCDF files from ASAP version 1.X.

The following functions are defined:

`Reader(filename)`:
An object for reading old-style (ASAP 1.x) NetCDF files.

`Reader.Read(frame = None)`:
Read a given frame into an ASAP ListOfAtoms object.  The
default frame is the one specifies by .Goto or in the
constructor.
   
`Read(filename, frame=-1)`:
    Reads a NetCDF file containing a simulation, returning a
    ListOfAtoms object.  Mainly for serial simulations.

"""

__docformat__ = "restructuredtext en"

import Asap as _Asap
from Scientific.IO.NetCDF import NetCDFFile
import Numeric

class Reader:
    """Reads configurations from a NetCDF file (ASAP 1.x format).

Read configuration number "N" from a NetCDF file:    
>>> atoms = Reader("stuff.nc").Read(N)

Read last configuration from a NetCDF file:
>>> atoms = Reader("stuff.nc").Read()

Read several configurations:
>>> r = Reader("things,.nc")
>>> atoms1 = r.Read(0)
>>> atoms2 = r.Read(1)
"""
    def __init__(self, filename, whatToRead = None):
        """Construct a reader from a filename."""        
        self.whatToRead = whatToRead
        try:
            self.nc = NetCDFFile(filename, 'r')
        except IOError:
            self.nc = NetCDFFile(filename + ".nc", 'r')
        
    def __del__(self):
        try:
            self.nc.close()
        except IOError:
            pass

    def Close(self):
        """Close the associated NetCDF file."""
        self.nc.close()

    def GetNetCDFFile(self):
        """Get the associated NetCDF file handle."""
        return self.nc

    def Read(self, frame=-1):
        """Reads a frame, returning a ListOfAtoms.

        Per default it reads the last frame, but other frames can be
        specified.
        """
        vars = self.nc.variables
        positions = vars["cartesianPositions"][frame]
        cell = vars["basisVectors"][frame]
        peri = vars["periodic"][:]
        atoms = _Asap.ListOfAtoms(positions=positions, cell=cell,
                                  periodic=peri)
        if self.whatToRead is None:
            self.whatToRead = []
            for stuff in ["CartesianMomenta", "Classes", "AtomicNumbers"]:
                if vars.has_key(stuff[0].lower() + stuff[1:]):
                    self.whatToRead.append(stuff)
        for stuff in self.whatToRead:
            if stuff == "CartesianMomenta":
                atoms.SetCartesianMomenta(vars["cartesianMomenta"][frame].astype("d"))
            elif stuff == "Classes":
                atoms.SetTags(vars["classes"][frame])
            elif stuff == "AtomicNumbers":
                if len(vars["atomicNumbers"]) == 1:
                    atoms.SetAtomicNumbers(vars["atomicNumbers"][0] +
                                           Numeric.zeros(len(atoms)))
                else:
                    atoms.SetAtomicNumbers(vars["atomicNumbers"][:])
            else:
                raise RuntimeError, "Don't know how to read " + stuff
        return atoms

def Read(filename, n = -1, optional = None):
    return Reader(filename, optional).Read(n)
