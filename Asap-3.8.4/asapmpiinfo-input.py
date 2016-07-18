# This file provides information about which MPI implementation was
# used to compile the parallel version of Asap.  The make process
# copies this file to Swig/`arch` and adds a line with the path to mpicc.

import os

def printmpiinfo():
    bindir = os.path.dirname(mpicc)
    mpidir = os.path.dirname(bindir)
    print mpidir


