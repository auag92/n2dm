"""Asap utility functions

This module defines the following functions:
PrintVersion
"""

__docformat__ = "restructuredtext en"


from asap3.Internal.Builtins import _asap, get_version, get_short_version
from asap3 import __file__ as _asapfile
import ase
import sys
import os

def print_version(level = 0):
    """Print the version number of the loaded version of Asap.

    If the optional argument is 1, also prints the pathnames of the
    most important files.
    """
    try:
        compiledfile = _asap.__file__
    except AttributeError:
        compiledfile = "<built-in>"
    print get_version()
    if level >= 1:
        print "  Python module:", _asapfile
        print "  C++ module:   ", compiledfile
        print "  ase module:   ", ase.__file__

def DebugOutput(filename, stdout=1, nomaster=False, sync=True):
    """Debugging output on each node goes to a different file.

    Redirect stderr to a different file on each node.  The filename should
    contain %d, which is replaced by the node number.  The file is opened
    with minimal buffering, and stderr (standard error) is redirected to
    it (also for C/C++ extensions).  If the optional argument stdout is
    true (the default), Python's sys.stdout is also redirected to the
    same file.  Standard output for C/C++ extensions is never touched.

    This is mainly useful for parallel simulations.
    """
    if stdout:
        sys.stdout = sys.stderr
    try:
        import asap3.mpi
        node = asap3.mpi.world.rank
    except (AttributeError, ImportError):
        node = 0
    if nomaster and node == 0:
        return
    flag = os.O_WRONLY|os.O_CREAT|os.O_TRUNC
    if sync:
        flag = flag|os.O_SYNC
    newerror = os.open((filename % (node,)), flag, 0660)
    os.dup2(newerror, sys.stderr.fileno())
    # This Python file must NOT go away.  Attach it to the sys module.
    sys._AsapStandardError = newerror

def print_memory(txt, a=None):
    import asap3.mpi
    procfile = open("/proc/self/status")
    vmsize = vmpeak = vmdata = vmrss = -1
    for line in procfile:
        words = line.split()
        if words[0] == "VmSize:":
            vmsize = int(words[1])
        elif words[0] == "VmPeak:":
            vmpeak = int(words[1])
        elif words[0] == "VmData:":
            vmdata = int(words[1])
        elif words[0] == "VmRSS:":
            vmrss = int(words[1])
    print >>sys.stderr, "Memory [proc %d '%s']: %d MB total (%d MB peak, %d MB data, %d MB rss)" % (
        asap3.mpi.world.rank, txt, (vmsize+512) / 1024,
        (vmpeak+512) / 1024, (vmdata+512) / 1024, (vmrss+512)/1024)        
    procfile.close()
    if a is not None:
        memory_usage(a)

def memory_usage(obj, total=True):
    """Print the memory usage of some kinds of objects.

    Supported objects are: atoms, EMT calculators and neighbor lists.
    """
    mem = 0
    if hasattr(obj, "arrays"):
        mem += _memory_usage_atoms(obj)
        try:
            calc = obj.get_calculator()
        except AttributeError:
            calc = None
        if calc is not None:
            mem += memory_usage(calc, total=False)
    elif hasattr(obj, "print_memory"):
        mem += obj.print_memory()
    else:
        print "*MEM* Memory usage of this object is not supported:", obj
        return 0
    if total:
        print "*MEM* Total %d MB." % (mem,)
    return mem

def _memory_usage_atoms(atoms):
    arr = atoms.arrays
    mem = 0
    nat = len(atoms)
    nvar = 0
    megabyte = 1024*1024
    for k in arr.keys():
        mem += arr[k].size * arr[k].itemsize
        nvar += 1
    gmem = 0
    gvar = 0
    if hasattr(atoms, "ghosts"):
        arr = atoms.ghosts
        for k in arr.keys():
            gmem += arr[k].size * arr[k].itemsize
            gvar += 1
    mem = (mem + gmem + megabyte/2) / megabyte
    gmem = (gmem + megabyte/2) / megabyte
    print "*MEM* Atoms %d MB.  [ %d atoms, %d arrays, %d gh_arr of %d MB ]" % (
        mem, nat, nvar, gvar, gmem)
    return mem
