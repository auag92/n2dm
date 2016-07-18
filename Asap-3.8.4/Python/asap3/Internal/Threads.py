"""Handles theads.

In particular checks if you are allowed to use them on a cluster.

This module also disables OpenMP threads initially unless the user
has requested them.
"""
import numpy 
import sys
import os
from asap3.Internal.Builtins import _asap

verbose = True

def AsapThreads(n = None, force=False):
    """Activates the use of threads in calculators that support it.

    The optional argument n can be used to specify the number of
    threads to use.  If not given, it will be set to the number of
    CPUs.

    The optional argument force can be set to true to override the
    security mechanisms on machines where it fails.  IT IS STRICTLY
    FORBIDDEN TO USE THE force PARAMETER ON NIFLHEIM!

    Return value: The number of threads
    """
    if not _asap.support_openmp():
        raise RuntimeError("AsapThreads: OpenMP threads not supported.")
    ncpu = _asap.get_num_procs()
    noriginal = n
    if n is None:
        if ncpu is not None:
            if verbose:
                print "AsapThreads: Choosing", ncpu, "threads."
            n = ncpu
        else:
            raise ValueError, "Cannot determine number of processors.  You must specify the number of threads."

    if n == 1:
        if noriginal == None:
            print >> sys.stderr, "AsapThreads warning: Only one CPU found, disabling threads."
        _asap.set_num_threads(1)
        return 1

    if force:
        print >> sys.stderr, "AsapThreads warning: force specified, hope you know what you are doing.  Enabling", n, "threads."
        _asap.set_num_threads(n)
        return n

    if n > ncpu:
        raise ValueError, ("Cannot start %d threads: only %d CPUs" % (n, ncpu))

    Check_MPI_Consistency(n, ncpu)
    Check_PBS_Allowance(n)
    
    _asap.set_num_threads(n)
    return n


def Check_MPI_Consistency(n, ncpu):
    """Check that we run the same number of threads on all nodes."""
    try:
        from asap3.mpi import world
        worldsize = world.size
        if verbose and worldsize > 1:
            print >>sys.stderr, "MPI world size is", worldsize
    except (ImportError, AttributeError):
        return   # No MPI support
    if worldsize == 1:
        return
    data = numpy.array([n, ncpu])
    world.broadcast(data, 0)
    if n != data[0] or ncpu != data[1]:
        raise RuntimeError, "AsapThreads: Inconsitent CPU layout across the MPI processor pool."

def Check_PBS_Allowance(n):
    "Check that we are not running more threads than allowed by the PBS batch system."
    try:
        from asap3.mpi import world
        worldsize = world.size
    except (ImportError, AttributeError):
        worldsize = 1
        
    try:
        nodefile = os.environ["PBS_NODEFILE"]
        nodelist = open(nodefile).readlines()   
    except (KeyError, IOError):
        return   # No PBS or not on master.

    nodes = {}
    for node in nodelist:
        try:
            nodes[node] += 1
        except KeyError:
            nodes[node] = 1
    
    if len(nodes) != worldsize:
        print >>sys.stderr, "MPI world size:", worldsize
        print >>sys.stderr, "Number of distinct nodes:", len(nodes)
        for node in nodelist:
            print sys.stderr, "   ", node
        raise RuntimeError, "Number of MPI nodes is different from number of distinct nodes in PBS_NODEFILE."
    for x in nodes.keys():
        if nodes[x] < n:
            raise RuntimeError, ("Attempting to run %d threads but %d jobs allowed on node %s" % (n, nodes[x], x))

# When this module is loaded as part of ASAP startup, it should
# disable OpenMP threads, unless the user has specified the
# OMP_NUM_THREADS environment variable.  Also disables processor
# affinity if compiled with a compiler that enables it per default
# (such as Open64 with OpenMP enabled), as affinity will force
# multiple MPI tasks to run on the same cpu while other cpus are
# empty.

if "OMP_NUM_THREADS" not in os.environ or os.environ["OMP_NUM_THREADS"] == "1":
    # We should turn of OpenMP threads until the user eventually enables them.
    if "O64_OMP_SET_AFFINITY" in os.environ:
        # The user controls affinity
        _asap.set_num_threads(1)
    else:
        # Also disable affinity.
        _asap.set_num_threads(1, noaffinity=True)
else:
    print "asap3...Threads: Leaving threads enabled as OMP_NUMTHREADS=%s" % (os.environ["OMP_NUM_THREADS"])
