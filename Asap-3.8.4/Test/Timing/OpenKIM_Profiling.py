from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.testtools import ReportTest
from asap3.md.velocitydistribution import MaxwellBoltzmannDistribution
from asap3 import mpi
from ase.data import reference_states
import numpy as np
import time
import sys

debug = 0
if debug == 1:
    DebugOutput("parallelopenkim%d.log", nomaster=True)
elif debug == 2:
    time.sleep(world.rank)
    print "PID:", os.getpid()
    time.sleep(20)

openkimmodel = "EMT_Asap_Standard_Jacobsen_Stoltze_Norskov_AlAgAuCuNiPdPt__MO_118428466217_002"

possibleLayouts = [(1,1,2), (1,2,2), (2,2,2), (2,2,3), (2,2,4), (2,3,3),
                   (3,3,3), (3,3,4), (3,4,4), (4,4,4), (4,4,5), (4,5,5),
                   (4,6,6), (6,6,6)]

nbl = "NEIGH_RVEC_H" 
useemt = False

ismaster = mpi.world.rank == 0
isparallel = mpi.world.size > 1

if ismaster:
    print_version(1)

T = 1500  # K   - initial T, real T will be the half.
presteps = 1
steps = 2
logsteps = 1

if isparallel:
    layouts = {}
    for l in possibleLayouts:
        layouts[ l[0]*l[1]*l[2] ] = l
    cpuLayout = layouts[mpi.world.size]
    
def MakeAtoms(elem1, elem2=None):
    if elem2 is None:
        elem2 = elem1
    a1 = reference_states[elem1]['a']
    a2 = reference_states[elem2]['a']
    a0 = (0.5 * a1**3 + 0.5 * a2**3)**(1.0/3.0) * 1.03
    if ismaster:
        # 50*50*50 would be big enough, but some vacancies are nice.
        print "Z1 = %i,  Z2 = %i,  a0 = %.5f" % (elem1, elem2, a0)
        atoms = FaceCenteredCubic(symbol='Cu', size=(51,51,51))
        nremove = len(atoms) - 500000
        assert nremove > 0
        remove = np.random.choice(len(atoms), nremove, replace=False)
        del atoms[remove]
        if isparallel:
            atoms = atoms.repeat(cpuLayout)
        if elem1 != elem2:
            z = atoms.get_atomic_numbers()
            z[np.random.choice(len(atoms), len(atoms)/2, replace=False)] = elem2
            atoms.set_atomic_numbers(z)
    else:
        atoms = None
    if isparallel:
        atoms = MakeParallelAtoms(atoms, cpuLayout)
    MaxwellBoltzmannDistribution(atoms, T * units.kB)
    return atoms
        
def RunTiming(atoms, label):
    dyn = VelocityVerlet(atoms, 5.0*units.fs)
    dyn.attach(MDLogger(dyn, atoms, sys.stderr, peratom=True), interval=logsteps)
    dyn.run(presteps)
    t = time.time()
    dyn.run(steps)
    t = time.time() - t
    t *= 2.0   # Only half a million atoms.
    if ismaster:
        print "%s: %.3f us/(atom*cpu)" % (label, t / steps)
    
def RunAll(name, *args):
    initial = MakeAtoms(*args)
    if isparallel:
        atoms = MakeParallelAtoms(initial, cpuLayout)
    else:
        atoms = Atoms(initial)
    if useemt:
        emt = EMT()
        emt.set_subtractE0(False)
        atoms.set_calculator(emt)
        RunTiming(atoms, name + ' Native EMT')
    else:
        atoms.set_calculator(OpenKIMcalculator(openkimmodel, allowed=nbl,
                                               stress=False, stresses=False))
        RunTiming(atoms, name + " " + nbl)
        
if ismaster:
    print "Running on %d processors" % (mpi.world.size,)    

RunAll("CuAu", 29, 79)
RunAll("Cu", 29)
RunAll("Au", 79)
