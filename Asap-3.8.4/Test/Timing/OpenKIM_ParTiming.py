from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.testtools import ReportTest
from asap3.md.velocitydistribution import MaxwellBoltzmannDistribution
from asap3 import mpi
from ase.data import reference_states
import numpy as np
import time
import sys
import os

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

nblist_candidates = ["NEIGH_RVEC_H", 
                     "NEIGH_RVEC_F", 
                     "NEIGH_PURE_H",
                     "NEIGH_PURE_F",
                     "MI_OPBC_H", 
                     "MI_OPBC_F",
                     ]

ismaster = mpi.world.rank == 0
isparallel = mpi.world.size > 1

if ismaster:
    print_version(1)
    os.system("cat /proc/cpuinfo | grep 'model name'")

T = 1500  # K   - initial T, real T will be the half.
presteps = 10
steps = 20
logsteps = 5

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
        print "Z1 = %i,  Z2 = %i,  a0 = %.5f" % (elem1, elem2, a0)
    # 50*50*50 would be big enough, but some vacancies are nice.
    atoms = FaceCenteredCubic(symbol='Cu', size=(51,51,51))
    nremove = len(atoms) - 500000
    assert nremove > 0
    remove = np.random.choice(len(atoms), nremove, replace=False)
    del atoms[remove]
    if elem1 != elem2:
        z = atoms.get_atomic_numbers()
        z[np.random.choice(len(atoms), len(atoms)/2, replace=False)] = elem2
        atoms.set_atomic_numbers(z)
    if isparallel:
        # Move this contribution into position
        uc = atoms.get_cell()
        x = mpi.world.rank % cpuLayout[0]
        y = (mpi.world.rank // cpuLayout[0]) % cpuLayout[1]
        z = mpi.world.rank // (cpuLayout[0] * cpuLayout[1])
        assert(0 <= x < cpuLayout[0])
        assert(0 <= y < cpuLayout[1])
        assert(0 <= z < cpuLayout[2])
        offset = x * uc[0] + y * uc[1] + z * uc[2]
        new_uc = cpuLayout[0] * uc[0] + cpuLayout[1] * uc[1] + cpuLayout[2] * uc[2]
        atoms.set_cell(new_uc, scale_atoms=False)
        atoms.set_positions(atoms.get_positions() + offset)
        # Distribute atoms. Maybe they are all on the wrong cpu, but that will
        # be taken care of.
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
        print "%s: %.3f us/(atom/cpu)" % (label, t / steps)
    
def RunAll(name, *args):
    initial = MakeAtoms(*args)
    if isparallel:
        atoms = MakeParallelAtoms(initial, cpuLayout)
    else:
        atoms = Atoms(initial)
    emt = EMT()
    emt.set_subtractE0(False)
    atoms.set_calculator(emt)
    RunTiming(atoms, name + ' Native EMT')
    del emt, atoms
    for nbl in nblist_candidates:
        if isparallel:
            atoms = MakeParallelAtoms(initial, cpuLayout)
        else:
            atoms = Atoms(initial)
        atoms.set_calculator(OpenKIMcalculator(openkimmodel, allowed=nbl,
                                               stress=False, stresses=False))
        RunTiming(atoms, name + " " + nbl)
        del atoms
        
if ismaster:
    print "Running on %d processors" % (mpi.world.size,)    

RunAll("CuAu", 29, 79)
RunAll("Cu", 29)
RunAll("Au", 79)
