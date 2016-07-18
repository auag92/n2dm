# Testing various potentials in parallel simulations

from asap3 import *
from asap3.mpi import world
from asap3.testtools import ReportTest
from asap3.GuptaParameters import PtY_parameters as PtY_Gupta
from asap3.EMT2013Parameters import PtY_parameters as PtY_EMT2013
from ase.lattice.cubic import FaceCenteredCubic
from ase.lattice.compounds import L1_2
from ase.visualize import view
from ase.data import atomic_numbers, reference_states

import os
import pickle
import time
import numpy as np

keep = 50   # Keep every 50 data point

ismaster = world.rank == 0
isparallel = world.size > 1
if world.size == 1:
    cpulayout = None
elif world.size == 2:
    cpulayout = [2,1,1]
elif world.size == 3:
    cpulayout = [1,3,1]
elif world.size == 4:
    cpulayout = [2,1,2]
else:
    raise RuntimeError("Not able to run on %i CPUs." % (world.size,))

# This may be run with Test or Test/parallel as working directory.
selfname = 'parallelPotentials.py'
resultname = 'parallelPotentials.pickle'
if not selfname in os.listdir('.'):
    if 'parallel' in os.listdir('.') and selfname in os.listdir('parallel'):
        resultname = os.path.join('parallel', resultname)
    else:
        raise RuntimeError("Cannot find " + selfname)
   
firsttime = (len(sys.argv) >= 2 and sys.argv[1] == '--first') 
if firsttime and os.path.exists(resultname):
    print "This will overwrite the result file '%s'." % (resultname,)
    print "If you don't want to do that press Ctrl-C within 10 seconds."
    time.sleep(10)
    print "Proceeding..."
elif not firsttime and not os.path.exists(resultname):
    print "'%s' not found!" % (resultname,)
    print "Maybe you need to create it with the --first option."
    sys.exit(1)

if firsttime:
    if isparallel:
        print "Cannot create result file (--first option) when running in parallel."
        sys.exit(1)
    result = None
    newresult = {}
else:
    result = pickle.load(open(resultname))
    
    
def myfcc(symbol, pbc):
    "Make an fcc lattice with standard lattice constant"
    if ismaster:
        a = FaceCenteredCubic(symbol=symbol, size=(10,10,10), pbc=pbc)    
        dx = 0.01 * np.sin(0.1 * np.arange(len(a) * 3))
        dx.shape = (len(a),3)
        a.set_positions(a.get_positions() + dx)
        a.set_momenta(np.zeros((len(a),3)))
        #view(a)
    else:
        a = None
    if isparallel:
        a = MakeParallelAtoms(a, cpulayout)
    return a

def myalloy(symbols, pbc):
    "Create an L1_2 alloy (AB3)."
    a1 = reference_states[atomic_numbers[symbols[0]]]['a']
    a2 = reference_states[atomic_numbers[symbols[1]]]['a']
    a0 = (a1 + 3*a2)/4
    if ismaster:
        a = L1_2(symbol=symbols, size=(10,10,10), latticeconstant=a0, pbc=pbc)    
        dx = 0.01 * np.sin(0.1 * np.arange(len(a) * 3))
        dx.shape = (len(a),3)
        a.set_positions(a.get_positions() + dx)
        a.set_momenta(np.zeros((len(a),3)))
        #view(a)
    else:
        a = None
    if isparallel:
        a = MakeParallelAtoms(a, cpulayout)
    return a
    
    
def test(atoms, name, dostress=True):
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    positions = atoms.get_positions()
    if dostress:
        stresses = atoms.get_stresses()
        stress = atoms.get_stress()
    if isparallel:
        ids = atoms.get_ids()
    else:
        ids = np.arange(len(atoms))
    if result:
        # Test that we still get the same
        data = result[name]
        ReportTest(name + ' (energy)', energy, data['energy'], 1e-10)
        compare(name, 'forces', forces, data['forces'], ids)
        compare(name, 'positions', positions, data['positions'], ids)
        if dostress:
            for i in range(6):
                ReportTest(name + '(stress %s)' % (i,), 
                           stress[i], data['stress'][i], 1e-10)
            compare(name, 'stresses', stresses, data['stresses'], ids)
    else:
        # Store new results
        data = {}
        data['energy'] = energy
        data['forces'] = forces[::keep]
        data['positions'] = positions[::keep]
        if dostress:
            data['stress'] = stress
            data['stresses'] = stresses[::keep]
        assert name not in data
        newresult[name] = data
        
def compare(name, quantity, a, b, ids):
    for i in range(len(a)):
        ii = ids[i]
        if ii % keep == 0:
            x = a[i]
            y = b[ii // keep]
            for j in range(len(x)):
                ReportTest('%i %s (%s[%i,%i])' % (world.rank, name, quantity, ii, j),
                           x[j], y[j], 1e-10, silent=True)

ReportTest.ExitAfterError()

# EMT
atoms = myfcc('Cu', pbc=True)
atoms.set_calculator(EMT())
test(atoms, 'Cu EMT periodic')

atoms = myfcc('Cu', pbc=False)
atoms.set_calculator(EMT())
test(atoms, 'Cu EMT free')

atoms = myfcc('Cu', pbc=(1,1,0))
atoms.set_calculator(EMT())
test(atoms, 'Cu EMT mixed')

atoms = myalloy(('Au', 'Cu'), pbc=True)
atoms.set_calculator(EMT())
test(atoms, 'AuCu3 EMT periodic')

atoms = myalloy(('Au', 'Cu'), pbc=False)
atoms.set_calculator(EMT())
test(atoms, 'AuCu3 EMT free')

# Lennard-Jones
atoms = myfcc('Cu', pbc=True)
atoms.set_calculator(LennardJones([29], [0.15], [2.7], 3*2.7, False))
test(atoms, 'LennardJones periodic')

atoms = myfcc('Cu', pbc=False)
atoms.set_calculator(LennardJones([29], [0.15], [2.7], 3*2.7, False))
test(atoms, 'LennardJones free')

atoms = myfcc('Cu', pbc=(1,1,0))
atoms.set_calculator(LennardJones([29], [0.15], [2.7], 3*2.7, False))
test(atoms, 'LennardJones mixed')

# Morse
elements = np.array([atomic_numbers['Ru'], atomic_numbers['Ar']])
epsilon = np.array([[5.720, 0.092], [0.092, 0.008]])
alpha = np.array([[1.475, 2.719], [2.719, 1.472]])
rmin = np.array([[2.110, 2.563], [2.563, 4.185]])
atoms = myfcc('Ar', pbc=True)
atoms.set_calculator(Morse(elements, epsilon, alpha, rmin))
test(atoms, 'Morse periodic', dostress=False)

atoms = myfcc('Ar', pbc=False)
atoms.set_calculator(Morse(elements, epsilon, alpha, rmin))
test(atoms, 'Morse free', dostress=False)

atoms = myfcc('Ar', pbc=(0,1,0))
atoms.set_calculator(Morse(elements, epsilon, alpha, rmin))
test(atoms, 'Morse mixed', dostress=False)

# Gupta
atoms = myalloy(('Y', 'Pt'), pbc=True)
atoms.set_calculator(Gupta(PtY_Gupta))
test(atoms, 'Gupta Pt3Y periodic')

atoms = myalloy(('Y', 'Pt'), pbc=False)
atoms.set_calculator(Gupta(PtY_Gupta))
test(atoms, 'Gupta Pt3Y free')

atoms = myalloy(('Y', 'Pt'), pbc=(0,0,1))
atoms.set_calculator(Gupta(PtY_Gupta))
test(atoms, 'Gupta Pt3Y mixed')

# EMT2013
atoms = myalloy(('Y', 'Pt'), pbc=True)
atoms.set_calculator(EMT2013(PtY_EMT2013))
test(atoms, 'EMT2013 Pt3Y periodic')

atoms = myalloy(('Y', 'Pt'), pbc=False)
atoms.set_calculator(EMT2013(PtY_EMT2013))
test(atoms, 'EMT2013 Pt3Y free')

atoms = myalloy(('Y', 'Pt'), pbc=(0,0,1))
atoms.set_calculator(EMT2013(PtY_EMT2013))
test(atoms, 'EMT2013 Pt3Y mixed')



if firsttime:
    print "Writing data to", resultname
    out = open(resultname, 'w')
    pickle.dump(newresult, out, -1)
    out.close()
    
ReportTest.Summary()
