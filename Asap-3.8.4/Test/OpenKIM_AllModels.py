from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.testtools import ReportTest
import numpy as np
from OpenKIM_modelname import openkimmodels
from numpy.random import RandomState
from ase.data import reference_states, atomic_numbers
from ase.lattice import bulk
import sys

brokenmodels = ['Dipole_Umeno_YSZ__MO_394669891912_001',  # Crash with Intel compiler, some consistency problems
                'IMD_EAM_Schopf_',   # Crash.  Also fails vc_forces_numer_deriv
                'Pair_Morse_Shifted_GirifalcoWeizer_HighCutoff_St__MO_497591319122_001',  # No such element.  St ????
                'Pair_Morse_Shifted_GirifalcoWeizer_LowCutoff_St__MO_801083489225_001',   # No such element.  St ????
                'Pair_Morse_Shifted_GirifalcoWeizer_MedCutoff_St__MO_964297938209_001',   # No such element.  St ????
                'kcc_meam_LiSi__MO_596436139350_001',  # Crash if simulator supports particleEnergies.
                'EMT_Asap_MetalGlass_CuMgZr__MO_655725647552_001',   # Crash - known to be broken.
                'EMT_Asap_Standard_Jacobsen_Stoltze_Norskov_AlAgAuCuNiPdPt__MO_118428466217_001', # do.
                'model_ArCHHeXe_BOP_AIREBO__MO_154399806462_001',  # Exception.  Too many particles.  FORTRAN
                'MEAM_2NN_Fe_to_Ga__MO_145522277939_001',  # Crash if simulator supports particleEnergies.
                'MEAM_2NN_GaInN__MO_117938381510_001',     # Crash if simulator supports particleEnergies.
                'Glue_Ercolessi_Adams_Al__MO_324507536345_001',  # Segmentation fault.  FORTRAN
               ]
    
tolerance = 1e-4
#set_verbose(1)

if not OpenKIMsupported:
    openkimmodels = []
    print "OpenKIM support is not compiled into Asap."
    
#rnd = RandomState(42)  # We want deterministic random numbers
rnd = RandomState()

known_states = ['bcc', 'fcc', 'hcp', 'diamond', 'sc']

delta = 0.01

openkimmodels = openkimmodels + ['EMT']

if len(sys.argv) > 1:
    openkimmodels = sys.argv[1:]
    brokenmodels = []

skipped = []
crashed = []
for model in openkimmodels:
    broken = False
    for brk in brokenmodels:
        if model.startswith(brk):
            broken = True
    if broken:
        print "\nSkipping broken KIM model:", model
        skipped.append(model)
        continue
    print "\nKIM model:", model
    if model == 'EMT':
        elements = ('Cu', 'Au', 'Ag')
        nblists = [None]
        access = [None]
    else:
        info = OpenKIMinfo(model)
        elements = info.get_supported_types()
        nblists = info.get_supported_nblist()
        access = info.get_supported_access()
    print "Supported elements:", elements
    print "Supported neighborlist methods:", nblists
    print "Supported access methods:", access
    if len(access) > 1:
        access = ['loca', 'iter']
    else:
        access = [None]
    if len(elements) == 1:
        main = elements[0]
        other = None
        state = reference_states[atomic_numbers[main]]
    else:
        elements = list(elements)
        rnd.shuffle(elements)
        for i in range(len(elements)):
            main = elements[i]
            other = elements[i-1]
            state = reference_states[atomic_numbers[main]]
            if state['symmetry'] in known_states:
                break
    if state['symmetry'] not in known_states:
        print "Cannot simulate %s, reference state '%s' not supported" % (main, state['symmetry'])
        print "SKIPPING MODEL!"
        continue
    
    init_atoms = bulk(main, orthorhombic=True).repeat((7,7,7))
    r = init_atoms.get_positions()
    r += rnd.normal(0.0, 0.1, r.shape)
    init_atoms.set_positions(r)
    z = init_atoms.get_atomic_numbers()
    if other:
        some_atoms = rnd.randint(0, 20, len(init_atoms)) == 0
        z[some_atoms] = atomic_numbers[other]
        init_atoms.set_atomic_numbers(z)
        z_other = atomic_numbers[other]
    else:
        z_other = 0
    print ("Generated a %s system with %i %s-atoms and %i %s-atoms"
           % (state['symmetry'], 
              np.equal(z, atomic_numbers[main]).sum(),
              main,
              np.equal(z, z_other).sum(),
              other))
    print "Lattice constant a =", state['a'] 
    old_energy = old_forces = None
    rndat = rnd.randint(len(init_atoms))
    for nbl in nblists:
        for ac in access:
            #print "Testing %s with %s" % (model, nbl)
            atoms = Atoms(init_atoms)
            if nbl == 'CLUSTER':
                atoms.set_pbc(False)
            try:
                if model == 'EMT':
                    atoms.set_calculator(EMT())
                else:
                    atoms.set_calculator(OpenKIMcalculator(model, allowed=nbl, access=ac))
                e = atoms.get_potential_energy()
                f = atoms.get_forces()
            except (RuntimeError, ValueError, AsapError) as exc:
                txt = ("%s with %s(%s) raised exception %s: %s" %
                       (model, nbl, str(ac), str(type(exc)), str(exc)))
                print txt
                crashed.append(txt)
                continue
            if nbl != 'CLUSTER':
                if old_energy == None:
                    old_energy = e
                    old_forces = f
                else:
                    ReportTest("%s with %s(%s); energy unchanged" % (model, nbl, str(ac)),
                               e, old_energy, 1e-6)
                    ReportTest("%s with %s(%s); force unchanged" % (model, nbl, str(ac)),
                               f[rndat,0], old_forces[rndat,0], 1e-6)                    
            atoms[rndat].position[0] += delta
            de = atoms.get_potential_energy() - e
            f = 0.5 * f + 0.5 * atoms.get_forces()
            exp_de = -delta * f[rndat,0]
            #print "Old energy: %.9f.   Change: %.9f    Expected: %.9f   Abs: %.9e   Relative: %.9f" % (e, de, exp_de, de-exp_de, (de-exp_de)/exp_de)
            ReportTest("%s with %s(%s) force consistent" % (model, nbl, str(ac)), de, exp_de, tolerance)
        
print "Skipped models (registerd as crashing):"
for m in skipped:
    print "   ", m
print "Crashing models:"
for m in crashed:
    print "   ", m
ReportTest.Summary()
