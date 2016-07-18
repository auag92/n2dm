import os
import sys
import time

from string import Template
from asap3.Tools.niflheim import submit_job
from ase.atoms import string2symbols
from asap3.Tools.MaterialProperties import MaterialPropertiesData

def use_template(target, template, args):
    f = open(template, 'r')
    temp = Template(f.read())
    f.close()

    f = open(target, 'w')
    f.write(temp.substitute(args))
    f.close()

def get_parameters(file, number=2):
    file = open(file)
    text = file.read()
    file.close()

    # Find elements
    s = -1
    for i in range(number):
        s = text.find('Optimization', s + 1)
    s = text.find('\n', s) + 1
    e = text.find('parameters', s)
    elements = tuple(string2symbols(text[s:e].strip()))

    # Find parameters
    s = text.find('\n', e) + 1
    e = text.find('Fitting', s) - 4
    parameters = []
    for line in text[s:e].split('\n'):
        rows = line.split('  ')
        parameters.append(float(rows[2]))
    
    return elements, parameters

# Global weighing of alloys - alloy properties are walloy times less important.
#walloy = 100
walloy = 1
dirprefix = 'onlyL12_CuAgAu_'+str(walloy)   # Set to none to get the date as prefix.

# Select parameter.  (A, B, C) means
# A is the processor number (i.e. fit-A.dat)
# B is the directory suffix.
# C is the parameter set within the file.
param_files = {'Ag': (10, '060314', 2),
               #'Al': (5, '151113', 2),
               'Au': (5, '060314', 2),
               'Cu': (13, '060314', 2),
               #'Ni': (0, '151113', 2),
               #'Pd': (6, '151113', 2),
               #'Pt': (1, '181113', 2),
               }
metals = param_files.keys()
metals.sort()

# Name of property in code, structure, name in data file, weight.
temp_metal_prop = [('lattice_constant_a', 'fcc', 'a', 0.001),
                   ('bulk_modulus', 'fcc', 'B', 0.01),
                   ('elastic_anisotropy', 'fcc', 'A', 0.05),
                   #('elastic_constant_C11', 'fcc', 'C11', 0.03),
                   #('elastic_constant_C12', 'fcc', 'C12', 0.03),
                   ('elastic_constant_C44', 'fcc', 'C44', 0.01),
                   ('cohesive_energy', 'fcc', 'Ecoh', 0.001),
                   ('surface_energy', 'fcc111', 'E111', 0.02),
                   #('surface_energy', 'fcc100', 'E100', 0.02),
                   ('surface_ratio', 'fcc111-fcc100', 'E111_100', 0.01),
                   ('stacking_fault', 'fcc', 'Esf', 0.02)
                   ]

# Name of property in code, structure, name in data file, weight, symmetric
# symmetric=True means that the alloy (A,B) = (B,A) and only the alphabetically
# ordered one should be used.
temp_alloy_prop = [('lattice_constant_a', 'l12', 'a', 0.01, False),
                   ('bulk_modulus', 'l12', 'B', 0.01, False),
                   ('heat_of_formation', 'l12', 'Eheat', 0.01, False),
                   #('lattice_constant_a', 'l10', 'a', 0.01, True),
                   #('lattice_ratio_ca', 'l10', 'c/a', 0.01, True),
                   #('heat_of_formation', 'l10', 'Eheat', 0.01, True),
                   ]

mp = MaterialPropertiesData(['properties_metals.dat', 'properties_alloys.dat'])

initparam = ''
varparam = ''
quantities = ''
latticeconst = ''
for i, m1 in enumerate(metals):
    if i > 0:
        initparam += 13 * ' '
        varparam += 12 * ' '
        latticeconst += 20 * ' '

    id, folder, number = param_files[m1]
    paramfile = '%s_%s/fit-%i.dat' % (m1, folder, id)
    e, p = get_parameters(paramfile, number)
    print e, p
    initparam += '%s: %s,\n' % (e, p)
    varparam += '%s: [False, False, True, False, True, False, True],\n' % (e,)

    latticeconst += "('fcc', '%s', %.4f),\n" % (m1, mp.get(m1, 'a'))

    for j, (name, struct, id, weight) in enumerate(temp_metal_prop):
        if id == 'E111_100':
            value = mp.get(m1, 'E111') / mp.get(m1, 'E100')
        else:
            value = mp.get(m1, id)
        if i == 0 and j == 0:
            pass
        else:
            quantities += 14*' '
        quantities += "('%s', '%s', '%s', %.4f, %.5f),\n" % (name, struct, m1,
                                                             value, weight)

    for m2 in metals:
        if m1 == m2:
            continue
        for name, struct, id, weight, symmetric in temp_alloy_prop:
            if symmetric and m1 > m2:
                continue
            value1 = mp.get((m1, m2), id, struct)
            args = (name, struct, (m1, m2), value1, weight * walloy)
            quantities += 14*' ' + "('%s', '%s', %s, %.4f, %.5f),\n" % args
            if name == 'lattice_constant_a' and struct == 'l12':
                a = mp.get((m1, m2), 'a', struct)
                latticeconst += 20*' ' + "('l12', %s, %.4f),\n" % ((m1, m2), a)
            if name == 'lattice_constant_a' and struct == 'l10':
                a = mp.get((m1, m2), 'a', struct)
                c = mp.get((m1, m2), 'c', struct)
                latticeconst += 20*' ' + "('l10', %s, (%.4f, %.4f)),\n" % ((m1, m2), a, c)

if not dirprefix:
    dirname = 'alloy_' + time.strftime('%d%m%y')
else:
    dirname = 'alloy_' + dirprefix
if not os.path.exists(dirname):
    os.mkdir(dirname)

args = {'initparam': initparam,
        'varparam': varparam,
        'quantities': quantities,
        'latticeconstants': latticeconst,
        'symbols': str(metals)}
use_template(dirname + '/fit.py', 'fit_alloy.py', args)

if 'start' in sys.argv:
    os.chdir(dirname)
    jobid =  submit_job(['gpaw-qsub', 'fit.py'])
    print 'Job id:', jobid
    os.chdir('..')

