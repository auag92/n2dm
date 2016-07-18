import sys

from ase.atoms import string2symbols
from asap3 import EMT
from asap3.Tools.ParameterOptimization import ParameterPerformance
from asap3.Tools.ParameterOptimization.EMT import *
from asap3.Tools.MaterialProperties import MaterialPropertiesData

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

param_files = {'Ag': (3, '151113'),
               'Al': (5, '151113'),
               'Au': (2, '151113'),
               'Cu': (5, '151113'),
               'Ni': (0, '151113'),
               'Pd': (6, '151113'),
               'Pt_1': (1, '181113'),
               #'Pt_2': (3, '181113'),
               #'Pt_3': (6, '151113'),
               }

temp_metal_prop = [('lattice_constant_a', 'fcc', 'a', 0.001),
                   ('bulk_modulus', 'fcc', 'B', 0.01),
                   ('elastic_anisotropy', 'fcc', 'A', 0.03),
                   ('elastic_constant_C11', 'fcc', 'C11', 0.03),
                   ('elastic_constant_C12', 'fcc', 'C12', 0.03),
                   ('elastic_constant_C44', 'fcc', 'C44', 0.01),
                   ('cohesive_energy', 'fcc', 'Ecoh', 0.001),
                   ('surface_energy', 'fcc111', 'E111', 0.02),
                   ('surface_energy', 'fcc100', 'E100', 0.02),
                   ('surface_ratio', 'fcc111-fcc100', 'E111_100', 0.01),
                   ('stacking_fault', 'fcc', 'Esf', 0.01),
                   ]

mp = MaterialPropertiesData(['properties_metals.dat', 'properties_alloys.dat'])

for i, (m, (id, folder)) in enumerate(param_files.items()):
    m = m.split('_')[0]
    paramfile = '%s_%s/fit-%i.dat' % (m, folder, id)
    e, p = get_parameters(paramfile)
    parameters = {e: p}
    #parameters = {e: EMTStdParameters(m, 'delta')}
    print paramfile, parameters

    latticeconstants = [('fcc', m, mp.get(m, 'a'))]

    quantities = []
    for j, (name, struct, id, weight) in enumerate(temp_metal_prop):
        if id == 'E111_100':
            value = mp.get(m, 'E111') / mp.get(m, 'E100')
        else:
            value = mp.get(m, id)
        quantities.append((name, struct, m, value, weight))

    #calculator = EMT()
    calculator = EMT2011Fit([m], parameters, 'delta')

    ParameterPerformance(calculator, quantities, latticeconstants, debug=False)

