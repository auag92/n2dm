import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt

from asap3.Tools.MaterialProperties import MaterialPropertiesData
from asap3.Tools.ParameterOptimization import ParameterPerformance
from asap3.Tools.ParameterOptimization.Optimization import ParameterOptimization
from asap3.Tools.ParameterOptimization.SearchParallel import ParameterSearch
from asap3.Tools.ParameterOptimization.EMT import EMT2011Fit, EMTStdParameters
from asap3.Tools.ParameterOptimization.GetParameters import get_parameters

paramdirs = sys.argv[1:]

#paramfile = sys.argv[1]


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
                   ('lattice_constant_a', 'l10', 'a', 0.01, True),
                   ('lattice_ratio_ca', 'l10', 'c/a', 0.01, True),
                   ('heat_of_formation', 'l10', 'Eheat', 0.01, True),
                   ]

mp = MaterialPropertiesData(['properties_metals.dat', 'properties_alloys.dat'])

paramfiles = []
for dir in paramdirs:
    for line in open(os.path.join(dir, 'fit.dat')):
        words = line.split()
        if words and words[0] == "0":
            paramfiles.append(os.path.join(dir, "fit-"+words[1]+".dat"))
            break

quantities = []
latticeconst = []

metals = []
initparam = {}
errfuncs = {}

for paramfile in paramfiles:
    param, errf = get_parameters(paramfile, 2)
    for k,v in param.iteritems():
        if k not in errfuncs or errf < errfuncs[k]:
            initparam[k] = v
            errfuncs[k] = errf
            print "Adding from %s: %s" % (paramfile, str(k))
    
for m1, m2 in initparam.keys():
    assert m1 == m2
    metals.append(m1)
metals.sort()
print metals

walloy = 1   # Without any effect.

for m1 in metals:
    latticeconst.append(('fcc', m1, mp.get(m1, 'a')))
    for name, struct, id, weight in temp_metal_prop:
        if id == 'E111_100':
            value = mp.get(m1, 'E111') / mp.get(m1, 'E100')
        else:
            value = mp.get(m1, id)
        quantities.append((name, struct, m1, value, weight))

    for m2 in metals:
        if m1 == m2:
            continue
        for name, struct, id, weight, symmetric in temp_alloy_prop:
            if symmetric and m1 > m2:
                continue
            value1 = mp.get((m1, m2), id, struct)
            args = (name, struct, (m1, m2), value1, weight * walloy)
            quantities.append(args)
            if name == 'lattice_constant_a' and struct == 'l12':
                a = mp.get((m1, m2), 'a', struct)
                latticeconst.append(('l12', (m1, m2), a))
            if name == 'lattice_constant_a' and struct == 'l10':
                a = mp.get((m1, m2), 'a', struct)
                c = mp.get((m1, m2), 'c', struct)
                latticeconst.append(('l10', (m1, m2), (a, c)))

#print latticeconst
#print quantities


calculator = EMT2011Fit(metals, initparam, 'delta')
result = ParameterPerformance(calculator, quantities, latticeconst, debug=False, noprint=True)

def makeplot(data, label, prop, structures=None, minimum=None, maximum=None):
    if structures is None:
        structures = {None: ('o', 'b')}
    xs = []
    ys = []
    fig = plt.figure()
    for wantedstruc, plottype in structures.iteritems():
        x = []
        y = []
        for name, struc, elements, target, value in result:
            #print name, struc, elements, target, value
            if name == prop:
                keep = (wantedstruc is None) or (struc == wantedstruc)
                if keep:
                    x.append(target)
                    y.append(value)
        m, c = plottype
        print "Found %i data points of type %s, %s" % (len(x), prop, wantedstruc)
        if x:
            plt.scatter(x, y, c=c, marker=m, label=wantedstruc)
            xs.extend(x)
            ys.extend(y)
    if minimum is None:
        minimum = min(xs + ys)
    if maximum is None:
        maximum = max(xs + ys)
    plt.plot([minimum, maximum], [minimum, maximum], color='black')
    plt.axes().set_aspect('equal')
    plt.title(label)
    plt.legend(loc='best')

style = {'fcc': ('o', 'r'),
         'l12': ('o', 'b'),
         'l10': ('o', 'y')}
makeplot(result, 'Lattice constants (a)', 'lattice_constant_a', minimum=3.4, maximum=4.3,
         structures=style)
makeplot(result, 'Bulk modulus (B)', 'bulk_modulus',
         structures=style)
makeplot(result, 'Shear modulus (C44)', 'elastic_constant_C44',
         structures=style)
makeplot(result, 'Heat of formation', 'heat_of_formation',
         structures=style)
plt.show()
