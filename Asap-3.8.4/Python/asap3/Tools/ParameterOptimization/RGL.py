import numpy as np

from asap3 import RGL

def RGLFit(elements, parameters, *args, **kwargs):
    n = len(elements)
    p = np.zeros((n, n))
    q = np.zeros((n, n))
    a = np.zeros((n, n))
    xi = np.zeros((n, n))
    r0 = np.zeros((n, n))

    elementmap = {}
    for i, symbol in enumerate(elements):
        elementmap[symbol] = i

    for key in parameters.keys():
        if len(key) == 1:
            i = j = elementmap[key[0]]
        elif len(key) == 2:
            i = elementmap[key[0]]
            j = elementmap[key[1]]
        else:
            raise KeyError('Number of elements in the key should be one or two.')

        p[i,j] = parameters[key][0]
        q[i,j] = parameters[key][1]
        a[i,j] = parameters[key][2]
        xi[i,j] = parameters[key][3]
        r0[i,j] = parameters[key][4]

        a1 = (np.log(np.sqrt(12) * a[i,j] * p[i,j] / (xi[i,j] * q[i,j])) /
              (p[i,j] - q[i,j]) + 1) * r0[i,j]
        if not 1.0 < a1 < 10.0:
            raise ValueError('Unreasonable parameters: a1 = %.3f' % (a1,))

    return RGL(elements, p, q, a, xi, r0, *args, **kwargs)

if __name__ == '__main__':
    # This script is an example of how to fit RGL parameters.
    import sys

    from asap3.Tools.ParameterOptimization import ParameterOptimization
    from asap3.Tools.OptimizationDatabase import get_data

    initparameters = {('Cu','Cu'): [10.55, 2.43, 0.0894, 1.2799, 2*1.28, 3.61*1.75, 3.61*0.1]}

    varparameters = {('Cu','Cu'): [True, True, True, True, True, False, False]}

    fitquantities = [('lattice_constant_a', 'fcc', 'Cu', get_data('Cu', 'a'), 0.005),
                     ('cohesive_energy', 'fcc', 'Cu', get_data('Cu', 'Ecoh'), 0.01),
                     ('bulk_modulus', 'fcc', 'Cu', get_data('Cu', 'B'), 0.01),
                     ('elastic_constant_C11', 'fcc', 'Cu', get_data('Cu', 'C11'), 0.01),
                     ('elastic_constant_C12', 'fcc', 'Cu', get_data('Cu', 'C12'), 0.01),
                     ('elastic_constant_C44', 'fcc', 'Cu', get_data('Cu', 'C44'), 0.01),
                     ]

    latticeconstants = [['fcc', 'Cu', 4.00]]

    opt = ParameterOptimization(['Cu'], RGLFit, initparameters, varparameters,
                                fitquantities, latticeconstants, debug=False)
    (optparam, optvalues) = opt.fit(log=sys.stdout)

    print "Initial parameters:", initparameters
    print "Optimal parameters:", optparam
    #print ""
    #print "Fitting values", fitquantities
    #print "Optimal values", optvalues
