from math import pi, sqrt

from ase.data import atomic_numbers

from asap3 import EMT, EMT2013
from asap3.Internal.EMTParameters import _std_emt_params

beta = (16.0 * pi / 3.0)**(1.0 / 3.0) / sqrt(2.0)

def EMTFit(elements, parameters):
    """Create an EMT object from a list (not dict!) of parameters."""
    #print "EMTFit called with:"
    #print "  elements =", elements
    #print "  parameters =", parameters
    parameternames = ['eta2', 'lambda', 'kappa', 'E0', 'V0', 'S0', 'n0']
    param = {}
    for k, v in parameters.items():
        if isinstance(k, tuple):
            k, k2 = k
            assert k == k2
        p = {}
        for name, value in zip(parameternames, v):
            p[name] = value
        # Check for resonable parameters
        #if p['eta2'] * 1.809 - p['kappa'] < 0.01:
        #    return None
        param[k] = p
    #print "  new parameters:", param
    return EMT(param)

def EMT2013Fit(elements, parameters, order='kappa'):
    """Create an EMT2013 object from a list (not dict!) of parameters."""
    #print "EMT2013Fit called with:"
    #print "  elements =", elements
    #print "  parameters =", parameters
    if order == 'kappa':
        parameternames = ['eta2', 'lambda', 'kappa', 'E0', 'V0', 'S0', 'n0']
    elif order == 'delta':
        parameternames = ['eta2', 'lambda', 'delta', 'E0', 'V0', 'S0', 'n0']
    else:
        raise ValueError('"order" must be either "kappa" or "delta"')

    V0max = 12.0
    param = {'Calculator': "EMT2013_V1"}
    for k, v in parameters.items():
        if isinstance(k, tuple):
            k, k2 = k
            assert k == k2

        # Copy parameters
        p = {}
        for name, value in zip(parameternames, v):
            if name != 'E0' and value <= 0.0:
                raise ValueError('Parameter %s <= 0.0' % (name,))
            elif name == 'E0' and value >= 0.0:
                raise ValueError('Parameter E0 >= 0.0')
            p[name] = value

        # Convert 'delta' to 'kappa'
        if order == 'delta':
            p['kappa'] = beta * p['eta2'] - p['delta']
            del p['delta']

        # Check for resonable parameters
        if beta * p['eta2'] - p['kappa'] < 0.01:
            raise ValueError('beta * eta2 - kappa < 0.01')
        if p['V0'] < 0.1:
            raise ValueError('V0 too small (less than 0.1)')
        if p['V0'] > V0max:
            raise ValueError('V0 too large (greater than %.1f)' % (V0max,))

        param[k] = p
    #print "  new parameters:", param
    return EMT2013(param)

EMT2011Fit = EMT2013Fit

def EMTStdParameters(z, order='kappa'):
    if isinstance(z, str):
        z = atomic_numbers[z]

    parameternames = ['eta2', 'lambda', 'kappa', 'E0', 'V0', 'S0', 'n0']
    parameters = [_std_emt_params[z][x] for x in parameternames]

    if order == 'delta':
        parameters[2] = beta * parameters[0] - parameters[2]

    return parameters

if __name__ == "__main__":
    # This script is an example of how to fit EMT parameters.
    import sys
    from asap3.Tools.ParameterOptimization import ParameterOptimization
    from asap3.Internal.EMTParameters import _std_emt_params
    from asap3.Tools.OptimizationDatabase import get_data

    initparameters = {('Cu','Cu'): EMTStdParameters('Cu')}

    varparameters = {('Cu','Cu'): [True, True, True, True, True, True, False]}

    # This should be generated semi-automatically
    fitquantities = [('lattice_constant_a', 'fcc', 'Cu', get_data('Cu', 'a'), 0.005),
                     ('cohesive_energy', 'fcc', 'Cu', get_data('Cu', 'Ecoh'), 0.01),
                     ('bulk_modulus', 'fcc', 'Cu', get_data('Cu', 'B'), 0.01),
                     ('elastic_constant_C11', 'fcc', 'Cu', get_data('Cu', 'C11'), 0.01),
                     ('elastic_constant_C12', 'fcc', 'Cu', get_data('Cu', 'C12'), 0.01),
                     ('elastic_constant_C44', 'fcc', 'Cu', get_data('Cu', 'C44'), 0.01),
                     ]

    # Initial guesses?
    latticeconstants = [['fcc', 'Cu', 4.00],  
                        ['bcc', 'Cu', 4.00]]

    opt = ParameterOptimization(['Cu'], EMT2013Fit, initparameters, varparameters,
                                fitquantities, latticeconstants, debug=False)
    (optparam, optvalues) = opt.fit(log=sys.stdout)

    print "Initial parameters:", initparameters
    print "Optimal parameters:", optparam
    #print ""
    #print "Fitting values", fitquantities
    #print "Optimal values", optvalues

