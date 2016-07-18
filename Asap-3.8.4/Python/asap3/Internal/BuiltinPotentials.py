"""The potentials built into Asap."""
# encoding: utf-8

__docformat__ = "restructuredtext en"

from asap3.Internal.Builtins import _asap
from asap3.Internal.CheckArray import check_parameter_array, _ChkLJarray

import numpy as np
import ase
import ase.calculators.emt
import ase.data
from copy import copy
import sys

EMTParameters = _asap.EMTParameters
EMTDefaultParameters = _asap.EMTDefaultParameters
EMTRasmussenParameters = _asap.EMTRasmussenParameters
#EMTVariableParameters = _asap.EMTVariableParameters
#EMThcpParameters = _asap.EMThcpParameters

BrennerPotential = _asap.BrennerPotential

#IntVector = asap.IntVector
#DoubleVector = asap.DoubleVector


# class FullNeighborList(asap.NeighborList2):
#     "A full neighbor list object."
#     def __init__(self, atoms, *args):
#         asap.NeighborList2.__init__(self, atoms, *args)
#         self.EnableFullNeighborLists()
#         self.atoms = atoms  # Keep the atoms alive.
#         self.CheckAndUpdateNeighborList()

# class NeighborCellLocator(asap.NeighborCellLocator):
#     "A Neighbor list object."
#     def __init__(self, atoms, *args):
#         asap.NeighborCellLocator.__init__(self, atoms, *args)
#         self.atoms = atoms  # Keep the atoms alive.
#         self.CheckAndUpdateNeighborList()
    
# class NeighborList(asap.NeighborList2):
#     "A Neighbor list object."
#     def __init__(self, atoms, *args):
#         asap.NeighborList2.__init__(self, atoms, *args)
#         self.atoms = atoms  # Keep the atoms alive.
#         self.CheckAndUpdateNeighborList()

def get_cell_heights(atoms):
    "Get the heights of the unit cells (perpendicular to the surfaces)"
    heights = np.zeros(3)
    cell = atoms.get_cell()
    for i in range(3):
        direction = np.cross(cell[(i+1) % 3], cell[(i+2) % 3])
        heights[i] = abs(np.dot(cell[i], direction) / np.sqrt(np.dot(direction,direction)))
    return heights

def smallest_pbc_direction(atoms):
    smallest = 1.0e20  # Very big
    pbc = atoms.get_pbc()
    h = get_cell_heights(atoms)
    for i in range(3):
        if pbc[i] and h[i] < smallest:
            smallest = h[i]
    return smallest

class EMT(_asap.EMT):
    """The Effective Medium Theory potential.

    Per default, the EMT potential uses standard parameters for the
    supported metals (Ni, Cu, Pd, Ag, Pt, Au).  If other parameters
    are desired, the relevant parameter provider can be provided as an
    optional argument.
    """
    def __init__(self, parameters=None, minimum_image=None):
        if parameters is None:
            parameters = EMTDefaultParameters()
        _asap.EMT.__init__(self, parameters)
        self.absolute_max_cutoff = parameters.get_max_cutoff_beforeinit()
        self._set_atoms_called = False
        self.use_minimum_image = minimum_image
        
    def set_atoms(self, atoms, *args):
        if not self._set_atoms_called:
            if self.use_minimum_image is None:
                self.use_minimum_image = smallest_pbc_direction(atoms) >= 3 * self.absolute_max_cutoff
            if not self.use_minimum_image:
                sys.stderr.write("Asap.EMT: Disabling minimum image convention.\n")
                self._use_imageatoms()
        self._set_atoms_called = True
        _asap.EMT.set_atoms(self, atoms, *args)

def MonteCarloEMT(parameters = None):
    """The Effective Medium Theory potential optimized for Monte Carlo.

    Per default, the EMT potential uses standard parameters for the
    supported metals (Ni, Cu, Pd, Ag, Pt, Au).  If other parameters
    are desired, the relevant parameter provider can be provided as an
    optional argument.

    This version is optimized for Monte Carlo simulations, at the
    price of somewhat slower ordinary energy and force calculations
    and a larger memory footprint.  In addition, this version cannot
    be parallelized.
    """
    if parameters is None:
        parameters = EMTDefaultParameters()
    return _asap.MonteCarloEMT(parameters)

def EMT2013(parameters, no_new_elements=False):
    """The Effective Medium Theory version 2011 potential.
    
    Parameters must be provided as a dictionary, the keys are
    elements (atomic numbers or strings), the values are dictionaries
    with parameters.  The names of the parameters are 'eta2', 'lambda',
    'kappa', 'E0', 'V0', 'S0' and 'n0'.
    
    If no_new_elements is set to True, the user promises that no
    new elements are introduced after the first time the calculator
    is attached to an Atoms object.  Only elements present at that
    instance will be initialized.
    
    If no_new_elements is False (the default), all elements present
    in the parameters dictionary will be initialized.  This may cause
    a performance penalty as neighbor lists will be large enough to
    accommodate the largest atoms supported.
    """
    try:
        par_ver = parameters['Calculator']
    except KeyError:
        raise ValueError("Dictionary with parameters should contain a" +
                         " Calculator item - be careful not to use an " +
                         "obsolete parameter definition !!")
    if par_ver != "EMT2013_V1":
        raise ValueError("The parameters appear to be for a calculator " +
                         "of type " + par_ver + 
                         " (expected EMT2013_V1).")
    params = {}
    for k in parameters.keys():
        if k == "Calculator":
            continue  # Skip the version specification.
        if isinstance(k, str):
            z = ase.data.atomic_numbers[k]
        else:
            z = k
        params[z] = copy(parameters[k])
        params[z]['mass'] = ase.data.atomic_masses[z]
    return _asap.EMT2013(params, no_new_elements)

def EMT2011(parameters):
    """EMT2011 has been renamed EMT2013, please use the new version."""
    sys.stderr.write("\nASAP Warning: EMT2011 is deprecated, use EMT2013 instead.\n")
    return EMT2013(parameters)

def RGL(elements, p=None, q=None, a=None, xi=None, r0=None,
        cutoff=1.73, delta=0.15, debug=False):
    """
    Calculator based on the RGL/Gupta semi empirical tight-binding potential.

    Parameters
    ----------
    elements: List of the elements that will be supported.

    The parameters can be given here as a dictionary with the form:
        {'Pt': [p, q, a, xi, r0], ('Pt','Y'): [...], ('Y','Y'): [...]}.

    p, q, a, xi, r0 (optional): The potential parameters listed in 2D arrays
    with each dimension equal to the number of elements. Only the upper
    triangular part is used. If there is only one element, then the
    parameters may be numbers.

    cutoff (optional): The distance at which the cutoff function should
    start given in nearest neighor distances. Default is 1.73.

    delta (optional): The length of the cutoff function given in nearest
    neighbor distances. Defauls is 0.15.

    The cutoff defaults are set to include the 3rd and not the 4th nearest
    neighors. The expected nearest neighbor distance is calculated based
    on the parameters.
    """

    # Check if parameters is given as a dictionary
    if isinstance(elements, dict):
        parameters = elements.copy()

        # Get elements
        elements = []
        elementmap = {}
        for key in parameters.keys():
            if isinstance(key, str):
                symbols = [key]
            elif isinstance(key, tuple) and len(key) == 2:
                symbols = list(key)
            else:
                raise KeyError('Key must either be a string or a tuple with ' +
                               'two elements.')
            for s in symbols:
                if not s in elements:
                    elementmap[s] = len(elements)
                    elements.append(s)
        n = len(elements)

        # Get parameters
        p = np.zeros((n, n)) 
        q = np.zeros((n, n)) 
        a = np.zeros((n, n)) 
        xi = np.zeros((n, n)) 
        r0 = np.zeros((n, n)) 

        visit = np.zeros(n)
        for key in parameters.keys():
            if isinstance(key, str):
                i = j = elementmap[key]
            elif isinstance(key, tuple) and len(key) == 2:
                i = elementmap[key[0]]
                j = elementmap[key[1]]
            else:
                raise KeyError('Key must either be a string or a tuple with ' +
                               'two elements.')

            visit[i] += 1
            visit[j] += 1

            p[i,j] = parameters[key][0]
            q[i,j] = parameters[key][1]
            a[i,j] = parameters[key][2]
            xi[i,j] = parameters[key][3]
            r0[i,j] = parameters[key][4]

        if np.any(visit != n + 1):
            raise ValueError('For some elements there are either too few or ' +
                             'too many parameters.')

    # Interpret elements
    if not isinstance(elements, (tuple, list, np.ndarray)):
        elements = [elements]
    n = len(elements)

    for i, symbol in enumerate(elements):
        if isinstance(symbol, str):
            elements[i] = ase.data.atomic_numbers[symbol]
    elements = np.array(elements)

         
    # Interpret parameters
    p = check_parameter_array(n, "p", p)
    q = check_parameter_array(n, "q", q)
    a = check_parameter_array(n, "a", a)
    xi = check_parameter_array(n, "xi", xi)
    r0 = check_parameter_array(n, "r0", r0)

    # Calculate cutoff
    a1 = (np.log(np.sqrt(12) * a * p / (xi * q)) / (p - q) + 1) * r0

    if not np.all((1.0 < a1) & (a1 < 10.0)):
        raise ValueError('Unreasonable parameters - the nearest neighbor ' +
                         'distance span [%.3f, %.3f]' % (a1.min(), a1.max()))

    rcs = cutoff * a1.max()
    rcd = delta * a1.max()
    rce = rcs + rcd
    #print "RGL cutoff: %.4f-%.4f" % (rcs, rcd)

    # Calculate parameters for the cutoff function
    qf = -xi * np.exp(-q * (rcs / r0 - 1)) / rcd**3
    qd = q / r0 * rcd
    q5 = (12.0*qf - 6.0*qf*qd + qf*qd*qd) / (2.0 * rcd**2)
    q4 = (15.0*qf - 7.0*qf*qd + qf*qd*qd) / rcd
    q3 = (20.0*qf - 8.0*qf*qd + qf*qd*qd) / 2.0

    pf = -a * np.exp(-p * (rcs / r0 - 1)) / rcd**3
    pd = p / r0 * rcd
    p5 = (12.0*pf - 6.0*pf*pd + pf*pd*pd) / (2.0 * rcd**2)
    p4 = (15.0*pf - 7.0*pf*pd + pf*pd*pd) / rcd
    p3 = (20.0*pf - 8.0*pf*pd + pf*pd*pd) / 2.0

    if debug:
        print "RGL potential parameters"
        print "p:", p
        print "q:", q
        print "a:", a
        print "xi:", xi
        print "r0:", r0
        print "rcs:", rcs
        print "rce:", rce
        print "q5:", q5
        print "q4:", q4
        print "q3:", q3
        print "p5:", p5
        print "p4:", p4
        print "p3:", p3

    return _asap.RGL(elements, p, q, a, xi * xi, r0, p3, p4, p5,
                     q3, q4, q5, rcs, rce)

Gupta = RGL

def LennardJones(elements, epsilon, sigma, rCut=-1.0, modified=True):
    """Lennard-Jones potential.

    Parameters:
    
    elements:  Lists the elements that will be supported.
    
    epsilon and sigma: The LJ parameters.  2D arrays with each
    dimension equal to the number of elements.  Only the lower
    triangular part is used.  Or a 1D array obtained by flattening the
    corresponding 2D array (this possibility is deprecated, and incurs
    less runtime testing).  If there is only one element, then epsilon
    and sigma may be numbers.

    rCut:  The cutoff distance.  Default: XXX

    modified: Should the potential be shifted so no jump occurs at the
    cutoff.  Default: True.
    """
    try:
        numelements = len(elements)
    except TypeError:
        numelements = 1
        elements = [elements]

    epsilon = _ChkLJarray(epsilon, numelements, "epsilon")
    sigma = _ChkLJarray(sigma, numelements, "sigma")
    masses = [ase.data.atomic_masses[z] for z in elements]
    
    return _asap.LennardJones(numelements, elements, epsilon, sigma,
                              masses, rCut, modified)

def Morse(elements, epsilon, alpha, rmin, rCut=-1.0, modified=True):
    """Lennard-Jones potential.

    Parameters:
    
    elements:  Lists the elements that will be supported.
    
    epsilon, alpha and rmin: The Morse potential parameters.
    2D arrays with each  dimension equal to the number of elements.
    Only the lower triangular part is used.  Or a 1D array obtained by
    flattening the corresponding 2D array (this possibility is deprecated,
    and incurs less runtime testing).  If there is only one element,
    then epsilon and sigma may be numbers.

    rCut:  The cutoff distance.  Default: XXX

    modified: Should the potential be shifted so no jump occurs at the
    cutoff.  Default: True.
    """
    try:
        numelements = len(elements)
    except TypeError:
        numelements = 1
        elements = [elements]
    epsilon = _ChkLJarray(epsilon, numelements, "epsilon", "Morse")
    alpha = _ChkLJarray(alpha, numelements, "alpha", "Morse")
    rmin = _ChkLJarray(rmin, numelements, "rmin", "Morse")

    # The cast below is a hack and should be done properly in C++ (ticket #44).
    return _asap.Morse(elements, epsilon, alpha,
                       rmin, rCut, modified)


# Disable ase.EMT to prevent catastrophic loss of performance.
ase.calculators.emt.EMT.disabled = "Disabled by loading asap3."

# Register if an OpenKIMcalculator is available
OpenKIMsupported = hasattr(_asap, 'OpenKIMcalculator')
