import numpy as np
# Delayed import:  scipy.optimize
# It is imported in the minimization function, so the main functionality
# of this module is preserved even if SciPy is not installed.

def ElasticConstants(atoms, symmetry, minimize=False, **kwargs):
    """Calculate the anisotropic elastic constants of a system.
    
    The atoms must be a periodic system with an orthorhombic cell,
    and the lattice constant(s) should be the equilibrium value(s).  
    If not, the function will crash, unless called with 
    minimize=True to perform an (inefficient) minimization first.
    
    The atoms are modified, but returned to the original state, 
    except if minimize=True in which case the unit cell is optimized.
    
    The symmetry of the system must be specified.  Currently, 
    the following symmetries are allowed::
    
    symmetry='cubic'
        Returns C11, C12, C44 and B.
        Note that B = (C11 + 2*C12)/3
        
    symmetry='hexagonal'
        Returns C11, C12, C13, C33, C44 and B.
        Note that B is given by the other constants (but the
        expression is not trivial)
        
    Note that the system must reflect the underlying symmetry 
    (not checked), and that the unit cell must be orthorombic 
    (may be relaxed later for hexagonal systems).
    """   
    ucell = atoms.get_cell()
    for i in range(3):
        for j in range(3):
            if i != j and abs(ucell[i,j]) > 1e-10:
                raise ValueError("Unit cell of atoms is not orthorhombic.")     
    if minimize:
        minimize_unit_cell(atoms, symmetry)
    if symmetry.lower() == 'cubic':
        return _elastic_constants_cubic(atoms, **kwargs)
    elif symmetry.lower() == 'hexagonal':
        return _elastic_constants_hexagonal(atoms, **kwargs)
    else:
        raise ValueError('Symmetry "%s" not supported' % (symmetry,))
    
def _elastic_constants_cubic(atoms, debug=False):
    """Calculate C11, C12, C44 and B for a system with cubic symmetry."""
    e0 = atoms.get_potential_energy()
    u0 = atoms.get_cell()
    C11 = elastic_constant('C11', atoms, e0, u0, (1,0,0,0,0,0),
                           debug=debug)
    B = elastic_constant('B', atoms, e0, u0, 
                         (1/3.0, 1/3.0, 1/3.0, 0, 0, 0), 
                         debug=debug)
    C44 = elastic_constant('C44', atoms, e0, u0, 
                           (0,0,0,1/2.0,0,0), debug=debug)
    C12 = (3*B - C11)/2.0
    return (C11, C12, C44, B)

def _elastic_constants_hexagonal(atoms, debug=False, sanity=True):
    """Calculate C11, C12, C13, C33, C44 and B for a hexagonal system."""
    e0 = atoms.get_potential_energy()
    u0 = atoms.get_cell()
    C11 = elastic_constant('C11', atoms, e0, u0, (1,0,0,0,0,0),
                           debug=debug)
    C33 = elastic_constant('C33', atoms, e0, u0, (0,0,1,0,0,0),
                           debug=debug)
    C11_C12 = elastic_constant('2C11+2C12', atoms, e0, u0,
                               (1,1,0,0,0,0), debug=debug) / 2
    C12 = C11_C12 - C11
    C11_C33_2C13 = elastic_constant('C11+C33+2C13', atoms, e0, u0,
                                    (1,0,1,0,0,0), debug=debug)
    C13 = 0.5 * (C11_C33_2C13 - C11 - C33)
    C44 = elastic_constant('C44', atoms, e0, u0, 
                           (0,0,0,1/2.0,0,0), debug=debug)
    C66 = elastic_constant('C66', atoms, e0, u0, 
                           (0,0,0,0,0,1/2.0), debug=debug)
    C66_alt = 0.5 * (C11 - C12)
    devC66 = 2 * abs(C66 - C66_alt) / (C66 + C66_alt) 
    if debug or devC66 > 0.02:
        print "Deviation in C66: %.1f%%" % (100 * devC66,)
    if devC66 > 0.02 and sanity:
        raise RuntimeError("Deviation in C66: %.1f%% > 2%%: Not a hexagonal system?"
                           % (100 * devC66,))
    alpha = (C11 + C12 - 2*C13) / (C33 - C13)
    if debug:
        print "alpha =", alpha
    B = (2 * C11 + 2 * C12 + 4 * alpha * C13 
         + C33 * alpha * alpha) / ((2 + alpha) * (2 + alpha))
    return (C11, C12, C13, C33, C44, B)

def elastic_constant(name, atoms, e0, u0, mode, strain=0.007, debug=False):
    """Calculate an elastic constant.
    
    Parameters::
    
    name
        Name of the constant (a string). 
        
    atoms
        The atoms.
        
    e0
        Energy in the equilibrium configuration.
        
    u0
        Unit cell in the equilibrium configuration.
        
    mode
        The deformation mode as six numbers, giving weights
        to the six strain components.
    """
    strainlist = (-1.0, 0.5, 0, 0.5, 1.0)
    energies = []
    order = np.array([[0, 5, 4],
                      [5, 1, 3],
                      [4, 3, 2]])
    straintensor = np.array(mode)[order]
    if debug:
        print "%s unit strain tensor:" % (name,)
        print straintensor
    for s in strainlist:
        if s == 0:
            energies.append(e0)
        else:
            # Apply strain
            s = s * strain * straintensor + np.identity(3)
            atoms.set_cell(np.dot(u0, s), scale_atoms=True)
            energies.append(atoms.get_potential_energy())
    atoms.set_cell(u0, scale_atoms=True)  # Reset atoms
    fit0 = np.poly1d(np.polyfit(strain * np.array(strainlist), energies, 3))
    fit1 = np.polyder(fit0, 1)
    fit2 = np.polyder(fit1, 1)
    x0 = None
    for x in np.roots(fit1):
        if fit2(x) > 0:
            if x0 is not None:
                raise RuntimeError("More than two roots found.")
            assert x0 is None
            x0 = x
    if x0 is None:
        raise RuntimeError("No roots found.")
    if np.abs(x0) > 0.5 * strain:
        raise RuntimeError("Unreasonable root (%f): " % (x0,) +
                           "Maybe the system was not at the equilibrium configuration")
    if debug:
        print "Root:", x0
    value = fit2(x0)
    return value / atoms.get_volume()

def minimize_unit_cell(atoms, symmetry):
    #from scipy.optimize import fmin_powell, fmin
    from asap3.Tools.ParameterOptimization.ScipyFmin import fmin
    if symmetry == 'cubic':
        var = [1.0]  # Scaling
    elif symmetry == 'hexagonal':
        var = [1.0, 1.0]
    else:
        raise ValueError('Symmetry "%s" not supported' % (symmetry,))
    u0 = atoms.get_cell()
    def energy(v, a=atoms, u=u0):
        newcell = v[0] * u
        if len(v) > 1:
            newcell[2] *= v[1]
        a.set_cell(newcell, scale_atoms=True)
        #if len(v) > 1:
        #    print "fmin: %.4f  %7.4f  %8.4f  %8.4f" % (v[0], v[1], a.get_positions()[1792,2], a.get_potential_energy())
        return a.get_potential_energy()
    xopt, fopt, iter, calls, flag = fmin(energy, var, delta=0.01, full_output=True, disp=False)
    #print "Minimize unit cell:", xopt, fopt, iter, calls

if __name__ == '__main__':
    from asap3 import EMT, units, EMThcpParameters
    from ase.lattice.cubic import FaceCenteredCubic
    from ase.lattice.hexagonal import HexagonalClosedPacked
    
    print "Calculating cubic constants for Cu"
    atoms = FaceCenteredCubic(size=(5,5,5), symbol='Cu', 
                              latticeconstant=3.59495722231)
    atoms.set_calculator(EMT())
    e = ElasticConstants(atoms, 'cubic')
    print np.array(e) / units.GPa
    
    print "Pretending that Cu is hexagonal"
    e = ElasticConstants(atoms, 'hexagonal', sanity=False)
    print np.array(e) / units.GPa
    
    print "Calculating elastic constants for Ru"
    atoms = HexagonalClosedPacked(size=(5,5,5), symbol='Ru',
                                  directions=[[1,-2,1,0], 
                                              [1,0,-1,0],
                                              [0,0,0,1]])
    atoms.set_calculator(EMT(EMThcpParameters()))
    e = ElasticConstants(atoms, 'hexagonal', minimize=True)
    print np.array(e) / units.GPa
    print "The EMT values are not even close to experimental values!"
    
    print "Pretending Ru is cubic"
    e = ElasticConstants(atoms, 'cubic')
    print np.array(e) / units.GPa
    
