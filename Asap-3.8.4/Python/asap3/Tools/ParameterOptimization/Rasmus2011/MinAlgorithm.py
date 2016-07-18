# A Script using the downhill simplex algorithm to optimize the seven parameters in the EMT calculator

# Written by: Rasmus E. Christiansen

import Optimize as Optim
from ErrorFunction import ErrorFunction as EF
import numpy


# Previously fittet Values for parameters:
parameters = {
#          E0     s0    V0     eta2    kappa   lambda  n0      
#          eV     bohr  eV     bohr^-1 bohr^-1 bohr^-1 bohr^-3 
    'H':  (-2.21, 0.71, 2.132, 1.652,  2.790,  1.892,  0.00547),
    'Al': (-3.28, 3.00, 1.493, 1.240,  2.000,  1.169,  0.00700),
    'Cu': (-3.51, 2.67, 2.476, 1.652,  2.740,  1.906,  0.00910),
    'Ag': (-2.96, 3.01, 2.132, 1.652,  2.790,  1.892,  0.00547),
    'Au': (-3.80, 3.00, 2.321, 1.674,  2.873,  2.182,  0.00703),
    'Ni': (-4.44, 2.60, 3.673, 1.669,  2.757,  1.948,  0.01030),
    'Pd': (-3.90, 2.87, 2.773, 1.818,  3.107,  2.155,  0.00688),
    'Pt': (-5.85, 2.90, 4.067, 1.812,  3.145,  2.192,  0.00802),
    'C':  (-1.97, 1.18, 0.132, 3.652,  5.790,  2.892,  0.01322),
    'N':  (-4.97, 1.18, 0.132, 2.652,  3.790,  2.892,  0.01222),
    'O':  (-2.97, 1.25, 2.132, 3.652,  5.790,  4.892,  0.00850)}



# The parameters which are to be optimized are specified. This is done by initializing the dictionary 
# VariableParameters indexed by the elements whose parameters are being optimized including len()=8 arrays 
# of the form [0.,0.,1.,1.,0.,1.,0.], where 0 means the parameter is to be kept constant while 1 means it 
# is to be optimized. The Initial values of the parameters to be optimized are then set in the list ParamInit.
VariableParameters = {'Cu': numpy.array([1.,1.,1.,1.,1.,1.,1.]),
                      'Au': numpy.array([1.,1.,1.,1.,1.,1.,1.])}
Elements = VariableParameters.keys()
# The alloys desired are given in ElementsAllyos
#ElementsAlloys = ['AuCu3']
ElementsAlloys = []


# Experimental values picked from tabel
ParamInit = []
for i in Elements:
    for j in range(7):
        if VariableParameters[i][j] == 1:
            ParamInit.append(parameters[i][j])

# The ErrorFunction object is created using the chosen elements and variable parameters
ErrorFunc = EF(Elements,ElementsAlloys,VariableParameters)

# The minimization simplex algorithm is run using the ErrorFunction.
print Optim.fmin(ErrorFunc.ErFu,ParamInit,xtol=0.001,ftol=0.001)



