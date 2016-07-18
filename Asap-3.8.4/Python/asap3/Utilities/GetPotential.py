import sys
import os
import string

from Asap.Internal.BuiltinPotentials import *


def GetParameterProvider(fileName):
    parameterFile = open(fileName,"r")
    paramList = []
    elementList = []
    for line in parameterFile.readlines():
        paramName, paramValueStr = string.split(line)
        paramParts = string.split(paramName,'_')
        elem = string.atoi(paramParts[1])
        paramValue = string.atof(paramValueStr)
        set_methodName = "Set_" + paramParts[0]
        paramList.append([set_methodName, elem, paramValue])

            
        if elem not in elementList:
            elementList.append(elem)

    parameterProvider = EMTVariableParameterProvider(elementList)
    
    for item in paramList:
        set_methodName, elem, paramValue = item
        set_method = getattr(parameterProvider, set_methodName)
        set_method(elem, paramValue)
    return parameterProvider


def GetPotential(fileName = None, subtractE0=1):
    """This function returns an EMT potential object given a filename and
    optionally an argument related to the definition of the zero of potential
    energy. The filename refers to a file that contains extra parameters
    for EMT. These can either define an element not already in EMT, or can
    override some or all of the default parameters for an existing element.
    Parameters for more than one element can be included. Lines in the file
    must have three entries, an example being

    S0_12	  1.766399

    which defines the S0 parameter of Mg (Z=12) to be 1.766399. The parameter
    names are S0, E0, kappa, V0, eta2, n0, lambda    
    """
    if fileName is None and len(sys.argv) == 2:
        fileName = sys.argv[1]
        
    if fileName is not None:
        parameterProvider = GetParameterProvider(fileName)
        potential = EMT(parameterProvider)
    else:
        print "WARNING: Using default EMT potential"
        potential = EMT()
    potential.SetSubtractE0(subtractE0)
    return potential
