from Numeric import array

__docformat__ = "restructuredtext en"


def MakeStrainVector(strainMatrix):
    """ Gives the strain expressed as a vector, the so-called engineering
    strain, where off diagonal components are multiplied by a factor of two.
    This ensures that when the dot-product is taken with the stress vector,
    the double counting is done correctly.    
    """
    return array([strainMatrix[0,0], strainMatrix[1,1],
                  strainMatrix[2,2], 2*strainMatrix[1,2],
                  2*strainMatrix[0,2], 2*strainMatrix[0,1]])

def MakeStrainMatrix(strainVector):
    """ Given the strain expressed as a vector, the so-called engineering
    strain, returns the matrix representation of the strain.
    """
    return array([[strainVector[0], 0.5*strainVector[5],0.5*strainVector[4]],
                  [0.5*strainVector[5],strainVector[1],0.5*strainVector[3]],
                  [0.5*strainVector[4],0.5*strainVector[3],strainVector[2]]])
