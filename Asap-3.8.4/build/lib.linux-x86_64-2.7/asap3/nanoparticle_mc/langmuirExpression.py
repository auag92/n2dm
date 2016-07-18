import numpy as np
from ase import units
import sys
sys.path.append('.')
from AdsorptionParameters import *
#Changes: Added dictionary to include multiple objects.
    
def getEnergies(species):
    
    p = adsorptionparameters(species)
    
    return p.get('Eads') 
    
    #return np.array([-0.25,-0.25,-0.25,-0.25,-0.25,-0.25,-0.2,-0.2,0.041,-0.107
    #,100,100,100,100,100,100])


def getCoverages(T=222.15,P=1E2,species="AuO"):
    #This script calculates the expression for the Langmuir isotherm(a function of Temperature and pressure):
    #It returns coverages ordered by coordinations such that CN=0 ~ index 0
    
    p = adsorptionparameters(species)
    
    #Constants:
    Pstd=1E5 #1E5 Pa is standard pressure
    Tstd = 298.15 #25 cecluis degrees
    
    Viben = np.array(p.get('viben')) #Vibration Energies
    Sgas = p.get('Sgas') #Entropy of gas (standard conditions)
    
    cp = 7./2*units.kB #Per atom diving by 2.heat capacity of diatomic ideal gas for constant pressure.
   
    #Adsorption Energies:
    dE = p.get('Eads')

    #Calculate entropy of vibrations:
    Svib = units.kB*((Viben/units.kB/T/(np.exp(Viben/units.kB/T)-1)-np.log(1-np.exp(-Viben/units.kB/T))).sum())
    
    #Rate Constants:
    dF = dE-T*(Svib-(Sgas+cp*np.log(T/Tstd)-units.kB*np.log(P/Pstd))/2.0)
    K = np.exp(-dF/units.kB/T)
    #Coverage Langmuir dissociative adsorption:
    theta = np.sqrt(K)/(1+np.sqrt(K))
    
    return theta


