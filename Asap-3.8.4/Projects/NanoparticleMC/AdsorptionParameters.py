import numpy as np
#Init Parameter dicts:
#Sgas is standard condition diatomic experimental entropy of gas in eV/K
#Eads is np.array of ads energies with [CN=0,CN=1,CN=....,CN=15]
#Viben is np.array of vibration energies(3 modes e.g. [10meV, 22.1meV,24meV]

#Au O
param_AuO = {'Sgas':0.002126,'Eads':[-0.35,-0.35,-0.35,-0.35,-0.35,-0.29,-0.45,-0.40,0.06,-0.08, 100,100,100,100,100,100] ,'viben':np.array([41.2E-3,43.3E-3,49.8E-3])}

#Cu H
param_CuH = {'Sgas':0.001354444194,'Eads':np.array([-0.25,-0.25,-0.25,-0.25,-0.25,-0.25,-0.2,-0.2,0.041,-0.107,100,100,100,100,100,100]),
'viben':np.array([89.1E-3,95.4E-3,120.7E-3]) }

#Au CO
param_AuCO = {'Sgas':0.00204860144,'Eads':np.array([-0.85,-0.85,-0.85,-0.85,-0.74,-0.66,-0.55,-0.30,-0.23,-0.12,100,100,100,100,100,100]),
'viben':np.array([13.5E-3,19.7E-3,21.9E-3,25.3E-3,43.6E-3,225.6E-3])}

#Now store them in dict:
params = {'AuO':param_AuO,'CuH':param_CuH,'AuCO':param_AuCO}


class adsorptionparameters:
    

    def __init__(self,sys_species="AuO"):
        self.spec = sys_species
        self.parameters = params[sys_species]
        #now this parameter holds the desired info.
        
    def get(self,param):
        return self.parameters[param]
        
        



