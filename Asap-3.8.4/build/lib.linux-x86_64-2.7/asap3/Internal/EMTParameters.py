"""Python definition of EMT parameters."""

from asap3.Internal.BuiltinPotentials import EMTParameters
import copy
import numpy as np
import ase

bohr = 0.5291772;  # Same number of digits as in EMTDefaultParameters
bohr3 = bohr * bohr * bohr

_std_emt_params = {
    13: {"name" :  "Al",
         "E0": -3.280,
         "S0": 3.000 * bohr,
         "V0": 1.493,
         "eta2": 1.240 / bohr,
         "kappa": 2.000 / bohr,
         "lambda": 1.169 / bohr,
         "n0": 0.007 / bohr3},
    29: {"name":  "Cu",
         "E0": -3.510,
         "S0": 2.67 * bohr,
         "V0": 2.476,
         "eta2": 1.652 / bohr,
         "kappa": 2.74 / bohr,
         "lambda": 1.906 / bohr,
         "n0": 0.0091 / bohr3},
    47: {"name": "Ag",
         "E0": -2.96,
         "S0": 3.01 * bohr,
         "V0": 2.132,
         "eta2": 1.652 / bohr,
         "kappa": 2.790 / bohr,
         "lambda": 1.892 / bohr,
         "n0": 0.00547 / bohr3},
    79: {"name":  "Au",
         "E0": -3.80,
         "S0": 3.00 * bohr,
         "V0": 2.321,
         "eta2": 1.674 / bohr,
         "kappa": 2.873 / bohr,
         "lambda": 2.182 / bohr,
         "n0": 0.00703 / bohr3},
    28: {"name":  "Ni",
         "E0": -4.44 ,
         "S0": 2.60  * bohr,
         "V0": 3.673,
         "eta2": 1.669 / bohr,
         "kappa": 2.757 / bohr,
         "lambda": 1.948 / bohr,
         "n0": 0.0103 / bohr3},
    46: {"name":  "Pd",
         "E0": -3.90 ,
         "S0": 2.87 * bohr,
         "V0": 2.773,
         "eta2": 1.818 / bohr,
         "kappa": 3.107 / bohr,
         "lambda": 2.155 / bohr,
         "n0": 0.00688 / bohr3},
    78: {"name":  "Pt",
         "E0": -5.85 ,
         "S0": 2.90 * bohr,
         "V0": 4.067,
         "eta2": 1.812 / bohr,
         "kappa": 3.145 / bohr,
         "lambda": 2.192 / bohr,
         "n0": 0.00802 / bohr3},
    }

class EMTStandardParameters(EMTParameters):
    """Standard EMT parameters.

    This is a Python implementation of EMTDefaultParameters,
    intended for use as a base class.
    """
    _defaultparameters = copy.deepcopy(_std_emt_params)
    beta = 1.809;  # Preserve same rounding as in ARTwork.
    
    def __init__(self):
        # Make a copy of parameters so we do not risk modifying the
        # orgininal ones in derived clases.
        EMTParameters.__init__(self)
        self.parameters = copy.deepcopy(self._defaultparameters)
        self.actual_params = []
        self.gammas_calculated = False
        self.shellpop = [12, 6, 24, 12, 24]
        self.shelldist = np.sqrt(np.arange(1,6))
        self.shell_cutoff = 3
        self.shell_gamma = 3
        self.listcutofffactor = 1.04500185048
        
    def get_maximal_cutoff(self):
        """The maximal cutoff value.  Defined before the elements is known."""
        maxS0 = 0.0
        for p in self.parameters.itervalues():
            S0 = p['S0']
            if S0 > maxS0:
                maxS0 = S0
        cutoff = 0.5 * maxS0 * self.beta * (np.sqrt(self.shell_cutoff) +
                                            np.sqrt(self.shell_cutoff + 1))
        return cutoff * self.listcutofffactor
    
    def get_parameters(self, element):
        assert self.gammas_calculated == False
        p = self.parameters[element]
        p["Z"] = element
        p["mass"] = ase.data.atomic_masses[element]
        p["name"] = ase.data.chemical_symbols[element]
        self.actual_params.append(p)
        return p

    def get_gammas_etc(self):
        assert self.gammas_calculated == False
        self.calc_cutoff()
        self.calc_gammas()
        self.calc_chi()
        self.gammas_calculated = True
        return (self.cutoff, self.gammas, self.chi)

    def calc_cutoff(self):
        maxS = 0.0
        for p in self.actual_params:
            if p["S0"] > maxS:
                maxS = p["S0"]
        cutoff = 0.5 * maxS * self.beta * (np.sqrt(self.shell_cutoff) +
                                           np.sqrt(self.shell_cutoff + 1))
        r = (cutoff * 2.0 * np.sqrt(self.shell_cutoff + 1)
             / (np.sqrt(self.shell_cutoff) + np.sqrt(self.shell_cutoff + 1)))
        cutslope = np.log(9999.0) / (r - cutoff);
        self.cutoff = (cutoff, cutslope, self.listcutofffactor)

    def calc_gammas(self):
        self.gammas = []
        for p in self.actual_params:
            gamma1 = gamma2 = 0.0
            for i in range(self.shell_gamma):
                d = self.shelldist[i] * self.beta * p["S0"]
                w = 1. / (1. + np.exp(self.cutoff[1] * (d - self.cutoff[0])))
                gamma1 += w * self.shellpop[i] * np.exp(-d * p["eta2"])
                gamma2 += w * self.shellpop[i] * np.exp(-d * p["kappa"]
                                                        / self.beta)
            gamma1 /= 12 * np.exp(-self.beta * p["S0"] *
                                                p["eta2"])
            gamma2 /= 12 * np.exp(-p["S0"] * p["kappa"])
            self.gammas.append((gamma1, gamma2))

    def calc_chi(self):
        n = len(self.actual_params)
        self.chi = np.zeros((n,n), float)
        for i in range(n):
            for j in range(n):
                self.chi[i,j] = (self.actual_params[j]["n0"]
                                 / self.actual_params[i]["n0"] )
        
class EMThcpParameters(EMTStandardParameters):
    """EMT parameters for elements in HCP structure.

    Currenlty, only Ruthenium is supported.
    """
    _hcp_parameters = {44: {"name": "Ru",
                            "E0": -6.350455,
                            "V0": 12.305218,
                            "S0": 1.498407,
                            "eta2": 3.945242,
                            "kappa": 6.561249,
                            "lambda": 6.701546,
                            "n0": 3.917210,
                            }
                       }
    
    def __init__(self):
        self.parameters = copy.deepcopy(self._hcp_parameters)
        self.actual_params = []
        self.gammas_calculated = False
        self.shellpop = (6, 6, 6, 2, 6, 12, 12, 6)
        c = 2*np.sqrt(2.0/3.0)
        self.shelldist = (1.0,
                          np.sqrt(c*c/4.0 + 1/3.0),
                          np.sqrt(c*c/4.0 + 4.0/3.0),
                          c,
                          np.sqrt(3.0),
                          np.sqrt(7.0/3.0 + c*c/4.0),
                          np.sqrt(1+c*c),
                          2.0)
        self.shell_gamma = 8
        self.shell_cutoff = 3  # Not really a shell number, but refers to FCC
        self.listcutofffactor = 1.04500185048
        
    def calc_cutoff(self):
        maxS = 0.0
        for p in self.actual_params:
            if p["S0"] > maxS:
                maxS = p["S0"]
        cutoff = maxS * self.beta * 1.824;
        r = (cutoff * 2.0 * np.sqrt(self.shell_cutoff + 1)
             / (np.sqrt(self.shell_cutoff) + np.sqrt(self.shell_cutoff + 1)))
        cutslope = np.log(9999.0) / (r - cutoff);
        listcutoffdistance = np.log(9999.0) / cutslope + cutoff;
        if listcutoffdistance > cutoff * self.listcutofffactor :
            listcutofffactor = listcutoffdistance / cutoff
        else:
            listcutofffactor = self.listcutofffactor
        self.cutoff = (cutoff, cutslope, listcutofffactor)

class EMTMetalGlassParameters(EMTStandardParameters):
    """EMT parameters for CuMg and CuZr metallic glasses.

    IMPORTANT: These copper parameters are different from the usual Cu
    parameters, and should ONLY be used to simulate metallic glasses.
    """
    _defaultparameters = {12: {"name": "Mg",
                               "E0":  -1.487000,
                               "eta2": 2.541137,
                               "kappa": 4.435425,
                               "lambda": 3.292725,
                               "n0": 0.035544,
                               "S0": 1.766399,
                               "V0": 2.229870},
                          29: {"name": "Cu",
                               "E0": -3.510000,
                               "eta2": 3.039871,
                               "kappa": 4.943848,
                               "lambda": 3.693666,
                               "n0": 0.063738,
                               "S0": 2.67 * bohr,  # Unchanged.
                               "V0": 1.993953},
                          40: {"name": "Zr",
                               "E0": -6.300000,
                               "eta2": 2.282059,
                               "kappa": 3.911040,
                               "lambda": 2.247348,
                               "n0": 0.031831,
                               "S0": 1.785939,
                               "V0": 2.320809}
                          }



del bohr, bohr3

                               
    
