# encoding: utf-8

# Written by: Rasmus E. Christiansen

# This Class implements an "error" function which gives a measure for how well the seven parameters used in the EMT() calculator are fittet for each element. This is done by comparing a set of quantities to their experimental values and weighting them with regards to the importance of their precision represented by some standard diviation.

# Imported packages
from QuantityCalculatorEMT import Calculator
from asap3.EMT2013 import EMT
import copy
import numpy

# Gobal parameters
parameters = {
#                      E0     s0    V0     eta2    kappa   lambda  n0     
#                      eV     bohr  eV     bohr^-1 bohr^-1 bohr^-1 bohr^-3 
    'H':  numpy.array([-2.21, 0.71, 2.132, 1.652,  2.790,  1.892,  0.00547]),
    'Al': numpy.array([-3.28, 3.00, 1.493, 1.240,  2.000,  1.169,  0.00700]),
    'Cu': numpy.array([-3.51, 2.67, 2.476, 1.652,  2.740,  1.906,  0.00910]),
    'Ag': numpy.array([-2.96, 3.01, 2.132, 1.652,  2.790,  1.892,  0.00547]),
    'Au': numpy.array([-3.80, 3.00, 2.321, 1.674,  2.873,  2.182,  0.00703]),
    'Ni': numpy.array([-4.44, 2.60, 3.673, 1.669,  2.757,  1.948,  0.01030]),
    'Pd': numpy.array([-3.90, 2.87, 2.773, 1.818,  3.107,  2.155,  0.00688]),
    'Pt': numpy.array([-5.85, 2.90, 4.067, 1.812,  3.145,  2.192,  0.00802]),
    'C':  numpy.array([-1.97, 1.18, 0.132, 3.652,  5.790,  2.892,  0.01322]),
    'N':  numpy.array([-4.97, 1.18, 0.132, 2.652,  3.790,  2.892,  0.01222]),
    'O':  numpy.array([-2.97, 1.25, 2.132, 3.652,  5.790,  4.892,  0.00850])}

beta = 1.809        # Calculated from the following equation
                    # beta = ((16 * Pi) / 3)^(1/3) / 2^(1/2)

# The quantities for pure systems that are being used for optimization of the seven parameters are given in the 
# dictionary: quantities.

# The structure of the quantities dictionary is given as follows 
# quantities = {Element: [[BM_exp,BM_sigma,Calc.BulkModulus],
#                    [a0_exp,a0_sigma,Calc.LatticeConstant],
#                    [Evoc_exp,Evoc_sigma,Calc.EvocCalculator],
#                    [EfccEhcp_exp,EfccEhcp_sigma,Calc.EfccEhcpCalculator],
#                    [EfccEbcc_exp,EfccEbcc_sigma,Calc.EfccEbccCalculator],
#                    [C44_exp,C44_sigma,Calc.C44Calculator],
#                    [C11_exp,C11_sigma,Calc.C11_Calculator],
#                    [E_coh_exp,E_coh_sigma,Calc.E_cohCalculator]
#                    [E_surf001_exp,E_surf001_sigma,Calc.Surface001]
#                    [E_surf111_exp,E_surf111_sigma,Calc.Surface001]]
#               ...                    

# The units for the experimental values are given as follows
# [BM_exp]: eV/Å^3 , [a0_exp]: Å , [C44_exp]: eV/Å^3 , [C11_exp]: eV/Å^3 ,
# [Evoc_exp]: eV ,[Ecoh_exp]: eV , [EfccEhcp_exp]: eV , [EfccEbcc_exp]: eV
# [E_surf001_exp]: eV / Å^2, [E_surf111_exp]: eV / Å^2.

# Here XX_exp is the experimental value of the quantity XX.
# XX_sigma is the accepted standard diviation of the quantity XX. 
# BM is the Bulk modulus, a0 the lattice constant, Evoc the vacation energy, EfccEhcp and EfccEbcc the energy 
# difference between a system of atoms of a given element in the FCC and HCP/BCC structures, C44 the Shear modulus,
# C11 another elastic constant, and lastly E_coh the Coherence energy of the system.

# The first column of the dictionary is filled (The experimental values for the parameters)
" FOR NOW THE Evoc_exp or surface energies of Al is not known! This should not be used!! they are all set to Val = 1 in order for the script to run "
# The experimental values are added to the dictionary first.
quantities = {'Al': [[0.474],
                     [4.050],
                     [1], # DONT USE THIS VALUE ITS WRONG!!!!
                     [0.0371],
                     [0.0925],
                     [0.177],
                     [0.666],
                     [3.39],
                     [1], # DONT USE THIS VALUE ITS WRONG!!!!
                     [1]], # DONT USE THIS VALUE ITS WRONG!!!!
              'Cu': [[0.874],
                     [3.615],
                     [1.3],
                     [0.00670],
                     [0.042],
                     [0.472],
                     [1.050],
                     [3.49],
                     [0.1352],
                     [0.1218]],
              'Ag': [[0.624],
                     [4.085],
                     [1.1],
                     [0.00567],
                     [0.0336],
                     [0.288],
                     [0.774],
                     [2.95],
                     [0.07489],
                     [0.07315]],
              'Au': [[1.373],
                     [4.078],
                     [0.9],
                     [0.02], # Uses estimate of EfccEhcp as DFT gives: -0.0201.
                     [0.0174],
                     [0.262],
                     [1.201],
                     [3.81],
                     [0.1015],
                     [0.08008]],
              'Ni': [[1.124],
                     [3.524],
                     [1.6],
                     [0.0325],
                     [0.0678],
                     [0.775],
                     [1.549],
                     [4.44],
                     [0.1514],
                     [0.1255]],
              'Pd': [[1.123],
                     [3.891],
                     [1.4],
                     [0.0363],
                     [0.0410],
                     [0.448],
                     [1.417],
                     [3.89],
                     [0.1452],
                     [0.1198]],
              'Pt': [[1.436],
                     [3.924],
                     [1.5],
                     [0.0825],
                     [0.0951],
                     [0.477],
                     [2.163],
                     [5.84],
                     [0.1706],
                     [0.1435]],
              }

# The procentual values for the standard diviation are converted to be procental values of the experimental values
for Elements,data in quantities.items():
    data[0].append(data[0][0] * 0.035)   # Bulk modulus
    data[1].append(data[1][0] * 0.005)  # Lattice constant
    data[2].append(data[2][0] * 0.05)   # Evoc
    if Elements == 'Au':
        data[3].append(data[3][0]) # STD for Au specifically 
    else:
        data[3].append(data[3][0] * 0.05)   # Hcp-Fcc
    data[4].append(data[4][0] * 0.05)   # Bcc-Fcc
    data[5].append(data[5][0] * 0.035)   # C44
    data[6].append(data[6][0] * 0.035)   # C11
    data[7].append(data[7][0] * 0.05)   # Ecoh
    data[8].append(data[8][0] * 0.1)   # E_surf001
    data[9].append(data[9][0] * 0.1)   # E_surf111

# The calculator objects from the QuantityCalculator package are created and appended with the relavant method to each 
for Elements,data in quantities.items():
    Calc = Calculator(Elements)
    data[0].append(Calc.BulkModulus)         # Bulk modulus
    data[1].append(Calc.LatticeConstant)     # Lattice constant
    data[2].append(Calc.EvocCalculator)      # Evoc
    data[3].append(Calc.EfccEhcpCalculator)  # Hcp-Fcc
    data[4].append(Calc.EfccEbccCalculator)  # Bcc-Fcc
    data[5].append(Calc.C44_Calculator)      # C44
    data[6].append(Calc.C11_Calculator)      # C11
    data[7].append(Calc.E_cohCalculator)     # Ecoh
    data[8].append(Calc.Surface001)          # E_surf001
    data[9].append(Calc.Surface111)          # E_surf111


# # The quantities for alloys used to fit the parameters are presented in quantitiesAlloys 
# quantitiesAlloys = {'AuCu3': [[0.878],       # Bulk modulus
#                               [3.790],       # Lattice constant
#                               [-0.0514]]}    # E_hof

# # The procentual values for the standard diviation are converted to be procental values of the experimental values
# for Elements,data in quantitiesAlloys.items():
#     data[0].append(data[0][0] * 0.035)  # L12_Bulk modulus
#     data[1].append(data[1][0] * 0.010)  # L12_Lattice constant
#     data[2].append(data[2][0] * 0.05)   # L12_HOFE


# # The calculator objects from the QuantityCalculator package are created and appended with the relavant method to each 
# for Elements,data in quantitiesAlloys.items():
#     Calc = Calculator([Elements[0:2],Elements[2:4]])
#     data[0].append(Calc.L12_BulkModulus)     # L12_BulkModulus
#     data[1].append(Calc.L12_LatticeConstant) # L12_LatticeConstant
#     data[2].append(Calc.L12_HOFE)            # L12_HOFE



class ErrorFunction:
    """ This class creates an ErFu object which is used by the MinAlgorithm script in order to optimize a number of  parameters used by the EMT() potential script. """

    def __init__(self,Elements,ElementsAlloys,VariableParameters):

        # The quantities for the elements used in the optimization are pulled from the quantities dictionary.
        self.q = {}
        self.qAlloys = {}
        for i in Elements:
            self.q[i] = quantities[i] 
        for i in ElementsAlloys:
            self.qAlloys[i] = quantitiesAlloys[i]
        # The parameters which are going to be allowed to vary during the optimization process are set.
        self.VP = VariableParameters
        # The initial parametervalues are set and the elements going through optimization are identified
        self.params = copy.deepcopy(parameters)
        self.Elements = Elements
        
    def ErFu(self,ParametersChanged):
        """ This method calculates and returns the value of an "Error" Function based on the current value of the 
            parameters in the input dictionary: ParametersChanged. """
        # The current value of the parameters undergoing optimization in this evaluation 
        # of the Error function are transfered to the used arrays.
        NewParamValues = copy.deepcopy(self.VP)
        ParamNumber = 0
        for i in self.Elements:
            for j in range(7):
                if self.VP[i][j] == 1:
                    NewParamValues[i][j] = ParametersChanged[ParamNumber]
                    ParamNumber += 1

        for i in self.Elements:
            self.params[i] = (self.params[i] - self.params[i] * self.VP[i]) + NewParamValues[i]

        # The calculator is initialized with the current values of the parameters.
        self.EMT = EMT(self.params)
        
        # The ErrorFunction is calculated using the objects given in self.q and self.qAlloys and it value is returned. 
        ErrFunc = 0
        for Elle,List in self.qAlloys.items():
            for i in range(len(List)):
                ErrFunc += (List[i][2](self.EMT,self.params) - List[i][0])**2 / List[i][1]**2
        for Elle,List in self.q.items():
            for i in range(len(List)):
                ErrFunc += (List[i][2](self.EMT,self.params) - List[i][0])**2 / List[i][1]**2
                print "Done something:", ErrFunc


        # The current parametervalues and Error function value are printed to the user
        for j in self.Elements:
            print self.params[j]
        print ErrFunc

        return ErrFunc










