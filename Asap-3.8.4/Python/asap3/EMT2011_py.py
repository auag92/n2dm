### Effective medium theory calculator. ###
# encoding: utf-8

# Written by: Rasmus E. Christiansen, s072162
# Version 1.02

# Imported Packages
import sys
import numpy
from math import sqrt,exp
from ase.data import atomic_numbers, chemical_symbols
from ase.units import Bohr
from asap3 import FullNeighborList

# Gobal parameters 
parameters = {
#          E0     s0    V0     eta2    kappa   lambda  n0      Lattice type
#          eV     bohr  eV     bohr^-1 bohr^-1 bohr^-1 bohr^-3 
    'H':  (-2.21, 0.71, 2.132, 1.652,  2.790,  1.892,  0.00547, 'dimer'),
    'Al': (-3.28, 3.00, 1.493, 1.240,  2.000,  1.169,  0.00700, 'fcc'),
    'Cu': (-3.51, 2.67, 2.476, 1.652,  2.740,  1.906,  0.00910, 'fcc'),
    'Ag': (-2.96, 3.01, 2.132, 1.652,  2.790,  1.892,  0.00547, 'fcc'),
    'Au': (-3.80, 3.00, 2.321, 1.674,  2.873,  2.182,  0.00703, 'fcc'),
    'Ni': (-4.44, 2.60, 3.673, 1.669,  2.757,  1.948,  0.01030, 'fcc'),
    'Pd': (-3.90, 2.87, 2.773, 1.818,  3.107,  2.155,  0.00688, 'fcc'),
    'Pt': (-5.85, 2.90, 4.067, 1.812,  3.145,  2.192,  0.00802, 'fcc'),
    'C':  (-1.97, 1.18, 0.132, 3.652,  5.790,  2.892,  0.01322, 'dimer'),
    'N':  (-4.97, 1.18, 0.132, 2.652,  3.790,  2.892,  0.01222, 'dimer'),
    'O':  (-2.97, 1.25, 2.132, 3.652,  5.790,  4.892,  0.00850, 'dimer')}

beta = 1.809    # Calculated from the following equation
                # beta = ((16 * Pi) / 3)^(1/3) / 2^(1/2)

# The "largest" element possibly supported by the calculator (this is determined by the list chemical_symbols).
NumElle = len(chemical_symbols)


class EMT:
    """ This class is an implementation of the revised edition of the Effective Medium Theory approach of calculating the energy of a given FCC crystal system. The functional form of the equations used can be found in ******* """
    def __init__(self,Params=None,ZeroPoint=1):
        """ Initializes the EMT object. The input Params is used to specify userdefined parameters and the input 
            Zeropoint is used to specify whether the potential energy measured realtive to E0 = Ecoh (ZeroPoint = 1) 
            or E0 = 0 (ZeroPoint = 0). """
        # Secures that the calculator is initialized correctly the first time it is used.
        self.energy = None
        self.ZeroPoint = ZeroPoint
        # If no parameters have been specified when creating the EMT object the default parameters are used.
        if Params is None:
            self.parameters = parameters
        else:
            self.parameters = Params

    def initialize(self, atoms):
        """ Method which initializes the EMT calculator by defining all needed values used in the
            calculations. """
        # A list, Z, of the element type for each atom in the system is defined:
        self.Z = atoms.get_atomic_numbers()
        # The number of atoms are calculated: 
        self.N = len(self.Z)        
        # Lists of the values for eta2, kappa, Seq, E0, V0, n0 and L (lambda) for the element types Z are 
        # defined:
        # The "largest" element number is corrected with regards to the method of counting in python.
        self.eta2 = numpy.zeros(NumElle)
        self.kappa = numpy.zeros(NumElle)
        self.Seq = numpy.zeros(NumElle)
        self.E0 = numpy.zeros(NumElle)
        self.V0 = numpy.zeros(NumElle)
        self.L = numpy.zeros(NumElle)
        self.n0 = numpy.zeros(NumElle)
        for i in range(NumElle):
            if chemical_symbols[i] in self.parameters:
                self.eta2[i] = self.parameters[chemical_symbols[i]][3] / Bohr
                self.kappa[i] = self.parameters[chemical_symbols[i]][4] / Bohr
                self.Seq[i] = self.parameters[chemical_symbols[i]][1] * Bohr
                self.E0[i] = self.parameters[chemical_symbols[i]][0]
                self.V0[i] = self.parameters[chemical_symbols[i]][2]
                self.L[i] = self.parameters[chemical_symbols[i]][5] / Bohr
                self.n0[i] = self.parameters[chemical_symbols[i]][6] / (Bohr**3)
        
        # Calculation of the X*X-arrays:    
        # r_cut; X*X-array of vaules for the cutoff length for the atompair of type (Z,Z') 
        # sigmaaRCUT; X*X-array of values for sigmaa evaluated in r_cut[Z,Z']
        # sigmabRCUT; X*X-array of values for sigmab evaluated in r_cut[Z,Z']
        # dsigmaadrRCUT; X*X-array of values for the deriviative of sigmaa evaluated in r_cut[Z,Z']
        # dsigmabdrRCUT; X*X-array of values for the deriviative of sigmab evaluated in r_cut[Z,Z']
        # chi; X*X-array of values for chi for the atompair of type (Z,Z')

        # The cutoff distances are calculated using that the lattice constant, a0 = sqrt(2)*beta*s0:
        self.r_cut = numpy.zeros([NumElle,NumElle])

        for i in range(NumElle):
            for j in range(NumElle):
                # Check to see if the two elements considered are defined or not. Only calculating an 
                # r_cut if both are. r_cut[Z,Z'] is calculated as the cutoff distance for the the one of 
                # elements Z,Z' which has the largest a0.
                if self.Seq[i] and self.Seq[j] != 0:
                    self.r_cut[i,j] = (1./2. * (sqrt(3. / 2.) + sqrt(4. / 2.)) * (sqrt(2) * beta) * 
                                      max(self.Seq[i],self.Seq[j]))

        
        ### Calculations for sigmaaRCUT, sigmabRCUT, d_sigmaadrRCUT and d_sigmabdrRCUT ###

        self.dsigmaadrRCUT = numpy.zeros([NumElle,NumElle])
        self.dsigmabdrRCUT = numpy.zeros([NumElle,NumElle])
        self.sigmaaRCUT = numpy.zeros([NumElle,NumElle])
        self.sigmabRCUT = numpy.zeros([NumElle,NumElle])

        for i in range(NumElle):
            for j in range(NumElle):
                # Check to see if r_cut[i,j] is defined for this pair of elements. r_cut[i,j] == 0 means 
                # that it is not defined.
                if self.r_cut[i,j] != 0:
                    self.sigmaaRCUT[i,j] = (numpy.exp(self.eta2[j] * 
                                           (-self.r_cut[i,j] + self.Seq[j] * beta)))
                    self.sigmabRCUT[i,j] = (numpy.exp(self.kappa[j] * 
                                           (-self.r_cut[i,j] / beta + self.Seq[j])))
                    self.dsigmaadrRCUT[i,j] = -self.eta2[j] * self.sigmaaRCUT[i,j]
                    self.dsigmabdrRCUT[i,j] = -self.kappa[j] / beta * self.sigmabRCUT[i,j]
        
        
        ### Calculations for chi[Z,Z'] ###

        self.chi = numpy.zeros([NumElle,NumElle])

        for i in range(NumElle):
            for j in range(NumElle):
                # Check to see if the elements i,j are defined.
                if self.n0[i] and self.n0[j] != 0:
                    self.chi[i,j] = self.n0[i] / self.n0[j]


        ### Calculations of gamma1 and gamma2 ###

	# Four (3 x NumElle)-arrays for lambda_1,2 (named L_1,2) and sigmaa_1,2 are calculated with the distance, 
        # r_ij = (beta * Seq, sqrt(2) * beta * Seq, sqrt(3) * beta * Seq) for all supportet elements. 
        # where Z = Z'.
        
        # The NumberNearestNeighbours variable is set to the number of nearest neighbors included in the model.
	NumberNearestNeighbours = 3
        # arrays for lambda and sigmaa are initialized 
	L_1_Z = numpy.zeros([NumberNearestNeighbours,NumElle])
	L_2_Z = numpy.zeros([NumberNearestNeighbours,NumElle])
	sigmaa_Z = numpy.zeros([NumberNearestNeighbours,NumElle])
	sigmab_Z = numpy.zeros([NumberNearestNeighbours,NumElle])
        # The values for each are calculated for each neighbour distance
        for i in range(NumberNearestNeighbours):
            L_1_Z[i] = (self.dsigmaadrRCUT[range(NumElle),range(NumElle)] *
                       ((sqrt(1 + i) * beta * self.Seq) - self.r_cut[range(NumElle),range(NumElle)]) +
                        self.sigmaaRCUT[range(NumElle),range(NumElle)])

            L_2_Z[i] = (self.dsigmabdrRCUT[range(NumElle),range(NumElle)] *
                       ((sqrt(1 + i) * beta * self.Seq) - self.r_cut[range(NumElle),range(NumElle)]) +
                        self.sigmabRCUT[range(NumElle),range(NumElle)])
              
            sigmaa_Z[i] = numpy.exp(self.eta2 * (-(sqrt(1 + i) * beta * self.Seq) + self.Seq * beta))
            sigmab_Z[i] = numpy.exp(self.kappa * (-(sqrt(1 + i) * self.Seq) + self.Seq))

        
        # The factor (self.Seq/self.Seq) is an array of zeros and ones and is only used to secure that only 
        # the elements which are actually defined in "parameters" gives a gamma_1,2 different from zero.
        self.gamma1 = ((self.Seq/self.Seq) * 
                      (12 * (sigmaa_Z[0] - L_1_Z[0]) + 
                       6 * (sigmaa_Z[1] - L_1_Z[1]) + 
                       24 * (sigmaa_Z[2] - L_1_Z[2]) ))
        self.gamma2 = ((self.Seq/self.Seq) * 
                      (12 * (sigmab_Z[0] - L_2_Z[0]) + 
                       6 * (sigmab_Z[1] - L_2_Z[1]) + 
                       24 * (sigmab_Z[2] - L_2_Z[2]) ))

        ### Construction of a Full Neighborlist for the system of atoms, 
        self.nbList = FullNeighborList(self.r_cut.max(),atoms)

        ### Initialization of the variables holding the forces on and energy of the atoms ###
        self.forces = None
        self.energy = None
        

    def NeighborList_rcutReduced(self,i):
        """ Method which makes sure that only the neighboratoms within the correct cutoff for the involved 
        element types Z,Z' are included in the calculations by modifying the output of the FullNeighborList
        function. """

        # Relavant data about the neighbor atom, j, for atom i which can possible give a contribution are 
        # selected
        (other_j,r_ij,rsq) = self.nbList.get_neighbors(i)

        # The neighbor atoms which will actually give a contribution to the energy, based on the individual 
        # cutoff distances between atom i of type Z and atom j of type Z', are selected.

        # The neighbor atoms which fullfill the condition are chosen
        keep = numpy.sqrt(rsq) <= self.r_cut[self.Z[i],self.Z[other_j]]
	
        # The lists of data about the neighbor atoms are updated
        if len(keep) != 0:
            return (other_j[keep],r_ij[keep],rsq[keep])
        else:
            # nbList returned empty lists, but we cannot index a shape (0,3) array (r_ij)
            # with an empty list (bug in numpy?).
            return (other_j,r_ij,rsq)

    
    def update(self, atoms):
        """ This method is called by the atoms object to which the calculator is attached, it secures that the 
            energy (and/or force) of the system is recalculated if this is required. """
        need_calc = False
        if (self.energy is None or len(self.Z) != len(atoms) or (self.Z != atoms.get_atomic_numbers()).any()):
            # The calculator is initialized with regards to the atoms object.
            self.initialize(atoms)
            need_calc = True
        elif (self.positions != atoms.get_positions()).any():
            # The atoms object has not changed enough for the calculator to need a reinitialization but a 
            # new calculation of the value for the energies are still needed.
            need_calc = True
        if need_calc:
            self.positions = atoms.get_positions()
            self.nbList.check_and_update(atoms)
            self.energy = self.calculate_Energy()
            self.forces = self.calculate_Force()


    # Returns the energy of the atoms (the method calculates the energy first if needed be)
    def get_potential_energy(self, atoms):
         self.update(atoms)
         return self.energy

    # Returns the forces on the atoms (the method calculates the forces first if needed be)
    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces.copy()

    def get_stress(self, atoms):
        raise NotImplementedError('No stresses implemented')


    ########## ENERGY Calculations ##########

    ### sigma_1,2 ###

    
    def calculate_sigma12(self):
        """ Calculates and returns sigma_1 and sigma_2. """
        # The N-arrays for sigma_1,2 are initialized
        sigma_1 = numpy.zeros(self.N)
        sigma_2 = numpy.zeros(self.N)

       
        for i in range(self.N):
            # The numbers of the neighbor atoms, the relative position vectors and length of the position 
            # vectors squared between atom i and the neighbor atoms j are defined in three arrays. 
            (other_j,r_ij,rsq) = self.NeighborList_rcutReduced(i)
            
            
            # The values for the linear subtracktion functions evaluated at norm(r_ij,2) for all the atom 
            # pairs, [i,other_j], are calculated.
            L_1_i = (self.dsigmaadrRCUT[self.Z[i],self.Z[other_j]] *
                    (numpy.sqrt(rsq) - self.r_cut[self.Z[i],self.Z[other_j]]) +
                     self.sigmaaRCUT[self.Z[i],self.Z[other_j]])
            L_2_i = (self.dsigmabdrRCUT[self.Z[i],self.Z[other_j]] *
                    (numpy.sqrt(rsq) - self.r_cut[self.Z[i],self.Z[other_j]]) +
                     self.sigmabRCUT[self.Z[i],self.Z[other_j]])
                       
            # sigmaa_i and sigmab_i are evaluated at norm(r_ij,2) for all the atom pairs, [i,other_j]. 
            sigmaa_i = (numpy.exp(self.eta2[self.Z[other_j]] * 
                       (-numpy.sqrt(rsq) + self.Seq[self.Z[other_j]] * beta)))
            sigmab_i = (numpy.exp(self.kappa[self.Z[other_j]] * 
                       (-numpy.sqrt(rsq) / beta + self.Seq[self.Z[other_j]])))

            #if (i == 1):
            #    print sigmaa_i

            # The values of sigma_1_i and sigma_2_i are calculated for
            # the atom i. Where max(a,b) is introduced in order to
            # secure a none zero minimumvalue for sigma_1 so the
            # following calculations wont result in an error.
            sigma_1[i] = max( pow(10,-9) , (self.chi[self.Z[i],self.Z[other_j]] * (sigmaa_i - L_1_i)).sum() )
            sigma_2[i] = (self.chi[self.Z[i],self.Z[other_j]] * (sigmab_i - L_2_i)).sum()

        #TESTER#
        #print sigma_1[:10]
        
        
        return (sigma_1,sigma_2)


    ### s_i ###

    def calculate_s(self,sigma_1):
        """ Calculates and returns an N-array containing the neutrality sphere radii, s, for the atoms of the 
        system is calculated."""    
        return self.Seq[self.Z] - numpy.log(sigma_1 / self.gamma1[self.Z]) / (beta * self.eta2[self.Z])
 

    ### E_tot ###

    
    def calculate_Energy_function(self,s,sigma_2):
        """ Calculates and returns the total energy of the system using s and sigma_2. """
        # Calculation of the N-array containing the cohesive energy for each of the N atoms.
        E_c = ( self.E0[self.Z] * (self.L[self.Z] * (s - self.Seq[self.Z]) + 1) * 
                numpy.exp(-self.L[self.Z] * (s - self.Seq[self.Z])) ) - self.E0[self.Z] * self.ZeroPoint

        # Calculation of the N-array containing the atomic sphere correction energy for each of the N atoms.
        E_as = (6 * ( self.V0[self.Z] * numpy.exp(-self.kappa[self.Z] * (s - self.Seq[self.Z])) 
                -self.V0[self.Z] * sigma_2 / self.gamma2[self.Z] ) )

        # Calculation of the total energy
        return (E_c + E_as).sum()


    ### Final Energy Calculator ###

    def calculate_Energy(self):
        """ Calculates and returnes the energy of the atoms in the atom object to which the 
            EMT calculator is attached. The calculations are done using the following methods, 
            also defined in EMT.py: calculate_sigma12(self), calculate_s(self,sigma_1), 
            calculate_Energy_function(self,s,sigma_2). """
        (sigma_1,sigma_2) = self.calculate_sigma12()
        s = self.calculate_s(sigma_1)
    
        # The total energy is calculated and returned
        return self.calculate_Energy_function(s,sigma_2)




    ########## FORCE Calculations ##########


    ### dsdsigma_1 ###

    
    def calculate_dsdsigma_1(self,sigma_1):
        """ Calculates and returns dsdsigma_1 using sigma_1. """
        # An N-array containing the the deriviative of neutrality sphere radii, s, with regards to sigma_1 for 
        # the atoms of the system is calculated.    
        dsdsigma_1 = -1 / (beta * self.eta2[self.Z] * sigma_1)

        return dsdsigma_1


    ### dE_cids, dE_asds, dE_asdsigma_2 ###

    def calculate_Deriviative_of_Energy(self,s):
        """ Calculates and returns the deriviatives of E_cs and E_as with regards to s and sigma_2. """
        # Calculation of the N-array containing the deriviative of the cohesive energy with regards to s for 
        # each of the N atoms.
        dE_cds = -( self.E0[self.Z] * self.L[self.Z] * self.L[self.Z] * 
                    numpy.exp(-self.L[self.Z] * (s - self.Seq[self.Z])) * (s - self.Seq[self.Z]) )

        # Calculation of the N-array containing the deriviative of the atomic sphere correction energy with 
        # regards to s for each of the N atoms.
        dE_asds = -6 * self.kappa[self.Z] * self.V0[self.Z] * numpy.exp(-self.kappa[self.Z] * (s - self.Seq[self.Z]))

        # Calculation of the N-array containing the deriviative of the atomic sphere correction energy with 
        # regards to sigma_2 for each of the N atoms.
        dE_asdsigma_2 = -6 * self.V0[self.Z] / (self.gamma2[self.Z])

        return (dE_cds,dE_asds,dE_asdsigma_2)


    ### F_kalpha ###
    def calculate_Force_function(self,dE_cds,dE_asds,dE_asdsigma_2,dsdsigma_1):
        """ Calculates the force on all k atoms in the three directions {x,y,z} representet by alpha. """ 
            
        # An array for the force is initialized
        F = numpy.zeros([self.N,3])

        for k in range(self.N):
            # The atoms interacting with atom k are selected.
            (other_i,r_ki,rsq) = self.NeighborList_rcutReduced(k)
            #print other_i
            #print r_ki
            #print numpy.sqrt(rsq)
            
            # The values for dr_ijdr_kalpha are calculated for the relevant atoms, k and other_i.
            dr_kidr_kalpha = r_ki / numpy.sqrt(rsq)[:,numpy.newaxis]
            
            ## The force on the k'th atom caused by the atoms, other_i's, interactions with k are calculated ##

            # The values for dsigmaa_idr_ij and dsigmab_idr_ij are calculated with regards to the k'th atom
            sigmaa_k = numpy.exp(self.eta2[self.Z[other_i]] * 
                       (-numpy.sqrt(rsq) + self.Seq[self.Z[other_i]] * beta))
            sigmab_k = numpy.exp(self.kappa[self.Z[other_i]] * 
                       (-numpy.sqrt(rsq) / beta + self.Seq[self.Z[other_i]]))
            
            dsigmaa_kdr_ki = (-self.eta2[self.Z[other_i]] * sigmaa_k)[:,numpy.newaxis]
            dsigmab_kdr_ki = (-self.kappa[self.Z[other_i]] / beta * sigmab_k)[:,numpy.newaxis]
            
            # Values for dL_1idr_ij and dL_2idr_ij are calculated with regards to the k'th atom
            dL_1kdr_ki = self.dsigmaadrRCUT[self.Z[k],self.Z[other_i]][:,numpy.newaxis]
            dL_2kdr_ki = self.dsigmabdrRCUT[self.Z[k],self.Z[other_i]][:,numpy.newaxis]
            
            # First the value of dsigma1_idr_kaplha and dsigma2_idr_kaplha are calculated for the k'th atom
            dsigma1_kdr_kalpha = ( self.chi[self.Z[k],self.Z[other_i]][:,numpy.newaxis] * 
                                   (dsigmaa_kdr_ki - dL_1kdr_ki) * dr_kidr_kalpha ).sum(axis=0)
            dsigma2_kdr_kalpha = ( self.chi[self.Z[k],self.Z[other_i]][:,numpy.newaxis] * 
                                   (dsigmab_kdr_ki - dL_2kdr_ki) * dr_kidr_kalpha ).sum(axis=0)
            
            """ TJEK DER SKAL FJERNES SENERE """
            assert len(dsigma1_kdr_kalpha) == 3
            """ TJEK DER SKAL FJERNES SENERE """

            # The contribution to the force on atom k from the k'th atoms interaction with the other_i atoms is 
            # calculated
            F[k] = (dE_cds[k] * dsdsigma_1[k] * dsigma1_kdr_kalpha +
                    dE_asds[k] * dsdsigma_1[k] * dsigma1_kdr_kalpha + 
                    dE_asdsigma_2[k] * dsigma2_kdr_kalpha)
            

            # The values for dsigmaa_idr_ij and dsigmab_idr_ij are calculated with regards to the atoms other_i 
            # where j = k for all other_i (thus we only need one value of dsigmaa_idr_ij).
            sigmaa_i = numpy.exp(self.eta2[self.Z[k]] * (-numpy.sqrt(rsq) + self.Seq[self.Z[k]] * beta))
            sigmab_i = numpy.exp(self.kappa[self.Z[k]] * (-numpy.sqrt(rsq) / beta + self.Seq[self.Z[k]]))
            
            dsigmaa_idr_ik = (-self.eta2[self.Z[k]] * sigmaa_i)[:,numpy.newaxis]
            dsigmab_idr_ik = (-self.kappa[self.Z[k]] / beta * sigmab_i)[:,numpy.newaxis]
            
            # Values for dL_1idr_ij and dL_2idr_ij are calculated with regards to the atoms other_i 
            # where j = k for all other_i.
            dL_1idr_ik = self.dsigmaadrRCUT[self.Z[other_i],self.Z[k]][:,numpy.newaxis]
            dL_2idr_ik = self.dsigmabdrRCUT[self.Z[other_i],self.Z[k]][:,numpy.newaxis]
            
            # First the value of dsigma1_idr_kaplha and dsigma2_idr_kaplha are calculated with regards to the atoms 
            # other_i where j are only the atom k for all other_i. (thus the sum only has one element for all other_i.
            # which results in the calculations leading to an [other_i,3]-array. 
            dsigma1_idr_kalpha = (self.chi[self.Z[other_i],self.Z[k]][:,numpy.newaxis] * 
                                  (dsigmaa_idr_ik - dL_1idr_ik) * (dr_kidr_kalpha) )
            dsigma2_idr_kalpha = (self.chi[self.Z[other_i],self.Z[k]][:,numpy.newaxis] * 
                                  (dsigmab_idr_ik - dL_2idr_ik) * (dr_kidr_kalpha) )
            
            # The contribution to the force on atom k from the other_i atoms interaction with the k'th atom is now 
            # calculated
            F[k] += (dE_cds[other_i][:,numpy.newaxis] * dsdsigma_1[other_i][:,numpy.newaxis] * dsigma1_idr_kalpha +
                     dE_asds[other_i][:,numpy.newaxis] * dsdsigma_1[other_i][:,numpy.newaxis] * dsigma1_idr_kalpha + 
                     dE_asdsigma_2[other_i][:,numpy.newaxis] * dsigma2_idr_kalpha).sum(axis=0)

            """ TJEK DER SKAL FJERNES SENERE """
            assert len(F[k]) == 3
            """ TJEK DER SKAL FJERNES SENERE """

        return F


    ### Final Force Calculator ###

    def calculate_Force(self):
        """ Calculates and returnes the force acting on each of the atoms in the atoms object to which the 
            EMT calculator is attached. These calculations are done using the following methods, 
            also defined in EMT.py: calculate_sigma12(self), calculate_s(self,sigma_1), calculate_dsdsigma_1
            (self,sigma_1), calculate_Deriviative_of_Energy(self,s) and calculate_Force_function
            (self,dE_cds,dE_asds,dE_asdsigma_2,dsdsigma_1) """
        (sigma_1,sigma_2) = self.calculate_sigma12()
        s = self.calculate_s(sigma_1)
        dsdsigma_1 = self.calculate_dsdsigma_1(sigma_1)
        (dE_cds,dE_asds,dE_asdsigma_2) = self.calculate_Deriviative_of_Energy(s)
    
        # The force is calculated and returned
        return self.calculate_Force_function(dE_cds,dE_asds,dE_asdsigma_2,dsdsigma_1)

    
 
