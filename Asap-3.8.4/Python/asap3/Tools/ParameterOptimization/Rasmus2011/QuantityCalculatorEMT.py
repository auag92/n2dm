# encoding: utf-8
# Class which creates "Calculator" objects that uses the EMT() model to calculate the quantities used in the ErrorFunction script in order to optimize the value of the seven parameters used in the EMT() model. 

# Written by: Rasmus E. Christiansen

# Imported Packages
import numpy
from asap3.EMT2013 import EMT
from ase.lattice.cubic import FaceCenteredCubic, BodyCenteredCubic
from ase.lattice.hexagonal import HexagonalClosedPacked
from ase.lattice.compounds import L1_0,L1_2
from ase.data import atomic_numbers, chemical_symbols 
from ase.units import Bohr
import numpy.lib.polynomial as numpyfit
from ase.utils.eos import EquationOfState
from ase.all import view
from ase.optimize import BFGS
from ase.io import read

# Global parameters
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

class Calculator:
    """ The class creates calculator objects which can be used to calculate the following 16 different quantities for 
        a crystal in the FCC structure using the new version of the EMT calculator. For pure crystals: Bulk modulus, 
        Lattice Constant, Vacancy formation energy, Energy difference between the element in question in an FCC 
        configuration and an HCP or BCC configuration, The elastic constants C44 and C11 along with the cohesive energy 
        and the surface energy for a 001 and 111 surface. 
        For the L10 and L12 alloys of two elements: The heat of formation energy, the bulk modulus and the lattice 
        constant for the L10 and L11 alloys.  """
    def __init__(self,Element,Size=3,Tol=0.002):
        # The element chosen for the calculator is set
        self.Element = Element
        self.Size = Size
        self.Tol = Tol

    def BM_LC_Calculator(self,EMT,PARAMETERS):
        # Return 0 is used when the method is not desired to be used
        # return 0
        """ The method BM_LC_Calculator uses the EMT calculator to find the lattice constant which gives 
            the lowest possible energy for a system of atoms and from a polynomial fit of the energy as a function of 
            the lattice constant it calculates and returns the Bulk modulus, the volume of the system at the lowest 
            possible energy and this energy """
        # Three values for the lattice constant, a0, are chosen. The highest and lowest value is chosen so that it     
        # is certain that these are to high/low respectively compared to the "correct" value. This is done by using the 
        # experimental value of the lattice constant, a0_exp, as the middle value and to high/low values are then: 
        # a0_exp +- a0_mod, a0_mod = a0_exp/10
        a0_exp = numpy.sqrt(2) * beta * PARAMETERS[self.Element][1] * Bohr
        a0_mod = a0_exp * 0.10 
        a0_guesses = numpy.array([a0_exp - a0_mod,a0_exp,a0_exp + a0_mod])
    
        # An atoms object consisting of atoms of the chosen element is initialized 
        atoms = FaceCenteredCubic(size=(self.Size,self.Size,self.Size), symbol=self.Element)
        atoms.set_calculator(EMT)

        # An identity matrix for the system is saved to a variable
        IdentityMatrix = numpy.array([[1,0,0],[0,1,0],[0,0,1]])

        # An array for the energy for the chosen guesses for a0 is initialized:
        E_guesses = numpy.zeros(5)

        # The energies are calculated
        for i in range(3):
        # Changes the lattice constant for the atoms object
            atoms.set_cell(a0_guesses[i] * self.Size * IdentityMatrix,scale_atoms=True)
            # Calculates the energy of the system for the new lattice constant
            E_guesses[i] = atoms.get_potential_energy()
           
        # Bisection is used in order to find a small interval of the lattice constant for the minimum of the energy (an 
        # abitrary interval length is chosen), This is possible because we are certian from theory that there is only 
        # one minimum of the energy function in the interval of interest and thus we wont "fall" into a local minimum 
        # and stay there by accident.
        while (a0_guesses[2]-a0_guesses[0]) >= self.Tol:
            if min([E_guesses[0],E_guesses[2]]) == E_guesses[0]:
                # A new guess for the lattice constant is introduced
                a0_new_guess = 0.67 * a0_guesses[1] + 0.33 * a0_guesses[2]
                # The energy for this new guess is calculated
                atoms.set_cell(a0_new_guess * self.Size * IdentityMatrix,scale_atoms=True)
                E_new_guess = atoms.get_potential_energy()
                # A check for changes in the energy minimum is made and the guesses for a0 and their corrosponding
                # energies are ajusted.
                if min(E_new_guess,min(E_guesses[0:3])) != E_new_guess:
                    a0_guesses[2] = a0_new_guess
                    E_guesses[2] = E_new_guess 
                else:
                    a0_guesses[0] = a0_guesses[1]
                    a0_guesses[1] = a0_new_guess
                    E_guesses[0] = E_guesses[1]
                    E_guesses[1] = E_new_guess
         

            elif min([E_guesses[0],E_guesses[2]]) == E_guesses[2]: 
                # A new guess for the lattice constant is introduced
                a0_new_guess = 0.33 * a0_guesses[0] + 0.67 * a0_guesses[1]
                # The energy for this new guess is calculated
                atoms.set_cell(a0_new_guess * self.Size * IdentityMatrix,scale_atoms=True)
                E_new_guess = atoms.get_potential_energy()
                # A check for changes in the energy minimum is made and the guesses for a0 and their corrosponding
                # energies are ajusted.
                if min(E_new_guess,min(E_guesses[0:3])) != E_new_guess:
                    a0_guesses[0] = a0_new_guess
                    E_guesses[0] = E_new_guess 
                else:
                    a0_guesses[2] = a0_guesses[1]
                    a0_guesses[1] = a0_new_guess
                    E_guesses[2] = E_guesses[1]
                    E_guesses[1] = E_new_guess

        
        # An estimate of the minimum energy can now be found from a second degree polynomial fit through the three 
        # current guesses for a0 and the corresponding values of the energy.
        Poly = numpyfit.polyfit(a0_guesses,E_guesses[0:3],2)

        # The lattice constant corresponding to the lowest energy from the Polynomiel fit is found
        a0 = - Poly[1] / (2 * Poly[0])

        # Now five guesses for a0 and the corresponding energy are evaluated from and around the current a0.
        a0_guesses = a0 * numpy.array([1 - 2 * self.Tol / 5,1 - self.Tol / 5,1,1 + self.Tol / 5,1 + 2 * self.Tol / 5])


        for i in range(5):
            # Changes the lattice constant for the atoms object
            atoms.set_cell(a0_guesses[i] * self.Size * IdentityMatrix,scale_atoms=True)
            # Calculates the energy of the system for the new lattice constant
            E_guesses[i] = atoms.get_potential_energy()

        # The method EquationOfState is now used to find the Bulk modulus and the minimum energy for the system.

        # The volume of the sample for the given a0_guesses
        Vol = (self.Size * a0_guesses )**3
        
        # The equilibrium volume, energy and bulk modulus are calculated
        (Vol0 , E0 , B ) = EquationOfState(Vol.tolist(),E_guesses.tolist()).fit()

        return (Vol0,E0,B,Vol0**(1./3) / self.Size)

    def BulkModulus(self,EMT,PARAMETERS):
        return self.BM_LC_Calculator(EMT,PARAMETERS)[2]

    def LatticeConstant(self,EMT,PARAMETERS):
        return self.BM_LC_Calculator(EMT,PARAMETERS)[3]

    def EvocCalculator(self,EMT,PARAMETERS):
        # Return 0 is used when the method is not desired to be used
        # return 0
        """ This method calculates the vacancy formation energy for the system of atoms. That is, this method calculates 
            the potential energy for a complete system of atoms, then removes an atom and calculates the potential energy
            again. Then the energy for the full system (scaled with the number of atoms in the reduced system) are  
            subtracted from that of the reduced system and the vacancy formation energy, Evoc, is returned. """
        # The atoms object is initialized for the chosen size and type of system
        atoms = FaceCenteredCubic(size=(self.Size,self.Size,self.Size), symbol=self.Element)
        # The EMT calculator given is attached to the atoms object
        atoms.set_calculator(EMT)
        
        # The energy of the full system is calculated
        E_FullSystem = atoms.get_potential_energy()
        
        # an atom is removed from the system
        atoms.pop()
        
        # The energy of the reduced system is calculated
        E_ReducedSystem = atoms.get_potential_energy()

        # The energy of a full system compared to the energy pr atom of the reduced system is calculated and returned
        return E_ReducedSystem - E_FullSystem * len(atoms) / (len(atoms) + 1)

    def EfccEhcpCalculator(self,EMT,PARAMETERS):
        # Return 0 is used when the method is not desired to be used
        # return 0
        """ This method uses the EMT calculator to calculate and return the difference in energy between a system of 
            atoms placed in the HCP and FCC structure. """            
        # The atoms objects are created using the input size and element and the energy calculator 
        # is set to the EMT calculator
        # The Lattice Constants, a,c, for the HCP lattice is here given by the nearest neighbor distance of 
        # the system in the FCC crystal structure, a = dnn, and the ideal relation between a and 
        # c: a/c = sqrt(8/3) => c = dnn / sqrt(8/3)
        a = beta * PARAMETERS[self.Element][1] * Bohr 
        c = a * numpy.sqrt(8./3.)
        
        # The HCP crystal is created, the size of the crystal is defined as 5,5,5, any smaller crystal will result in 
        # Neighborlist errors.
        atoms1 = HexagonalClosedPacked(size=(5,5,5),directions=[[2,-1,-1,0],[0,1,-1,0],[0,0,0,1]], 
                 symbol=self.Element,latticeconstant={'a':a,'c':c})
        atoms1.set_calculator(EMT)
        
        # The FCC crystal is created
        atoms2 = FaceCenteredCubic(size=(self.Size,self.Size,self.Size),symbol=self.Element)
        atoms2.set_calculator(EMT)

        # The energy difference pr atom is calculated and returned
        return atoms1.get_potential_energy() / len(atoms1) - atoms2.get_potential_energy() / len(atoms2)

    def EfccEbccCalculator(self,EMT,PARAMETERS):
        # Return 0 is used when the method is not desired to be used
        # return 0
        """ This method uses the given EMT calculator to calculate and return the difference in energy between a system 
            of atoms placed in the BCC and FCC structure. """ 
        # The atoms objects are created using the input size and element and the energy calculator 
        # is set to the EMT calculator
        # The Lattice Constant for the BCC lattice is here given so that the volume pr atom of the system is 
        # held fixed for the FCC and BCC lattice, that is Vol_pr_atom_BCC = Vol_pr_atom_FCC.
        atoms1 = FaceCenteredCubic(size=(self.Size,self.Size,self.Size),symbol=self.Element)
        atoms1.set_calculator(EMT)
        # The Lattice constant which produces the same volume pr atom for an BCC crystal is calculated
        LCBcc = 1./(2.**(1./3.)) * beta * PARAMETERS[self.Element][1] * Bohr * numpy.sqrt(2)
        atoms2 = BodyCenteredCubic(size=(self.Size,self.Size,self.Size),symbol=self.Element,latticeconstant=LCBcc)
        atoms2.set_calculator(EMT)
        
        # The energy difference pr atom is calculated and returned
        return atoms2.get_potential_energy() / len(atoms2) - atoms1.get_potential_energy() / len(atoms1)

    def C44_Calculator(self,EMT,PARAMETERS):
        # Return 0 is used when the method is not desired to be used
        # return 0
        """ This method uses the given EMT calculator to calculate and return the value of the matrix element C44 for 
            a system of atoms of a given element type. The calculation is done by using that:
            C44 = 1 / Volume * d^2/depsilon^2 (E_system) where epsilon is the displacement in one direction of the 
            system along one axis diveded by the highth of the system. """

        # An atom object is created and the calculator attached
        atoms = FaceCenteredCubic(size=(self.Size,self.Size,self.Size),symbol=self.Element)
        atoms.set_calculator(EMT)

        # The volume of the sample is calculated
        Vol = atoms.get_volume()    

        # The value of the relative displacement, epsilon, is set 
        epsilon = 1. / 1000

        # The matrix used to change the unitcell by n*epsilon is initialized
        LMM = numpy.array([[1,0,-10*epsilon],[0,1,0],[0,0,1]])
        
        # The original unit cell is conserved
        OCell = atoms.get_cell()

        # The array for storing the energies is initialized
        E_calc = numpy.zeros(20)

        # The numerical value of C44 is calculated
        for i in range(20):
            # The new system cell based on the pertubation epsilon is set
            atoms.set_cell(numpy.dot(OCell, LMM), scale_atoms=True)
            # The energy of the system is calculated
            E_calc[i] = atoms.get_potential_energy()
            # The value of LMM is updated
            LMM[0,2] += epsilon
            
        # A polynomial fit is made for the energy as a function of epsilon 

        # The displaced axis is defined
        da = numpy.arange(-10, 10) * epsilon
        
        # The fit is made
        Poly = numpyfit.polyfit(da,E_calc,2)

        # Now C44 can be estimated from this fit by the second derivative of Poly = a * x^2 + b * x + c , multiplied
        # with 1 / (2 * Volume) of system
        C44 = 2. / Vol * Poly[0]
        
        return C44

    def C11_Calculator(self,EMT,PARAMETERS):
        # Return 0 is used when the method is not desired to be used
        # return 0
        """ This method uses the given EMT calculator to calculate and return the value of the matrix element C11 for 
            a system of atoms of a given element type. The calculation is done by using that:
            C11 = 1 / Volume * d^2/depsilon^2 (E_system). """
        # An atom object is created and the calculator attached
        atoms = FaceCenteredCubic(size=(self.Size,self.Size,self.Size),symbol=self.Element)
        atoms.set_calculator(EMT)

        # The volume of the sample is calculated
        Vol = atoms.get_volume()    

        # The value of the relative displacement, epsilon, is set 
        epsilon = 1. / 1000

        # The matrix used to change the unitcell by n*epsilon is initialized
        LMM = numpy.array([[1-10*epsilon,0,0],[0,1,0],[0,0,1]])
        
        # The original unit cell is conserved
        OCell = atoms.get_cell()

        # The array for storing the energies is initialized
        E_calc = numpy.zeros(20)

        # The numerical value of C11 is calculated
        for i in range(20):
            # The new system cell based on the pertubation epsilon is set
            atoms.set_cell(numpy.dot(OCell, LMM), scale_atoms=True)
            # The energy of the system is calculated
            E_calc[i] = atoms.get_potential_energy()
            # The value of LMM is updated
            LMM[0,0] += epsilon
        

        # A second degree polynomial fit is made for the energy as a function of epsilon 

        # The displaced axis is defined
        da = numpy.arange(-10, 10) * epsilon
        
        # The fit is made
        Poly = numpyfit.polyfit(da,E_calc,2)
        
        # Now C11 can be estimated from this fit by the second derivative of Poly = a * x^2 + b * x + c, multiplied
        # with 1 / (2 * Volume) of system
        C11 = 2. / Vol * Poly[0]
        
        return C11

    def E_cohCalculator(self,EMT,PARAMETERS):
        # Return 0 is used when the method is not desired to be used
        # return 0
        """ Calculates the Cohesive energy of a system of atoms using the EMT calculator specified. """

        # As the EMT calculator calculates the energy of the system such that the energy of the individual atoms in 
        # their equilibrium distance in the crystal is zero it is only needed to calculate the energy of a single atom
        # in an empty system.
        
        # The crystal is created
        atoms = FaceCenteredCubic(size=(self.Size,self.Size,self.Size),symbol=self.Element)
        
        # a single atom is taken out of the crystal
        atoms2 = atoms[[0,]]

        # The calculator is attached to the atoms objects
        atoms.set_calculator(EMT)
        atoms2.set_calculator(EMT)

        # The energy difference between the atom alone in vacuum and in the crystal structure is calculated and returned
        return atoms2.get_potential_energy() - atoms.get_potential_energy() / len(atoms)
        

    def SurfaceEnergy_Calculator(self,EMT,PARAMETERS):
        """ The Method calculates and returns the surface energy for the given element along the [0,0,1] and [1,1,1] 
            directions in the FCC crystal structure. """
        # The size of the crystals are set " NEED TO THINK ABOUT THESE LATER! " 
        S001 = 3,3,5
        S111 = 5,5,5
        # The surfaces (slabs) are created (pbc=(1,1,0) creates periodic boudry conditions 
        # in two of three directions and thus leaves the last direction as two surfaces.
        Surface001 = FaceCenteredCubic(size=S001,symbol=self.Element,pbc=(1,1,0))
        Surface111 = FaceCenteredCubic(size=S111,directions=[[1,-1,0],[1,1,-2],[1,1,1]],
                                       symbol=self.Element,pbc=(1,1,0))
        Surface001.set_calculator(EMT)
        Surface111.set_calculator(EMT)

        # A structural relaxsation is run for the surface crystal in order to secure 
        # the correct structure of the crystal.
        dyn001 = BFGS(Surface001, trajectory='relaxedsurface001.traj')
        dyn111 = BFGS(Surface111, trajectory='relaxedsurface111.traj')
        dyn001.run(fmax=0.01)
        dyn111.run(fmax=0.01)        

        # The referance bulk crystals are created
        Bulk001 = FaceCenteredCubic(size=S001,symbol=self.Element)
        Bulk111 = FaceCenteredCubic(size=S111,directions=[[1,-1,0],[1,1,-2],[1,1,1]],symbol=self.Element)

        # The calculator is assigned
        Bulk001.set_calculator(EMT)
        Bulk111.set_calculator(EMT)

        # The surface area is calculated
        # The cross product between the x and y axis in the crystal is determined
        Cross001 =  numpy.cross(Bulk001.get_cell()[:,0],Bulk001.get_cell()[:,1])
        Cross111 =  numpy.cross(Bulk111.get_cell()[:,0],Bulk111.get_cell()[:,1])
        # The area of the surface is determined from the formular A = |X x Y|. 
        area001 = numpy.sqrt(numpy.dot(Cross001,Cross001))
        area111 = numpy.sqrt(numpy.dot(Cross111,Cross111))

        # The surface energy is calculated and returned (two surfaces are present in 
        # SurfaceRelaxed)
        return ( (Surface001.get_potential_energy() - Bulk001.get_potential_energy()) / 2 / area001,
                 (Surface111.get_potential_energy() - Bulk111.get_potential_energy()) / 2 / area111)

    def Surface001(self,EMT,PARAMETERS):
        """ The Method calculates and returns the surface energy for the given element along the [0,0,1] direction in   
            the FCC crystal structure. """
        # The size of the crystals are set: 
        S001 = 3,3,5
         # The surfaces (slabs) are created (pbc=(1,1,0) creates periodic boudry conditions 
        # in two of three directions and thus leaves the last direction as two surfaces.
        Surface001 = FaceCenteredCubic(size=S001,symbol=self.Element,pbc=(1,1,0))
        Surface001.set_calculator(EMT)

        # A structural relaxsation is run for the surface crystal in order to secure 
        # the correct structure of the crystal.
        dyn001 = BFGS(Surface001, logfile=None)
        dyn001.run(fmax=0.01)

        # The referance bulk crystals are created
        Bulk001 = FaceCenteredCubic(size=S001,symbol=self.Element)

        # The calculator is assigned
        Bulk001.set_calculator(EMT)
  
        # The surface area is calculated
        # The cross product between the x and y axis in the crystal is determined
        Cross001 =  numpy.cross(Bulk001.get_cell()[:,0],Bulk001.get_cell()[:,1])
        # The area of the surface is determined from the formular A = |X x Y|. 
        area001 = numpy.sqrt(numpy.dot(Cross001,Cross001))
    
        # The surface energy is calculated and returned (two surfaces are present in 
        # SurfaceRelaxed)
        return ( (Surface001.get_potential_energy() - Bulk001.get_potential_energy()) / 2 / area001)
   
    def Surface110(self,EMT,PARAMETERS):
        """ The Method calculates and returns the surface energy for the given element along the [1,1,0] direction in   
            the FCC crystal structure. """
        # The size of the crystals are set: 
        S110 = 3,5,5
         # The surfaces (slabs) are created (pbc=(1,1,0) creates periodic boudry conditions 
        # in two of three directions and thus leaves the last direction as two surfaces.
        Surface110 = FaceCenteredCubic(size=S110,directions=[[1,-1,0],[0,0,1],[1,1,0]],
                                       symbol=self.Element,pbc=(1,1,0))
        Surface110.set_calculator(EMT)

        # A structural relaxsation is run for the surface crystal in order to secure 
        # the correct structure of the crystal.
        dyn001 = BFGS(Surface110, logfile=None)
        dyn001.run(fmax=0.01)

        # The referance bulk crystals are created
        Bulk110 = FaceCenteredCubic(size=S110,directions=[[1,-1,0],[0,0,1],[1,1,0]],
                                    symbol=self.Element)

        # The calculator is assigned
        Bulk110.set_calculator(EMT)
  
        # The surface area is calculated
        # The cross product between the x and y axis in the crystal is determined
        Cross110 =  numpy.cross(Bulk110.get_cell()[:,0],Bulk110.get_cell()[:,1])
        # The area of the surface is determined from the formular A = |X x Y|. 
        area110 = numpy.sqrt(numpy.dot(Cross110,Cross110))
    
        # The surface energy is calculated and returned (two surfaces are present in 
        # SurfaceRelaxed)
        return ( (Surface110.get_potential_energy() - Bulk110.get_potential_energy()) / 2 / area110)
  
    def Surface111(self,EMT,PARAMETERS):
        """ The Method calculates and returns the surface energy for the given element along the [1,1,1] 
            direction in the FCC crystal structure. """
        # The size of the crystals are set " NEED TO THINK ABOUT THESE LATER! " 
        S111 = 4,4,4
        # The surfaces (slabs) are created (pbc=(1,1,0) creates periodic boudry conditions 
        # in two of three directions and thus leaves the last direction as two surfaces.
        Surface111 = FaceCenteredCubic(size=S111,directions=[[1,-1,0],[1,1,-2],[1,1,1]],
                                       symbol=self.Element,pbc=(1,1,0))
        Surface111.set_calculator(EMT)

        # A structural relaxsation is run for the surface crystal in order to secure 
        # the correct structure of the crystal.
        dyn111 = BFGS(Surface111, logfile=None)
        dyn111.run(fmax=0.01)        

        # The referance bulk crystals are created
        Bulk111 = FaceCenteredCubic(size=S111,directions=[[1,-1,0],[1,1,-2],[1,1,1]],symbol=self.Element)

        # The calculator is assigned
        Bulk111.set_calculator(EMT)

        # The surface area is calculated
        # The cross product between the x and y axis in the crystal is determined
        Cross111 =  numpy.cross(Bulk111.get_cell()[:,0],Bulk111.get_cell()[:,1])
        # The area of the surface is determined from the formular A = |X x Y|. 
        area111 = numpy.sqrt(numpy.dot(Cross111,Cross111))

        # The surface energy is calculated and returned (two surfaces are present in 
        # SurfaceRelaxed)
        return ( (Surface111.get_potential_energy() - Bulk111.get_potential_energy()) / 2 / area111)


    def L12_BM_LC_Calculator(self,EMT,PARAMETERS):
        # Return 0 is used when the method is not desired to be used
        # return 0
        """ The method L12_BM_LC_Calculator uses the EMT calculator to find the lattice constant which gives 
            the lowest possible energy for an alloy of two elements and from a polynomial fit of the energy as a 
            function of the lattice constant it calculates and returns the Bulk modulus, the volume of the system at 
            the lowest possible energy and this energy """
        # Three values for the lattice constant, a0, are chosen. The highest and lowest value is chosen so that it     
        # is certain that these are to high/low respectively compared to the "correct" value. This is done by using the 
        # experimental value of the lattice constant, a0_exp, as the middle value and to high/low values are then: 
        # a0_exp +- a0_mod, a0_mod = a0_exp/10
        a0_exp = ( ( PARAMETERS[self.Element[0]][1] + 3 * PARAMETERS[self.Element[1]][1] ) 
                     * numpy.sqrt(2) * beta * Bohr / 4 )
        a0_mod = a0_exp * 0.20 
        a0_guesses = numpy.array([a0_exp - a0_mod,a0_exp,a0_exp + a0_mod])
    
        # An atoms object consisting of atoms of the chosen element is initialized 
        atoms = L1_2(size=(self.Size,self.Size,self.Size),symbol=self.Element,latticeconstant=a0_exp)
        atoms.set_calculator(EMT)

        # An identity matrix for the system is saved to a variable
        IdentityMatrix = numpy.array([[1,0,0],[0,1,0],[0,0,1]])

        # An array for the energy for the chosen guesses for a0 is initialized:
        E_guesses = numpy.zeros(5)

        # The energies are calculated
        for i in range(3):
        # Changes the lattice constant for the atoms object
            atoms.set_cell(a0_guesses[i] * self.Size * IdentityMatrix,scale_atoms=True)
            # Calculates the energy of the system for the new lattice constant
            E_guesses[i] = atoms.get_potential_energy()
           
        # Bisection is used in order to find a small interval of the lattice constant for the minimum of the energy (an 
        # abitrary interval length is chosen), This is possible because we are certian from theory that there is only 
        # one minimum of the energy function in the interval of interest and thus we wont "fall" into a local minimum 
        # and stay there by accident.
        while (a0_guesses[2]-a0_guesses[0]) >= self.Tol:
            if min([E_guesses[0],E_guesses[2]]) == E_guesses[0]:
                # A new guess for the lattice constant is introduced
                a0_new_guess = 0.67 * a0_guesses[1] + 0.33 * a0_guesses[2]
                # The energy for this new guess is calculated
                atoms.set_cell(a0_new_guess * self.Size * IdentityMatrix,scale_atoms=True)
                E_new_guess = atoms.get_potential_energy()
                # A check for changes in the energy minimum is made and the guesses for a0 and their corrosponding
                # energies are ajusted.
                if min(E_new_guess,min(E_guesses[0:3])) != E_new_guess:
                    a0_guesses[2] = a0_new_guess
                    E_guesses[2] = E_new_guess 
                else:
                    a0_guesses[0] = a0_guesses[1]
                    a0_guesses[1] = a0_new_guess
                    E_guesses[0] = E_guesses[1]
                    E_guesses[1] = E_new_guess
         

            elif min([E_guesses[0],E_guesses[2]]) == E_guesses[2]: 
                # A new guess for the lattice constant is introduced
                a0_new_guess = 0.33 * a0_guesses[0] + 0.67 * a0_guesses[1]
                # The energy for this new guess is calculated
                atoms.set_cell(a0_new_guess * self.Size * IdentityMatrix,scale_atoms=True)
                E_new_guess = atoms.get_potential_energy()
                # A check for changes in the energy minimum is made and the guesses for a0 and their corrosponding
                # energies are ajusted.
                if min(E_new_guess,min(E_guesses[0:3])) != E_new_guess:
                    a0_guesses[0] = a0_new_guess
                    E_guesses[0] = E_new_guess 
                else:
                    a0_guesses[2] = a0_guesses[1]
                    a0_guesses[1] = a0_new_guess
                    E_guesses[2] = E_guesses[1]
                    E_guesses[1] = E_new_guess

        
        # An estimate of the minimum energy can now be found from a second degree polynomial fit through the three 
        # current guesses for a0 and the corresponding values of the energy.
        Poly = numpyfit.polyfit(a0_guesses,E_guesses[0:3],2)

        # The lattice constant corresponding to the lowest energy from the Polynomiel fit is found
        a0 = - Poly[1] / (2 * Poly[0])

        # Now five guesses for a0 and the corresponding energy are evaluated from and around the current a0.
        a0_guesses = a0 * numpy.array([1 - 2 * self.Tol / 5,1 - self.Tol / 5,1,1 + self.Tol / 5,1 + 2 * self.Tol / 5])


        for i in range(5):
            # Changes the lattice constant for the atoms object
            atoms.set_cell(a0_guesses[i] * self.Size * IdentityMatrix,scale_atoms=True)
            # Calculates the energy of the system for the new lattice constant
            E_guesses[i] = atoms.get_potential_energy()

        # The method EquationOfState is now used to find the Bulk modulus and the minimum energy for the system.

        # The volume of the sample for the given a0_guesses
        Vol = (self.Size * a0_guesses )**3
        
        # The equilibrium volume, energy and bulk modulus are calculated
        (Vol0 , E0 , B ) = EquationOfState(Vol.tolist(),E_guesses.tolist()).fit()



        return (Vol0,E0,B,Vol0**(1./3) / self.Size)

    def L12_BulkModulus(self,EMT,PARAMETERS):
        return self.L12_BM_LC_Calculator(EMT,PARAMETERS)[2]

    def L12_LatticeConstant(self,EMT,PARAMETERS):
        return self.L12_BM_LC_Calculator(EMT,PARAMETERS)[3]

  




