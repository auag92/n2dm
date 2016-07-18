import Asap
from ASE.Dynamics.Langevin import Langevin
from ASE.Dynamics.MDMin import MDMin

#from Structures.IonDynamics import RobustLangevin
from Asap.Filters.HomogeneousStrain import HomogeneousStrain
import Numeric as num
import LinearAlgebra as LA
from Asap.Utilities.GetPotential import GetPotential
from Asap.Utilities.StrainVectorMatrix import *
import Scientific.MPI

if Scientific.MPI.world.size > 1:
    parallel = 1
else:
    parallel = 0


def ApplyStrain(atoms, strain):
    #parallel = hasattr(atoms, "parallel")

    if strain.shape == (6,):
        strainMatrix = MakeStrainMatrix(strain)
    elif strain.shape == (3,3):
        strainMatrix = strain
    else:
        print "bad array shape"
        raise StandardError

    defGrad = num.identity(3).astype(num.Float) + strainMatrix
    superCellVectors = atoms.GetUnitCell().GetBasis()
    strainedSVectors =  matrixmultiply(superCellVectors, transpose(defGrad))

    if getattr(atoms, "parallel ", 0):
        self.comm = Scientific.MPI.world
        # Be absolutely sure to keep thing exactly identical everywhere
        self.comm.broadcast(strainedSVectors, 0)

    atoms.GetUnitCell().SetBasis(strainedSVectors)

    # make the original vector available in order to undo the strain
    return superCellVectors


def ReducePressure(atoms, numIterations, bulkModulus):
    stress = atoms.GetStress()
    pressure = -(stress[1] + stress[2] + stress[2])/3.
    
    for idx in xrange(numIterations):
        volStrain = pressure / bulkModulus
        strainMtx = num.identity(3).astype(num.Float) *  volStrain/3.
    
        ApplyStrain(atoms, strainMtx)
        stress = atoms.GetStress()
        pressure = -(stress[1] + stress[2] + stress[2])/3.


class BoxMinimizer:
    """Note: This should work in parallel as is. If convergence checks are
    added, it will be important to make them consistent for parallel
    simulations.
    """
    def __init__(self, qmSteps, boxSteps, majSteps, tol, boxTimeStep = 0.001):
        self.atoms = None
        self.qmSteps = qmSteps
        self.boxSteps = boxSteps
        self.majSteps = majSteps
        self.tol = tol
        self.debug = 0
        self.lastEnergy = None
        self.boxTimeStep = boxTimeStep
        self.selector = None
        self.temperature = None
        self.targetStress = None
        self.timestep = 0.2
        self.oneBoxStep = 0
        self.bulkModulus = None
        self.shearModulus = None
        self.thermAveStress = None

    def SetQMSteps(self, qmSteps):
        self.qmSteps = qmSteps

    def SetTemperature(self, temp):
        self.temperature = temp
        
    def SetDebug(self, debug):
        self.debug = debug

    def SetSelector(self,selector):
        self.selector = selector

    def SetTargetStress(self, targetStress):
        self.targetStress = targetStress

    def SetOneBoxStep(self, oneBS, bm=None, sm=None):
        self.oneBoxStep = oneBS
        self.bulkModulus = bm
        self.shearModulus  = sm        
        self.ComputeComplianceMatrix()

    def ComputeComplianceMatrix(self):
        print "In BoxMinimizer.ComputeComplianceMatrix"
        bm = self.bulkModulus
        mu = self.shearModulus
        Lambda = bm - (2./3) * mu

        ecMatrix = num.zeros((6,6),num.Float)

        for idx in xrange(0,3):
            ecMatrix[idx, idx] = Lambda + 2*mu

        ecMatrix[0,1] = Lambda
        ecMatrix[1,0] = Lambda
        ecMatrix[0,2] = Lambda
        ecMatrix[2,0] = Lambda
        ecMatrix[1,2] = Lambda
        ecMatrix[2,1] = Lambda
        for idx in xrange(3,6):
            ecMatrix[idx, idx] = mu


        self.compliance = LA.inverse(ecMatrix)


    def GetConvergence(self):
        energy, forces, stress = self.atoms.GetPotentialEnergy(), self.atoms.GetCartesianForces(), self.atoms.GetStress()
        pressure = -(stress[0] + stress[1] + stress[2])

        print "ave. force:",num.sqrt(num.dot(forces.flat, forces.flat) / len(forces))
        print "pressure",pressure
        #print stress[0],stress[1],stress[2],stress[3],stress[4],stress[5]
        print "energy",energy,
        if self.lastEnergy is not None:
            print "decrease by",self.lastEnergy - energy
        else:
            print

        self.lastEnergy = energy
        
    def EnergyFromBoxShape(self, strain):
        if self.atoms == None:
            return

        if self.debug:
            unstrainedEnergy = self.atoms.GetPotentialEnergy()

        unstrainedSCVs = ApplyStrain(self.atoms, strain)

        if self.debug >= 2:
            print "volume:",LA.determinant(self.atoms.GetUnitCell())
    
        energy = self.atoms.GetPotentialEnergy()

        # restore original state
        self.atoms.GetUnitCell().SetBasis(unstrainedSCVs)

        if self.debug:
            if abs(unstrainedEnergy - self.atoms.GetPotentialEnergy()) > 1.e-10:
                print unstrainedEnergy,self.atoms.GetPotentialEnergy()
                raise StandardError

        
        return energy

    def OneStepZeroStress(self):
        if self.thermAveStress is None:
            stress = self.atoms.GetStress()
        else:
            stress = self.thermAveStress

        print "OneStepZeroStress:stress",stress
        elasticStrain = num.matrixmultiply(self.compliance, stress)

        strainChange = - elasticStrain
        print "strainChange",strainChange
        # want to subtract this much strain
        ApplyStrain(self.atoms, strainChange)

    def Minimize(self, atoms):
        self.atoms = atoms
        nAtoms = self.atoms.GetNumberOfAtoms()

        if self.debug:
            self.GetConvergence()

        if self.temperature is None:
            self.quickmin = MDMin(self.atoms, self.timestep)
            # pre-minimization
            self.MinimizeBoxQM(2 * self.boxSteps)
            self.MinimizeAtoms(2 * self.qmSteps)
            if self.debug:
                self.GetConvergence()
        else:
            friction = 0.1
            if parallel:
                self.langevin = RobustLangevin(atoms, self.timestep, self.temperature, friction)
            else:
                self.langevin = Langevin(atoms, self.timestep, self.temperature, friction)


        for idx in xrange(self.majSteps):
            if self.oneBoxStep:
                self.OneStepZeroStress()
            else:
                self.MinimizeBoxQM(self.boxSteps)

            if self.temperature is None:
                self.MinimizeAtoms(self.qmSteps)
            else:
                self.ThermalizeAtoms(self.qmSteps)
            if self.debug:
                self.GetConvergence()

        if self.temperature is None:
            self.MinimizeAtoms(2 * self.qmSteps)
        minEnergy = self.atoms.GetPotentialEnergy()

        if self.debug:
            print "final energy",minEnergy
        return minEnergy
        
    def MinimizeAtoms(self, nSteps):
        self.quickmin.Run(nSteps)
        
    def MinimizeBoxQM(self, nSteps):
        supercell = HomogeneousStrain(self.atoms, self.selector, self.targetStress)
        if self.debug:
            stressBefore = self.atoms.GetStress()
            energyBefore =  self.atoms.GetPotentialEnergy()
        boxMin = MDMin(supercell, dt=self.boxTimeStep)
        boxMin.Run(nSteps)
        if self.debug:
            print "stress change from box-min",self.atoms.GetStress() - stressBefore
            print "energy change from box-min", self.atoms.GetPotentialEnergy() - energyBefore

    def ThermalizeAtoms(self, nSteps):
        aveStress = num.zeros(6, num.Float)
        if self.oneBoxStep:
            nMinor = 10
            nMaj = nSteps / nMinor
            for idx in xrange(nMaj):
                self.langevin.Run(nMinor)
                aveStress += self.atoms.GetStress()
            aveStress /= nMaj
            self.thermAveStress = aveStress
        else:
            self.langevin.Run(nSteps)
    
