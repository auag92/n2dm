from Numeric import *
from Simulations.Asap import *
#from IO.AtomsInCmFile import *
from Setup.Lattice.FCC111Ortho import *
from Structures.ChemicalElements import Copper, Silver
import Structures.ChemicalElements.AtomicWeight
import Structures.ChemicalElements.CrystalStructure
from Structures.IonDynamics import VelocityVerlet, Langevin
from Scientific.Functions.LeastSquares import *
import sys

def makefits(atoms):
    energies = []
    pressures = []
    for epsilon in (-1, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0):
        scaling = 1 + 0.01 * epsilon
        atoms.GetUnitCell().SetBasis(scaling * basis)
        vol = atoms.GetUnitCell().GetVolume()
        atoms.CalculateEnergiesForcesAndStresses()
        energy = atoms.GetEnergy()
        stress = atoms.GetStress()
        #print ""
        #print epsilon, energy
        #print stress
        energies.append(((vol-originalvolume)/originalvolume, energy/vol))
        pressures.append(((vol-originalvolume)/originalvolume,
                          -(stress[0]+stress[1]+stress[2])/3))

    energyfit = polynomialLeastSquaresFit((0.0, 0.0, 0.0), energies)
    pressurefit = polynomialLeastSquaresFit((0.0, 0.0), pressures)
    return (energyfit, pressurefit)

element = Copper
bulkmodulus = 134.3
expansioncoeff = 50e-6
assert element.GetProperty('CrystalStructure')['Symmetry'] == "FCC"
latticeconstant = element.GetProperty('CrystalStructure')['LatticeConstant']
initial = FCC111Ortho((21,20,24), element, latticeconstant)
basis = initial.GetCoordinateBasis()
supercell = SuperCell(basis, [1,1,1])
originalvolume = supercell.GetVolume()

atoms = Atoms(initial, supercell, potential=EMT())
atoms.SetCartesianMomenta(zeros((len(atoms),3), Float))
atoms.CalculateEnergies()  # Forces creation of a neigborlist

# Find correct lattice constant
for i in (0,1):
    (energyfit, pressurefit) = makefits(atoms)

    dilation = -pressurefit[0][0]/pressurefit[0][1]
    print "Optimizing lattice constant:", latticeconstant, "->", latticeconstant*(1+dilation/3)
    latticeconstant = latticeconstant*(1+dilation/3)
    basis = (1+dilation/3) * basis
    atoms.GetUnitCell().SetBasis(basis)
    originalvolume = atoms.GetUnitCell().GetVolume()

# Now do the calculation.
(energyfit, pressurefit) = makefits(atoms)
kB = 8.61734e-5  # Boltzmann's constant in eV/K
eV_per_cubicaangstrom = 1.602e-19/1e-30

print ""
print energyfit
print pressurefit

print ""
bmod1 = 2*energyfit[0][2] * eV_per_cubicaangstrom
print "Bulk modulus (from energies): ", bmod1 * 1e-9, "GPa"
bmod2 = -pressurefit[0][1] * eV_per_cubicaangstrom
print "Bulk modulus (from pressures):",  bmod2 * 1e-9, "GPa"
print "Textbook value:               ", bulkmodulus, "GPa"
bmod = (bmod1+bmod2)/2

ReportTiming()
