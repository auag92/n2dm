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
        newsupercell = SuperCell(scaling * basis, [1,1,1])
        atoms.DeformAtoms(newsupercell)
        vol = newsupercell.GetVolume()
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
atoms.GetEnergy()  # Forces creation of a neigborlist

# Find correct lattice constant
for i in (0,1):
    (energyfit, pressurefit) = makefits(atoms)

    dilation = -pressurefit[0][0]/pressurefit[0][1]
    print "Optimizing lattice constant:", latticeconstant, "->", latticeconstant*(1+dilation/3)
    latticeconstant = latticeconstant*(1+dilation/3)
    basis = (1+dilation/3) * basis
    supercell = SuperCell(basis, [1,1,1])
    originalvolume = supercell.GetVolume()
    atoms.DeformAtoms(supercell)

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

newsupercell = SuperCell(basis, [1,1,1])
atoms.DeformAtoms(newsupercell)
stress = atoms.GetStress()

p = -(stress[0]+stress[1]+stress[2])/3
p0 = p
print ""
print "Pressure at T=0 K:", p * eV_per_cubicaangstrom * 1e-6, "MPa"

def getTemp(atoms):
    p = atoms.GetCartesianMomenta()
    m = atoms.GetMasses()
    ekin = 0.5*sum(sum(transpose(p*p))/m)/len(atoms)
    return 2.0/(3.0*kB) * ekin

def pressureattemperature(atoms, temp):
    mover = Langevin(atoms, lambda:None, 0.4, 3*temp*kB, 0.01)
    while getTemp(atoms) < temp:
        mover(5)
        sys.stderr.write("+")
    sys.stderr.write("\n")
    print "Reached temperature:", getTemp(atoms),"K"
    mover.SetTemperature(kB*temp)
    for i in range(25):
        mover(5)
        sys.stderr.write("-")
    pressures = []
    temps = []
    for i in range(50):
        mover(5)
        sys.stderr.write(".")
        stress = atoms.GetStress()
        temps.append(getTemp(atoms))
        pressures.append(-(stress[0]+stress[1]+stress[2])/3.0)
    sys.stderr.write("\n")
    p = sum(array(pressures))/len(pressures)
    t = sum(array(temps))/len(temps)
    return (t, p)

tp = []
for T in range(0,601,100):
    (t, p) = pressureattemperature(atoms, T)
    print "T = %f (%d): p = %f\n" % (t, T, p)
    tp.append((t,p))

pressurefit = polynomialLeastSquaresFit((0.0, 0.0), tp)
print pressurefit
const = pressurefit[0][1] * eV_per_cubicaangstrom
coeff = const / bmod
print "Thermal expansion coefficient:", coeff
print "Textbook value:               ", expansioncoeff
    


