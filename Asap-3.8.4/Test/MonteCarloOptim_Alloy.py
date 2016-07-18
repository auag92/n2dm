from numpy import *
from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from ase import data
from asap3.testtools import *
import time

random.seed([114543, 8487])

listcutoff = 4.98409   # The EMT value for Cu
element = "Cu"
element2 = "Ag"
latconst = data.reference_states[data.atomic_numbers[element]]['a']
elementnumber = data.atomic_numbers[element]
element2number = data.atomic_numbers[element2]

def checkenergies(a, a2):
    err = 0
    if 0:
        print "    Checking Cartesian positions"
        err = abs((a.get_positions() - a2.get_positions()).ravel())
        idx = argmax(err)
        at, co = idx/3, idx%3
        print "      err =", err[idx], "at", idx, ":", at, co
        ReportTest("      Worst Cartesian position (%d,%d)" % (at, co),
                   a.get_positions()[at,co],
                   a2.get_positions()[at,co], 1e-10)
    print "    Checking energies"
    e1 = a.get_potential_energies()
    e2 = a2.get_potential_energies()
    err = e1 - e2
    idx = argmax(err)
    ReportTest("      Worst energy (%s)" % (idx,), e1[idx], e2[idx], 1e-5)

for i in range(10):
    print
print_version(1)
#print "RandomArray seeds:", RandomArray.get_seed()
if len(sys.argv) >= 2 and sys.argv[1] == '--pbc':
    if len(sys.argv) < 3:
        raise ValueError("No argument to --pbc")
    pbcch = sys.argv[2]
    if len(pbcch) != 3:
        raise ValueError("Argument to pbc should be a string of three digits.")
    pbc = []
    for i in (0,1,2):
        pbc.append(int(pbcch[i]))
    print "Periodic boundary conditions forced!"
else:
    pbc = tuple(random.randint(0,2,3))
print "Periodic boundaries:", pbc

#Verbose(1)

atoms = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]], size=(10,10,10),
                          symbol=element, pbc=pbc)
atoms2 = atoms.copy()
atoms = MonteCarloAtoms(atoms)
print "Number of atoms:", len(atoms)
# Replace 10% of the atoms with element2
newz = where(greater(random.random((len(atoms),)), 0.9),
             element2number, elementnumber)
atoms.set_atomic_numbers(newz)
atoms2.set_atomic_numbers(newz)

atoms.set_calculator(MonteCarloEMT())
atoms2.set_calculator(EMT())
atoms.get_potential_energy()

for e in ((element,elementnumber), (element2,element2number)):
    print ("Number of %s atoms (Z=%d): %d"
           % (e[0], e[1], sum(equal(atoms.get_atomic_numbers(), e[1]))))

print 
print "Testing perturbations of single atoms"

magnarr = array((0.0, 0.01, 0.1, 1.0, 3.0, 10.0))
numarr = array((1, 3, 10, len(atoms)/2, -3))
pick1 = argsort(random.random((len(magnarr),)))
for magn in take(magnarr, pick1):
    pick2 = argsort(random.random((len(numarr),)))
    for number in take(numarr, pick2):
        # Pick number random atoms.
        if number < 0:
            # Find N neighboring atoms and perturb them
            number = -number
            nblist = atoms.get_calculator().get_neighborlist()
            n = 0
            while len(nblist[n]) < number:
                n += 1
            pick = concatenate([array((n,)), nblist[n][:number-1]])
            del nblist
        else:
            pick = argsort(random.random((len(atoms),)))[:number]
        if number > 15:
            s = "<<< %s atoms >>>" % (number,)
        else:
            s = str(list(pick))
        print "  dr = %.3f: %d atom(s): %s" % (magn, number, s)
        for i in pick:
            dr = random.standard_normal(3)
            dr *= magn/sqrt(dot(dr,dr))
            atom = atoms[i]
            atom.position = atom.position + dr
            atom = atoms2[i]
            atom.position = atom.position + dr

        checkenergies(atoms, atoms2)

print
print "Testing alchemy"

pick2 = argsort(random.random((len(numarr),)))
for number in take(numarr, pick2):
    # Pick number random atoms.
    if number < 0:
        # Find N neighboring atoms and perturb them
        number = -number
        nblist = atoms.get_calculator().get_neighborlist()
        n = 0
        while len(nblist[n]) < number:
            n += 1
        pick = concatenate([array((n,)), nblist[n][:number-1]])
        del nblist
    else:
        pick = argsort(random.random((len(atoms),)))[:number]
    if number > 15:
        s = "<<< %s atoms >>>" % (number,)
    else:
        s = str(list(pick))
        print "  Alchemy, %d atoms: %s" % (number, s)
        for i in pick:
            atom = atoms[i]
            atom2 = atoms2[i]
            z = atom.number
            assert z == atom2.number
            if (z == elementnumber):
                atom.number = element2number
                atom2.number = element2number
            else:
                assert z == element2number
                atom.number = elementnumber
                atom2.number = elementnumber
        checkenergies(atoms, atoms2)
        
ReportTest.Summary()

