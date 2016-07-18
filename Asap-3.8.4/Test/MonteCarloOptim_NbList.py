from numpy import *
from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from ase import data
from asap3.testtools import *
import time

#random.seed([114173, 98658])

listcutoff = 4.98409   # The EMT value for Cu
element = "Cu"
latconst = data.reference_states[data.atomic_numbers[element]]['a']

def checklist(nb, a, cut):
    error = 0
    #nb2 = FullNeighborList(a, cut)
    a2 = Atoms(a)
    #nb2 = NeighborList(cut, a2, 0, full=True)
    nb2 = FullNeighborList(cut, a2, 0)
    print "    Checking Cartesian positions"
    err = abs((a.get_positions() - a2.get_positions()).ravel())
    idx = argmax(err)
    at, co = idx/3, idx%3
    print "      err =", err[idx], "at", idx, ":", at, co
    ReportTest("      Worst Cartesian position (%d,%d)" % (at, co),
               a.get_positions()[at,co],
               a2.get_positions()[at,co], 1e-10)
    #print "    Checking internal positions"
    #err = abs((a.GetInternalPositions() - a2.GetInternalPositions()).flat)
    #idx = argmax(err)
    #at, co = idx/3, idx%3
    #print "      err =", err[idx], "at", idx, ":", at, co
    #ReportTest("      Worst internal position (%d,%d)" % (at, co),
    #           a.GetInternalPositions()[at,co],
    #           a2.GetInternalPositions()[at,co], 1e-10)
    print "    Checking full list"
    for i in range(len(a)):
        l1 = list(nb[i])
        l1.sort()
        l2 = list(nb2[i])
        l2.sort()
        if not l1 == l2:
            ReportTest.BoolTest("NB lists for atom %d should be identical"
                                % (i,), False)
            print "nb1[%d] = %s" % (i, str(l1))
            print "nb2[%d] = %s" % (i, str(l2))
            print "pos1[%d] = %s" % (i, str(a.get_positions()[i]))
            print "pos2[%d] = %s" % (i, str(a2.get_positions()[i]))
            error += 1
            if error > 3:
                print "Too many errors, giving up!"
                raise RuntimeError, "Too many errors, giving up!"

for i in range(10):
    print
print_version(1)
#print "RandomArray seeds:", RandomArray.get_seed()
pbc = tuple(random.random_integers(0,1,3))
print "Periodic boundaries:", pbc

atoms = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]], size=(10,10,10),
                          symbol=element, pbc=pbc, debug=0)

nblist = NeighborList(listcutoff, atoms, 0, full=True)

# This would only work with full PBC
#ReportTest("Total neighborlist length", totallength, 42*len(atoms), 0)

if 1:
    print
    print "Checking number of affected atoms."
    #Verbose(1)

    np_atoms = Atoms(atoms, pbc=False)
    np_nblist = NeighborList(listcutoff, np_atoms, 0, full=True)

    i = random.random_integers(0, len(np_atoms)-1)
    dr = random.standard_normal(3)
    dr *= 10.0/sqrt(dot(dr,dr))
    atom = np_atoms[i]
    oldpos = array(atom.position)
    newpos = oldpos + dr
    atom.position = newpos
    nr = np_nblist.test_partial_update(array([i,], int32), np_atoms)
    # Count who should be affected.
    r = np_atoms.get_positions() - oldpos
    r2 = sum(r * r, 1)
    a = less_equal(r2, listcutoff*listcutoff)
    r = np_atoms.get_positions() - newpos
    r2 = sum(r * r, 1)
    b = less_equal(r2, listcutoff*listcutoff)

    nr2 = sum(logical_or(a,b))
    print nr, "atoms were affected"
    print sum(a), "near old position and", sum(b), "near new position."
    print nr2, "atoms should be affected"
    ReportTest("Number of affected atoms", nr, nr2, 0)
    
    del np_atoms, np_nblist
    
if 1:
    print
    print "Testing perturbations of all the atoms"

    magnarr = array((0.0, 0.01, 0.1, 1.0, 3.0, 10.0))
    pick = argsort(random.random((len(magnarr),)))
    for magn in take(magnarr, pick):
        dx = magn * random.uniform(-1, 1, (len(atoms), 3))
        print "  Perturbation:", magn
        atoms.set_positions(dx + atoms.get_positions())
        nblist.check_and_update(atoms)
        checklist(nblist, atoms, listcutoff)



print 
print "Testing perturbations of single atoms"

magnarr = array((0.0, 0.01, 0.1, 1.0, 3.0, 10.0))
pick1 = argsort(random.random((len(magnarr),)))
for magn in take(magnarr, pick1):
    numarr = array((1, 3, 10, len(atoms)/2, -3))
    pick2 = argsort(random.random((len(numarr),)))
    for number in take(numarr, pick2):
        # Pick number random atoms.
        if number < 0:
            # Find N neighboring atoms and perturb them
            number = -number
            n = 0
            while len(nblist[n]) < number:
                n += 1
            pick = concatenate([array((n,)), nblist[n][:number-1]])
        else:
            pick = argsort(random.random((len(atoms),)))[:number]
        if number > 15:
            s = "<<< %s atoms >>>" % (number,)
        else:
            s = str(list(pick))
        print "  dr = %.3f: %d atoms: %s" % (magn, number, s)
        for i in pick:
            atom = atoms[i]
            dr = random.standard_normal(3)
            dr *= magn/sqrt(dot(dr,dr))
            atom.position += dr

        nblist.test_partial_update(pick.astype(int32), atoms)
        checklist(nblist, atoms, listcutoff)


ReportTest.Summary()

