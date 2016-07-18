from numpy import *
from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from ase.constraints import Filter
from asap3.testtools import *

atoms = FaceCenteredCubic(size=(5, 5, 5), symbol="Cu")
#p = RasMol(atoms)

ReportTest("Total number of atoms", len(atoms), 500, 0)

ss = Filter(atoms, indices=[2,5,7,10,25,100])
ReportTest("number of atoms in subset",len(ss), 6, 0)

a = atoms[25]
b = ss[4]
assert a.index == b.index
assert a.atoms is b.atoms

a2 = Atoms(ss)
ReportTest("number of atoms in copy of subset",len(a2), 6, 0)

a = Filter(atoms, indices=(5,6,7,8))
ReportTest("Number of atoms in subset", len(a), 4, 0)

a.set_atomic_numbers(47*ones(4, int32))
ReportTest("subset: z[0]", a.get_atomic_numbers()[0], 47, 0)
ReportTest("subset: z[1]", a.get_atomic_numbers()[1], 47, 0)
ReportTest("subset: z[2]", a.get_atomic_numbers()[2], 47, 0)
ReportTest("subset: z[3]", a.get_atomic_numbers()[3], 47, 0)
z = atoms.get_atomic_numbers()
ReportTest("full: z[0]", z[0], 29, 0)
ReportTest("full: z[4]", z[4], 29, 0)
ReportTest("full: z[5]", z[5], 47, 0)
ReportTest("full: z[8]", z[8], 47, 0)
ReportTest("full: z[9]", z[9], 29, 0)
ReportTest("full: z[100]", z[100], 29, 0)

r = a.get_positions()
r[2] = array([10.0, 20.0, 30.0])
a.set_positions(r)

r = atoms.get_positions()
ReportTest("full: r[7,0]", r[7,0], 10.0, 0)
ReportTest("full: r[7,1]", r[7,1], 20.0, 0)
ReportTest("full: r[7,2]", r[7,2], 30.0, 0)

ra = atoms[7].position
rb = a[2].position
ReportTest("x coordinate of number 2/7", ra[0], rb[0], 0)
ReportTest("y coordinate of number 2/7", ra[1], rb[1], 0)
ReportTest("z coordinate of number 2/7", ra[2], rb[2], 0)

#print a.get_positions()

ReportTest.Summary()
