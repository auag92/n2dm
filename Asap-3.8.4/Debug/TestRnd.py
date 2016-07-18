from asap3 import *
from ase import units
from ase.lattice.cubic import FaceCenteredCubic
import numpy as np
import matplotlib.pyplot as plt
import time

repeat = 25

AsapThreads()

#set_verbose(1)

atoms = FaceCenteredCubic(size=(50,50,50), symbol='Cu')
atoms.set_calculator(EMT())
lgv = Langevin(atoms, 5*units.fs, 0.001, 300*units.kB)
lgv.run(1)

print len(atoms)
x1, x2 = lgv.get_random(False)
x1.shape = (-1,)
x2.shape = (-1,)
x = np.concatenate((x1, x2))
h1 = np.histogram(x, 1000, range=(-1.0,2.0))
h2 = np.histogram(np.random.random(len(x)), 1000, range=(-1.0,2.0))

print h1[0].shape, h1[1].shape

plt.figure()
plt.plot(h1[1][:-1], h1[0], 'b-')
plt.plot(h2[1][:-1], h2[0], 'r-')

t1 = time.time()
for i in range(repeat):
    x1, x2 = lgv.get_random(True)
t1 = time.time() - t1
x1.shape = (-1,)
x2.shape = (-1,)
x = np.concatenate((x1, x2))
t2 = time.time()
for i in range(repeat):
    y = np.random.standard_normal(len(x))
t2 = time.time() - t2
h1 = np.histogram(x, 1000, range=(-3.0,3.0))
h2 = np.histogram(y, 1000, range=(-3.0,3.0))

print h1[0].shape, h1[1].shape

plt.figure()
plt.plot(h1[1][:-1], h1[0], 'b-')
plt.plot(h2[1][:-1], h2[0], 'r-')

print "T1 =", t1, "   T2 =", t2

plt.show()
