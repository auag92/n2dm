import matplotlib.pyplot as plt
import cPickle as pickle
import numpy as np
import sys

fn = sys.argv[1]#Filename

f = open(fn,"r")
p = pickle.load(f)

E_mc = p['E_mc']
E_urel = p['E_urel']
E_rel = p['E_rel']

E = E_mc + E_rel - E_urel
I = range(0,len(E_mc),1)
print E_mc[0]
plt.figure(1)

plt.subplot(211)
titstr = fn.split('_')[-1].split('.')[0]
plt.title(titstr)
plt.plot(I,E[I], 'k',label="Ereal")
plt.legend( loc=2, borderaxespad=0.)



plt.subplot(212)

plt.plot(I,E_mc[I],'b',label = "Emc")
plt.legend( loc=2, borderaxespad=0.)


plt.show()
