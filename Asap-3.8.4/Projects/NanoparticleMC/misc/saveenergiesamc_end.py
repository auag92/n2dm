#PBS -l nodes=1:ppn=1:xeon8
#PBS -q verylong
#PBS -N E_plot_data_gen
#PBS -m ae
import sys
import os
from atommontecarlodata import AtomMonteCarloData
import numpy as np
from ase import units #used to get kB in eV/K
import cPickle as pickle

fdir = sys.argv[1] #Directory of .amc files to join
#e_filter = float(sys.argv[2])

os.chdir(fdir)
files = []

for f in os.listdir("./"):
	if f.endswith(".amc.gz"):
		files.append(f)
#Now we know the file names. Read merge the informations.
dtemp = AtomMonteCarloData()
#Init arrays
#relaxed_energies = np.array([])
#unrelaxed_energies = np.array([])

for f in files:
	dtemp.read(f)
	
	unrelaxed_energies = np.array(dtemp.lists['unrelaxed_energies'])
	relaxed_energies = np.array(dtemp.lists['relaxed_energies'])
	energy = np.array(dtemp.lists['energy'])
	print energy

	fout = open("Energies_for_plot_"+str(f)+".pickle","w")
	pickle.dump({'E_mc':energy,'E_urel':unrelaxed_energies,'E_rel':relaxed_energies},fout)
	fout.close()
