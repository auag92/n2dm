import os
import math
import cPickle as pickle
from numpy import random
import numpy as np
import sys
from ase.optimize import FIRE
from ase.io import write
from ase.io.trajectory import PickleTrajectory
import subprocess
import gzip
#Changelog:
#v1.1 : Added last atoms configuration to write out for visualizing
#v1.2 : Save only counts in totals not atomwise
#v1.3 : Adscalc or EMT is determined in relax.
#v1.4 : Save SMC energy initially
#v1.5 : Add multiplicity
#v1.6 : Only relax if n_step >= n_step/2 ;  Corrected multiplicities loading.

class AtomMonteCarloData:
    version=1.6
    def sumcoordinations(self,coordinations):
		#Takes coordinations listed atomwise and returns coordinations listed 
		#by counts(coordinationnumber)
		maxcoor=15
		out = np.zeros(maxcoor, int)
		for i in range(maxcoor):
			out[i] = np.equal(coordinations,i).sum() #growing array
		return out


    def __init__(self, filename=None, atoms=None,
                 surfmove=None, symbol=None, symmetry=None,
                 latticeconstant=None, size=None, temp=None,
                 mode=None,Esmc = None, multiplicity=1,total_steps=None):
        #First find atomic number if symbol is a string, else user must have             specified Z
        if isinstance(symbol,str):
            self.atomic_number = data.atomic_numbers[symbol]
        else:
            self.atomic_number = symbol
        #Set attributes:
        self.filename = filename
        self.symmetry = symmetry
        self.lattice_constant = latticeconstant
        self.size = size
        self.temp = temp #Simulation Temperature really.
        self.atoms = atoms
        self.surfmove = surfmove #This holds info from class SurfaceMove from             Surface.py in ase such as coordinations, types etc.
        #Maybe:
        if not self.atoms==None:
        	self.natoms = len(atoms)
        else:
        	self.natoms = 0
        self.ground_state_energy = 1e10
        self.multiplicity = multiplicity

        #The lists the object shall store:
        self.lists = {}
        self.lists['energy'] = []
        self.lists['counts'] = []
        self.lists['types'] = [] #This holds the atom site type e.g. hcp 
        self.lists['coordinations']=[]
        self.lists['relaxed_energies'] = []
        self.lists['unrelaxed_energies'] =[]
        #The vacant sites:
        self.lists['vac_coordinations']=[]
        self.lists['vac_types']=[]
        
        #Store smc energy for future analysis of E_filter
        self.Esmc = Esmc
        
        #Store MC step number : MUST be increased in accept move:
        self.mc_step = 0
        self.total_steps = total_steps
        self.id_relax = None #used to determine useable data.
        self.n_accept = 0
    
    def get(self,list2get):
        if self.lists[str(list2get)] in self.lists:
            return lists
        else:
            raise ValueError("No data named "+list2get)
        
        
    def read(self,filename=None):
        if filename is None:
            filename = self.filename  
        
        if os.path.isfile(filename) and filename.endswith(".amc.gz"):
            f = gzip.open(filename,'r')
        elif os.path.isfile(filename) and filename.endswith(".amc"):
            f = open(filename)
        else: #File does not exist else probably
            print >> sys.stderr
            sys.exit(1)
        #load pickleand store all lists:
        d = pickle.load(f)
        assert(d['filetype'] == "AtomMonteCarloData") #File check.
        assert(d['version'] == self.version)
        self.lists = d['lists']
        #store the non list properties:
        self.atomic_number = d['atomic_number']
        self.lattice_constant = d['lattice_constant']
        #self.size = d['size']
        self.temp = d['temp']
        #MaybeStore these:
        self.natoms = d['natoms']
        self.Esmc = d['Esmc']
        self.multiplicity = d['multiplicity']
        self.id_relax = d['id_relax']
        
    def write(self,filename=None):
        if filename is None:
            filename = self.filename
        if not filename.endswith('.gz'):
        	filename = filename + '.gz'    
        if os.path.isfile(filename):
            os.rename(filename,filename + '.bak')
        #Data is stored as dict.
        d = {'filetype' : 'AtomMonteCarloData',
             'version' : self.version, 
             'atomic_number':          self.atomic_number,
             'symmetry': self.symmetry,
             'lattice_constant':                   self.lattice_constant,
             'temp': self.temp,
             'Esmc':self.Esmc,
             'natoms': self.natoms,
             'lists':          self.lists,
             'atomspositions':self.atoms.get_positions(),
             'multiplicity': self.multiplicity,
             'id_relax':self.id_relax
             }
        
        #f = open(filename,'w')
        f = subprocess.Popen("gzip > "+filename, shell=True,
        					 stdin=subprocess.PIPE).stdin
        pickle.dump(d,f, protocol=-1)
        f.close()
        
    def relax(self):
        """Relax the atoms return the relaxed energy.
        This function leaves the atoms unchanged (it changes them
        but restores the positions).
        
        If a single atoms moves "too much" the relaxation is rejected.
        """
        if self.atoms.get_calculator().__class__.__name__ == 'adscalc':
            self.atoms.get_calculator().use_adsorbates(False)
        energy_before = self.atoms.get_potential_energy()
        #print "unrelaxed energy: "+str(energy)
        old_pos = self.atoms.get_positions()
        
        for i in range(3):
            dyn = FIRE(atoms=self.atoms,logfile=None)
            dyn.run(fmax=0.05, steps=2000)
        
        reten=self.atoms.get_potential_energy()
        
        # Sanity check.
        nblist = self.atoms.get_calculator().get_neighbor_list()
        r = self.atoms.get_positions() - old_pos
        dr = np.zeros((len(self.atoms),3))
        for i in range(len(self.atoms)):
            nblist_i = nblist[i]
            dr[i] = r - r[nblist_i].sum(axis=0)/len(nblist[i])
        dr = np.sqrt( (dr*dr).sum(axis=1).max() )
        # dr is now how far the most mobile atom has moved
        # relative to its neighbors.
        # Now, identify how far it is reasonable that it has moved.
        x = old_pos[0] - old_pos[nblist[0]]
        acceptable = np.sqrt( (x*x).sum(axis=1).min() )
        assert acceptable > 0.01
        if dr > 0.5 * acceptable:
            reten = energy_before  # Reject
        # Restore the atoms
        self.atoms.set_positions(old_pos)
        if self.atoms.get_calculator().__class__.__name__ == 'adscalc':
            self.atoms.get_calculator().use_adsorbates(True)
        return energy_before, reten 
    
    def accept_move(self):
        e = self.atoms.get_potential_energy()
        self.lists['energy'].append(e)
        self.lists['counts'].append(1)
        self.lists['coordinations'].append(self.sumcoordinations(self.surfmove.coordinations))
        self.lists['types'].append(self.sumcoordinations(self.surfmove.types)) #See sumcoordinations
        self.lists['vac_coordinations'].append(self.sumcoordinations(self.surfmove.vacant_coordinations))
        self.lists['vac_types'].append(self.sumcoordinations(self.surfmove.vacant_types))
        
        
        if self.mc_step>=self.total_steps/2:
            if self.id_relax == None:
        	    self.id_relax = self.n_accept
            unrelax_energy, relax_energy = self.relax()
            self.lists['unrelaxed_energies'].append(unrelax_energy)
            self.lists['relaxed_energies'].append(relax_energy)
            e += relax_energy - unrelax_energy
        
            if e < self.ground_state_energy:
                #now store .amc in tmpfn[1] and 0000* in tmpfn[0]
                tmpfn = self.filename.split(".")
                traj = PickleTrajectory(tmpfn[0]+".traj", 'w', self.atoms, master=True, backup=False)
                traj.write()
                traj.close()
                #write(filename=tmpfn[0]+".traj",images=self.atoms)
            self.ground_state_energy = e
        else:
            self.lists['unrelaxed_energies'].append(np.nan)
            self.lists['relaxed_energies'].append(np.nan)
        self.mc_step += 1
        self.n_accept += 1
    
    def accept_move_save_all(self):
        e = self.atoms.get_potential_energy()
        self.lists['energy'].append(e)
        self.lists['counts'].append(1)
        self.lists['coordinations'].append(self.sumcoordinations(self.surfmove.coordinations))
        self.lists['types'].append(self.sumcoordinations(self.surfmove.types)) #See sumcoordinations
        self.lists['vac_coordinations'].append(self.sumcoordinations(self.surfmove.vacant_coordinations))
        self.lists['vac_types'].append(self.sumcoordinations(self.surfmove.vacant_types))
        unrelax_energy, relax_energy = self.relax()
        self.lists['unrelaxed_energies'].append(unrelax_energy)
        self.lists['relaxed_energies'].append(relax_energy)
        e += relax_energy - unrelax_energy
        #now store .amc in tmpfn[1] and 0000* in tmpfn[0]
        i=0
        tmpfn = self.filename.split(".")
        outfn = tmpfn[0]+"_step"+str(i)+"_"+".traj"
        while os.path.exists(outfn):
        	i+=1
        	outfn = tmpfn[0]+"_step"+str(i)+"_"+".traj"
        
        write(filename=outfn,images=self.atoms)
        self.ground_state_energy = e
    
        
    def reject_move(self):
    	self.mc_step+=1
    	l = self.lists['counts']
    	if l:
    		l[-1] += 1
    
	
		
