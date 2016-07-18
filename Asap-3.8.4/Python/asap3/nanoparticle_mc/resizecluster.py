from asap3 import CoordinationNumbers
import numpy as np
from asap3.MonteCarlo.Moves import SurfaceMove
from ase import Atom
import random as rd

def mychoice(a, n):
	"Replaces numpy.random.choice(a, n, False) as our numpy is ancient."
	a = list(a)
	while len(a) > n:
		# Remove a random element
		del a[np.random.randint(len(a))]
	return np.array(a)

def resizecluster(atoms, nwanted):
	nactual = len(atoms)
	if nactual == nwanted:
		return
	elif nactual > nwanted:
		removeatoms(atoms, nactual - nwanted)		
	else:
		addatoms(atoms, nwanted - nactual)
		
def removeatoms(atoms, n):
	"Remove n atoms from the cluster."
	while n > 0:
		# We still have atoms to remove
		# Find the ones with lowest coordination number
		coords = CoordinationNumbers(atoms)
		coordination = coords.min()
		idx = np.arange(len(atoms))
		candidates = idx[np.less_equal(coords, coordination)]
		# candidates now contains the indices of the atoms with
		# low coordination number
		if len(candidates) > n:
			# We have too many candidates, must choose.
			candidates = mychoice(candidates, n)
		del atoms[candidates]
		n -= len(candidates)

def addatoms(atoms,n):
	element = atoms[0].symbol
	SM = SurfaceMove()
	SM.set_atoms(atoms)
	idx = SM.vacant_indexes()
	#Find the n highest coordinated sites(this is in ids of SM)
	coords = SM.vacant_coordinations[idx]
	candidates = idx[np.greater_equal(coords,coords.max())]
	vacPos = SM.vacant_positions
	if len(candidates)>=n: #If have sufficient no. candidates
		chosenIDS = rd.sample(range(0,len(candidates)),n) #Random ids
		id1 = candidates[chosenIDS] #the random ids
	else:
		id1 = candidates
	for i in id1:
		#print "Adding atom at site", id1
		atoms.append(Atom(element, position=vacPos[i]))
	if len(id1) < n:
		# We did not add enough atoms
		addatoms(atoms, n - len(id1))
		
		
def objarraysize(arr): #Outputs the size of array objects of arrays
	out = 0
	for i in range(0,arr.size):
		out+=arr[i].size
	return out
