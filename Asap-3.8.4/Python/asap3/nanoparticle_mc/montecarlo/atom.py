import os
import math
import pickle
from numpy import random
import numpy as np

from datetime import datetime
from ase.io.trajectory import PickleTrajectory
from asap3.nanoparticle_mc.cluster import Cluster
from asap3.nanoparticle_mc.clusteratom import ClusterAtom
from asap3.nanoparticle_mc import data

class AtomMonteCarlo:
    def __init__(self, cluster=None, debug=0, asap=None):
        self.debug = debug
        if asap is not None:
            self.use_asap = True
            self.asap_calc = asap
        else:
            self.use_asap = False

        if not isinstance(cluster, Cluster):
            raise Warning('The cluster is not a valid Cluster instance.')

        if cluster.symmetry != 'fcc':
            raise Warning('Can only work with fcc structures.')

        #Set the structual parameters
        self.cluster = cluster
        self.size = len(self.cluster)

        self.neighbor_positions = (data.lattice[cluster.symmetry]['neighbor_positions'] * 
                                   self.cluster.lattice_constant)
        self.neighbor_mapping = data.lattice[cluster.symmetry]['neighbor_mapping']
        self.neighbor_count = data.lattice[cluster.symmetry]['neighbor_count']
        self.type_count = data.lattice[cluster.symmetry]['type_count']
        self.type_names = data.lattice[cluster.symmetry]['type_names']
        self.type_numbers = data.lattice[cluster.symmetry]['type_numbers']
        self.type_data = data.lattice[cluster.symmetry]['type_data']

        self.arrays = {}

        #Generate a list with vacant positions and their nearest neighbors
        self.cluster.make_neighborlist()

        positions = self.get_positions()
        neighbors = self.get_neighbors()

        for i in self.surface_indexes():
            self.add_vacancies(positions[i], neighbors[i], i)

        if self.debug > 2:
            if self.check_vacancies():
                raise Warning('Something is wrong, read the message above')

    def run(self, steps=1000, temp=1000, txt=None, amc=None, trajectory=None):
        #Open the output textfile
        if txt is not None:
            if not isinstance(txt, str):
                raise Warning('You must specify a valid filename.')
            f = open(txt, 'w')
        else:
            f = None

        #Initialize the dataobject
        if amc is not None:
            amcdata = AtomMonteCarloData(amc,
                                         self.cluster.atomic_number,
                                         self.cluster.symmetry,
                                         self.cluster.lattice_constant,
                                         self.size,
                                         temp)
            amcdata.set_dimensions(self.type_count + 1, self.neighbor_count + 1)
        else:
            amcdata = None

        #Open the trajectoryfile
        if trajectory is not None:
            traj = PickleTrajectory(trajectory, 'w', self.cluster, None, False)
            #traj.write_energy = False
            traj.write_forces = False
            traj.write_stress = False
            traj.write_momenta = False
            traj.write()

        #Print header to outputfile
        if isinstance(f, file):
            print >> f, ('Monte Carlo Sumulation on Cluster Surface Atoms \n' +
                         'Date: %s\n' % datetime.now().strftime("%A, %d %b %Y %H:%M:%S") +
                         'Material: %s\n' % data.chemical_symbols[self.cluster.atomic_number] +
                         'Lattice constant: %s\n' % self.cluster.lattice_constant +
                         'Cluster size: %i\n' % self.size +
                         'Multiplicity: %i\n' % self.cluster.multiplicity +
                         'Temperature: %i K\n' % temp +
                         'Thermal energy: %.5f eV\n' % (data.kB * temp) +
                         'Iteration  Atom          Position                ' +
                         'Desination           Energy    Accept  Probability')

        #Set initial values for the simulation
        old_energy = self.cluster.get_potential_energy()
        acceptances = 0

        if amcdata is not None:
            amcdata.append((old_energy,
                            0.0, #Releaxed energy not found yet
                            self.cluster.get_diameter(),
                            self.cluster.multiplicity,
                            1,
                            np.bincount(self.get_types()),
                            np.bincount(self.get_coordinations()),
                            np.bincount(self.get_vacant_types()),
                            np.bincount(self.get_vacant_coordinations()),
                           ))

        #Print initial energy to outputfile
        if isinstance(f, file):
            print >> f, '%6i %70.4f' % (0, old_energy)

        # Check for ASAP calculator, and enable optimization
        #use_asap = 'asap' in str(self.cluster.get_calculator())
        if self.use_asap:
            auxcluster = self.cluster.copy()
            auxcluster.set_calculator(self.asap_calc)
            auxcluster.get_potential_energy()
            
        #Monte Carlo simulering
        for i in range(1, steps + 1):
            #Chose a move
            surface_i, vacant_i = self.chose_move_old()

            if self.use_asap:
                new_atom = auxcluster[surface_i]
                new_atom.set_position(self.get_vacant_positions()[vacant_i])
                new_energy = auxcluster.get_potential_energy()
                # Move atom back.  It is moved forward again if accepted.
                #new_atom.set_position(self.cluster[surface_i].position)
            else:
                #Make new cluster with change and calculate energy
                new_cluster = self.cluster.copy()
                new_cluster.calc = self.cluster.calc
                new_atom = ClusterAtom(new_cluster.atomic_number)
                new_atom.position = self.get_vacant_positions()[vacant_i]
                new_atom.neighbors = self.get_vacant_neighbors()[vacant_i]
                new_atom.type = self.get_vacant_types()[vacant_i]
                new_atom.coordination = self.get_vacant_coordinations()[vacant_i]
                new_cluster.move_atom(surface_i, new_atom)
                new_energy = new_cluster.get_potential_energy()

            #Calculate the energy change and accept or reject the change
            energy_change = new_energy - old_energy
            accept = False

            if energy_change < 0:
                p = 1.0
                accept = True
            else:
                p = math.exp(- energy_change / (data.kB * temp))
                pc = random.random()
                if pc < p: accept = True

            #Print info on MC move til outputfile
            if isinstance(f, file):
                old_position = self.cluster[surface_i].position
                new_position = new_atom.position

                print >> f, ('%6i %8i   (%6.2f%7.2f%7.2f)   (%6.2f%7.2f%7.2f) %11.4f' %  
                             (i, surface_i, 
                              old_position[0], old_position[1], old_position[2], 
                              new_position[0], new_position[1], new_position[2], 
                              new_energy)),

            if accept:
                old_energy = new_energy
                acceptances += 1
                self.move_atom(surface_i, vacant_i)

                if isinstance(f, file):
                    print >> f, '%7s %10.2f' % ('Yes', p)

                if amcdata is not None:
                    amcdata.append((old_energy,
                                    0.0, #Releaxed energy not found yet
                                    self.cluster.get_diameter(),
                                    self.cluster.multiplicity,
                                    1,
                                    np.bincount(self.get_types()),
                                    np.bincount(self.get_coordinations()),
                                    np.bincount(self.get_vacant_types()),
                                    np.bincount(self.get_vacant_coordinations()),
                                   ))

                if trajectory is not None:
                    traj.write(self.cluster)
            else:
                if self.use_asap:
                    # Move atom back
                    new_atom.set_position(self.cluster[surface_i].position)
                if isinstance(f, file):
                    print >> f, '%7s %10.2f' % ('No', p)

                if amcdata is not None:
                    amcdata.add_count(-1, 1)

        #Print footer to the outputfile
        if isinstance(f, file):
            print >> f, '\nDate: %s' % datetime.now().strftime("%A, %d %b %Y %H:%M:%S")
            print >> f, 'Acceptances ratio: %.2f' % (acceptances * 1.0 / steps)

        #Save the data object
        if amcdata is not None:
            amcdata.write()

        #Close the trajectory
        if trajectory is not None:
            traj.close()

    def chose_move(self):
        #Choose a surface atom to move and where to move it
        #coordinations = np.arange(1, 12)
        coordinations = np.array([1,1,2,2,3,3,4,4,5,5,6,7,8,9,10,11])
        indexes = np.empty(0)
        while len(indexes) == 0:
            coordination = self.choice(coordinations)
            mask = (self.get_coordinations() == coordination)
            indexes = np.arange(len(mask))[mask]
        surface_i = self.choice(indexes)

        #coordinations = np.arange(1, 12)
        coordinations = np.array([1,2,3,4,5,6,7,8,9,10,11])
        indexes = np.empty(0)
        self_neighbor = True
        while True:
            while len(indexes) == 0:
                coordination = self.choice(coordinations)
                mask = (self.get_vacant_coordinations() == coordination)
                indexes = np.arange(len(mask))[mask]
            vacant_i = self.choice(indexes)
            #The atom must not move to a site where its only neighbor is it self
            #If just one other neighbor is found then its fine => break
            neighbors = self.get_vacant_neighbors()[vacant_i]
            if ((neighbors != -1) * (neighbors != surface_i)).any():
                break

        return surface_i, vacant_i

    def chose_move_old(self):
        #Choose a surface atom to move and where to move it
        surface_i = self.choice(self.surface_indexes())

        self_neighbor = True
        while self_neighbor:
            vacant_i = self.choice(self.vacant_indexes())
            #The atom must not move to a site where its only neighbor is it self
            #If just one other neighbor is found then its fine => break
            for n in self.get_vacant_neighbors()[vacant_i]:
                if n != -1 and n != surface_i:
                    self_neighbor = False
                    break

        return surface_i, vacant_i

    def resize(self, size=None, method="random"):
        if method == "random":
            if size > self.size: #Add
                adding_atoms = size - self.size 
                added_atoms = 0
                coordination = 12

                while True:
                    coordination -= 1
                    coordinations = self.get_vacant_coordinations()
                    indexes = np.arange(len(coordinations))[coordinations >= coordination]
                    difference = adding_atoms - added_atoms

                    if self.debug:
                        print ("Found %i vacancies with coordination %i." %
                               (len(indexes), coordination))

                    if len(indexes) > difference:
                        for i in range(difference):
                            n = self.choice(range(len(indexes)))
                            self.add_atom(indexes[n])
                            added_atoms += 1
                            indexes = indexes - (indexes > indexes[n]).astype(int)
                            indexes = np.delete(indexes, n)
                    else:
                        for n in range(len(indexes)):
                            self.add_atom(indexes[n])
                            added_atoms += 1
                            indexes = indexes - (indexes > indexes[n]).astype(int)

                    if added_atoms == adding_atoms:
                        break

            elif size < self.size: #Remove
                removing_atoms = self.size - size
                removed_atoms = 0
                coordination = 0

                while True:
                    coordination += 1
                    coordinations = self.get_coordinations()
                    indexes = np.arange(len(coordinations))[coordinations <= coordination]
                    difference = removing_atoms - removed_atoms

                    if self.debug:
                        print ("Found %i atoms with coordination %i." %
                               (len(indexes), coordination))

                    if len(indexes) > difference:
                        for i in range(difference):
                            n = self.choice(range(len(indexes)))
                            self.remove_atom(indexes[n])
                            removed_atoms += 1
                            indexes = indexes - (indexes > indexes[n]).astype(int)
                            indexes = np.delete(indexes, n)
                    else:
                        for n in range(len(indexes)):
                            self.remove_atom(indexes[n])
                            removed_atoms += 1
                            indexes = indexes - (indexes > indexes[n]).astype(int)

                    if removed_atoms == removing_atoms:
                        break

        elif method == "cluster":
            raise Warning('The method is not supported yet.')

        self.size = len(self.cluster)

        if self.debug > 2:
            check_vacancies = self.check_vacancies()
            check_neighbors = self.check_neighbors()

            if check_vacancies or check_neighbors:
                raise Warning('Something is wrong, read message above')

    #Manipulation functions for atoms
    def add_atom(self, index):
        """ Adds an atom to the vacant position given by "index". The
        vacant position is removed and new vacancies are added around
        the atom. The atoms neighbors neighborlists are updated with
        the atoms new index.
        """

        #Set up the new atom with position, neighbors, type and coordination
        atom = ClusterAtom(self.cluster.atomic_number)
        atom.position = self.get_vacant_positions()[index]
        atom.neighbors = self.get_vacant_neighbors()[index]
        atom.coordination = self.get_vacant_coordinations()[index]
        atom.type = self.get_vacant_types()[index]

        #Add the atom
        self.cluster.add_atom(atom)

        #Print that an atom has been added
        if self.debug:
            position = atom.position
            print (">> Adding atom %i at vacancy %i (%2.2f %2.2f %2.2f):" %
                   (atom.index, index, position[0], position[1], position[2]))

        #Remove old vacacy and add new vacancies
        self.remove_vacancy(index)
        self.add_vacancies(atom.position, atom.neighbors, atom.index)

    def remove_atom(self, index):
        """ Removes the atom given by "index". Its position is added as
        a vacancy and vacancies only associated with the atom is removed.
        Neighbor references are updated, because the atoms gets new 
        indexes.
        """

        #Remove and save the atom (neighborlist indexes are shifted due to the delete)
        atom = self.cluster.remove_atom(index)

        #Print that an atom is removed
        if self.debug:
            position = atom.position
            print (">> Removing atom %i (%2.2f %2.2f %2.2f):" % 
                   (index, position[0], position[1], position[2]))
            
        #Remove old vacancies
        self.remove_vacancies(index)

        #Subtract one from all neighbor references that is greater than "index"
        vacant_neighbors = self.get_vacant_neighbors()
        vacant_neighbors = vacant_neighbors - (vacant_neighbors > index).astype(int)
        self.set_vacant_neighbors(vacant_neighbors)

        #Add new vacancy
        self.add_vacancy(atom.position, atom.neighbors, atom.coordination, atom.type)

    def move_atom(self, i, vacant_i):
        """ Moves an atom given by "i" to the vacant positon given by
        "vacant_i". The old position is added as vacancy and the new 
        position is removed. Atoms and vacancies around the old and 
        new position are updated.
        """

        #Save the old position, neighbors, type and coordination
        old_atom = self.cluster[i]
        old_atom.cut_reference_to_atoms()

        #Set up the new atom with position, neighbors, type and coordination
        new_atom = ClusterAtom(self.cluster.atomic_number)
        new_atom.position = self.get_vacant_positions()[vacant_i]
        new_atom.neighbors = self.get_vacant_neighbors()[vacant_i]
        new_atom.coordination = self.get_vacant_coordinations()[vacant_i]
        new_atom.type = self.get_vacant_types()[vacant_i]

        self.cluster.move_atom(i, new_atom)

        #Print that an atom is moved
        if self.debug:
            old_position = old_atom.position
            new_position = self.cluster[i].position
            print (">> Moving atom %i (%2.2f %2.2f %2.2f) " %
                   (i, old_position[0], old_position[1], old_position[2]) +
                   "to vacancy %i (%2.2f %2.2f %2.2f):" % 
                   (vacant_i, new_position[0], new_position[1], new_position[2]))
            
        #Add the old and remove the new position
        self.remove_vacancy(vacant_i)
        self.add_vacancy(old_atom.position,
                         old_atom.neighbors,
                         old_atom.coordination,
                         old_atom.type)

        #Update vacancies around the moved atom old and new position
        self.remove_vacancies(i)
        self.add_vacancies(self.cluster[i].position, self.cluster[i].neighbors, i)

        if self.debug > 2:
            check_vacancies = self.check_vacancies()
            check_neighbors = self.check_neighbors()

            if check_vacancies or check_neighbors:
                raise Warning('Something is wrong, read message above')

    #Manipulation functions for vacancies
    def add_vacancy(self, position, neighbors, coordination, type):
        """Adds a single vacancy, without updating any lists."""
        vacant_positions = self.get_vacant_positions()
        vacant_neighbors = self.get_vacant_neighbors()
        vacant_coordinations = self.get_vacant_coordinations()
        vacant_types = self.get_vacant_types()

        if self.debug:
            print ("Adding vacacy %i (%2.2f %2.2f %2.2f)" %
                   (len(vacant_positions), position[0], position[1], position[2]))

        vacant_positions = np.append(vacant_positions,
                                     position.reshape(1, 3),
                                     axis=0)
        vacant_neighbors = np.append(vacant_neighbors,
                                     neighbors.reshape(1, self.neighbor_count),
                                     axis=0)
        vacant_coordinations = np.append(vacant_coordinations, 
                                         coordination)
        vacant_types = np.append(vacant_types, 
                                 type)

        self.set_vacant_positions(vacant_positions)
        self.set_vacant_neighbors(vacant_neighbors)
        self.set_vacant_coordinations(vacant_coordinations)
        self.set_vacant_types(vacant_types)

    def remove_vacancy(self, index):
        """Removes a single vacancy, without updating any lists."""
        vacant_positions = self.get_vacant_positions()
        vacant_neighbors = self.get_vacant_neighbors()
        vacant_coordinations = self.get_vacant_coordinations()
        vacant_types = self.get_vacant_types()

        if self.debug:
            position = vacant_positions[index]
            print ("Removing vacacy %i (%2.2f %2.2f %2.2f)" %
                   (index, position[0], position[1], position[2]))

        vacant_positions = np.delete(vacant_positions, index, axis=0)
        vacant_neighbors = np.delete(vacant_neighbors, index, axis=0)
        vacant_coordinations = np.delete(vacant_coordinations, index)
        vacant_types = np.delete(vacant_types, index)

        self.set_vacant_positions(vacant_positions)
        self.set_vacant_neighbors(vacant_neighbors)
        self.set_vacant_coordinations(vacant_coordinations)
        self.set_vacant_types(vacant_types)

    def add_vacancies(self, position, neighbors, index):
        """Adds and updates vacancies associated with the atom 
        at "position" indexed with "index". It will only add a
        vacancy if it is not pressent in the list of vacancies
        already."""

        vacant_positions = self.get_vacant_positions()
        vacant_neighbors = self.get_vacant_neighbors()
        vacant_coordinations = self.get_vacant_coordinations()
        vacant_types = self.get_vacant_types()

        for i, n in enumerate(neighbors):
            if n < 0:
                vacant_position = position + self.neighbor_positions[i] 
                vacant_append = True

                for j, vp in enumerate(vacant_positions):
                    if (np.abs(vp - vacant_position) < 1e-3).all():
                        vacant_neighbors[j][self.neighbor_mapping[i]] = index
                        vacant_coordinations[j] += 1
                        vacant_types[j] = self.get_atom_type(vacant_neighbors[j])
                        vacant_append = False

                if vacant_append:
                    if self.debug:
                        print ("Adding vacacy %i (%2.2f %2.2f %2.2f)" %
                               (len(vacant_positions),
                                vacant_position[0], 
                                vacant_position[1], 
                                vacant_position[2]))

                    vacant_positions = np.append(vacant_positions,
                                                 vacant_position.reshape(1,3),
                                                 axis=0)
                    vacant_neighbors = np.append(vacant_neighbors,
                                                 -np.ones((1,self.neighbor_count), dtype=int),
                                                 axis=0)
                    vacant_neighbors[-1][self.neighbor_mapping[i]] = index
                    vacant_coordinations = np.append(vacant_coordinations, 1)
                    vacant_types = np.append(vacant_types, 
                                             self.get_atom_type(vacant_neighbors[-1]))

        self.set_vacant_positions(vacant_positions)
        self.set_vacant_neighbors(vacant_neighbors)
        self.set_vacant_coordinations(vacant_coordinations)
        self.set_vacant_types(vacant_types)

    def remove_vacancies(self, index):
        """Removes and updates vacancies associated with the 
        atom indexed "index". It will only remove vacancies
        if they are not associated with other atoms."""

        vacant_positions = self.get_vacant_positions()
        vacant_neighbors = self.get_vacant_neighbors()
        vacant_coordinations = self.get_vacant_coordinations()
        vacant_types = self.get_vacant_types()

        i = 0
        end = len(vacant_neighbors)
        while i < end:
            vacant_remove = True

            for j, n in enumerate(vacant_neighbors[i]):
                if n >= 0 and n != index:
                    vacant_remove = False
                if n == index:
                    vacant_neighbors[i][j] = -1
                    vacant_coordinations[i] -= 1
                    vacant_types[i] = self.get_atom_type(vacant_neighbors[i])

            if vacant_remove:
                if self.debug:
                    position = vacant_positions[i]
                    print ("Removing vacacy %i (%2.2f %2.2f %2.2f)" %
                           (i, position[0], position[1], position[2]))

                end -= 1
                vacant_positions = np.delete(vacant_positions, i, axis=0)
                vacant_neighbors = np.delete(vacant_neighbors, i, axis=0)
                vacant_coordinations = np.delete(vacant_coordinations, i)
                vacant_types = np.delete(vacant_types, i)
            else:
                i += 1

        self.set_vacant_positions(vacant_positions)
        self.set_vacant_neighbors(vacant_neighbors)
        self.set_vacant_coordinations(vacant_coordinations)
        self.set_vacant_types(vacant_types)

    #Functions to check for faults in the neighbor lists
    def check_neighbors(self):
        positions = self.get_positions()
        neighbors = self.get_neighbors()
        coordinations = self.get_coordinations()
        types = self.get_types()
        trouble = False

        for i, n in enumerate(neighbors):
            for j, m in enumerate(n):
                found = False
                p1 = positions[i] + self.neighbor_positions[j]

                for k, p2 in enumerate(positions):
                    if (np.abs(p1 - p2) < 1e-3).all():
                        found = True
                        if m < 0:
                            print 'Atom %i is missing neighbor atom %i at %s' % (i, k, p2)
                            trouble = True
                        elif m != k:
                            print ('Atom %i is having atom %i as ' +
                                   'neighbor instead of atom %i' % (i, m, k))
                            trouble = True

                if m >= 0 and found is False:
                    trouble = True
                    print 'Atom %i does not have atom %i as neighbor' % (i, m)

            if self.get_atom_coordination(n) != coordinations[i]:
                trouble = True
                print 'Atom %i coordination does not match' % i

            if self.get_atom_type(n) != types[i]:
                trouble = True
                print 'Atom %i type does not match' % i

        return trouble

    def check_vacancies(self):
        vacant_positions = self.get_vacant_positions()
        vacant_neighbors = self.get_vacant_neighbors()
        positions = self.get_positions()
        trouble = False
        
        for i, p1 in enumerate(vacant_positions):
            for j, p2 in enumerate(positions):
                if (np.abs(p1 - p2) < 1e-3).all():
                    trouble = True
                    print "Vacancy %i is identical to %i at %s" % (i, j, p1)

        for i, n in enumerate(vacant_neighbors):
            for j, m in enumerate(n):
                p1 = vacant_positions[i] + self.neighbor_positions[j]

                for k, p2 in enumerate(positions):
                    if (np.abs(p1 - p2) < 1e-3).all():
                        if m < 0:
                            print 'Vacancy %i is missing neighbor atom %i at %s' % (i, k, p2)
                            trouble = True
                        elif m != k:
                            print ('Vacancy %i is having atom %i as ' % (i, m) +
                                   'neighbor instead of atom %i' % k)
                            trouble = True

        return trouble

    #Some useful functions
    def choice(self, a):
        return a[random.randint(len(a))]

    def surface_indexes(self):
        """Returns the indexes for all the surface atoms (coordination < 12)."""
        mask = self.get_coordinations() < 12
        return np.arange(len(self.cluster))[mask]

    def vacant_indexes(self):
        """Returns the indexes for the vacancies."""
        return np.arange(len(self.get_vacant_positions()))

    def get_atom_coordination(self, neighbors):
        neighbors = np.array(neighbors)
        neighbors = neighbors[neighbors >= 0]
        return len(neighbors)

    def get_atom_type(self, neighbors):
        neighbors = np.array(neighbors)
        typename = tuple((neighbors >= 0).astype(int))
        if typename in self.type_names:
            return self.type_data[self.type_numbers[typename]]['type']
        else:
            return 0

    #Positions (edited on self.cluster)
    def get_positions(self, index=None):
        return self.cluster.get_positions()

    def set_positions(self, positions):
        return self.cluster.set_positions(positions)

    #Neighbors
    def get_neighbors(self):
        return self.cluster.get_neighbors()

    def set_neighbors(self, neighbors):
        self.cluster.set_neighbors(neighbors)

    #Coordination
    def get_coordinations(self):
        return self.cluster.get_coordinations()

    def set_coordinations(self, coordinations):
        self.cluster.set_coordinations(coordinations)

    #Types
    def get_types(self):
        return self.cluster.get_types()

    def set_types(self, types):
        self.cluster.set_types(types)

    #Vacant positions
    def get_vacant_positions(self):
        if 'vacant_positions' not in self.arrays:
            return np.empty((0, 3), dtype=float)
        else:
            return self.arrays['vacant_positions'].copy()

    def set_vacant_positions(self, positions):
        self.arrays['vacant_positions'] = np.array(positions, dtype=float)

    #Vacant neighbors
    def get_vacant_neighbors(self):
        if 'vacant_neighbors' not in self.arrays:
            return np.empty((0, self.neighbor_count), dtype=int)
        else:
            return self.arrays['vacant_neighbors'].copy()

    def set_vacant_neighbors(self, neighbors):
        self.arrays['vacant_neighbors'] = np.array(neighbors, dtype=int)

    #Vacant coordinations
    def get_vacant_coordinations(self):
        if 'vacant_coordinations' not in self.arrays:
            return np.empty(0, dtype=int)
        else:
            return self.arrays['vacant_coordinations'].copy()

    def set_vacant_coordinations(self, coordinations):
        self.arrays['vacant_coordinations'] = np.array(coordinations, dtype=int)

    #Vacant types
    def get_vacant_types(self):
        if 'vacant_types' not in self.arrays:
            return np.empty(0, dtype=int)
        else:
            return self.arrays['vacant_types'].copy()

    def set_vacant_types(self, types):
        self.arrays['vacant_types'] = np.array(types, dtype=int)

class AtomMonteCarloData:
    def __init__(self, filename=None, symbol=None, symmetry=None, latticeconstant=None,
                 size=None, temp=None):
        if isinstance(symbol, str):
            self.atomic_number = data.atomic_numbers[symbol]
        else:
            self.atomic_number = symbol

        self.filename = filename
        self.symmetry = symmetry
        self.lattice_constant = latticeconstant
        self.size = size
        self.temp = temp

        self.arrays = {}
        self.arrays['energy'] = np.empty(0, dtype=float)
        self.arrays['energyrelax'] = np.empty(0, dtype=float)
        self.arrays['size'] = np.empty(0, dtype=float)
        self.arrays['multi'] = np.empty(0, dtype=int)
        self.arrays['counts'] = np.empty(0, dtype=int)
        self.arrays['types'] = np.empty(0, dtype=int)
        self.arrays['coordinations'] = np.empty(0, dtype=int)
        self.arrays['vtypes'] = np.empty(0, dtype=int)
        self.arrays['vcoordinations'] = np.empty(0, dtype=int)

    def set_dimensions(self, type_count, coord_count):
        self.arrays['types'] = self.arrays['types'].reshape((0, type_count))
        self.arrays['coordinations'] = self.arrays['coordinations'].reshape((0, coord_count))
        self.arrays['vtypes'] = self.arrays['vtypes'].reshape((0, type_count))
        self.arrays['vcoordinations'] = self.arrays['vcoordinations'].reshape((0, coord_count))

    def __len__(self):
        return len(self.arrays['energy'])

    def __getitem__(self, index=-1):
        energy = self.arrays['energy'][index]
        energyrelax = self.arrays['energyrelax'][index]
        size = self.arrays['size'][index]
        counts = self.arrays['counts'][index]
        multi = self.arrays['multi'][index]
        types = self.arrays['types'][index]
        coordinations = self.arrays['coordinations'][index]
        vtypes = self.arrays['vtypes'][index]
        vcoordinations = self.arrays['vcoordinations'][index]
        
        return (energy, energyrelax, size, multi, counts,
                types, coordinations, vtypes, vcoordinations)

    def __delitem__(self, index=-1):
        mask = np.ones(len(self), bool)
        mask[index] = False

        for name, a in self.arrays.items():
            self.arrays[name] = a[mask]

    def __setitem__(self, index, values):
        try:
            if values[0] is not None:
                self.arrays['energy'][index] = values[0]

            if values[1] is not None:
                self.arrays['energyrelax'][index] = values[1]

            if values[2] is not None:
                self.arrays['size'][index] = values[2]

            if values[3] is not None:
                self.arrays['multi'][index] = values[3]

            if values[4] is not None:
                self.arrays['counts'][index] = values[4]

            if values[5] is not None:
                a = self.arrays['types']
                b = np.zeros(a.shape[1], a.dtype)
                b[:len(values[5])] = np.array(values[5], a.dtype)
                self.arrays['types'][index] = b

            if values[6] is not None:
                a = self.arrays['coordinations']
                b = np.zeros(a.shape[1], a.dtype)
                b[:len(values[6])] = np.array(values[6], a.dtype)
                self.arrays['coordinations'][index] = b

            if values[7] is not None:
                a = self.arrays['vtypes']
                b = np.zeros(a.shape[1], a.dtype)
                b[:len(values[7])] = np.array(values[7], a.dtype)
                self.arrays['vtypes'][index] = b

            if values[8] is not None:
                a = self.arrays['vcoordinations']
                b = np.zeros(a.shape[1], a.dtype)
                b[:len(values[8])] = np.array(values[8], a.dtype)
                self.arrays['vcoordinations'][index] = b
        except:
            print values
            raise

    def append(self, values):
        for name, a in self.arrays.items():
            b = np.zeros((a.shape[0] + 1,) + a.shape[1:], a.dtype)
            b[:a.shape[0]] = a
            self.arrays[name] = b

        self[-1] = values

    def extend(self, other):
        if not isinstance(other, AtomMonteCarloData):
            raise Warning('The added object is not an AtomMonteCarloData instance!')

        if len(self) == 0:
            self.set_dimensions(other.get_types().shape[1],
                                other.get_coordinations().shape[1])

        for name, a in self.arrays.items():
            b = other.arrays[name]
            if a.shape[1:] != b.shape[1:]:
                raise Warning('Array %s has different shapes: %s != %s' % (name,
                                                                           a.shape[1:],
                                                                           b.shape[1:]))
            c = np.zeros((a.shape[0] + b.shape[0],) + a.shape[1:], a.dtype)
            c[:a.shape[0]] = a
            c[a.shape[0]:] = b
            self.arrays[name] = c

    def add_count(self, index, value):
        self.arrays['counts'][index] += value

    def filter_mask(self, mask):
        """Function that filters the data with the given mask."""

        if len(mask) != len(self):
            raise Warning('The mask has a wrong size: %s != %s' % (len(mask), len(self)))

        for name, a in self.arrays.items():
            self.arrays[name] = a[mask]

    def filter_energy(self, excess):
        """Function that filters the data leaving only entries with
        energy lower than 'excess'."""

        e = self.get_energies()
        de = (e - e.min())
        mask = de < excess
        self.filter_mask(mask)

    def filter_relaxed(self):
        """Function that filters the data leaving only relaxed entries."""

        mask = (np.abs(self.get_relaxed_energies()) > 1e-3)
        self.filter_mask(mask)

    def get_energies(self):
        return self.arrays['energy'].copy()

    def get_relaxed_energies(self):
        return self.arrays['energyrelax'].copy()

    def get_sizes(self):
        return self.arrays['size'].copy()

    def get_multiplicities(self):
        return self.arrays['multi'].copy()

    def get_counts(self):
        return self.arrays['counts'].copy()

    def get_counts_mask(self):
        mask = np.empty(0, int)
        for i, c in enumerate(self.get_counts()):
            mask_add = np.ones(c, int) * i
            mask_old = mask.copy()
            mask = np.zeros(len(mask_add) + len(mask_old), int)
            mask[:len(mask_old)] = mask_old
            mask[len(mask_old):] = mask_add
        return mask

    def get_types(self):
        return self.arrays['types'].copy()

    def get_coordinations(self):
        return self.arrays['coordinations'].copy()

    def get_vacant_types(self):
        return self.arrays['vtypes'].copy()

    def get_vacant_coordinations(self):
        return self.arrays['vcoordinations'].copy()

    #Functions to write and read the data on the disk
    def write(self, filename=None, zipped=True, backup=True):
        if filename is None:
            filename = self.filename

        if not isinstance(filename, str):
            raise Warning('You must specify a valid filename.')

        if backup:
            if os.path.isfile(filename):
                os.rename(filename, filename + '.bak')
            elif os.path.isfile(filename + ".gz"):
                os.rename(filename + ".gz", filename + '.gz.bak')

        d = {'atomicnumber': self.atomic_number,
             'symmetry': self.symmetry,
             'latticeconstant': self.lattice_constant,
             'size': self.size,
             'temp': self.temp}

        if zipped:
            f = os.popen('gzip > ' + filename + '.gz', 'w')
        else:
            f = open(filename, 'w')

        f.write('AtomMonteCarlo')
        pickle.dump(d, f)
        pickle.dump(self.arrays, f)
        f.close()

    def read(self, filename=None):
        if filename is None:
            filename = self.filename

        if not isinstance(filename, str):
            raise Warning('You must specify a valid filename.')

        if os.path.isfile(filename):
            if filename[-3:] == '.gz':
                filename = filename[:-3]
                f = os.popen("zcat " + filename + ".gz")
            else:
                f = open(filename, 'r')
        elif os.path.isfile(filename + '.gz'):
            f = os.popen("zcat " + filename + '.gz')
        else:
            raise Warning('The file %s do not exist.' % (filename,))

        self.filename = filename

        try:
            if f.read(len('AtomMonteCarlo')) != 'AtomMonteCarlo':
                raise Warning('This is not a compatible file.')
            d = pickle.load(f)
            self.arrays = pickle.load(f)
        except EOFError:
            raise Warinig('Bad file.')

        f.close()

        if 'multiplicity' in d:
            self.multiplicity = d['multiplicity']
            self.arrays['multi'] = np.ones(len(self), float) * d['multiplicity']

        if 'energyrelax' not in self.arrays:
            self.arrays['energyrelax'] = np.zeros(len(self))

        if 'size' not in self.arrays:
            self.arrays['size'] = np.zeros(len(self))

        self.atomc_number = d['atomicnumber']
        self.lattice_constant = d['latticeconstant']
        self.size = d['size']
        self.temp = d['temp']

