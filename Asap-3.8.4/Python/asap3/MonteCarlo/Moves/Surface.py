import numpy as np

from asap3 import FullNeighborList
from asap3.MonteCarlo.Moves.Data import lattice as latticedata
from asap3.MonteCarlo.Moves.Base import Move, RandChoice
from ase.cluster.cluster import Cluster

class SurfaceMove(Move):
    def __init__(self, debug=0):
        self.debug = debug
        self.surface_i = -1
        self.vacant_i = -1
        Move.__init__(self)

    def set_atoms(self, atoms):
        Move.set_atoms(self, atoms)

        if not isinstance(atoms, Cluster):
            raise Warning('The cluster is not a valid Cluster instance.')

        if atoms.symmetry != 'fcc':
            raise Warning('Can only work with fcc structures.')

        self.natoms = len(atoms)

        # Get the structual parameters
        sym = atoms.symmetry
        self.neighbor_positions = (latticedata[sym]['neighbor_positions'] *
                                   atoms.lattice_basis[0, 0])
        self.neighbor_cutoff = (latticedata[sym]['neighbor_cutoff'] *
                                atoms.lattice_basis[0, 0])
        self.neighbor_numbers = latticedata[sym]['neighbor_numbers']
        self.neighbor_mapping = latticedata[sym]['neighbor_mapping']
        self.neighbor_count = latticedata[sym]['neighbor_count']
        self.type_count = latticedata[sym]['type_count']
        self.type_names = latticedata[sym]['type_names']
        self.type_numbers = latticedata[sym]['type_numbers']
        self.type_data = latticedata[sym]['type_data']

        # Set neighbors, coordination and type for all atoms
        get_neighbors = FullNeighborList(self.neighbor_cutoff, atoms)

        positions = atoms.get_positions()
        neighbors = []
        coordinations = []
        types = []

        for i, pos in enumerate(positions):
            nl = get_neighbors[i]
            dl = (positions[nl] - pos) / atoms.lattice_basis[0, 0]

            neighbors.append([-1] * self.neighbor_count)
            for n, d in zip(nl, dl):
                name = tuple(d.round(1))
                if name in self.neighbor_numbers:
                    neighbors[i][self.neighbor_numbers[name]] = n

            coordinations.append(self.get_atom_coordination(neighbors[i]))
            types.append(self.get_atom_type(neighbors[i]))

        self.positions = positions
        self.neighbors = np.array(neighbors, dtype=int)
        self.coordinations = np.array(coordinations, dtype=int)
        self.types = np.array(types, dtype=int)

        # Generate vacancies (position, neighbors, coordination and type)
        self.vacant_positions = np.zeros((0, 3), dtype=float)
        self.vacant_neighbors = np.zeros((0, self.neighbor_count), dtype=int)
        self.vacant_coordinations = np.zeros(0, dtype=int)
        self.vacant_types = np.zeros(0, dtype=int)

        for i in self.surface_indexes():
            self.add_vacancies(self.positions[i], self.neighbors[i], i)

        if self.debug > 2:
            if self.check_vacancies():
                raise Warning('Something is wrong, read the message above')

    def __call__(self, atoms):
        self.chose()
        self.backup(self.surface_i)
        atoms[self.surface_i].position = self.vacant_positions[self.vacant_i].copy()

    def chose(self):
        # Choose a surface atom to move and where to move it
        surface_i = RandChoice(self.surface_indexes())

        self_neighbor = True
        while self_neighbor:
            vacant_i = RandChoice(self.vacant_indexes())
            # The atom must not move to a site where its only neighbor is it self
            # If just one other neighbor is found then its fine => break
            for n in self.vacant_neighbors[vacant_i]:
                if n != -1 and n != surface_i:
                    self_neighbor = False
                    break

        self.surface_i = surface_i
        self.vacant_i = vacant_i

    def accept(self):
        # Save the old position, neighbors, type and coordination
        old_position = self.positions[self.surface_i].copy()
        old_neighbors = self.neighbors[self.surface_i].copy()
        old_coordination = self.coordinations[self.surface_i].copy()
        old_type = self.types[self.surface_i].copy()

        # Set the position, neighbors, type and coordination
        new_position = self.vacant_positions[self.vacant_i].copy()
        new_neighbors = self.vacant_neighbors[self.vacant_i].copy()
        new_coordination = self.vacant_coordinations[self.vacant_i].copy()
        new_type = self.vacant_types[self.vacant_i].copy()

        self.positions[self.surface_i] = new_position
        self.neighbors[self.surface_i] = new_neighbors
        self.coordinations[self.surface_i] = new_coordination
        self.types[self.surface_i] = new_type

        # Add "surface_i" in the new neighbors neighborlists
        for j, n in enumerate(new_neighbors):
            if n >= 0 and n != self.surface_i:
                self.neighbors[n][self.neighbor_mapping[j]] = self.surface_i
                self.coordinations[n] += 1
                self.types[n] = self.get_atom_type(self.neighbors[n])
            elif n == self.surface_i:
                self.neighbors[n][j] = -1
                self.coordinations[n] -= 1
                self.types[n] = self.get_atom_type(self.neighbors[n])

        # Set "-1" in the old neighbors neighborlists
        for j, n in enumerate(old_neighbors):
            if n >= 0:
                self.neighbors[n][self.neighbor_mapping[j]] = -1
                self.coordinations[n] -= 1
                self.types[n] = self.get_atom_type(self.neighbors[n])

        # Add the old and remove the new position
        self.remove_vacancy(self.vacant_i)
        self.add_vacancy(old_position,
                         old_neighbors,
                         old_coordination,
                         old_type)

        # Update vacancies around the moved atom old and new position
        self.remove_vacancies(self.surface_i)
        self.add_vacancies(self.positions[self.surface_i],
                           self.neighbors[self.surface_i],
                           self.surface_i)

        #Print that an atom is moved
        if self.debug:
            print (">> Moving atom %i (%2.2f %2.2f %2.2f) " %
                   (self.surface_i, old_position[0],
                    old_position[1], old_position[2]) +
                   "to vacancy %i (%2.2f %2.2f %2.2f):" %
                   (self.vacant_i, new_position[0],
                    new_position[1], new_position[2]))

        if self.debug > 2:
            check_vacancies = self.check_vacancies()
            check_neighbors = self.check_neighbors()

            if check_vacancies or check_neighbors:
                raise Warning('Something is wrong, read messages above')

    # Manipulation functions for vacancies
    def add_vacancy(self, position, neighbors, coordination, type):
        """Adds a single vacancy, without updating any lists."""
        if self.debug:
            print ("Adding vacancy %i (%2.2f %2.2f %2.2f)" %
                   (len(self.vacant_positions), position[0], position[1], position[2]))

        self.vacant_positions = np.append(self.vacant_positions,
                                          position.reshape(1, 3),
                                          axis=0)
        self.vacant_neighbors = np.append(self.vacant_neighbors,
                                          neighbors.reshape(1, self.neighbor_count),
                                          axis=0)
        self.vacant_coordinations = np.append(self.vacant_coordinations, 
                                              coordination)
        self.vacant_types = np.append(self.vacant_types, 
                                      type)

    def remove_vacancy(self, index):
        """Removes a single vacancy, without updating any lists."""
        if self.debug:
            position = self.vacant_positions[index]
            print ("Removing vacacy %i (%2.2f %2.2f %2.2f)" %
                   (index, position[0], position[1], position[2]))

        self.vacant_positions = np.delete(self.vacant_positions, index, axis=0)
        self.vacant_neighbors = np.delete(self.vacant_neighbors, index, axis=0)
        self.vacant_coordinations = np.delete(self.vacant_coordinations, index)
        self.vacant_types = np.delete(self.vacant_types, index)

    def add_vacancies(self, position, neighbors, index):
        """Adds and updates vacancies associated with the atom 
        at "position" indexed with "index". It will only add a
        vacancy if it is not pressent in the list of vacancies
        already."""

        for i, n in enumerate(neighbors):
            if n < 0:
                vacant_position = position + self.neighbor_positions[i] 
                vacant_append = True

                for j, vp in enumerate(self.vacant_positions):
                    if (np.abs(vp - vacant_position) < 1e-3).all():
                        self.vacant_neighbors[j][self.neighbor_mapping[i]] = index
                        self.vacant_coordinations[j] += 1
                        self.vacant_types[j] = self.get_atom_type(self.vacant_neighbors[j])
                        vacant_append = False

                if vacant_append:
                    if self.debug:
                        print ("Adding vacancy %i (%2.2f %2.2f %2.2f)" %
                               (len(self.vacant_positions),
                                vacant_position[0], 
                                vacant_position[1], 
                                vacant_position[2]))

                    self.vacant_positions = np.append(self.vacant_positions,
                                                      vacant_position.reshape(1,3),
                                                      axis=0)
                    self.vacant_neighbors = np.append(self.vacant_neighbors,
                                                      -np.ones((1,self.neighbor_count), dtype=int),
                                                      axis=0)
                    self.vacant_neighbors[-1][self.neighbor_mapping[i]] = index
                    self.vacant_coordinations = np.append(self.vacant_coordinations, 1)
                    self.vacant_types = np.append(self.vacant_types, 
                                             self.get_atom_type(self.vacant_neighbors[-1]))

    def remove_vacancies(self, index):
        """Removes and updates vacancies associated with the 
        atom indexed "index". It will only remove vacancies
        if they are not associated with other atoms."""

        i = 0
        end = len(self.vacant_neighbors)
        while i < end:
            vacant_remove = True

            for j, n in enumerate(self.vacant_neighbors[i]):
                if n >= 0 and n != index:
                    vacant_remove = False
                if n == index:
                    self.vacant_neighbors[i][j] = -1
                    self.vacant_coordinations[i] -= 1
                    self.vacant_types[i] = self.get_atom_type(self.vacant_neighbors[i])

            if vacant_remove:
                if self.debug:
                    position = self.vacant_positions[i]
                    print ("Removing vacacy %i (%2.2f %2.2f %2.2f)" %
                           (i, position[0], position[1], position[2]))

                end -= 1
                self.vacant_positions = np.delete(self.vacant_positions, i, axis=0)
                self.vacant_neighbors = np.delete(self.vacant_neighbors, i, axis=0)
                self.vacant_coordinations = np.delete(self.vacant_coordinations, i)
                self.vacant_types = np.delete(self.vacant_types, i)
            else:
                i += 1

    # Some useful functions
    def surface_indexes(self):
        """Returns the indexes for all the surface atoms (coordination < 12)."""
        mask = self.coordinations < 12
        return np.arange(self.natoms)[mask]

    def vacant_indexes(self):
        """Returns the indexes for the vacancies."""
        return np.arange(len(self.vacant_positions))

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

    # Functions to check the neighbor lists and vacancies
    def check_neighbors(self):
        print ">> Checking the neighbor list..."
        positions = self.positions
        neighbors = self.neighbors
        coordinations = self.coordinations
        types = self.types
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
        print ">> Checking the vacancies..."
        vacant_positions = self.vacant_positions
        vacant_neighbors = self.vacant_neighbors
        positions = self.positions
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

