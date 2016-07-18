import os
import math
import pickle
import numpy as np
import time

try:
    import matplotlib.pyplot as plt
    import pylab as pl
    _graphics = True
except:
    _graphics = False

from datetime import datetime
from ase.visualize.primiplotter import PrimiPlotter, GifFile
from asap3.nanoparticle_mc.cluster import Cluster
from asap3.nanoparticle_mc import data

class SurfaceMonteCarlo:
    def __init__(self, symbol=None, layers=None, size=None,
                 latticeconstant=None, symmetry=None, debug=0):

        self.debug = debug
        self.set_calculator()

        #Find the atomic number
        if symbol is not None:
            if isinstance(symbol, str):
                self.atomic_number = data.atomic_numbers[symbol]
            else:
                self.atomic_number = symbol
        else:
            raise Warning('You must specify a atomic symbol or number!')

        #Find the crystal structure
        if symmetry is not None:
            if symmetry.lower() in ['bcc', 'fcc', 'hcp']:
                self.symmetry = symmetry.lower()
            else:
                raise Warning('The %s symmetry does not exist!' % symmetry.lower())
        else:
            self.symmetry = data.reference_states[self.atomic_number]['symmetry'].lower()

        if self.debug:
            print 'Crystal structure:', self.symmetry

        #Find the lattice constant
        if latticeconstant is None:
            if self.symmetry == 'fcc':
                self.lattice_constant = data.reference_states[self.atomic_number]['a']
            else:
                raise Warning(('Cannot find the lattice constant ' +
                               'for a %s structure!' % self.symmetry))
        else:
            self.lattice_constant = latticeconstant

        if self.debug:
            print 'Lattice constant(s):', self.lattice_constant

        #Set the initial layers
        if layers is None:
            raise Warning('You must specify a set of surfaces.')
        else:
            self.init_layers = list(layers)
            
        #Find the size that the clusters will be fitted to
        if size is None:
            self.fitsize = self.get_size_from_layers(layers)
        else:
            self.fitsize = size

    def run(self, steps=1000, temp=1000, txt=None, smc=None):
        #Open the output textfile
        if txt is not None:
            if not isinstance(txt, str):
                raise Warning('You must specify a valid filename.')
            f = open(txt, 'w', 1)

        #Initialize the dataobject
        if smc is not None:
            smc = SurfaceMonteCarloData(smc,
                                        self.atomic_number,
                                        self.symmetry,
                                        self.lattice_constant,
                                        self.fitsize)
            smc.set_dimensions(data.lattice[self.symmetry]['surface_count'])
            self.data = smc

        #Initialize the MC simulation
        oldset = self.fit_to_size(self.init_layers)
        oldsize = self.get_size_from_layers(oldset)
        oldenergy = self.get_energy_from_layers(oldset, True)
        diameter = self.get_diameter_from_layers(oldset)

        if smc is not None:
            smc.append((oldenergy, oldset, oldsize, diameter, None))
        
        #Print header to outputfile
        if isinstance(f, file):
            print >> f, 'Monte Carlo Simulation on Cluster Surface Parameters \n'
            print >> f, 'Date: %s' % datetime.now().strftime("%A, %d %b %Y %H:%M:%S")
            print >> f, 'Material: %s' % data.chemical_symbols[self.atomic_number]
            print >> f, 'Lattice constant: %s' % self.lattice_constant
            print >> f, 'Fitting size: %s' % self.fitsize
            print >> f, 'Temperature: %i K' % temp
            print >> f, 'Thermal energy: %.3f eV\n' % (data.kB * temp)
            print >> f, ('Iteration  Size   Energy   Energy/Atom  ' +
                         'Surface no.  Change  Accept  Probability')
            print >> f, ('%7i %7i %9.3f %10.4f' %  
                         (0, oldsize, oldenergy * oldsize, oldenergy))

        #Start simulation
        acceptances = 0
        accepted = False
        for i in range(1, steps + 1):
            #Choose change to make and do it
            id = self.choice(range(len(oldset)))
            change = self.choice([-1,1])

            newset = oldset[:]
            newset[id] += change
            newsize = self.get_size_from_layers(newset)

            #If layers = 0 there cannot be removed one and the change has to be real
            while (newset[id] < 1 and change < 0) or newsize == oldsize:
                id = self.choice(range(len(oldset)))
                #change = self.choice([-1,1])

                newset = oldset[:]
                newset[id] += change
                newsize = self.get_size_from_layers(newset)

            #Remove excess layers there have been removed
            c = self.get_cluster_from_layers(newset)
            newset = list(c.get_layers())
            del c

            #Fit the size and calculate energy and size
            newset = self.fit_to_size(newset, id)
            newsize = self.get_size_from_layers(newset)
            newenergy = self.get_energy_from_layers(newset, True)

            #Print MC step info to outputfile
            if isinstance(f, file):
                print >> f, ('%7i %7i %9.3f %10.4f %9i %10i' %  
                             (i, newsize, newenergy * newsize, newenergy, id, change)),

            #Accept or reject the change
            de = (newenergy - oldenergy) * (newsize + oldsize) / 2

            if de < 0:
                p = 1.0
                accepted = True
            else:
                p = math.exp(- de / (data.kB * temp))
                pc = self.random()
                if p > pc: accepted = True

            #If accepted do...
            if accepted:
                accepted = False
                acceptances += 1

                oldset = newset[:]
                oldsize = newsize
                oldenergy = newenergy
                diameter = self.get_diameter_from_layers(oldset)

                if smc is not None:
                    smc.append((oldenergy, oldset, oldsize, diameter, None))

                if isinstance(f, file):
                    print >> f, '%8s %10.2f ' % ('Yes', p) + time.ctime()
            else:
                if isinstance(f, file):
                    print >> f, '%8s %10.2f ' % ('No', p) + time.ctime()

        if isinstance(f, file):
            print >> f, '\nDate: %s' % datetime.now().strftime("%A, %d %b %Y %H:%M:%S")
            print >> f, 'Acceptance ratio: %.2f' % (acceptances * 1.0 / steps)

        if smc is not None:
            smc.write()

    def fit_to_size(self, surfaces, restriction=None):
        oldset = list(surfaces)
        oldsize = self.get_size_from_layers(oldset)

        if self.debug > 1:
            print 42 * '-'
            print 'Fitting surfaces to match the size %i.' % self.fitsize
            print 42 * '-'
            print 'Start: ' + str(oldset) + '\n'
            print 'Iter  Size  Surface no.  Change'
            print '  0 %6i  Time is now: ' % oldsize + time.ctime()

        i = 0
        while True:
            i += 1

            #Chose a change to make in the right direction
            id = self.choice(range(len(oldset)))

            if oldsize < self.fitsize:
                change = 1
            elif oldsize > self.fitsize:
                change = -1
            else:
                newsize = None
                break

            newset = oldset[:]
            newset[id] += change
            newsize = self.get_size_from_layers(newset)

            #If layers = 0 there cannot be removed one, the change has to be real and
            #the restricted parameter is keeped fixed
            while (newset[id] < 1 and change < 0) or newsize == oldsize or id == restriction:
                id = self.choice(range(len(newset)))

                newset = oldset[:]
                newset[id] += change
                newsize = self.get_size_from_layers(newset)

            #Remove excess layers there have been removed
            c = self.get_cluster_from_layers(newset)
            newset = list(c.get_layers())
            del c

            if self.debug > 1:
                print '%3i %6i %7i %10i Time: ' % (i, newsize, id, change) + time.ctime()

            #If the wanted size is passed or reached then break
            if oldsize < self.fitsize and newsize > self.fitsize:
                break
            elif oldsize > self.fitsize and newsize < self.fitsize:
                break
            elif newsize == self.fitsize:
                break

            #Save the new paremeter set and start again
            oldset = newset[:]
            oldsize = newsize

        #Chosing between the last two configurations the one with the lowest energy
        if newsize is not None:
            oldenergy = self.get_energy_from_layers(oldset, True)
            newenergy = self.get_energy_from_layers(newset, True)

            if self.debug > 1:
                print '\nIter   Energy   Energy/Atom'
                print '%3i %9.2f %10.3f' % (i - 1, oldenergy * oldsize, oldenergy)
                print '%3i %9.2f %10.3f' % (i, newenergy * newsize, newenergy)
                print 'Time is now: ' + time.ctime()

            if newenergy < oldenergy:
                chosenset = newset
            else:
                chosenset = oldset
        else:
            chosenset = oldset

        if self.symmetry == 'fcc':
            chosenset = list(data.fcc.surface_centering(chosenset))

        #Check for different flaws in the chosen set of layers
        if self.debug > 1:
            print '\nFinal: ' + str(chosenset) + '\n'
            print 'Time is now: ' + time.ctime()

            if newsize is not None and newenergy < oldenergy:
                oldcluster = self.get_cluster_from_layers(newset)
            else:
                oldcluster = self.get_cluster_from_layers(oldset)

            if (np.array(chosenset) < 0).any():
                raise Warning('Negative number of layers at a surface!')

            newcluster = self.get_cluster_from_layers(chosenset)

            if len(oldcluster) != len(newcluster):
                raise Warning(('Size does not match: %i != %i' % 
                               (len(oldcluster), len(newcluster))))

            if (np.array(chosenset) != newcluster.get_layers()).any():
                raise Warning('None existing layers in chosen set')

            if self.debug > 2:
                oldcluster.make_neighborlist()
                newcluster.make_neighborlist()

                oldtypes = np.bincount(oldcluster.get_types())
                newtypes = np.bincount(newcluster.get_types())

                if (oldtypes != newtypes).any():
                    raise Warning('Types does not match: %s != %s' % (oldtypes, newtypes))

        return chosenset

    def set_calculator(self, calc=None):
        self.calc = calc

    def get_calculator(self):
        return self.calc

    #Some useful functions
    def choice(self, data):
        return data[np.random.randint(len(data))]

    def random(self):
        return np.random.random()

    def get_cluster_from_layers(self, layers):
        return Cluster(symbol=self.atomic_number,
                       layers=layers,
                       symmetry=self.symmetry,
                       latticeconstant=self.lattice_constant)

    def get_size_from_layers(self, layers):
        c = self.get_cluster_from_layers(layers)
        return len(c)
        
    def get_energy_from_layers(self, layers, atoms=False):
        if self.get_calculator() is None:
            raise Warning('You must specify a calculator to calculate the energy.')

        c = self.get_cluster_from_layers(layers)
        c.set_calculator(self.get_calculator())
        
        if atoms:
            return c.get_potential_energy() / len(c)
        else:
            return c.get_potential_energy()

    def get_diameter_from_layers(self, layers):
        """Makes an estimate of the cluster diameter based on the average
        distance between opposite layers"""
        surface_mapping = data.lattice[self.symmetry]['surface_mapping']
        surface_data = data.lattice[self.symmetry]['surface_data']
        surface_count = data.lattice[self.symmetry]['surface_count']

        d = 0.0
        for i, l1 in enumerate(layers):
            l2 = layers[surface_mapping[i]]
            dl = surface_data[i]['d'] * self.lattice_constant
            d += (l1 + l2) * dl / surface_count

        return d

class SurfaceMonteCarloData:
    def __init__(self, filename=None, symbol=None, symmetry=None,
                 latticeconstant=None, fitsize=None, debug=0):
        if isinstance(symbol, str):
            self.atomic_number = data.atomic_numbers[symbol]
        else:
            self.atomic_number = symbol

        self.debug = debug
        self.filename = filename
        self.symmetry = symmetry
        self.lattice_constant = latticeconstant
        self.fitsize = fitsize

        self.arrays = {}
        self.arrays['surfaces'] = np.empty(0, dtype=int)
        self.arrays['energy'] = np.empty(0, dtype=float)
        self.arrays['size'] = np.empty(0, dtype=int)
        self.arrays['diameter'] = np.empty(0, dtype=float)
        self.arrays['counts'] = np.empty(0, dtype=int)

    def set_dimensions(self, surface_count):
        self.arrays['surfaces'] = self.arrays['surfaces'].reshape((0, surface_count))

    def __len__(self):
        return len(self.arrays['energy'])

    def __getitem__(self, index=-1):
        energy = self.arrays['energy'][index]
        surfaces = self.arrays['surfaces'][index]
        size = self.arrays['size'][index]
        diameter = self.arrays['diameter'][index]
        counts = self.arrays['counts'][index]

        return (energy, surfaces, size, diameter, counts)

    def __delitem__(self, index=-1):
        mask = np.ones(len(self), bool)
        mask[index] = False

        for name, a in self.arrays.items():
            self.arrays[name] = a[mask]

    def __setitem__(self, index, values):
        if values[0] is not None:
            self.arrays['energy'][index] = values[0]

        if values[1] is not None:
            self.arrays['surfaces'][index] = np.array(values[1], dtype=int)

        if values[2] is not None:
            self.arrays['size'][index] = values[2]

        if values[3] is not None:
            self.arrays['diameter'][index] = values[3]

        if values[4] is not None:
            self.arrays['counts'][index] = values[4]

    def append(self, values):
        for name, a in self.arrays.items():
            if name == 'counts':
                b = np.ones((a.shape[0] + 1,) + a.shape[1:], a.dtype)
            else:
                b = np.zeros((a.shape[0] + 1,) + a.shape[1:], a.dtype)

            b[:a.shape[0]] = a
            self.arrays[name] = b

        self[-1] = values

    def extend(self, other):
        if not isinstance(other, SurfaceMonteCarloData):
            raise Warning('The added object is not a SurfaceMonteCarloData instance!')

        if len(self) == 0:
            self.set_dimensions(other.arrays['surfaces'].shape[1])

        for name, a in self.arrays.items():
            b = other.arrays[name]
            if a.shape[1:] != b.shape[1:]:
                raise Warning('Different shapes: %s != %s' % (a.shape[1:], b.shape[1:]))
            c = np.zeros((a.shape[0] + b.shape[0],) + a.shape[1:], a.dtype)
            c[:a.shape[0]] = a
            c[a.shape[0]:] = b
            self.arrays[name] = c

    def write(self, filename=None):
        if filename is None:
            filename = self.filename

        if not isinstance(filename, str):
            raise Warning('You must specify a valid filename.')

        if os.path.isfile(filename):
            os.rename(filename, filename + '.bak')
        
        d = {'atomicnumber': self.atomic_number,
             'symmetry': self.symmetry,
             'latticeconstant': self.lattice_constant,
             'fitsize': self.fitsize}
             
        f = open(filename, 'w')
        f.write('ParameterMonteCarlo')
        pickle.dump(d, f)
        pickle.dump(self.arrays, f)
        f.close()

    def read(self, filename=None):
        if filename is None:
            filename = self.filename

        if not isinstance(filename, str):
            raise Warning('You must specify a valid filename.')

        if not os.path.isfile(filename):
            raise Warning('The file specified do not exist.')
        
        f = open(filename, 'r')

        try:
            if f.read(len('ParameterMonteCarlo')) != 'ParameterMonteCarlo':
                raise Warning('This is not a compatible file.')
            d = pickle.load(f)
            self.arrays = pickle.load(f)
        except EOFError:
            raise Warinig('Bad file.')
        
        f.close()

        if 'energies' in self.arrays:
            self.arrays['energy'] = self.arrays['energies']
            self.arrays['size'] = self.arrays['sizes']
            self.arrays['counts'] = np.zeros(len(self))
            del self.arrays['energies']
            del self.arrays['sizes']

        if 'multi' in self.arrays:
            self.arrays['counts'] = self.arrays['multi']
            del self.arrays['multi']

        self.atomic_number = d['atomicnumber']
        self.symmetry = d['symmetry']
        self.lattice_constant = d['latticeconstant']
        self.fitsize = d['fitsize']

    #Function to filter the data
    def filter_energy(self, excess):
        """Function that filters the data leaving only entries with
        energy lower than 'excess'."""

        e = self.get_energies()
        s = self.get_sizes()
        de = (e - e.min()) * s.mean()

        mask = de < excess
        for name, a in self.arrays.items():
            self.arrays[name] = a[mask]

    def filter_number(self, number):
        """Function that filters the data leaving only entries with
        the lowest energy in number less than 'number'."""

        if len(self) > number:
            e = self.get_energies()
            ind = e.argsort().argsort()
            mask = ind < number
            for name, a in self.arrays.items():
                self.arrays[name] = a[mask]

    def filter_surfaces(self):
        """Function that filters the data leaving only geometric
        different clusters."""

        symmetries = data.lattice[self.symmetry]['surface_symmetries']

        arrays_new = {}
        for name, a in self.arrays.items():
            arrays_new[name] = np.empty((0,) + a.shape[1:], dtype=a.dtype)
        lookup_table = {}
            
        for i, s in enumerate(self.arrays['surfaces']):
            f = -1

            for j, sym in enumerate(symmetries):
                if len(sym) != len(s):
                    raise Warning('Lenght of surfaces does not match with the symmetry.')

                s1 = s[sym]

                if self.symmetry == 'fcc':
                    s1 = data.fcc.surface_centering(s1)

                k = lookup_table.get(tuple(s1), None)

                if k is not None:
                    if f < 0 or f == k:
                        if self.debug:
                            print 'Surfaces %i found as %i with symmetry %i' % (i, k, j)
                        f = k
                    else:
                        print 'Surfaces #%i:\n%s' % (i, s)
                        print 'Symmetry %i' % j
                        print 'Found as filtered surfaces %i and ealier as %i:' % (k, f)
                        print '%s\n%s' % (list(s1), list(arrays_new['surfaces'][f]))
                        raise Warning(('Surfaces with applied symmetry ' +
                                       'found for different filtered surfaces!'))

            if f < 0:
                for name, a in arrays_new.items():
                    a1 = np.zeros((a.shape[0] + 1,) + a.shape[1:], a.dtype)
                    a1[:a.shape[0]] = a
                    a1[a.shape[0]:] = self.arrays[name][i]
                    arrays_new[name] = a1
                    length = len(a1)
                s_as_tuple = tuple(data.fcc.surface_centering(s))
                lookup_table[s_as_tuple] = length-1
            else:
                arrays_new['counts'][f] += 1

        self.arrays = arrays_new

    def get_multiplicity(self, i):
        symmetries = data.lattice[self.symmetry]['surface_symmetries']
        s = self.get_surfaces()[i]

        m = 0
        for i, sym in enumerate(symmetries):
            s1 = data.fcc.surface_centering(s[sym])
            if (s == s1).all(): m += 1

        return m

    def get_multiplicities(self):
        m = []
        for i in range(len(self)):
            m.append(self.get_multiplicity(i))

        return np.array(m, int)

    def get_surfaces(self):
        return self.arrays['surfaces'].copy()

    def get_energies(self):
        return self.arrays['energy'].copy()

    def get_sizes(self):
        return self.arrays['size'].copy()

    def get_diameters(self):
        return self.arrays['diameter'].copy()

    def get_counts(self):
        return self.arrays['counts'].copy()

    def get_cluster(self, i):
        return Cluster(symbol=self.atomic_number,
                       layers=self.get_surfaces()[i],
                       symmetry=self.symmetry,
                       latticeconstant=self.lattice_constant,
                       multiplicity=self.get_multiplicity(i))

    #Function to make a movie
    def make_movie(self, filename=None, frames=10, interval=None):
        if not isinstance(filename, str):
            raise Warning('The specified filename is not valid.')

        if isinstance(interval, list):
            if (interval[0] < 0 or interval[0] >= len(self) or
                interval[1] < 0 or interval[1] >= len(self)):
                raise Warning('Index out of range.')

            if interval[0] < interval[1]:
                start = interval[0]
                end = interval[1]
            else:
                start = interval[1]
                end = interval[0]
        elif interval is None:
            start = 0
            end = len(self) - 1
        else:
            raise Warning('Interval is not a list.')

        if frames < (end - start):
            delta = int(round((end - start) * 1.0 / frames))
        else:
            delta = 1

        images = range(start, end + 1, delta)

        if (end - images[-1]) > delta / 2.0:
            images.append(end - 1)
        else:
            images[-1] = end

        p = PrimiPlotter(self.get_cluster(images[0]))
        p.set_output(GifFile('SMCimg'))
        p.set_rotation([45,45,45])
        p.autoscale('now')
        p.autoscale('off')
        p.set_scale(0.8 * p.get_scale())
        p.update()
        #os.popen(('mogrify -draw \'text 10,10 "Frame %i of %i"\' SMCimg0000.gif' % 
        #          (images[0], len(self))))
        del images[0]
        
        for i, j in enumerate(images):
            p.update(self.get_cluster(j))
            #os.popen(('mogrify -draw \'text 10,10 "Frame %i of %i"\' SMCimg%04i.gif' %
            #          (j, len(self) - 1, i + 1)))

        os.popen('convert -delay 70 -loop 10 SMCimg*.gif ' + filename)
        os.popen('rm SMCimg*.gif -f')

    #OLD STUFF
    def filter_energy_old(self, number=None, energy=None):
        """Function that filters the data leaving only entries with
        energy lower the 'energy' and in a number that is lower or
        equal to 'number'."""

        if number is None:
            number = len(self)

        if energy is None:
            energy = self.arrays['energy'].max()

        arrays_new = {}
        for name, a in self.arrays.items():
            arrays_new[name] = np.empty((0,) + a.shape[1:], dtype=a.dtype)

        appended = 0
        for i, e in enumerate(self.arrays['energy']):
            if e <= energy and appended < number:
                appended += 1

                for name, a in arrays_new.items():
                    a1 = np.zeros((a.shape[0] + 1,) + a.shape[1:], a.dtype)
                    a1[:a.shape[0]] = a
                    a1[a.shape[0]:] = self.arrays[name][i]
                    arrays_new[name] = a1

            elif e <= energy and appended >= number:
                emax = arrays_new['energy'].max()
                imax = arrays_new['energy'].argmax()

                if e < emax:
                    for name, a in arrays_new.items():
                        a1 = np.zeros((a.shape[0] + 1,) + a.shape[1:], a.dtype)
                        a1[:a.shape[0]] = a
                        a1[a.shape[0]:] = self.arrays[name][i]
                        arrays_new[name] = a1

                    mask = np.ones(len(arrays_new['energy']), dtype=bool)
                    mask[imax] = False

                    for name, a in arrays_new.items():
                        arrays_new[name] = a[mask]

        self.arrays = arrays_new
