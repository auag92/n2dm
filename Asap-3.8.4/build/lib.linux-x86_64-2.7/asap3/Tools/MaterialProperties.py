"""MaterialPropertiesData - module to handle material property data

The module loads the metarial data form one or more files with the format:

    Element(s);Structure;Property;Value;Temp;Reference
    Al;fcc;Ecoh;3.39;0;C. Kittel, Introduction to solid state physics, ...
    Al,Ag;l12;Eheat;-0.08708;0;Jakob Schi\otz 2013

Property names and units::

a, c
    Lattice constants

B
    Bulk modulus (GPa)

A
    Cubic anisotopy ratio = 2*C44/(C11-C12) (GPa)

C11, C12, C44
    Cubic elastic constants (GPa)
    
C11, C12, C13, C33, C44
    Hexagonal elastic constants (GPa)
    
Ecoh
    Cohesive energy (eV)  
    
Evac
    Vacancy formation energy (eV)

Eheat
    Alloy heat of formation (eV)

E100, E111
    Surface energies of 100, 111 (etc) surfaces (eV/atom)

Eisf
    Intrinsic stacking fault energy (J/m^2)

Ehcpfcc
    Energy difference between HCP and FCC structure, both with
    their optimal lattice constants

Ebccfcc
    Similar, but BCC minus FCC energy.

"""

from ase.data import atomic_numbers, reference_states

class MaterialPropertiesData:
    def __init__(self, filenames):
        self.data = {}
        if isinstance(filenames, str):
            filenames = [filenames,]
        for file in filenames:
            self.read(file)

    def get(self, elements, property, structure=None, temperature=None,
            calculate=False, reference=False):
        if isinstance(elements, str):
            elements = (elements,)
        elif isinstance(elements, list):
            elements = tuple(elements)

        if structure == None:
            if len(elements) == 1:
                Z = atomic_numbers[elements[0]]
                structure = reference_states[Z]['symmetry']
            else:
                raise ValueError('Structure must be specified for an alloy.')

        if not self.data.has_key(elements):
            if len(elements) == 1 and property in ['a', 'c']:
                Z = atomic_numbers[elements[0]]
                assert structure == reference_states[Z]['symmetry']
                a = reference_states[Z]['a']
                if property == 'c':
                    return a * reference_states[Z]['c/a']
                else:
                    return a
            else:
                raise ValueError('No data for the element(s): %s.' % (elements,))

        key = structure.lower() + '_' + property
        data = self.data[elements]
        if not data.has_key(key) or calculate:
            if property == 'B':
                value = self.calc_B(elements, structure, temperature)
            elif property == 'A':
                value = self.calc_A(elements, structure, temperature)
            else:
                raise ValueError('Cannot find or calculate the property %s %s.' % (key, property))
        else:
            found = False
            if temperature != None:
                for (val, temp, ref) in data[key]:
                    if temp == temperature:
                        found = True
                        break
            if not found:
                (val, temp, ref) = data[key][0]

            if reference:
                value = ref
            else:
                value = val

        return value

    def calc_B(self, elements, structure=None, temperature=None):
        C11 = self.get(elements, 'C11', structure, temperature)
        C12 = self.get(elements, 'C12', structure, temperature)
        return 1.0/3.0*(C11 + 2*C12)

    def calc_A(self, elements, structure=None, temperature=None):
        C11 = self.get(elements, 'C11', structure, temperature)
        C12 = self.get(elements, 'C12', structure, temperature)
        C44 = self.get(elements, 'C44', structure, temperature)
        return 2.0*C44/(C11 - C12)

    def read(self, filename):
        f = open(filename)
        for line in f.readlines():
            line = line[:-1].split(';') #'\t')

            elements = tuple(line[0].split(','))
            if not atomic_numbers.has_key(elements[0]):
                continue

            if not self.data.has_key(elements):
                self.data[elements] = {}

            property = '_'.join(line[1:3])
            if not self.data[elements].has_key(property):
                self.data[elements][property] = []

            data = (float(line[3]), int(line[4]), line[5])
            self.data[elements][property].append(data)
        f.close()

    def write(self, filename):
        f = open(filename, 'w')
        f.write('Element(s);Structure;Property;Value;Temperature;Reference\n')
        for elements in self.data.keys():
            for property in self.data[elements].keys():
                for data in self.data[elements][property]:
                    data = elements + tuple(property.split('_')) + data
                    if len(elements) > 1:
                        format = '%s,%s;%s;%s;%f;%i;%s\n'
                    else:
                        format = '%s;%s;%s;%f;%i;%s\n'
                    f.write(format % tuple(data))
        f.close()

if __name__ == '__main__':
    print 'Testing...'
    mp_data = MaterialPropertiesData('material_properties.csv')
    print mp_data.get('Pt', 'Ecoh')

