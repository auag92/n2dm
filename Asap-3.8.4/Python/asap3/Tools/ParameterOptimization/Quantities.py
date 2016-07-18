import sys
import time
import signal
import numpy as np

from asap3.Tools.ParameterOptimization.ScipyFmin import fmin, fmin_bfgs
from asap3.Tools.ElasticConstants import ElasticConstants, minimize_unit_cell
from asap3.Tools.SurfaceEnergies import SurfaceEnergy
from asap3.Tools.AtomicEnergies import CohesiveEnergy, VacancyEnergy

from ase import Atoms
from ase.io import read
from ase.data import reference_states, atomic_numbers, chemical_symbols
from ase.units import kJ, GPa
from ase.optimize import polyfit, BFGS, QuasiNewton
from ase.constraints import StrainFilter
from ase.lattice.cubic import BodyCenteredCubic, FaceCenteredCubic
from ase.lattice.hexagonal import HexagonalClosedPacked
from ase.lattice.surface import fcc100, fcc110, fcc111, hcp0001
from ase.lattice.compounds import B2, L1_2, L1_0
from ase.cluster.octahedron import Octahedron

class QuantityCalc:
    properties = {'fcc_size': (4, 4, 4), # Size of the FCC and HCP unit cell
                  'bcc_size': (6, 6, 6),
                  'hcp_size': (4, 4, 3),
                  'l12_size': (5, 5, 5),
                  'l10_size': (5, 5, 5),
                  'B2_size': (8, 8, 8),
                  'fcc100_size': (6, 6), # FCC and HCP surface size
                  'fcc111_size': (6, 6),
                  'hcp0001_size': (4, 6),
                  'surface_layers': [7, 9, 11, 13],
                  'sf_repeat' : (5,5,2),
                  }

    def __init__(self, latticeconstants, calc=None, timer=None, debug=False):
        self.debug = debug
        self.values = {}
        self.oldparams = None

        self.set_calculator(calc)

        if timer != None:
            self.timer = timer
        else:
            self.timer = None

        # Interpet lattice constants
        self.latticeconstants = {}
        for structure, elements, lc in latticeconstants:
            if isinstance(elements, (str, int)):
                elements = tuple([elements])
            else:
                elements = tuple(elements)
            self.latticeconstants[(structure, elements)] = lc

    def set_calculator(self, calc=None):
        if calc != None and not (hasattr(calc, 'get_potential_energy')
                                 or callable(calc)):
            raise ValueError('The calculator is not valid')
        self.calc = calc
        self.reset_values()
        
    def reset_values(self):
        """Reset all the calculated values that may have changed"""
        fullreset = (self.calc == None or self.oldparams == None)
        params = None
        try:
            params = self.calc.get_extra('parameters')
        except (ValueError, AttributeError):
            fullreset = True

        if fullreset:
            self.values = {}
            print >>sys.stderr, "Full reset performed."
        else:
            # Find which parameters have changed
            nbefore = len(self.values)
            changed = []
            for k,v in params.iteritems():
                assert isinstance(k, int)
                kstr = chemical_symbols[k]
                assert isinstance(kstr, str)
                if k not in self.oldparams or v != self.oldparams[k]:
                    # Delete any calculated value that includes this element
                    changed.append(k)
                    for prop in self.values.keys():  # NOT iterkeys !
                        if kstr in prop[2]:
                            del self.values[prop]
            print >>sys.stderr, ("Resetting %d values of %d: %s" %
                                 (nbefore - len(self.values), nbefore, str(changed)))
        self.oldparams = params
            
    ### Get quantity functions ###
    def get_quantity(self, name, structure, elements):
        func = getattr(self, 'get_' + name)
        return func(structure, elements)

    def get_cohesive_energy(self, structure, elements):
        return self.get_value('cohesive_energy', structure, elements,
                              self.calc_cohesive_energy)

    def get_vacancy_energy(self, structure, elements):
        return self.get_value('vacancy_energy', structure, elements,
                              self.calc_vacancy_energy)

    def get_phase_energy(self, structure, elements):
        return self.get_value('phase_energy', structure, elements,
                              self.calc_phase_energy)

    def get_lattice_constant_a(self, structure, elements):
        return self.get_value('lattice_constant_a', structure, elements,
                              self.calc_lattice_constants)

    def get_lattice_constant_c(self, structure, elements):
        return self.get_value('lattice_constant_c', structure, elements,
                              self.calc_lattice_constants)

    def get_lattice_ratio_ca(self, structure, elements):
        return self.get_value('lattice_ratio_ca', structure, elements,
                              self.calc_lattice_ratio_ca)

    def get_volume_per_atom(self, structure, elements):
        return self.get_value('volume_per_atom', structure, elements,
                              self.calc_volume_per_atom)

    def get_bulk_modulus(self, structure, elements):
        return self.get_value('bulk_modulus', structure, elements,
                              self.calc_elastic_constants)

    def get_elastic_anisotropy(self, structure, elements):
        return self.get_value('elastic_anisotropy', structure, elements,
                              self.calc_elastic_anisotropy)

    def get_elastic_constant_C11(self, structure, elements):
        return self.get_value('elastic_constant_C11', structure, elements,
                              self.calc_elastic_constants)

    def get_elastic_constant_C12(self, structure, elements):
        return self.get_value('elastic_constant_C12', structure, elements,
                              self.calc_elastic_constants)

    def get_elastic_constant_C13(self, structure, elements):
        return self.get_value('elastic_constant_C13', structure, elements,
                              self.calc_elastic_constants)

    def get_elastic_constant_C33(self, structure, elements):
        return self.get_value('elastic_constant_C33', structure, elements,
                              self.calc_elastic_constants)

    def get_elastic_constant_C44(self, structure, elements):
        return self.get_value('elastic_constant_C44', structure, elements,
                              self.calc_elastic_constants)

    def get_surface_energy(self, structure, elements):
        return self.get_value('surface_energy', structure, elements,
                              self.calc_surface_energy)

    def get_surface_ratio(self, structure, elements):
        return self.get_value('surface_ratio', structure, elements,
                              self.calc_surface_ratio)

    def get_heat_of_formation(self, structure, elements):
        return self.get_value('heat_of_formation', structure, elements,
                              self.calc_heat_of_formation)

    #def get_heat_of_sulution(self, structure, elements):
    #    return self.get_value('heat_of_solution', structure, elements,
    #                          self.calc_heat_of_solution)

    def get_impurity_energy(self, structure, elements):
        return self.get_value('impurity_energy', structure, elements,
                              self.calc_impurity_energy)

    def get_impurity_ratio(self, structure, elements):
        return self.get_value('impurity_ratio', structure, elements,
                              self.calc_impurity_ratio)

    def get_cutoff_energy(self, structure, elements):
        return self.get_value('cutoff_energy', structure, elements,
                              self.calc_cutoff_energy)

    def get_scaling_energy(self, structure, elements):
        return self.get_value('scaling_energy', structure, elements,
                              self.calc_scaling_energy)

    def get_force_match(self, structure, elements):
        return self.get_value('force_match', structure, elements,
                              self.calc_force_match)
        
    def get_stacking_fault(self, structure, elements):
        return self.get_value('stacking_fault', structure, elements,
                              self.calc_stacking_fault)
        
    ### Get and set function ###
    def get_value(self, name, structure, elements, function):
        if isinstance(elements, str):
            elements = (elements,)
        key = (name, structure, elements)
        if key not in self.values:
            if self.calc != None:
                if self.timer != None and self.debug:
                    self.timer.start(function.func_name)
                try:
                    function(structure, elements)
                except:
                    if self.debug:
                        sys.stderr.write("Calculation Error..." +
                                         "\n Name: %s" % (name,) + 
                                         "\n Structure: %s" % (structure,) +
                                         "\n Elements: %s\n\n" % (elements,))
                        sys.stderr.flush()
                    raise
                if self.timer != None and self.debug:
                    self.timer.stop(function.func_name)
            else:
                raise RuntimeError('No calculator set')

        return self.values[key]

    def set_value(self, name, structure, elements, value):
        key = (name, structure, elements)
        self.values[key] = value

    ### Calculation functions ###
    def calc_cohesive_energy(self, structure, elements):
        if structure not in ['fcc', 'hcp', 'bcc']:
            raise ValueError("Cannot calculate cohesive energy for '%s'" % (structure,))

        if self.debug:
            print "Calculating cohesive energy..."
            print 40 * "*"
            print "Structure: %s" % (structure,)
            print "Elements: %s\n" % (elements,)

        e_full = self.get_value('energy', structure, elements,
                                 self.calc_lattice_constants)

        atoms = Atoms([elements[0]], [(0.0, 0.0, 0.0)])
        atoms.set_calculator(self.get_calc())
        e_reduced = atoms.get_potential_energy()

        e = e_reduced - e_full

        self.values[('cohesive_energy', structure, elements)] = e

        if self.debug:
            print "Cohesive energy: %.5f eV/atom" % (e,)
            print 40 * "-" + "\n"

    def calc_vacancy_energy(self, structure, elements):
        if structure not in ['fcc', 'hcp', 'bcc']:
            raise ValueError("Cannot calculate vacancy energy for '%s'" % (structure,))

        if self.debug:
            print "Calculating vacancy energy..."
            print 40 * "*"
            print "Structure: %s" % (structure,)
            print "Elements: %s\n" % (elements,)

        atoms = self.get_structure(structure, elements)
        atoms.set_calculator(self.get_calc())
        e = VacancyEnergy(atoms)

        self.values[('vacancy_energy', structure, elements)] = e

        if self.debug:
            print "Vacancy energy: %.5f eV/atom" % (e,)
            print 40 * "-" + "\n"

    def calc_phase_energy(self, structure, elements):
        if structure not in ['fcc-hcp', 'fcc-bcc', 'bcc-hcp']:
            raise ValueError("Cannot calculate phase energy for '%s'" % (structure,))

        if self.debug:
            print "Calculating phase energy..."
            print 40 * "*"
            print "Structure: %s" % (structure,)
            print "Elements: %s\n" % (elements,)

        if structure[:3] == 'fcc':
            e1 = self.get_cohesive_energy('fcc', elements)
        else:
            e1 = self.get_cohesive_energy('bcc', elements)

        if structure[-3:] == 'hcp':
            e2 = self.get_cohesive_energy('hcp', elements)
        else:
            e2 = self.get_cohesive_energy('bcc', elements)

        self.values[('phase_energy', structure, elements)] = e1 - e2

        if self.debug:
            print "Phase energy: %.5f eV/atom" % (e1 - e2,)
            print 40 * "-" + "\n"

    def calc_elastic_constants(self, structure, elements):
        if structure not in ['fcc', 'hcp', 'bcc', 'l12', 'B2']:
            raise ValueError("Cannot calculate elasitc constants for '%s'" % (structure,))

        if self.debug:
            print "Calculating elastic constants..."
            print 40 * "*"
            print "Structure: %s" % (structure,)
            print "Elements: %s\n" % (elements,)

        atoms = self.get_structure(structure, elements)
        atoms.set_calculator(self.get_calc())

        if structure in ['fcc', 'bcc', 'l12', 'B2']:
            C11, C12, C44, B = ElasticConstants(atoms, 'cubic')
            self.values[('elastic_constant_C11', structure, elements)] = C11 / GPa
            self.values[('elastic_constant_C12', structure, elements)] = C12 / GPa
            self.values[('elastic_constant_C44', structure, elements)] = C44 / GPa
            self.values[('bulk_modulus', structure, elements)] = B / GPa
        else:
            C11, C12, C13, C33, C44, B = ElasticConstants(atoms, 'hexagonal', sanity=False)
            self.values[('elastic_constant_C11', structure, elements)] = C11 / GPa
            self.values[('elastic_constant_C12', structure, elements)] = C12 / GPa
            self.values[('elastic_constant_C13', structure, elements)] = C13 / GPa
            self.values[('elastic_constant_C33', structure, elements)] = C33 / GPa
            self.values[('elastic_constant_C44', structure, elements)] = C44 / GPa
            self.values[('bulk_modulus', structure, elements)] = B / GPa

        if self.debug:
            print "C11: %.5f GPa\nC12: %.5f GPa" % (C11 / GPa, C12 / GPa)
            if structure == 'hcp':
                print "C13: %.5f GPa\nC33: %.5f GPa" % (C13 / GPa, C33 / GPa)
            print "C44: %.5f GPa\nB: %.5f GPa" % (C44 / GPa, B / GPa)
            print 40 * "-" + "\n"

    def calc_elastic_anisotropy(self, structure, elements):
        if structure in ['fcc', 'bcc', 'l12', 'B2']:
            C11 = self.get_elastic_constant_C11(structure, elements)
            C12 = self.get_elastic_constant_C12(structure, elements)
            C44 = self.get_elastic_constant_C44(structure, elements)
            if abs(C11-C12) < 1e-5:
                raise ValueError('Cannot calculate anisotropy: C11=%g  C12=%g  C44=%g'
                                  % (C11, C12, C44))
            self.values[('elastic_anisotropy', structure, elements)] = 2*C44/(C11 - C12)
        else:
            raise ValueError("Cannot calculate elastic anisotropy for '%s'" % (structure,))

    def calc_lattice_constants(self, structure, elements):
        if structure not in ['fcc', 'hcp', 'bcc', 'l12', 'l10', 'B2']:
            raise ValueError("Cannot calculate lattice constants for '%s'" % (structure,))

        if self.debug:
            print "Calculating lattice constants..."
            print 40 * "*"
            print "Structure: %s" % (structure,)
            print "Elements: %s" % (elements,)

        if structure in ['fcc', 'bcc', 'l12', 'B2']:
            a = self.latticeconstants[(structure, elements)]
            c = None
        else:
            lc = self.latticeconstants[(structure, elements)]
            a = lc[0]
            c = lc[1]

        if self.debug:
            print "Guess: %s, %s\n" % (a,c)

        atoms = self.get_structure(structure, elements, a=a, c=c)
        atoms.set_calculator(self.get_calc())
        #filter = StrainFilter(atoms)
        #dyn = BFGS(filter, logfile=None, trajectory=None)
        #dyn.run(fmax=0.1)

        if structure in ['fcc', 'bcc', 'l12', 'B2']:
            MinimizeUnitCell(atoms, 'cubic', 0.0001)
            a = atoms.get_cell()[0][0] / self.properties[structure + '_size'][0]
            self.values[('lattice_constant_a', structure, elements)] = a
            #self.latticeconstants[(structure, elements)] = a
        else:
            MinimizeUnitCell(atoms, 'hexagonal', 0.0001)
            a = atoms.get_cell()[0][0] / self.properties[structure + '_size'][0]
            c = atoms.get_cell()[2][2] / self.properties[structure + '_size'][2]
            self.values[('lattice_constant_a', structure, elements)] = a
            self.values[('lattice_constant_c', structure, elements)] = c
            #self.latticeconstants[(structure, elements)] = [a, c]

        E = atoms.get_potential_energy() / len(atoms)
        self.values[('energy', structure, elements)] = E

        if self.debug:
            print "a: %.5f A" % (a,)
            if structure in ['hcp', 'l10']:
                print "c: %.5f A" % (c,)
            print "E: %.5f eV" % (E,)
            print 40 * "-" + "\n"

    def calc_lattice_ratio_ca(self, structure, elements):
        a = self.get_lattice_constant_a(structure, elements)
        c = self.get_lattice_constant_c(structure, elements)

        self.values[('lattice_ratio_ca', structure, elements)] = c / a

    def calc_volume_per_atom(self, structure, elements):
        if structure == 'fcc':
            a = self.get_lattice_constant_a(structure, elements)

            self.values[('volume_per_atom', structure, elements)] = a**3 / 4.0
        elif structure == 'hcp':
            a = self.get_lattice_constant_a(structure, elements)
            c = self.get_lattice_constant_c(structure, elements)

            self.values[('volume_per_atom', structure, elements)] = np.sqrt(3) * a**2 * c / 4.0
        else:
            raise ValueError("Cannot calculate lattice constants for '%s'" % (structure,))

    def calc_surface_energy(self, structure, elements):
        if structure not in ['fcc100', 'fcc111', 'hcp0001', 'l12100',
                             'l12111']:#, 'hcp1010A', 'hcp1010B']:
            raise ValueError("Cannot calculate surface energy for '%s'" % (structure,))

        if self.debug:
            print "Calculating surface energy..."
            print 40 * "*"
            print "Structure: %s" % (structure,)
            print "Elements: %s\n" % (elements,)

        images = []
        for n in self.properties['surface_layers']:
            images.append(self.get_structure(structure, elements, l=n))

        natoms = len(images[-1]) / self.properties['surface_layers'][-1]
        e = SurfaceEnergy(images, natoms, self.get_calc(), fmax=0.01, debug=self.debug)

        self.values[('surface_energy', structure, elements)] = e

        if self.debug:
            print "\nSurface energy: %.5f eV/atom" % (e,)
            print 40 * "-" + "\n"

    def calc_surface_ratio(self, structure, elements):
        structures = structure.split('-')
        if len(structures) != 2:
            raise ValueError("You must specify exactly two surfaces")

        e = []
        for struct in structures:
            e.append(self.get_surface_energy(struct, elements))

        self.values[('surface_ratio', structure, elements)] = e[0] / e[1]

    def calc_heat_of_formation(self, structure, elements):
        if structure not in ['l12', 'l10', 'B2']:
            raise ValueError("Cannot calculate heat of formation for '%s'" % (structure,))

        if self.debug:
            print "Calculating heat of formation..."
            print 40 * "*"
            print "Structure: %s" % (structure,)
            print "Elements: %s\n" % (elements,)

        e_single = {}
        for symbol in elements:
            sym = reference_states[atomic_numbers[symbol]]['symmetry']
            #atoms = self.get_structure(sym, (symbol,))
            #atoms.set_calculator(self.get_calc())
            #e_single[symbol] = atoms.get_potential_energy() / len(atoms)
            e_single[symbol] = self.get_value('energy', sym, (symbol,),
                                              self.calc_lattice_constants)


        #atoms = self.get_structure(structure, elements)
        #atoms.set_calculator(self.get_calc())
        #e_alloy = atoms.get_potential_energy() / len(atoms)
        e_alloy = self.get_value('energy', structure, elements,
                                 self.calc_lattice_constants)

        if structure == 'l12':
            e = e_alloy - 0.25 * e_single[elements[0]] \
                        - 0.75 * e_single[elements[1]]
        else:
            e = e_alloy - 0.5 * e_single[elements[0]] \
                        - 0.5 * e_single[elements[1]]

        self.values[('heat_of_formation', structure, elements)] = e

        if self.debug:
            print "Heat of formation: %.5f eV/atom" % (e,)
            print 40 * "-" + "\n"

    #def calc_heat_of_solution(self, structure, elements):
    # lim(N->inf) {E(Pt_N Y) - N*E(Pt) - E(Y)}

    def calc_impurity_energy(self, structure, elements):
        if structure not in ['oct38-center', 'oct38-face', 'oct38-edge']:
            raise ValueError("Cannot calculate impurity energy for '%s'" % (structure,))

        if len(elements) != 2:
            raise ValueError("Tuple of elements must be of length two")

        if self.debug:
            print "Calculating impurity energy..."
            print 40 * "*"
            print "Structure: %s" % (structure,)
            print "Elements: %s\n" % (elements,)

        name, impurity = structure.split('-')
        sym = reference_states[atomic_numbers[elements[1]]]['symmetry']
        latticeconstant = self.get_lattice_constant_a(sym, (elements[1],))

        if name == 'oct38':
            sites = {'center': 10, 'face': 9, 'edge': 0}

            atoms = Octahedron(elements[1], 4, 1, latticeconstant)
            atoms.set_calculator(self.get_calc())

            dyn = BFGS(atoms, logfile=None, trajectory=None)
            dyn.run(fmax=0.001, steps=100)
            assert dyn.converged()
            s_clean = dyn.get_number_of_steps()
            e_clean = atoms.get_potential_energy()

            atoms[sites[impurity]].symbol = elements[0]

            dyn.run(fmax=0.001, steps=100)
            assert dyn.converged()
            s_impurity = dyn.get_number_of_steps()
            e_impurity = atoms.get_potential_energy()

        self.values[('impurity_energy', structure, elements)] = e_impurity - e_clean

        if self.debug:
            print "BFGS steps: %i    %i" % (s_clean, s_impurity)
            print "Impurity energy: %.5f - %.5f = %.5f eV" % (e_impurity, e_clean,
                                                              e_impurity - e_clean)
            print 40 * "-" + "\n"

    def calc_impurity_ratio(self, structure, elements):
        structures = structure.split('/')
        if len(structures) != 2:
            raise ValueError("You must specify exactly two structures")

        e = []
        for struct in structures:
            e.append(self.get_impurity_energy(struct, elements))

        self.values[('impurity_ratio', structure, elements)] = e[0] / e[1]

    def calc_cutoff_energy(self, structure, elements):
        if structure not in ['fcc', 'hcp', 'l12']:
            raise ValueError("Cannot calculate cutoff energy for '%s'" % (structure,))

        if self.debug:
            print "Calculating cutoff energy..."
            print 40 * "*"
            print "Structure: %s" % (structure,)
            print "Elements: %s\n" % (elements,)

        if structure in ['fcc', 'l12']:
            a = self.latticeconstants[(structure, elements)]
            c = None
            f = 1.73
        else:
            lc = self.latticeconstants[(structure, elements)]
            a = lc[0]
            c = lc[1]
            f = 1.73

        atoms = self.get_structure(structure, elements, a=a, c=c)
        atoms.set_calculator(self.get_calc())

        e_equilibrium = atoms.get_potential_energy()

        new_cell = f * atoms.get_cell()
        atoms.set_cell(new_cell, True)

        e_cutoff = atoms.get_potential_energy()

        self.values[('cutoff_energy', structure, elements)] = e_cutoff / e_equilibrium

        if self.debug:
            print "Equilibrium: %.5f eV" % (e_equilibrium,)
            print "Cutoff: %.5f eV" % (e_cutoff,)
            print "Ratio: %.5f" % (e_cutoff / e_equilibrium,)
            print 40 * "-" + "\n"

    def calc_scaling_energy(self, structure, elements):
        t = structure.split('-')
        structure = t[0]
        scaling = float(t[1])

        if structure not in ['fcc', 'hcp', 'l12']:
            raise ValueError("Cannot calculate scaling energy for '%s'" % (structure,))

        if self.debug:
            print "Calculating scaling energy..."
            print 40 * "*"
            print "Structure: %s" % (structure,)
            print "Elements: %s\n" % (elements,)

        if structure in ['fcc', 'l12']:
            a = self.latticeconstants[(structure, elements)]
            c = None
            f = 1.73
        else:
            lc = self.latticeconstants[(structure, elements)]
            a = lc[0]
            c = lc[1]
            f = 1.73

        atoms = self.get_structure(structure, elements, a=a, c=c)
        atoms.set_calculator(self.get_calc())
        new_cell = scaling * atoms.get_cell()
        atoms.set_cell(new_cell, True)
        energy = atoms.get_potential_energy() / len(atoms)

        structure = '%s-%.2f' % (structure, scaling)
        self.values[('scaling_energy', structure, elements)] = energy

        if self.debug:
            print "Energy: %.5f eV" % (energy,)
            print 40 * "-" + "\n"

    def calc_force_match(self, structure, elements):
        if self.debug:
            print "Calculating force matching..."
            print 40 * "*"
            print "Structure: %s" % (structure,)
            print "Elements: %s\n" % (elements,)
        atoms = read(structure)
        f_dft = atoms.get_forces()
        atoms = atoms.repeat((2,2,2))
        atoms.set_calculator(self.get_calc())
        f = atoms.get_forces()[:len(f_dft)]
        df = f_dft - f
        err = np.sqrt((df*df).sum())
        err /= np.sqrt((f_dft*f_dft).sum())
        self.values[('force_match', structure, elements)] = err
        if self.debug:
            print "Mean force error: %.5f eV/A" % (err,)
            print 40 * "-" + "\n"


    def calc_stacking_fault(self, structure, elements):
        if structure != 'fcc':
            raise ValueError("Cannot calculate stacking fault energy of structure "+structure)
        a = self.latticeconstants[(structure, elements)]
        atoms = fcc111(elements[0], (1,2,5), orthogonal=True, a=a)
        atoms.set_pbc(True)
        atoms = atoms.repeat(self.properties['sf_repeat'])
        atoms.set_calculator(self.get_calc())
        dyn = QuasiNewton(atoms, logfile=None, trajectory=None)
        dyn.run(fmax=0.02)
        e_sf = atoms.get_potential_energy()
        n_sf = len(atoms)
        atoms = fcc111(elements[0], (1,2,6), orthogonal=True, a=a)
        atoms.set_pbc(True)
        atoms = atoms.repeat(self.properties['sf_repeat'])
        atoms.set_calculator(self.get_calc())
        dyn = QuasiNewton(atoms, logfile=None, trajectory=None)
        dyn.run(fmax=0.02)
        e_bulk = atoms.get_potential_energy()
        n_bulk = len(atoms)
        layers = self.properties['sf_repeat'][2]
        uc = atoms.get_cell()
        area = uc[0,0] * uc[1,1]
        result = (e_sf - e_bulk * n_sf / n_bulk) / layers / area
        result /= 1e-3 * kJ * 1e-20
        self.values[('stacking_fault', structure, elements)] = result
        
    ### Bulk structures ###
    def get_structure(self, name, elements, a=None, c= None, l=None):
        # Check number of elements
        if name[:3] in ['fcc', 'hcp']:
            if len(elements) != 1:
                raise ValueError("Tuple of elements must be of length one")
        if name[:3] in ['l12', 'l10'] or name[:2] == 'B2':
            if len(elements) != 2:
                raise ValueError("Tuple of elements must be of length two")

        # Get lattice constants
        if a is None:
            if name[:2] == 'B2':
                a = self.get_lattice_constant_a(name[:2], elements)
            elif name[:3] in ['fcc', 'hcp', 'bcc', 'l12', 'l10']:
                a = self.get_lattice_constant_a(name[:3], elements)

        if c is None:
            if name[:3] in ['hcp', 'l10']:
                c = self.get_lattice_constant_c(name[:3], elements)

        # Get size
        if name in ['fcc', 'hcp', 'bcc', 'l12', 'l10', 'B2']:
            size = self.properties[name + '_size']
        elif name in ['fcc100', 'fcc111', 'hcp0001']:
            size = self.properties[name + '_size'][:2] + (l,)

        # Make structure
        if name == 'fcc':
            atoms = FaceCenteredCubic(symbol=elements[0],
                                      size=size,
                                      latticeconstant=a)
        elif name == 'hcp':
            atoms = HexagonalClosedPacked(symbol=elements[0], size=size,
                                          directions=[[2,-1,-1,0],[0,1,-1,0],[0,0,0,1]],
                                          latticeconstant=(a, c))
        elif name == 'bcc':
            atoms = BodyCenteredCubic(symbol=elements[0],
                                      size=size,
                                      latticeconstant=a)
        elif name == 'B2':
            atoms = B2(symbol=elements, size=size, latticeconstant=a)
        elif name == 'l12':
            atoms = L1_2(symbol=elements, size=size, latticeconstant=a)
        elif name == 'l10':
            atoms = L1_0(symbol=elements, size=size, latticeconstant=(a, c))
        elif name == 'fcc100':
            atoms = fcc100(symbol=elements[0], size=size, a=a, vacuum=10.0)
        elif name == 'fcc111':
            atoms = fcc111(symbol=elements[0], size=size, a=a,
                           vacuum=10.0, orthogonal=True)
        elif name == 'hcp0001':
            atoms = hcp0001(symbol=elements[0], size=size, a=a, c=c,
                            vacuum=10.0, orthogonal=True)
        elif name == 'hcp1010A':
            raise ValueError("Structure '%s' not supported" % (name,))
            atoms = None
        elif name == 'hcp1010B':
            raise ValueError("Structure '%s' not supported" % (name,))
            atoms = None
        elif name == 'l12100':
            n = (l + 1) / 2
            atoms = L1_2(symbol=elements, size=(8, 8, n), latticeconstant=a)
            atoms.set_pbc([True, True, False])
            # Remove layers
            atoms = atoms[atoms.get_positions()[:,2] > 0.1 * a]
            # Set vacuum
            atoms.center(axis=2, vacuum=10.0)
        elif name == 'l12111':
            if l % 3 == 0:
                n = l / 3
                c = 0
            else:
                n = l / 3 + 1
                c = 3 - l % 3
            atoms = L1_2(symbol=elements, size=(8, 4, n),
                         #directions=[[1,-1,0],[1,0,-1],[1,1,1]], latticeconstant=a)
                         directions=[[1,-1,0],[1,1,-2],[1,1,1]], latticeconstant=a)
            atoms.set_pbc([True, True, False])
            # Wrap positions
            scpos = atoms.get_scaled_positions()
            scpos[scpos > (1.0 - 1e-12)] = 0.0
            atoms.set_scaled_positions(scpos)
            # Remove layers
            if c > 0:
                atoms = atoms[atoms.get_positions()[:,2] > (c - 0.5) * a / np.sqrt(3.0)]
            # Set vacuum
            atoms.center(axis=2, vacuum=10.0)
        else:
            raise ValueError("Structure '%s' not supported" % (name,))
        return atoms
    
    def get_calc(self):
        if hasattr(self.calc, 'get_potential_energy'):
            return self.calc
        else:
            return self.calc()
        

def MinimizeUnitCell(atoms, symmetry, ftol):
    from asap3.Tools.ParameterOptimization.ScipyFmin import fmin
    if symmetry == 'cubic':
        var = [1.0]  # Scaling
    elif symmetry == 'hexagonal':
        var = [1.0, 1.0]
    else:
        raise ValueError('Symmetry "%s" not supported' % (symmetry,))
    u0 = atoms.get_cell()
    def energy(v, a=atoms, u=u0):
        newcell = v[0] * u
        if len(v) > 1:
            newcell[2] *= v[1]
        a.set_cell(newcell, scale_atoms=True)
        return a.get_potential_energy()
    xopt, fopt, iter, calls, flag = fmin(energy, var, ftol=ftol, delta=0.01,
                                         full_output=True, disp=False)
    newcell = xopt[0] * u0
    if len(xopt) > 1:
        newcell[2] *= xopt[1]
    atoms.set_cell(newcell, scale_atoms=True)

if __name__ == '__main__':
    from asap3 import EMT, EMT2013
    from asap3.Tools.Timing import Timer
    from asap3.Tools.ParameterOptimization.EMT import EMT2011Fit

    timer = Timer()
    quantities = [('lattice_constant_a', 'fcc', 'Pt', 3.92),
                  ('lattice_constant_a', 'hcp', 'Pt', 2.77),
                  ('lattice_ratio_ca', 'hcp', 'Pt', 4.78 / 2.77),
                  ('bulk_modulus', 'fcc', 'Pt', 278.3),
                  ('elastic_anisotropy', 'fcc', 'Pt', 1.594),
                  #('elastic_constant_C11', 'fcc', 'Pt', 346.7),
                  #('elastic_constant_C12', 'fcc', 'Pt', 250.7),
                  ('elastic_constant_C44', 'fcc', 'Pt', 76.5),
                  ('cohesive_energy', 'fcc', 'Pt', 5.84),
                  #('phase_energy', 'fcc-hcp', 'Pt', -0.05),
                  ('surface_energy', 'fcc111', 'Pt', 0.631),
                  ('surface_ratio', 'fcc111-fcc100', 'Pt', 0.631 / 0.892),
                  #('lattice_constant_a', 'l12', ('Pt', 'Cu'), 3.67860),
                  #('heat_of_formation', 'l12', ('Pt', 'Cu'), -0.12331),
                  ]

    latticeconstants = [('fcc', 'Pt', 3.92),
                        ('hcp', 'Pt', (2.77, 4.78)),
                        ('fcc', 'Cu', 3.5),
                        ('l12', ('Pt', 'Cu'), 3.8)]

    
    calc = QuantityCalc(latticeconstants, EMT(), timer, True)
    for name, structure, elements, value in quantities:
        value = calc.get_quantity(name, structure, elements)

    timer.write()

