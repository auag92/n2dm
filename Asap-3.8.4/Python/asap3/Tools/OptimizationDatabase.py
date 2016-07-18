# encoding: utf-8

"""OptimizationDatabase - target values for parameter optimization.

This file contains a database of expermental and DFT-derived values
useful for optimizing classical potentials.

There are two dictionaries, elements and alloys, defining data for 
pure elements and for alloys, respectively.  The keys are chemical 
symbols.  The values are dictionaries of data.

The value dictionaries:  The keys are standard property names,
defined below.  If it makes sense, the string '_0' or '_300' is 
appended to the name to indicate the approximate temperature.
The values are tuples of (value, type, source).

The values should be given in standard units.  The type is either 
EXP (experimental) or DFT-PBE (calculated with DFT-PBE).  The 
source is a string (define a variable, so they can be reused!) 
giving the reference or a descriptive text.

Lattice constants are NOT included in the database, they are available
in ASE.

Data is accessed throught the access function XXXXX.

Property names and units::

B
    Bulk modulus (GPa)

C11, C12, C44
    Cubic elastic constants (GPa)
    
C11, C12, C13, C33, C44
    Hexagonal elastic constants (GPa)
    
Ecoh
    Cohesive energy (eV)  
    
Evac
    Vacancy formation energy (eV)
    
E100, E111
    Surface energies of 100, 111 (etc) surfaces (J/m^2)
    
Eisf
    Intrinsic stacking fault energy (J/m^2)
    
Ehcpfcc
    Energy difference between HCP and FCC structure, both with
    their optimal lattice constants
    
Ebccfcc
    Similar, but BCC minus FCC energy.
    
"""

import ase.data

##################################################
###                                            ###
### List of variables defining the references  ###
###                                            ###
##################################################

Frederikse = 'H.P.R Frederikse: "Elastic constants of single crystals", in CRC Handbook of Chemistry and Physics, 91st edition (2010).'
webelements = "www.webelements.com"
Kittel = 'C. Kittel: "Introduction to Solid State Physics", 8th ed.'
Triftshauser = "TRIFTSHAUSER, W., & MCGERVEY, J. (1975). Monovacancy Formation Energy in Copper, Silver, and Gold by Positron-Annihilation. Applied Physics, 6(2), 177�180."
JS_Surf = "Calculated by J. Schi�tz with GPAW (PBE), see /home/camp/schiotz/simulations/PotentialFit/SurfaceEnergies"

#####################################
###                               ###
### Properties for pure elements  ###
###                               ###
#####################################

elements = \
{
    'Cu': {'B_300': (140, 'EXP', webelements),
           'C11_300': (168.3, 'EXP', Frederikse),
           'C12_300': (122.1, 'EXP', Frederikse),
           'C44_300': (75.7, 'EXP', Frederikse),
           'C11_0': (176.2, 'EXP', Kittel),
           'C12_0': (124.9, 'EXP', Kittel),
           'C44_0': (81.8, 'EXP', Kittel),
           'B_0' : True,  # Can be calculated
           'Ecoh': (3.49, 'EXP', Kittel),
           'Evac': (1.29, 'EXP', Triftshauser),
           'E100': (2.888, 'PBE', JS_Surf),
           'E111': (2.635, 'PBE', JS_Surf),
           },
    'Ag': {
           'Evac': (1.29, 'EXP', Triftshauser),
           },
    'Au': {
           'Evac': (0.97, 'EXP', Triftshauser),
           },
}


#########################
###                   ###
###  Access function  ###
###                   ###
#########################

def get_data(element, property, T=0):
    """Get a data item from the database."""
    elemdata = elements[element]
    data = None
    for key in [property, "%s_%i" % (property, T), property+'_0', property+'_300']:
        if key in elemdata:
            data = elemdata[key]
            break
    if data is True:
        return calculate_data(element, key)
    elif data is None:
        # Not found.
        if property in ['a', 'c/a']:
            z = ase.data.atomic_numbers[element]
            return ase.data.reference_states[z][property]
        else:
            raise KeyError("Property %s for %s not found in database." 
                           % (property, element))
    else:
        return data[0]
    
def calculate_data(element, key):
    """Calculate data derived from other data points."""
    if '_' in key:
        key, suffix = key.split('_')
        suffix = '_' + suffix
    else:
        suffix = ""
    if key == 'B':
        return calculate_B(element, suffix)
    else:
        raise RuntimeError("Cannot calculate quantity " + key + suffix)
    
def calculate_B(element, suffix):
    z = ase.data.atomic_numbers[element]
    struct = ase.data.reference_states[z]['symmetry']
    elemdata = elements[element]
    if struct == 'fcc':
        C11 = elemdata['C11'+suffix][0]
        C12 = elemdata['C12'+suffix][0]
        return (C11 + 2*C12) / 3.0
    else:
        raise RuntimeError("Cannot calculate B for %s in structure %s" % (element, struct))
