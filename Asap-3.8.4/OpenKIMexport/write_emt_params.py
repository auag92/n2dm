
from asap3.Internal.EMTParameters import _std_emt_params as std_emt_params
import ase.data

outfile = 'asap_emt/asap_emt.params'
makefileinclude = 'asap_emt/Species.mk'

#parameternames = ("E0", "S0", "V0", "eta2", 
#         "kappa": 2.000 / bohr,
#         "lambda": 1.169 / bohr,
#         "n0": 0.007 / bohr3},)

comment = '''
(anything below the STOP line is ignored)

These are the default EMT parameters from ASAP.  
Units are based on eV, Angstrom and atomic mass units.
The mass parameter will not actually be used by the KIM calculator.
'''

out = open(outfile, 'w')
spec = open(makefileinclude, 'w')

fakez = 0
for z, params in std_emt_params.iteritems():
    fakez += 1
    out.write('NEWELEMENT %i %s\n' % (fakez, params['name']))
    spec.write('SPECIES%d_NAME := %s\n' % (fakez, params['name']))
    params["mass"] = ase.data.atomic_masses[z] 
    for name, value in params.iteritems():
        if name != 'name':
            out.write('  %s %.15f\n' % (name, value))
    out.write('ENDELEMENT -1\n')
    
out.write('STOP -1\n\n') 
out.write(comment)       
out.close()
spec.close()


        