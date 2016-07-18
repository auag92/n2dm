"""Experimental parameters for EMT2013."""

from ase.units import Bohr

EMT2013_experimental_parameters = {
    'Cu': {'E0': -3.51, 'S0': 2.67*Bohr, 'V0': 2.476, 'eta2': 1.652/Bohr,  'kappa': 2.740/Bohr,  
           'lambda': 1.906/Bohr,  'n0': 0.00910/(Bohr**3)}
    }

PtY_parameters = {
    'Calculator': "EMT2013_V1",
    'Pt': {'E0': -5.82286,
           'S0': 1.55073,
           'V0': 2.69717,
           'eta2': 2.41957,
           'kappa': 3.86730,
           'lambda': 4.02350,
           'n0': 0.05412},

    'Y': {'E0': -4.22967,
          'S0': 2.12796,
          'V0': 5.93342,
          'eta2': 1.65912,
          'kappa': 2.78458,
          'lambda': 2.14360,
          'n0': 0.02626882869955157},  # Absurdly many digits to maintain result after a bug fix.
}
