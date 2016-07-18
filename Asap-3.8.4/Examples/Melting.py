"Melting a copper cluster."

from numpy import *
from asap3 import Atoms, EMT, units
from ase.visualize.primiplotter import *
from ase.lattice.cubic import FaceCenteredCubic
from asap3.md.langevin import Langevin

# Create the atoms
atoms = FaceCenteredCubic(size=(5,5,5), symbol="Cu", pbc=False)


# Associate the EMT potential with the atoms
atoms.set_calculator(EMT())

# Temperature profile
temperatures = (250, 500, 750, 1000, 1250, 1500, 1750)

# How many steps at each temperature
nsteps = 5000

# Interval between plots
plotinterval = 500

# Make the Langevin dynamics module
dyn = Langevin(atoms, 5*units.fs, units.kB*temperatures[0], 0.002)

# The plotter
plotter = PrimiPlotter(atoms)
plotter.set_output(X11Window())
plotter.set_rotation((10.0, 5.0, 0))
dyn.attach(plotter.plot, interval=plotinterval)

# Some functions for calculating the actual temperature, energy, ...

def temperature(a):
    return 2.0/3.0 * a.get_kinetic_energy() / (len(a) * units.kB)

def etotal(a):
    return (a.get_kinetic_energy() + a.get_potential_energy()) / len(a)

# The main loop

for t in temperatures:
    dyn.set_temperature(units.kB*t)
    for i in range(nsteps/100):
        dyn.run(100)
        print "E_total = %-10.5f  T = %.0f K  (goal: %.0f K, step %d of %d)" %\
              (etotal(atoms), temperature(atoms), t, i, nsteps/100)
        
