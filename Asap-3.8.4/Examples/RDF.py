from asap3 import *
from asap3.analysis.rdf import RadialDistributionFunction
from ase.lattice.cubic import *
import matplotlib.pyplot as plt

# We want RDF out to a distance of 15 Angstrom, with 200 bins
rng=15.0
bins = 200

# Create a block of Cu atoms, expand the crystal 10% to leave space
# for thermal expansion and to facilitate melting.  Describe
# interactions with EMT, and set a very high initial temperature so it
# melts.
atoms = FaceCenteredCubic(size=(20,20,20), symbol="Cu", pbc=True)
atoms.set_cell(atoms.get_cell() * 1.1, scale_atoms=True)

atoms.set_calculator(EMT())
MaxwellBoltzmannDistribution(atoms, 3500*units.kB)

# First, run 25 time steps before taking data.
dyn = VelocityVerlet(atoms, 4*units.fs, logfile='-')
dyn.run(25)

# Make the RDF object, attach it to the dynamics
RDFobj = RadialDistributionFunction(atoms, rng, bins, verbose=True)
dyn.attach(RDFobj.update, interval=5)

# Run the simulation
dyn.run(100)

# Get the RDF and plot it.
rdf = RDFobj.get_rdf()
x = np.arange(bins) * rng / bins
plt.plot(x, rdf)
plt.show()
