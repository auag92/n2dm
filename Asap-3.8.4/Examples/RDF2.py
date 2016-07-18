from asap3 import *
from asap3.analysis.rdf import RadialDistributionFunction
from ase.lattice.compounds import *
import matplotlib.pyplot as plt
import numpy as np

# We want RDF out to a distance of 15 Angstrom, with 200 bins
rng=10.0
bins = 200

# Helper function calculating coordination number from RDF
def getcoord(rdf, r, dens):
    "Calculate coordination number from RDF."
    # Find first peak by assuming it is the largest
    imax = rdf.argmax()
    # Find first minimum, assuming it is the deepest after the maximum.
    imin = imax + rdf[imax:].argmin()
    print "Position of first peak and following minimum:", imax, imin
    dr = r[1] - r[0]
    integral = 0.0
    for i in range(imin):
        integral += 4 * np.pi * r[i] * r[i] * dr * rdf[i]
    return integral * dens

# Create a block of Cu atoms, expand the crystal 10% to leave space
# for thermal expansion and to facilitate melting.  Describe
# interactions with EMT, and set a very high initial temperature so it
# melts.
atoms = L1_2(size=(20,20,20), symbol=("Au", "Cu"), pbc=True,
             latticeconstant=3.75)

atoms.set_calculator(EMT())
MaxwellBoltzmannDistribution(atoms, 600*units.kB)

# First, run 25 time steps before taking data.
dyn = VelocityVerlet(atoms, 5*units.fs, logfile='-')
dyn.run(25)

# Make the RDF object, attach it to the dynamics
RDFobj = RadialDistributionFunction(atoms, rng, bins, verbose=True)
dyn.attach(RDFobj.update, interval=5)

# Run the simulation
dyn.run(100)

# Get the RDF and plot it.
rdf = RDFobj.get_rdf()
x = np.arange(bins) * rng / bins
plt.plot(x, rdf, 'k')  # Black

# Get the partial RDFs and plot them
Au = data.atomic_numbers['Au']
Cu = data.atomic_numbers['Cu']
rdfAuAu = RDFobj.get_rdf(elements=(Au, Au))
plt.plot(x, rdfAuAu, 'r')
rdfCuAu = RDFobj.get_rdf(elements=(Cu, Au))
plt.plot(x, rdfCuAu, 'b--')
rdfAuCu = RDFobj.get_rdf(elements=(Au, Cu))
plt.plot(x, rdfAuCu, 'c')
rdfCuCu = RDFobj.get_rdf(elements=(Cu, Cu))
plt.plot(x, rdfCuCu, 'g')

rho = len(atoms)/atoms.get_volume()
c = getcoord(rdf, x, rho)
cAuAu = getcoord(rdfAuAu, x, rho)
cAuAu = 0.0  # RDF found next-nearest neighbors as no nearest.
cCuAu = getcoord(rdfCuAu, x, rho) / 0.75
cAuCu = getcoord(rdfAuCu, x, rho) / 0.25
cCuCu = getcoord(rdfCuCu, x, rho) / 0.75

print "Total coodination number:", c
print "Cu atoms around an Cu atom:", cCuCu
print "Au atoms around an Au atom:", cAuAu
print "Cu atoms around an Au atom:", cAuCu
print "Au atoms around an Cu atom:", cCuAu

plt.show()

