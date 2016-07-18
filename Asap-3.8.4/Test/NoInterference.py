# Test that various observers and analysis tools do not interfere with
# the energy calculator e.g. by messing up the neighbor list.

from asap3 import *
from asap3.md.verlet import VelocityVerlet
from asap3.md.langevin import Langevin
from ase.lattice.cubic import *
from asap3.analysis.localstructure import RestrictedCNA
from asap3.analysis.rdf import RadialDistributionFunction
from numpy import *
from asap3.testtools import ReportTest

def Compare(name, atoms, observer, ref, interval=4):
    dyn = VelocityVerlet(atoms, 5*units.fs)
    dyn.attach(observer, interval=interval)
    r = []
    for i in range(50):
        dyn.run(5)
        epot = atoms.get_potential_energy() / len(atoms)
        ekin = atoms.get_kinetic_energy() / len(atoms)
        r.append([epot, ekin, epot+ekin])
    r = array(r)
    diff = r - ref
    maxdiff = diff.max()
    print "Maximal difference is", maxdiff
    ReportTest(name, maxdiff, 0.0, 1e-6)
    
print_version(1)

print "Making initial system"
iniatoms = FaceCenteredCubic(directions=[[1,0,0],[0,1,0],[0,0,1]],
                             size=(10,5,5), symbol="Cu", pbc=(1,1,0))
iniatoms.set_calculator(EMT())
inidyn = Langevin(iniatoms, 5*units.fs, 450*units.kB, 0.05)
inidyn.run(100)
print "Temperature is now", iniatoms.get_kinetic_energy() / (1.5*units.kB*len(iniatoms)), "K"

print "Making reference simulation"
refatoms = Atoms(iniatoms)
refatoms.set_calculator(EMT())
dyn = VelocityVerlet(refatoms, 5*units.fs)
ref = []
for i in range(50):
    dyn.run(5)
    epot = refatoms.get_potential_energy() / len(refatoms)
    ekin = refatoms.get_kinetic_energy() / len(refatoms)
    ref.append([epot, ekin, epot+ekin])

ref = array(ref)

print "Testing RestrictedCNA"
atoms = Atoms(iniatoms)
atoms.set_calculator(EMT())
cna = RestrictedCNA(atoms)
Compare("RestrictedCNA", atoms, cna.analyze, ref)

print "Testing RadialDistributionFunction"
atoms = Atoms(iniatoms)
atoms.set_calculator(EMT())
rdf = RadialDistributionFunction(atoms, 4.0, 25)
Compare("RadialDistributionFunction", atoms, rdf.update, ref)

