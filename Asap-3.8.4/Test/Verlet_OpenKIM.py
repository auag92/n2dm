from asap3 import *
from cPickle import *
from numpy import *
from asap3.testtools import ReportTest
from OpenKIM_modelname import openkimmodel

timeunit = 1.018047e-14             # Seconds
femtosecond = 1e-15 / timeunit      # Marginally different from units.fs

print_version(1)

if OpenKIMsupported:
    data = load(file("testVerlet.pickle"))
    init_pos = array(data["initial"])
    init_pos.shape = (-1,3)
    init_box = array(data["box"])
    init_box.shape = (3,3)
    atoms = Atoms(positions=init_pos, cell=init_box)
    atoms.set_atomic_numbers(47*ones((len(atoms),)))
    atoms.set_calculator(OpenKIMcalculator(openkimmodel))
    
    dyn = VelocityVerlet(atoms, 2 * femtosecond)
    dyn.attach(MDLogger(dyn, atoms, '-', peratom=True), interval=5)
    
    etot = None
    
    for i in range(10):
        dyn.run(20)
        epot = atoms.get_potential_energy() / len(atoms)
        ekin = atoms.get_kinetic_energy() / len(atoms)
        if etot is None:
            etot = epot + ekin
        else:
            ReportTest("Energy conservation:", epot + ekin, etot, 0.001)
            
    final_pos = array(data["final"])
    diff = max(abs(atoms.get_positions().flat - final_pos))
    print "Maximal deviation of positions:", diff
    ReportTest("Maximal deviation of positions", diff, 0, 1e-9)
    
    #diff = max(abs(atoms.get_stresses().flat - array(data["stress"])))
    #print "Maximal deviation of stresses:", diff
    #ReportTest("Maximal deviation of stresses", diff, 0, 1e-9)
    
    ReportTest.Summary()

else:
    print "OpenKIM support is not compiled into Asap."
    
