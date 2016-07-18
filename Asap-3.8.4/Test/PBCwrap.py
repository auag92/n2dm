"""Check that energy is correct even after wrapping through periodic boundary conditions.
"""

from ase.lattice.cubic import FaceCenteredCubic
from asap3 import *
from asap3.testtools import *
import random

ref_atoms = FaceCenteredCubic(size=(7,7,7), symbol="Cu", pbc=(True, False, True))
ref_atoms.set_calculator(EMT())

ref_energy = ref_atoms.get_potential_energy()
ref_energies = ref_atoms.get_potential_energies()
ref_forces = ref_atoms.get_forces()

passes = 5
for ps in range(passes):
    print "Pass", ps, "of", passes

    atoms = ref_atoms.copy()
    atoms.set_calculator(EMT())
    nat = random.randint(0, len(atoms))
    assert nat < len(atoms)
    pos0 = atoms[nat].position
    cell = atoms.get_cell()
    for d in range(1,4):
        for dx in (-d, 0, d):
            #for dy in (-d, 0, d):
            for dy in (0,):
                for dz in (-d, 0 ,d):
                    deltar = dx * cell[0] + dy * cell[1] + dz * cell[2]
                    atoms[nat].position = pos0 + deltar
                    label = "(%2d, %2d, %2d)" % (dx, dy, dz)
                    ReportTest("Pot. energy   "+label,
                               atoms.get_potential_energy(), ref_energy, 1e-6)
                    de = (atoms.get_potential_energies() - ref_energies)
                    ReportTest("Pot. energies "+label, de.max(), 0.0, 1e-6)
                    df = (atoms.get_forces() - ref_forces)
                    ReportTest("Forces        "+label, df.max(), 0.0, 1e-6)
                
ReportTest.Summary()

