import numpy as np

def CohesiveEnergy(atoms):
    e_full = atoms.get_potential_energy() / len(atoms)

    atoms2 = atoms[[0,]]
    atoms2.set_calculator(atoms.get_calculator())
    e_reduced = atoms2.get_potential_energy()

    return e_reduced - e_full

def VacancyEnergy(atoms):
    e_full = atoms.get_potential_energy()

    atoms2 = atoms.copy()
    atoms2.set_calculator(atoms.get_calculator())
    atoms2.pop()
    e_reduced = atoms.get_potential_energy()

    return e_reduced - e_full * len(atoms) / (len(atoms) + 1)

def HeatOfFormation(atoms):
    return None

def HeatOfSolution(atoms):
    return None
