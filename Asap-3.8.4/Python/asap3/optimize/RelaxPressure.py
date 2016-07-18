"Relax the pressure in a simulation."

from ASE.Dynamics.Dynamics import Dynamics

class RelaxPressure(Dynamics):
    """This class implements Behrendsen's method to relax the pressure."""

    def __init__(self, atoms, p=0.0, coeff=0.1):
        Dynamics.__init__(self, atoms)
        self.p = p
        self.coeff = coeff

    def _Step(self):
        stress = self.atoms.GetStress()
        ucell = self.atoms.GetUnitCell()
        p = - sum(stress[0:3]) / 3.0
        delV = self.coeff * (p - self.p)
        ucell *= (1.0 + delV/3.0)
        self.atoms.SetUnitCell(ucell)
        
