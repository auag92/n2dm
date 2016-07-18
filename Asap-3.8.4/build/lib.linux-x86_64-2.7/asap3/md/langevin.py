"""Langevin dynamics class."""

import numpy as np
from numpy.random import standard_normal, randint
import asap3
import ase
from asap3.md.md import MolecularDynamics, ParallelMolDynMixin
from ase.md.langevin import Langevin as Langevin_ASE
import asap3.constraints
import sys

class Langevin_New(object, Langevin_ASE, ParallelMolDynMixin):
    def __init__(self, atoms, timestep, temperature, friction, fixcm=True,
                 trajectory=None, logfile=None, loginterval=1):
        assert(getattr(Langevin_ASE, "_lgv_version", 1) >= 2)
        ParallelMolDynMixin.__init__(self, "Langevin", atoms)
        Langevin_ASE.__init__(self, atoms, timestep, temperature, friction, fixcm,
                              trajectory, logfile, loginterval, communicator=None)
        
    sdpos = property(lambda s: s.get("sdpos"), lambda s, x: s.set("sdpos", x))
    sdmom = property(lambda s: s.get("sdmom"), lambda s, x: s.set("sdmom", x))
    c1 = property(lambda s: s.get("c1"), lambda s, x: s.set("c1", x))
    c2 = property(lambda s: s.get("c2"), lambda s, x: s.set("c2", x))
    act0 = property(lambda s: s.get("act0"), lambda s, x: s.set("act0", x))
    c3 = property(lambda s: s.get("c3"), lambda s, x: s.set("c3", x))
    c4 = property(lambda s: s.get("c4"), lambda s, x: s.set("c4", x))
    pmcor = property(lambda s: s.get("pmcor"), lambda s, x: s.set("pmcor", x))
    cnst = property(lambda s: s.get("cnst"), lambda s, x: s.set("cnst", x))
    masses = property(lambda s: s.get("masses"), lambda s, x: s.set("masses", x))

class Langevin_Fast(object, Langevin_ASE, ParallelMolDynMixin):
    def __init__(self, atoms, timestep, temperature, friction, fixcm=True,
                 trajectory=None, logfile=None, loginterval=1):
        assert(getattr(Langevin_ASE, "_lgv_version", 1) >= 2)
        ParallelMolDynMixin.__init__(self, "Langevin", atoms)
        self._uselocaldata = False # Need to store on atoms for serial simul too.
        self.calculator = atoms.get_calculator()            
        self.asap_md = asap3._asap.Langevin(atoms, self.calculator, timestep,
                                            self.prefix+"sdpos", self.prefix+"sdmom",
                                            self.prefix+"c1", self.prefix+"c2",
                                            fixcm, randint(1 << 30))
        if not atoms.has('momenta'):
            atoms.set_momenta(np.zeros((len(atoms), 3), float))
        if atoms.constraints:
            assert len(atoms.constraints) == 1
            constraint = atoms.constraints[0]
            assert isinstance(constraint, asap3.constraints.FixAtoms)
            constraint.prepare_for_asap(atoms)            
            # Make all constants arrays by making friction an array
            friction = friction * np.zeros(len(atoms))
        Langevin_ASE.__init__(self, atoms, timestep, temperature, friction, fixcm,
                              trajectory, logfile, loginterval, communicator=None)
        
    def updatevars(self):
        Langevin_ASE.updatevars(self)
        if len(self.atoms.constraints) == 1:
            # Process the FixAtoms constraint
            constr = self.atoms.constraints[0].index
            self.sdpos[constr] = 0.0
            self.sdmom[constr] = 0.0
            self.c1[constr] = 0.0
            self.c2[constr] = 0.0
            self.c3[constr] = 0.0
            self.c4[constr] = 0.0
            self.act0[constr] = 0.0
        if self._localfrict:
            self.asap_md.set_vector_constants(self.prefix+"act0", self.prefix+"c3", 
                                              self.prefix+"c4", self.prefix+"pmcor", 
                                              self.prefix+"cnst")
        else:
            self.asap_md.set_scalar_constants(self.act0, self.c3, self.c4,
                                              self.pmcor, self.cnst)
            
    def run(self, steps):
        assert(self.calculator is self.atoms.get_calculator())
        self.asap_md.run(steps, self.observers, self)
        
    def get_random(self, gaussian):
        return self.asap_md.get_random(gaussian)
    
    # Properties are not inherited, need to repeat them
    sdpos = property(lambda s: s.get("sdpos"), lambda s, x: s.set("sdpos", x))
    sdmom = property(lambda s: s.get("sdmom"), lambda s, x: s.set("sdmom", x))
    c1 = property(lambda s: s.get("c1"), lambda s, x: s.set("c1", x))
    c2 = property(lambda s: s.get("c2"), lambda s, x: s.set("c2", x))
    act0 = property(lambda s: s.get("act0"), lambda s, x: s.set("act0", x))
    c3 = property(lambda s: s.get("c3"), lambda s, x: s.set("c3", x))
    c4 = property(lambda s: s.get("c4"), lambda s, x: s.set("c4", x))
    pmcor = property(lambda s: s.get("pmcor"), lambda s, x: s.set("pmcor", x))
    cnst = property(lambda s: s.get("cnst"), lambda s, x: s.set("cnst", x))
    masses = property(lambda s: s.get("masses"), lambda s, x: s.set("masses", x))

        
            
class Langevin_Old(MolecularDynamics):
    """Langevin (constant N, V, T) molecular dynamics.

    Usage: Langevin(atoms, dt, temperature, friction)

    atoms
        The list of atoms.
        
    dt
        The time step.

    temperature
        The desired temperature, in energy units.

    friction
        A friction coefficient, typically 1e-4 to 1e-2.

    fixcm
        If True, the position and momentum of the center of mass is
        kept unperturbed.  Default: True.

    The temperature and friction are normally scalars, but in principle one
    quantity per atom could be specified by giving an array.

    This dynamics accesses the atoms using Cartesian coordinates."""
    
    def __init__(self, atoms, timestep, temperature, friction, fixcm=True,
                 trajectory=None, logfile=None, loginterval=1):
        MolecularDynamics.__init__(self, atoms, timestep, trajectory,
                                   logfile, loginterval)
        self.temp = temperature
        self.frict = friction
        self.fixcm = fixcm  # will the center of mass be held fixed?

        # We need to store some data on the atoms, so they can migrate
        # along with the atoms in parallel asap simulations.  Since
        # there could be more than one Langevin object acting on the
        # same atoms, a unique prefix to the variable names is needed.
        i = 0
        while True:
            self.prefix = "langevin%d_" % (i,)
            if not self.atoms.has(self.prefix+"sdmom"):
                break  # Found unique prefix
            i += 1
        self.updatevars()
        
    def set_temperature(self, temperature):
        self.set("temp", temperature)
        self.updatevars()

    def set_friction(self, friction):
        self.set("frict", friction)
        self.updatevars()

    def set_timestep(self, timestep):
        self.dt = timestep
        self.updatevars()

    def updatevars(self):
        dt = self.dt
        # If the friction is an array some other constants must be arrays too.
        frict = self.get("frict")
        temp = self.get("temp")
        self._localfrict = hasattr(frict, 'shape')
        lt = frict * dt
        masses = self.atoms.get_masses()
        sdpos = dt * np.sqrt(temp / masses * (2.0/3.0 - 0.5 * lt) * lt)
        sdpos.shape = (-1, 1)
        sdmom = np.sqrt(temp * masses * 2.0 * (1.0 - lt) * lt)
        sdmom.shape = (-1, 1)
        pmcor = np.sqrt(3.0)/2.0 * (1.0 - 0.125 * lt)
        cnst = np.sqrt((1.0 - pmcor) * (1.0 + pmcor))

        act0 = 1.0 - lt + 0.5 * lt * lt
        act1 = (1.0 - 0.5 * lt + (1.0/6.0) * lt * lt)
        act2 = 0.5 - (1.0/6.0) * lt + (1.0/24.0) * lt * lt
        c1 = act1 * dt / masses
        c1.shape = (-1, 1)
        c2 = act2 * dt * dt / masses
        c2.shape = (-1, 1)
        c3 = (act1 - act2) * dt
        c4 = act2 * dt
        del act1, act2
        if self._localfrict:
            # If the friction is an array, so are these
            act0.shape = (-1, 1)
            c3.shape = (-1, 1)
            c4.shape = (-1, 1)
            pmcor.shape = (-1, 1)
            cnst.shape = (-1, 1)
        self.set("sdpos", sdpos)
        self.set("sdmom", sdmom)
        self.set("c1", c1)
        self.set("c2", c2)
        self.set("act0", act0)
        self.set("c3", c3)
        self.set("c4", c4)
        self.set("pmcor", pmcor)
        self.set("cnst", cnst)
        self.natoms = self.atoms.get_number_of_atoms()

    def step(self, f):
        atoms = self.atoms
        p = self.atoms.get_momenta()

        random1 = standard_normal(size=(len(atoms), 3))
        random2 = standard_normal(size=(len(atoms), 3))
        
        rrnd = self.get("sdpos") * random1
        prnd = (self.get("sdmom") * self.get("pmcor") * random1 +
                self.get("sdmom") * self.get("cnst") * random2)

        if self.fixcm:
            rrnd = rrnd - np.sum(rrnd, 0) / len(atoms)
            prnd = prnd - np.sum(prnd, 0) / len(atoms)
            factor = np.sqrt(self.natoms / (self.natoms - 1.0)) 
            rrnd *= factor
            prnd *= factor

        atoms.set_positions(atoms.get_positions() +
                            self.get("c1") * p +
                            self.get("c2") * f + rrnd)
        p *= self.get("act0")
        p += self.get("c3") * f + prnd
        atoms.set_momenta(p)
                      
        f = atoms.get_forces()
        atoms.set_momenta(atoms.get_momenta() + self.get("c4") * f)
        return f

def Langevin(atoms, *args, **kwargs):
    if (isinstance(atoms, ase.Atoms) 
        and asap3.constraints.check_asap_constraints(atoms)
        and getattr(Langevin_ASE, "_lgv_version", 1) >= 2
       ):
        sys.stderr.write("Using Asap-optimized C++-Langevin algorithm\n")
        return Langevin_Fast(atoms, *args, **kwargs)
    elif getattr(Langevin_New, "_lgv_version", 1) >= 2:
        sys.stderr.write("Using ASE-based Langevin algorithm\n")
        return Langevin_New(atoms, *args, **kwargs)
    else:
        sys.stderr.write("Using reimplemented Python-Langevin algorithm\n")
        return Langevin_Old(atoms, *args, **kwargs)
    
