from ASE.Filters.HomogeneousStrain import HomogeneousStrain
_HomogeneousStrain = HomogeneousStrain

__docformat__ = "restructuredtext en"

class HomogeneousStrain(_HomogeneousStrain):
    """Modify the supercell while keeping the scaled positions fixed.

    Presents the strain of the supercell as the generalized positions,
    and the global stress tensor (times the volume) as the generalized
    force.

    This filter can be used to relax the unit cell until the stress is zero.

    The stress and strain are presented as 6-vectors, the order of the
    components follow the standard engingeering practice: xx, yy, zz,
    yz, xz, xy.  
    
    This special version supports also supports parallel ASAP atoms.
    """

    def __init__(self, atoms, mask = None, appliedStress=None):
        """Create a filter applying a homogeneous strain to a list of atoms.

        The first argument, atoms, is the list of atoms.

        The optional second argument, mask, is a list of six booleans,
        indicating which of the six independent components of the
        strain that are allowed to become non-zero.  It defaults to
        [1,1,1,1,1,1].  Please note that in any case there are 6
        components of the generalized force, position and momenta
        (this may be considered a bug and may change in the future).
        
        """
        _HomogeneousStrain.__init__(self, atoms, mask, appliedStress)
        # Now try to identify a parallel simulation
        if getattr(atoms, "parallel ", 0):
            import Scientific.MPI
            self.comm = Scientific.MPI.world
            # Be absolutely sure to keep thing exactly identical everywhere
            self.comm.broadcast(self.origbasis, 0)
        else:
            self.comm = None

    def SetGeneralizedPositions(self, new):
        if self.comm:
            self.comm.broadcast(new, 0)
        _HomogeneousStrain.SetGeneralizedPositions(self, new)

    def SetGeneralizedMomenta(self, mom):
        if self.comm:
            self.comm.broadcast(mom, 0)
        _HomogeneousStrain.SetGeneralizedMomenta(self, mom)

            
