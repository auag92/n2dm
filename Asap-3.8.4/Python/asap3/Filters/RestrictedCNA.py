"""CNA (common neighbor analysis) Filter looking for FCC and HCP structure."""

__docformat__ = "restructuredtext en"

import ASE.ChemicalElements
import ASE.Utilities.Subject
import Numeric
import Asap

class RestrictedCNA(ASE.Utilities.Subject.Subject):
    """CNA (common neighbor analysis) Filter looking for FCC and HCP structure.

    Restricted CNA analysis, looking for FCC and HCP crystal
    structures.  The crystal structure is reported as an integer per
    atom, where 0 signifies FCC, 1 signifies HCP and 2 signifies any
    other structure.  These numbers are either placed in the Tags of
    the list of atoms, or made available through methods of this
    filter, depending on the mode string (see below).

    The RestrictedCNA filter also acts like an observer, and observe
    the dynamics moving the atoms (or a `TimeAveragedPositions` filter).
    when the dynamics notifies its observes, RestrictedCNA does its
    calculation (or every Nth time, if interval is set to N>1).  The
    RestrictedCNA filter also acts like a subject, notifying any
    observers every time it has done a calculation.  These
    subject-observer patterns can be used to collaborate with e.g. a
    Trajectory, see the example below.

    The filter is created with the following arguments:

        atoms: The list of atoms being analysed.
        
        cutoff (default None): The cutoff for the CNA analysis.
            Should be between the first and second neighbor distances.
            If unspecified, an attempt is made to extract it from the
            ChemicalElements database.  If a chemical element is
            given, a reasonable value for that element is used.

        interval (default 1): How often the analysis should be done.

        mode (default 'modify'): Specifies how the result is reported.
            'filter' means that it will only be available through this
            filters GetTags and GetCNA methods.  'special' means that
            it will only be available though the special GetCNA
            method.  'modify' (the default) means that in addition,
            the Tags of the original atoms are modified.

    The Common Neighbor Analysis can be confused by thermal
    vibrations, and it can therefore be advantageous to combine this
    filter with the `TimeAveragedPositions` filter.
    
    Example:

    A simulation with output going to a trajectory could be done like this

    dyn = VelocityVerlet(atoms, timestep)
    trajectory = NetCDFTrajectory(outfile, atoms, interval=1000)
    dyn.Attach(trajectory)
    dyn.Run(nsteps)

    With CNA analysis, this could be done like this (doing a CNA every
    1000 steps, using the positions averaged over the last 100 steps).

    dyn = VelocityVerlet(atoms, timestep)
    avg = TimeAveragedPositions(atoms, 100, 1000)
    dyn.Attach(avg)
    cna = RestrictedCNA(atoms)
    avg.Attach(cna)
    trajectory = NetCDFTrajectory(outfile, atoms, interval=1)
    cna.Attach(trajectory)
    """
    def __init__(self, atoms, cutoff=None, interval=1, mode='modify',
                 verbose=0):
        """Create the RestrictedCNA filter.

        atoms: The list of atoms being analysed.
        
        cutoff (default None): The cutoff for the CNA analysis.
            Should be between the first and second neighbor distances.
            If unspecified, an attempt is made to extract it from the
            ChemicalElements database.  If a chemical element is
            given, a reasonable value for that element is used.

        interval (default 1): How often the analysis should be done.

        mode (default 'modify'): Specifies how the result is reported.
            'filter' means that it will only be available through this
            filters GetTags and GetCNA methods.  'special' means that
            it will only be available though the special GetCNA
            method.  'modify' (the default) means that in addition,
            the Tags of the original atoms are modified.
        """
        ASE.Utilities.Subject.Subject.__init__(self)
        self.atoms = atoms
        self.interval = interval
        assert mode in ['modify', 'filter', 'special']
        self.mode = mode
        self.verbose = verbose
        self.count = 0
        # The cutoff may be a
        #   None - the cutoff is extracted from the atoms.
        #   number - the cutoff
        #   string - the element name
        #   ChemicalElements.Element object - the element
        if cutoff is None:
            z = atoms.GetAtomicNumbers()
            if min(z) != max(z):
                raise ValueError, "RestrictedCNA: cutoff must be given when multiple elements are present."
            cutoff = ASE.ChemicalElements.Element(min(z))
        elif isinstance(cutoff, str):
            cutoff = ASE.ChemicalElements.Element(cutoff)
        # Now cutoff should be a number or an element instance.
        try:
            struct = cutoff.crystal_structure
        except AttributeError:
            # stuct better be a number :-)
            pass
        else:
            if struct['symmetry'].lower() != 'fcc':
                raise RuntimeError, "Can only extract cutoffs from FCC elements"
            cutoff = struct['a'] * 0.5 * (1 + 1/Numeric.sqrt(2))
        if not (isinstance(cutoff, int) or isinstance(cutoff, float)):
            raise RuntimeError, "Cutoff is not a number, it is a "+str(type(cutoff))+": "+ str(cutoff)
        self.cutoff = cutoff

    def Update(self):
        "Do the CNA analysis"
        self.count += 1
        if self.count >= self.interval:
            if self.verbose:
                print "Doing CNA analysis."
            self.count = 0
            self.cna = self.do_cna()
            if self.mode == 'modify':
                self.atoms.SetTags(self.cna)
            self.Notify()

    def GetCNA(self):
        "Get the last CNA"
        return self.cna

    def GetTags(self):
        "Depending on the mode, get the last CNA or the tags of the atoms."
        if mode == 'special':
            return self.atoms.GetTags()
        else:
            return self.GetCNA()

    # do_cna can be replaced by another CNA function, if needed.
    def do_cna(self):
        "Internal method doing the actual work."
        if isinstance(self.atoms, Asap.CommonListOfAtoms):
            a = self.atoms
            if self.verbose:
                print "CNA: Called with Asap atoms."
        else:
            if self.verbose:
                print "CNA: Copying atoms to Asap atoms."
            a = Asap.ListOfAtoms(self.atoms)
        return Asap.asap.CommonNeighborAnalysis(a, self.cutoff)

    # All other calls are passed onto the atoms.
    # Is this a good idea, or should there be a "whitelist" ??    XXXXX

    def __getattr__(self, attr):
        return getattr(self.atoms, attr)

