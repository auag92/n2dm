"""A radial distribution function observer."""

import numpy
import asap3
from asap3.Internal.Subject import Subject
from asap3.Internal.ListOfElements import ListOfElements
import cPickle
import new
import sys

class RadialDistributionFunction(Subject):
    def __init__(self, atoms, rMax, nBins, groups=None, interval=1, average=1,
                 autoclear=False, verbose=False):
        """Create a RadialDistributionFunction observer.

        Arguments:
        
        atoms: The atoms being observed.
        
        rMax: The maximum distance in the RDF.
        
        nBins: The number of bins in the histograms.
        
        groups (optional): A non-negative integer per atom, used to
        classify atoms into groups, the RDF is calculated for each
        group.

        interval (optional, default=1): How often should the RDF be
        calculated.  DEPRECATED: Use interval argument when attaching
        to dynamics instead.

        average (optional, default=1): How many times should the RDF
        be calculated and averaged before any observer of this object
        is notified and/or the RDF is saved to a file.

        autoclear (optional, default=False): Should the RDF be cleared
        after the RDF has been processed by observers and/or written
        to a file?  The default is to continue to accumulate data.
        """
        Subject.__init__(self)
        self.atoms = atoms
        self.natoms = len(atoms)
        self.globalrdf = None
        self.rMax = rMax * 1.0  # Guard against integer divisions.
        self.nBins = nBins
        self.dr = self.rMax / nBins
        self.interval = interval
        self.average = average
        self.autoclear = autoclear
        self.verbose = verbose
        self.autosave = False
        self.n1 = 0  # Counter associated with interval
        self.n2 = 0  # Counter associated with average
        self.clearnext = False
        self.countRDF = 0
        self.listOfElements = ListOfElements(atoms)
        
        if groups is None:
            self.groups = numpy.zeros(len(atoms), numpy.int32)
            self.ngroups = 1
        else:
            if groups.shape != (len(atoms),):
                raise ValueError, "groups must be an integer per atom"
            if min(groups) < 0 or max(groups) >= len(atoms):
                raise ValueError, "groups array is unreasonable"
            self.groups = groups.astype(numpy.int32)
            self.ngroups = max(self.groups) + 1

    def update(self, atoms=None):
        """Calculate the RDF of the atoms.

        Make an RDF calculation now (or if interval=n was specified,
        every n time this method is called.

        If average = 1, then any observers of this object is notified.

        If average > 1, the calculated RDFs are accumulated, and once
        a sufficient number of RDFs have been accumulated, the
        observers are notified and the averaged RDFs are ready.
        """
        if atoms is not None:
            self.atoms = atoms

        self.n1 += 1
        if self.n1 >= self.interval:
            if self.clearnext:
                self.clear()
            self.n1 = 0
            self.do_rdf()
            self.n2 += 1
            if self.n2 >= self.average:
                self.clearnext = self.autoclear
                self.call_observers()  # I have data ready
                if self.autosave:
                    self.save()
                self.n2 = 0

    def do_rdf(self):
        """Do the actual RDF calculation.  Do not call directly."""
        if self.verbose:
            print >>sys.stderr, "Calculating RDF"
        assert self.natoms == len(self.atoms)
        grdf, rdfs, cnts = asap3._asap.RawRDF(self.atoms, self.rMax, self.nBins,
                                              self.groups, self.ngroups,
                                              self.listOfElements)

        ###
        print 'Natoms:', self.natoms
        print 'Global RDF count:', grdf.sum()
        for g in range(len(rdfs)):
            print 'Group:', g
            for e in rdfs[g].keys():
                print '  RDF count %s: %i' % (e, rdfs[g][e].sum())
            for e in cnts[g].keys():
                print '  Count (%i): %i' % (e, cnts[g][e])
        ###

        # Sanity check on counts
        c = 0
        for cnt in cnts: # Loop over groups
            for t in cnt.keys(): # Loop over elements
                c += cnt[t]
        assert c == self.natoms
        # Now store the data
        if self.globalrdf is None:
            self.globalrdf = grdf
            self.rdfs = rdfs
            self.atomcounts = cnts
            # Ensure that all elements are represented
            for i in range(len(rdfs)): # Loop over groups
                for e in self.listOfElements: # Loop over elements
                    if not self.atomcounts[i].has_key(e):
                        self.atomcounts[i][e] = 0
        else:
            self.globalrdf += grdf
            for i in range(len(rdfs)): # Loop over groups
                for t in rdfs[i].keys(): # Loop over pairs
                    self.rdfs[i][t] += rdfs[i][t]
                for t in cnts[i].keys(): # Loop over elements
                    self.atomcounts[i][t] += cnts[i][t]
        self.countRDF += 1
        self.volume = self.atoms.get_volume()
        if self.verbose:
            print >>sys.stderr, "RDF done!"

    def clear(self):
        """Clear the accumulated RDFs.  Called automatically."""
        if self.verbose:
            print >>sys.stderr, "Clearing RDF"
        self.globalrdf = self.rdfs = self.atomcounts = None
        self.countRDF = 0
        
    def get_rdf(self, groups=None, elements=None, normalize='volume'):
        """Get an RDF.

        Arguments (both optional):
        
        groups: Only get the RDF for atoms in these groups. Can either be an
        integer or a list integer refering to the groups (default: all groups)
        
        elements: Get the partial RDF for these two elements.
        elements must be a tuple of two atomic numbers (a, b), the
        returned RDF tells how many b neighbors an a atom has.

        If get_rdf is called on a newly created
        RadialDistributionFunction object, and this object has been
        created without specifying interval and average parameters, it
        is assumed that the user is not using the object as an
        observer, but wants an immediate calculation of the RDF.  In
        that case, calling get_rdf triggers the calculation of the
        RDFs.  In all other cases previously stored RDFs are returned.
        """
        if self.globalrdf is None and self.interval == 1 and self.average == 1:
            self.update()
            
        if groups is None and elements is None:
            # Return global RDF
            return self.normalize(self.globalrdf,
                                  self.countRDF * self.natoms,
                                  normalize)

        # Either a group or a pair of elements have been specified.
        # Sum over selected groups
        if groups is None:
            groups = range(len(self.rdfs))
        else:
            if not isinstance(groups, list):
                groups = [groups]

        rdfs = self.rdfs[groups[0]].copy()  # Shallow copy
        atomcounts = self.atomcounts[groups[0]].copy()  # Shallow copy
        for i in groups[1:]: # Loop over groups
            if i >= len(self.rdfs):
                continue
            for t in self.rdfs[i].keys(): # Loop over pairs
                rdfs[t] += self.rdfs[i][t]
            for t in self.atomcounts[i].keys(): # Loop over elements
                atomcounts[t] += self.atomcounts[i][t]

        # Sum over selected element pairs
        elementpairs = elements
        if elementpairs is None:
            elements = self.listOfElements
            elementpairs = rdfs.keys()
        else:
            if isinstance(elementpairs, tuple):
                elements = [elementpairs[0]]
                elementpairs = [elementpairs]
            elif isinstance(elementpairs, int):
                elements = [elementpairs]
                elementpairs = []
                for e1 in self.listOfElements:
                    elementpairs.append((elements[0], e1))
            else:
                raise ValueError('Elements must be either an interger or ' +
                                 'a tuple of two integers.')

        rdf = None
        for pair in elementpairs:
            if rdf is None:
                rdf = numpy.array(rdfs[pair])
            else:
                rdf += rdfs[pair]

        # The atomcounts should always be summed over elements!
        atomcount = None
        for t in atomcounts.keys():
            if atomcount is None:
                atomcount = atomcounts[t]
            else:
                atomcount += atomcounts[t]

        ###
        atomcount = None
        for e in elements:
            if atomcount is None:
                atomcount = atomcounts[e]
            else:
                atomcount += atomcounts[e]
        print 'Number of selected atoms:', atomcount
        ###
        return self.normalize(rdf, atomcount, normalize)

    def normalize(self, rdf, ncount, type):
        """Normalize the raw RDF returned by the C++ module."""
        if type == 'volume':
            factor = (4 * numpy.pi / 3.0) * ncount * self.natoms / self.volume
            r = numpy.arange(self.nBins) * self.dr
            r3low = r * r * r
            r += self.dr
            r3high = r * r * r
            normalization = factor * (r3high - r3low)
        elif type == 'atoms':
            normalization = 1.0 * ncount
        else:
            normalization = 1.0
        return rdf / normalization
    
    # Saving and loading RDF's
    def output_file(self, prefix):
        """Give the file prefix for saving RDFs, and turn on saving."""
        self.autosave = True
        self.savefilecounter = 0
        if "%" in prefix:
            # Assume the user knows what (s)he is doing
            self.filenames = prefix
        else:
            self.filenames = prefix + "%04d.rdf"

    def save(self, filename=None):
        if self.verbose:
            print >>sys.stderr, "Saving RDF"
        if filename == None:
            filename = (self.filenames % (self.savefilecounter,))
            self.savefilecounter += 1
        data = {"globalrdf": self.globalrdf,
                "rdfs" : self.rdfs,
                "atomcounts": self.atomcounts,
                "countRDF": self.countRDF,
                "rMax": self.rMax,
                "dr": self.dr,
                "nBins": self.nBins,
                "natoms": self.natoms,
                "volume": self.volume}
        f = open(filename, "w")
        cPickle.dump(data, f, 2) # Protocol 2 (efficient bin from Python 2.3)
        f.close()
        
    # Not supported in python 2.3:  @classmethod
    #@classmethod
    def load(cls, filename):
        f = open(filename)
        data = cPickle.load(f)
        f.close()
        obj = new.instance(cls, data)
        return obj
    load=classmethod(load)
    
