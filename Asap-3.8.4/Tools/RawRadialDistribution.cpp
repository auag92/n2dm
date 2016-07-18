// Make Radial Distribution Functions (RDFs) grouped by elements and
// possibly by another integer tag, for example encoding the location.
// Only the counting is done here.  The normalization etc is done in Python.

// Copyright (C) 2008-2011 Jakob Schiotz and Center for Individual
// Nanoparticle Functionality, Department of Physics, Technical
// University of Denmark.  Email: schiotz@fysik.dtu.dk
//
// This file is part of Asap version 3.
//
// This program is free software: you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// version 3 as published by the Free Software Foundation.  Permission
// to use other versions of the GNU Lesser General Public License may
// granted by Jakob Schiotz or the head of department of the
// Department of Physics, Technical University of Denmark, as
// described in section 14 of the GNU General Public License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// and the GNU Lesser Public License along with this program.  If not,
// see <http://www.gnu.org/licenses/>.


#include "AsapPython.h"
#include "RawRadialDistribution.h"
#include "Exception.h"
#include "Atoms.h"
#include "NeighborCellLocator.h"
#include "Potential.h"
#include "GetNeighborList.h"
#include <math.h>
#include <iostream>
using std::cerr;
using std::endl;

namespace ASAPSPACE {

void RawRadialDistribution(PyObject *pyatoms, int nGroups, const int *group,
			   double rCut, int nBins,
			   RDFtype &result, RDFcountType &counts,
			   long *globalresult)
{
#if 0
  // Check that the group array is meaningful  XXX THIS LOOKS WRONG! Loop over atoms!
  for (int i = 0; i < nGroups; i++)
    if (group[i] < 0 || group[i] >= nGroups)
      throw AsapError("Invalid group for atom ") << i << ": " << group[i];
#endif
  
  assert(result.size() == nGroups);
  assert(counts.size() == nGroups);

#if 0
  // Check that the output data is zero.
  typedef map< pair<int,int>, long* >::const_iterator iter1;
  typedef map<int,long>::const_iterator iter2;
  for (int i = 0; i < nGroups; i++)
    {
      for (iter1 j = result[i].begin(); j != result[i].end(); ++j)
	{
	  cerr << "Checking group " << i << " " << j->first.first << ", "
	       << j->first.second << endl;
	  for (int k = 0; k < nBins; k++)
	    assert(j->second[k] == 0);
	}
      for (iter2 j = counts[i].begin(); j != counts[i].end(); ++j)
	assert(j->second == 0);
    }    
#endif
  
  // Get a neighbor list that can access the atoms.
  PyObject *py_nblist = NULL;
  Atoms *atoms = NULL;
  GetNbList_FromAtoms(pyatoms, rCut, &py_nblist, &atoms);
  NeighborLocator *nblist = ((PyAsap_NeighborLocatorObject*)py_nblist)->cobj;
  
  // Get some info about the atoms
  const asap_z_int *z = atoms->GetAtomicNumbers();
  int nAtoms = atoms->GetNumberOfAtoms();

#if 0
  // Do we have ghost atoms?
  int nAtomsInclGhosts = nAtoms;
  GhostAtoms *ghostAtoms = dynamic_cast<GhostAtoms *>(atoms);
  if (ghostAtoms)
    nAtomsInclGhosts += ghostAtoms->GetNumberOfGhosts();
#endif
  

  // Now make the RDFs
  double deltaR = rCut/nBins;  // Take care later if r == rCut!
  int maxnb = nblist->MaxNeighborListLength() + 1;
  vector<int> other(maxnb);
  vector<double> sqdist(maxnb);
  vector<Vec> rnb(maxnb);
  pair<int,int> forward;
  pair<int,int> reverse;
  for (int atom = 0; atom < nAtoms; atom++)
    {
      counts[group[atom]][z[atom]] += 1;        // Count the atom
      forward.first = reverse.second = z[atom];
      int maxnb2 = maxnb;  // Will be changed by the next line
      int nnb = nblist->GetNeighbors(atom, &other[0], &rnb[0], &sqdist[0],
				     maxnb2, rCut);
      assert(nnb < maxnb);
      // Loop over neighbors
      for (int i = 0; i < nnb; i++)
	{
	  int nb = other[i];
	  double r = sqrt(sqdist[i]);
	  int bin = int(r/deltaR);
	  if (bin >= nBins)
	    {
	      // If r == rCut, bin==nBins.  Discard this case.
	      assert(bin == nBins);  
	      break;
	    }
	  forward.second = reverse.first = z[nb];
	  result[group[atom]][forward][bin] += 1;
	  // Check if the neighbor is a ghost
	  if (nb < nAtoms)
	    {
	      // No
	      result[group[nb]][reverse][bin] += 1;
	      globalresult[bin] += 2;
	    }
	  else
	    {
	      // Yes
	      globalresult[bin] += 1;
	    }
	}
    }

  atoms->End();
  CHECKREF(py_nblist);
  Py_DECREF(py_nblist);
  AsapAtoms_DECREF(atoms);
}

} // end namespace
