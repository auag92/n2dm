/*
 * NeighborList2013.cpp
 *
 *  Created on: May 5, 2011
 *      Author: s072162
 */

#include "NeighborList2013.h"
#include "Debug.h"
#include "Atoms.h"
#include <iostream>

NeighborList2013::NeighborList2013(Atoms *a, double rCut, double driftfactor,
    const TinyMatrix<double> &rcut2_NB)
: NeighborList(a, rCut, driftfactor)
{
  rcut2_byz.CopyFrom(rcut2_NB);
};

//// retrieve neighbors based on rcut between the involved atoms (TinyMatrix version)
int NeighborList2013::GetNeighbors(int a1, int *neighbors, Vec *diffs,
    double *diffs2, int& size, double r) const
{
  if (r > 0.0)
    return NeighborList::GetNeighbors(a1, neighbors, diffs, diffs2, size, r);

  if (invalid)
    {
      DEBUGPRINT;
      throw AsapError("NeighborList has been invalidated, possibly by another NeighborList using the same atoms.");
    }

  if (size < nbList[a1].size())
    {
      DEBUGPRINT;
      throw AsapError("NeighborList::GetNeighbors: list overflow.");
    }

  const vector<Vec> &positions = cells->GetWrappedPositions();
  const asap_z_int *z = atoms->GetAtomicNumbers();

  // Need to use GET_CELL instead of GetCell as the atoms are not open
  // when called from the Python interface.
  const Vec *superCell = atoms->GET_CELL();
  Vec pos1 = positions[a1];
  int nNeighbors = 0;
  asap_z_int a1Element = z[a1];

  typedef vector< pair<int,translationsidx_t> >::const_iterator iterat;
  iterat terminate = nbList[a1].end();
  for (iterat a2 = nbList[a1].begin(); a2 < terminate; ++a2)
    {
	  // Check to see if the potential neighboring atom should be
	  // on the neighbor list of atom a1

      assert(a2->second >= 0 && a2->second < translationTable.size());
      assert(a2->first >= 0 && a2->first < nAllAtoms);
      assert(nNeighbors < size);
      const IVec &celltranslation = translationTable[a2->second];
      diffs[nNeighbors] = positions[a2->first] - pos1
	- celltranslation[0] * superCell[0]
	- celltranslation[1] * superCell[1]
	- celltranslation[2] * superCell[2];
      double d2 = Length2(diffs[nNeighbors]);
      asap_z_int a2Element = z[a2->first];

      /* Assuming that the problem about accessing the element
      * type has been solved a1/2Element should be substituted
      * with the right command! */
      if (d2 < rcut2_byz[a1Element][a2Element])
	{
	  diffs2[nNeighbors] = d2;
	  neighbors[nNeighbors] = a2->first;
	  nNeighbors++;
	}
    }
  size -= nNeighbors;
  assert(size >= 0);
  return nNeighbors;
}

/// PYTHON VERSION!
void NeighborList2013::GetNeighbors(int a1, vector<int> &neighbors) const
{
  if (invalid)
    throw AsapError("NeighborList has been invalidated, possibly by another NeighborList using the same atoms.");

  neighbors.clear();
  const vector<Vec> &positions = cells->GetWrappedPositions();
  const asap_z_int *z = atoms->GetAtomicNumbers();

  const Vec *superCell = atoms->GET_CELL();  // atoms may not be open
  Vec pos1 = positions[a1];
  asap_z_int a1Element = z[a1];

  typedef vector< pair<int,translationsidx_t> >::const_iterator iterat;
  iterat terminate = nbList[a1].end();
  for (iterat a2 = nbList[a1].begin(); a2 < terminate; ++a2)
    {
      const IVec &celltranslation = translationTable[a2->second];
      Vec diff = positions[a2->first] - pos1
	- celltranslation[0] * superCell[0]
	- celltranslation[1] * superCell[1]
	- celltranslation[2] * superCell[2];
      /* Assuming that the problem about accessing the element
       * type has been solved a1/2Element should be substituted
       * with the right command! */
      double d2 = Length2(diff);
      asap_z_int a2Element = z[a2->first];
      if (d2 < rcut2_byz[a1Element][a2Element])
	    neighbors.push_back(a2->first);
    }
}



