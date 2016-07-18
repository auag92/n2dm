/*
 * NeighborList2013.h
 *
 *  Created on: May 5, 2011
 *      Author: s072162
 */

#ifndef _NEIGHBORLIST2011_H
#define _NEIGHBORLIST2011_H

#include "NeighborList.h"
#include "TinyMatrix.h"

namespace ASAPSPACE {

PyAsap_NeighborLocatorObject *PyAsap_NewNeighborList2013(Atoms *atoms, double rCut,
                                                       double driftfactor,
                                                       const TinyMatrix<double> &rcut2);

class NeighborList2013 : public NeighborList
{
protected:
  /// Generate a modified neighbor list for atoms a with cutoff rCut, chosen from a TinyMatrix RCUT.

  /* New Constructor */
  /// The neighbor list will contain all neighbors within the distance
  /// rCut, based on the individual elements involved in the interactions.
  /// The neighborlist can be reused until an atom has moved
  /// more than rCut*driftfactor.
  NeighborList2013(Atoms *a, double rCut, double driftfactor,
      const TinyMatrix<double> &rcut2);

  friend PyAsap_NeighborLocatorObject *PyAsap_NewNeighborList2013(Atoms *atoms,
      double rCut,
      double driftfactor,
      const TinyMatrix<double> &rcut2);

  /* Changed functions */
  /// Return value: The number of neighbors.
  virtual int GetNeighbors(int n, int *neighbors, Vec *diffs, double *diffs2,
      int& size, double r = -1.0) const;

  /// Get information about the neighbors of atom n ("half" neighbor list)
  ///
  /// This version of GetNeighbors only returns the numbers of the neighbors.
  /// It is intended for the Python interface.
  virtual void GetNeighbors(int n, vector<int> &neighbors) const;

  /* New data */
  TinyMatrix<double> rcut2_byz;
};

} // end namespace

#endif /* NEIGHBORLIST2011_H_ */
