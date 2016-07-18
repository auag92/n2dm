// -*- C++ -*-
// NeighborList.h:  The main neighbor list object.
//
// Copyright (C) 2001-2011 Jakob Schiotz and Center for Individual
// Nanoparticle Functionality, Department of Physics, Technical
// University of Denmark.  Email: schiotz@fysik.dtu.dk
//
// This file is part of Asap version 3.
// Asap is released under the GNU Lesser Public License (LGPL) version 3.
// However, the parts of Asap distributed within the OpenKIM project
// (including this file) are also released under the Common Development
// and Distribution License (CDDL) version 1.0.
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


#ifndef NEIGHBORLIST2
#define NEIGHBORLIST2

#include "AsapPython.h"
#include "Asap.h"
#include "NeighborLocator.h"
#include "NeighborCellLocator.h"
#include "IVec.h"
#include "Templates.h"
#include <vector>
using std::vector;
#include <utility>
using std::pair;
#include <math.h>
#include <string.h>

namespace ASAPSPACE {

PyAsap_NeighborLocatorObject *PyAsap_NewNeighborList(Atoms *atoms, double rCut,
						     double driftfactor);

class NeighborList : public NeighborLocator
{
protected:
  /// Generate a neighbor list for atoms a with cutoff rCut.

  /// The neighbor list will contain all neighbors within the distance
  /// rCut.  The neighborlist can be reused until an atom has moved
  /// more than rCut*driftfactor. 
  NeighborList(Atoms *a, double rCut, double driftfactor);
  virtual ~NeighborList();

  friend ASAPSPACE::PyAsap_NeighborLocatorObject *PyAsap_NewNeighborList(Atoms *atoms,
     double rCut, double driftfactor);

  friend void PyAsap_Dealloc<PyAsap_NeighborLocatorObject>(PyObject *self);
  
public:
  /// Enable full neighbor lists by calling this just after the constructor
  void EnableFullNeighborLists();

  /// Check if full lists are enabled
  bool HasFullNeighborLists() const {return fulllists;}

  /// Get wrapped positions of all the atoms
  const vector<Vec> &GetWrappedPositions() const
  {return cells->GetWrappedPositions();}

  void GetWrappedPositions(vector<Vec> &wp) const
  {cells->GetWrappedPositions(wp);}

  /// Get scaled positions of all the atoms
  const vector<Vec> &GetScaledPositions() const
  {return cells->GetScaledPositions();} 
  
  /// Check the neighbor list.
  ///
  /// Check if the neighbor list can still be reused, return true if
  /// it should be updated.
  virtual bool CheckNeighborList();

  /// Update neighbor list
  virtual void UpdateNeighborList();
  
  /// Check if the neighbor list can still be reused, update if not.
  bool CheckAndUpdateNeighborList();

  /// Check if the neighbor list can still be reused, update if not.
  ///
  /// This version is used when called from Python
  virtual bool CheckAndUpdateNeighborList(PyObject *atoms);

  /// Get information about the neighbors of atom n ("half" neighbor list)
  ///
  /// Input values: n is the number of the atom.  r (optional) is a
  /// cutoff, must be less than rCut in the constructor (not
  /// checked!).
  ///
  /// In-out values: size contains the maximum space in the arrays.
  /// It is decremented by the number of neighbors placed in the
  /// arrays.  It is an error to call GetNeighbors with too small a
  /// value of size.
  /// 
  /// Out values: neighbors[] contains the numbers of the atoms,
  /// diffs[] contains the \em relative positions of the atoms,
  /// diffs2[] contains the norms of the diffs vectors.
  ///
  /// Return value: The number of neighbors.
  virtual int GetNeighbors(int n, int *neighbors, Vec *diffs, double *diffs2,
		   int& size, double r = -1.0) const;

  /// Get information about the neighbors of atom n ("half" neighbor list)
  ///
  /// This version of GetNeighbors only returns the numbers of the neighbors.
  /// It is intended for the Python interface.
  virtual void GetNeighbors(int n, vector<int> &neighbors) const;
  
  /// GetFullNeighbors is as GetNeighbors, but return a full list
  int GetFullNeighbors(int n, int *neighbors, Vec *diffs, double *diffs2,
		       int& size, double r = -1.0) const;

  /// Get information about the neighbors of atom n (full neighbor list)
  ///
  /// This version of GetNeighbors only returns the numbers of the neighbors.
  /// It is intended for the Python interface.
  void GetFullNeighbors(int n, vector<int> &neighbors) const;
  
  /// Return the guaranteed maximal length of a single atom's NB list.

  /// Call this before using GetNeighbors() to make sure the arrays
  /// are big enough.  The value may change when the neighbor list is
  /// updated. 
  int MaxNeighborListLength() const {return maxLength;}

  /// Get the number of atoms in the corresponding list of atoms.
  int GetNumberOfAtoms() const {return nAtoms;}  // Used by interface.

  /// Return the cutoff distance (rCut) specified when creating this nblist.
  double GetCutoffRadius() const {return rCut;}

  /// Return the cutoff distance including twice the drift.
  double GetCutoffRadiusWithDrift() const {return rCut + 2*drift;}

  /// Remake the list for one or more atoms that have moved.
  ///
  /// Their neighbor's lists will also be remade, their identities
  /// will be reported in the set 'affected'.
  void RemakeLists(const set<int> &modified, set<int> &affected);

  /// Test the partial remaking of lists.  NEVER CALL ON IN-USE NB LIST!
  ///
  /// This function is the Python interface to RemakeLists.  It should
  /// only be called directly for testing purposes.  Calling it on a
  /// neighbor list used by a potential will lead to INCORRECT
  /// energies and forces!
  int TestPartialUpdate(set<int> modified, PyObject *pyatoms);

  /// Normalize the positions and calculate scaled space version
  ///
  /// This is used when a neighbor list is updated
  void ScaleAndNormalizePositions();

  /// Normalize some positions and calculate scaled space version
  ///
  /// The first argument is a set of atoms to be normalized, the
  /// corresponding scaled positions are placed in scaledpos.

  void ScaleAndNormalizePositions(const set<int> &modified,
                                 vector<Vec> &scaledpos);

  /// Return the atoms access object.  Used by a few tool functions.
  virtual Atoms *GetAtoms() const {return atoms;}

  string GetName() const {return "NeighborList";}
  
  /// Print internal info about an atom
  virtual void print_info(int n);

  /// Print memory usage
  virtual long PrintMemory() const;

protected:
  /// Generate a new neighbor list.
  virtual void MakeList();

  /// Make the lists of neighboring cells.
  void MakeNeighboringCellLists();

  void CheckFullListConsistency(const string where, bool chkdst = true);

  void printlist(int n) const;

  double GetMaxStrainDisplacement();

protected:
  Atoms *atoms;   ///< A pointer to the atoms.
  int nAtoms;     ///< The number of atoms excluding ghosts.
  int nAllAtoms;  ///< The number of atoms including ghosts.
  double rCut;    ///< The cutoff radius.
  double rCut2;   ///< The square of the cutoff radius.
  double drift;   ///< The maximally allowed drift of an atom.
  double drift2;  ///< The square of the maximally allowed drift of an atom.
  int maxLength;  ///< The lenght of the longest neighbor list.
  bool firsttime; ///< True during the very first update.
  bool fulllists; ///< True if full neighbor lists are supported.
  bool pbc[3];    ///< Boundary conditions at last update.
  Vec storedSuperCell[3]; ///< So full neighbor list can be accessed.
  Vec referenceSuperCell[3];  ///< For detecting shape changes.
  
  /// Cell locator used when constructing the neighbor list.
  NeighborCellLocator *cells;
  PyObject *cells_obj;
  
  /// Table of possible translation vectors
  vector<IVec> translationTable;

  /// The actual neigbor list (half list)
  vector< vector< pair<int, translationsidx_t> > > nbList;

  /// The complementary neighbor list if full lists are enabled.
  vector< vector< pair<int, translationsidx_t> > > complNbList;
};

} // end namespace

#endif //  NEIGHBORLIST2

