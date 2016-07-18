// -*- C++ -*-
// NeighborCellLocator.h:  Cell-based algorithm for finding neighbors.
//
// Copyright (C) 2008-2011 Jakob Schiotz and Center for Individual
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


#ifndef NEIGHBORCELLLOCATOR
#define NEIGHBORCELLLOCATOR

#include "AsapPython.h"
#include "Asap.h"
#include "Vec.h"
#include "NeighborLocator.h"
#include "IVec.h"
#include "Templates.h"
#include <vector>
using std::vector;
#include <utility>
using std::pair;
#include <math.h>

namespace ASAPSPACE {

//typedef unsigned char translationsidx_t;   // Could also be int
typedef int translationsidx_t;   // Could also be unsigned char
//typedef long translationsidx_t;   // Silly!

PyAsap_NeighborLocatorObject *PyAsap_NewNeighborCellLocator(Atoms *a,
     double rCut, double driftfactor = 0.05, bool slave = false);

class NeighborCellLocator : public NeighborLocator
{
protected:
  /// Generate a neighbor list for atoms a with cutoff rCut.

  /// The neighbor list will contain all neighbors within the distance
  /// rCut.  The neighborlist can be reused until an atom has moved
  /// more than rCut*driftfactor.  The slave parameter is set to true
  /// if this NeighborCellLocator is used to help building another
  /// NeighborLocator, it suppresses the code that would otherwise
  /// trigger migrations in parallel simulations.
  NeighborCellLocator(Atoms *a, double rCut, double driftfactor = 0.05,
		      bool slave = false);
  virtual ~NeighborCellLocator();

  friend PyAsap_NeighborLocatorObject *PyAsap_NewNeighborCellLocator(
       Atoms *a, double rCut, double driftfactor, bool slave);

  friend void PyAsap_Dealloc<PyAsap_NeighborLocatorObject>(PyObject *self);

public:
  /// Check if the neighbor list can still be reused, update if not.
  bool CheckAndUpdateNeighborList();

  /// Check if the neighbor list can still be reused, update if not.
  ///
  /// This version is used when called from Python
  virtual bool CheckAndUpdateNeighborList(PyObject *atoms);

  /// Check the neighbor list.
  ///
  /// Check if the neighbor list can still be reused, return true if
  /// it should be updated.
  virtual bool CheckNeighborList();

  /// Update neighbor list
  virtual void UpdateNeighborList();

  /// Get wrapped positions of all the atoms
  const vector<Vec> &GetWrappedPositions() const {assert(wrappedPositionsValid); 
    return wrappedPositions;}

  void GetWrappedPositions(vector<Vec> &wp) const {assert(wrappedPositionsValid); 
    wp.insert(wp.begin(), wrappedPositions.begin(), wrappedPositions.end());}
  
  /// Get scaled positions of all the atoms
  const vector<Vec> &GetScaledPositions() const;

  /// Get reference positions (used by NeighborList)
  const Vec *GetReferencePositions() const {return &referencePositions[0];}
  
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
  int GetNeighbors(int n, int *neighbors, Vec *diffs, double *diffs2,
		   int& size, double r = -1.0) const;

  
  /// Get information about the neighbors of atom n ("half" neighbor list)
  ///
  /// This version of GetNeighbors only returns the numbers of the neighbors.
  /// It is intended for the Python interface.
  void GetNeighbors(int n, vector<int> &neighbors) const;

  int GetFullNeighbors(int n, int *neighbors, Vec *diffs, double *diffs2,
		       int& size, double r = -1.0) const;

  /// Get information about the neighbors of atom n (full neighbor list)
  ///
  /// This version of GetNeighbors only returns the numbers of the neighbors.
  /// It is intended for the Python interface.
  void GetFullNeighbors(int n, vector<int> &neighbors) const;


  /// Return the neighbors and the corresponding translations.
  ///
  /// The vectors are cleared before data is put into them.
  ///
  /// Return value: The number of neighbors.
#if 0
  int GetListAndTranslations(int n, vector<int> &neighbors,
			     vector<translationsidx_t> &translations) const;
#endif
  int GetListAndTranslations(int n, vector< pair<int,translationsidx_t> >
			     &neighbors) const;
  
  int GetComplementaryListAndTranslations(int n,
					  vector< pair<int,translationsidx_t> >
					  &neighbors) const;

  /// Remake neighbor lists when a few atoms have been modified.
  ///
  /// This version, unlike NeighborList::RemakeList does
  /// not report back which other atoms have been affected.
  void RemakeLists_Simple(const set<int> &modified);
  
  /// Return the guaranteed maximal length of a single atom's NB list.

  /// Call this before using GetNeighbors() to make sure the arrays
  /// are big enough.  The value may change when the neighbor list is
  /// updated. 
  int MaxNeighborListLength() const {return maxLength;}

  /// Get the number of atoms in the corresponding list of atoms.
  int GetNumberOfAtoms() const {return nAtoms;}  // Used from swig.

  /// Return the cutoff distance (rCut) specified when creating this nblist.
  double GetCutoffRadius() const {return rCut;}

  /// Return the cutoff distance including twice the drift.
  double GetCutoffRadiusWithDrift() const {
    throw AsapError("NeighborCellLocator cannot predict its drift radius.");}

  /// Get a copy of the table of translations (27 entries)
  void GetTranslationTable(vector<IVec> &table) const;

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

  /// Normalizing new positions using the normalization already calculated.
  ///
  /// This is used when checking a neighborlist.
  void RenormalizePositions();

#if 0
  /// Renormalize the old positions using old and new positions of the atoms.
  ///
  /// This is called by a master neighbor list when it has normalized
  /// the positions of the atoms.  The argument is the positions of
  /// the atoms *before* renormalization, this method uses them
  /// together with the current positions to update its own list of
  /// old positions.
  void RenormalizeReferencePositions(const vector<Vec> &oldpos);
#endif

  /// Update reference positions of some atoms
  void UpdateReferencePositions(const set<int> &modified);
  
  /// Return the atoms access object.  Used by a few tool functions.
  virtual Atoms *GetAtoms() const {return atoms;}

  string GetName() const {return "NeighborCellLocator";}

    /// Print internal info about an atom
  virtual void print_info(int n);

  /// Print memory usage
  virtual long PrintMemory() const;

protected:
  /// Generate a new neighbor list.
  virtual void MakeList();

  /// Make the lists of neighboring cells.
  void MakeNeighboringCellLists();

  /// Make translation table
  void MakeTranslationTable();

  typedef vector<pair<int, translationsidx_t> > nbcell_t;

  const nbcell_t &makeNbCells(int thiscell, nbcell_t &NbCells) const;

  int CommonGetNeighbors(int n, int *neighbors, Vec *diffs, double *diffs2,
			 int& size, double r, bool wantfull) const;

  void CommonGetNeighbors(int a1, vector<int> &neighbors,
			  bool wantfull) const;

  double get_drift() const;
  
protected:
  Atoms *atoms;   ///< A pointer to the atoms.
  int nAtoms;     ///< The number of atoms excluding ghosts.
  int nAllAtoms;  ///< The number of atoms including ghosts.
  double rCut;    ///< The cutoff radius.
  double rCut2;   ///< The square of the cutoff radius.
  double minboxsize;  ///< The minimal box size
  bool periodic[3];   ///< The boundary conditions.
  bool oldperiodic[3];    ///< The boundary conditions at the last update.
  int maxLength;  ///< The guaranteed max length of a neighbor list.
  int nCells[3];  ///< Number of cells
  int nTotalCells[4];
  int nCellsTrue[3];
  int nCellsGapStart[3];
  int nCellsGapSize[3];
  bool slave;     ///< This locator is slave of another one.
  double size[3];     ///< Size of system in scaled space
  double minimum[3];  ///< Offset of system in scaled space
  vector<Vec> referencePositions;  ///< Positions at last update.
  vector<Vec> wrappedPositions;    ///< Wrapped positions.
  vector<Vec> scaledPositions;     ///< Scaled positions.
  vector<Vec> offsetPositions;     ///< wrappedPositions - positions.
  vector<Vec> scaledOffsetPositions;  
  bool scaledPositionsValid;
  bool wrappedPositionsValid;

  Vec old_inverse_cell[3];  ///< Inverse unit cell of last renormalization.
  int supercell_counter; ///< When was old_inverse_cell last updated?
  
  /// The list of cells containing the atoms
  vector< vector<int> > cells;

  /// The number of the cell to which an atom belongs
  vector<int> cellIndices;

  /// For each cell, a list of the neighboring cells, and their offset
  /// across the periodic boundaries.
  vector<IVec> neighborCellOffsets;

  /// List of neighboring cells, valid for a cell not touching the boundary.
  nbcell_t nbCells_inside;
  nbcell_t nbCells_left;
  nbcell_t nbCells_right;
  nbcell_t nbCells_top;
  nbcell_t nbCells_bottom;
  nbcell_t nbCells_front;
  nbcell_t nbCells_back;
  
  /// Table of possible translation vectors
  vector<IVec> translationTable;
};

} // end namespace

#endif //  NEIGHBORCELLLOCATOR
