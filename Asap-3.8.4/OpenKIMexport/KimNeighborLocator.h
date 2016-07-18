// -*- C++ -*-
//
// KimNeighborLocator.h: Common base class for KIM interface neighbor locators.
//
// Copyright (C) 2012-2013 Jakob Schiotz and the Department of Physics,
// Technical University of Denmark.  Email: schiotz@fysik.dtu.dk
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




#ifndef KIMNEIGHBORLOCATOR_H
#define KIMNEIGHBORLOCATOR_H

#include "NeighborLocator.h"

#define INVALIDMETHOD {throw AsapError("Invalid KimNeighborLocator method called.");}

namespace ASAPSPACE {

class KimNeighborLocator : public NeighborLocator
{
public:
  KimNeighborLocator(intptr_t *pkim, KimAtoms *atoms);

  virtual ~KimNeighborLocator();

  virtual string GetName() const INVALIDMETHOD;

  /// Check if the neighbor list can still be reused, update if not.
  ///
  /// KIM version: Call UpdateNeighborList and then return false.
  virtual bool CheckAndUpdateNeighborList() {return CheckNeighborList();}

  /// Check if the neighbor list can still be reused, update if not.
  ///
  /// This version is used when called from Python
  virtual bool CheckAndUpdateNeighborList(PyObject *atoms) {return CheckNeighborList();}

  /// Check the neighbor list.
  ///
  /// Check if the neighbor list can still be reused, return true if
  /// it should be updated.
  ///
  /// KIM version: Call UpdateNeigborList and then return false.

  virtual bool CheckNeighborList();

  /// Update neighbor list
  ///
  /// KIM version: Extract any necessary info from the API object.
  virtual void UpdateNeighborList() = 0;

  /// Get wrapped positions of all the atoms
  virtual const vector<Vec> &GetWrappedPositions() const INVALIDMETHOD

  virtual void GetWrappedPositions(vector<Vec> &wp) const INVALIDMETHOD;

  /// Get scaled positions of all the atoms
  virtual const vector<Vec> &GetScaledPositions() const INVALIDMETHOD;

  /// Get info about the neighbors of atom n.  The most important method :-)
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
                           int& size, double r = -1.0) const = 0;

  /// Get the neighbors of atom n (half neighbor list).
  ///
  /// This version of GetNeighbors only returns the numbers of the neighbors.
  /// It is intended for the Python interface.
  virtual void GetNeighbors(int n, vector<int> &neighbors) const INVALIDMETHOD;

  /// Return the guaranteed maximal length of a single atom's NB list.

  /// Call this before using GetNeighbors() to make sure the arrays
  /// are big enough.  The value may change when the neighbor list is
  /// updated.
  virtual int MaxNeighborListLength() const {return 512;} // Compiled into API.

  /// Return the cutoff distance of this neighbor locator.
  virtual double GetCutoffRadius() const INVALIDMETHOD;

  /// Return the cutoff distance including twice the drift.
  virtual double GetCutoffRadiusWithDrift() const INVALIDMETHOD;

  /// Get the number of atoms in the corresponding list of atoms.
  virtual int GetNumberOfAtoms() const INVALIDMETHOD;  // Used by the Python interface

  /// Return the atoms access object.  Used by a few tool functions.
  virtual Atoms *GetAtoms() const INVALIDMETHOD;

  /// Print internal info about an atom
  virtual void print_info(int n) INVALIDMETHOD;

  /// Print memory usage
  virtual long PrintMemory() const INVALIDMETHOD;

protected:
  int check_iterator(int n) const;

protected:  // Data
  intptr_t *pkim;
  KimAtoms *atoms;
  bool nbmode;
  int nAtoms;
  int nGhosts;
};

} // end namespace

#undef INVALIDMETHOD
#endif // !KIMNEIGHBORLOCATOR_H

