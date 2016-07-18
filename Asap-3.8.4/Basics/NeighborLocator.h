// -*- C++ -*-
//
// NeighborLocator.h: Abstract base class for all neighbor lists and
// similar.
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


#ifndef NEIGHBORLOCATOR_H
#define NEIGHBORLOCATOR_H

#include "AsapPython.h"
#include "AsapObject.h"
#include "Templates.h"
#include "Exception.h"
#include "Atoms.h"
#include <set>
#include <vector>
using std::set;
using std::vector;

namespace ASAPSPACE {

class Vec;

class NeighborLocator;  // Defined later in this file.

/// The Python object corresponding to a NeighborLocator object.
typedef struct {
  PyObject_HEAD
  NeighborLocator *cobj;
  PyObject *weakrefs;
  bool fulllist;
} PyAsap_NeighborLocatorObject;

/// Mostly abstract base class for neighbor locators and neighbor lists.

/// A little common functionality is included in this otherwise
/// abstract base class, in particular a number of flags best handled
/// with inline member functions, and unlikely to need redefinition.

class NeighborLocator : public AsapObject
{
protected:
  NeighborLocator() {invalid=true;}

public:
  // Destructor should be protected, but friend template functions do
  // not work with GNU C++.
  virtual ~NeighborLocator() {};

  /// Check if the neighbor list can still be reused, update if not.
  ///
  /// If an implementation does not support multithreading, it must
  /// throw an exception if called with non-default arguments.
  virtual bool CheckAndUpdateNeighborList() = 0;

  /// Check if the neighbor list can still be reused, update if not.
  ///
  /// This version is used when called from Python
  virtual bool CheckAndUpdateNeighborList(PyObject *atoms) = 0;

  /// Check the neighbor list.
  ///
  /// Check if the neighbor list can still be reused, return true if
  /// it should be updated.
  virtual bool CheckNeighborList() = 0;

  /// Update neighbor list
  virtual void UpdateNeighborList() = 0;
  
  /// Get wrapped positions of all the atoms
  virtual const vector<Vec> &GetWrappedPositions() const = 0;

  virtual void GetWrappedPositions(vector<Vec> &wp) const = 0;
  
  /// Get scaled positions of all the atoms
  virtual const vector<Vec> &GetScaledPositions() const = 0;

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
  virtual void GetNeighbors(int n, vector<int> &neighbors) const = 0;

  /// Return the guaranteed maximal length of a single atom's NB list.

  /// Call this before using GetNeighbors() to make sure the arrays
  /// are big enough.  The value may change when the neighbor list is
  /// updated. 
  virtual int MaxNeighborListLength() const = 0;

  /// Return the cutoff distance of this neighbor locator.
  virtual double GetCutoffRadius() const = 0;

  /// Return the cutoff distance including twice the drift.
  virtual double GetCutoffRadiusWithDrift() const = 0;
    
  /// Get the number of atoms in the corresponding list of atoms.
  virtual int GetNumberOfAtoms() const = 0;  // Used by the Python interface

  /// Invalidate a neighbor list.
  void Invalidate() {invalid = true;}

  /// Return true if neighbor list has been invalidated.
  bool IsInvalid() {return invalid;}

  /// Return the atoms access object.  Used by a few tool functions.
  virtual Atoms *GetAtoms() const = 0;

  // The following methods are only implemented in a "Full" neighbor locator.

  virtual int GetFullNeighbors(int n, int *neighbors, Vec *diffs,
			       double *diffs2, int& size,
			       double r = -1.0) const
  {throw
      AsapError("Internal error: Calling half neighbor locator as a full one.");}

  virtual void GetFullNeighbors(int n, vector<int> &neighbors) const
  {throw
      AsapError("Internal error: Calling half neighbor locator as a full one.");}

  /// Print internal info about an atom
  virtual void print_info(int n) = 0;

  /// Print memory usage
  virtual long PrintMemory() const = 0;
  
protected:
  bool invalid;   ///< True if a list has been invalidated.

};

} // end namespace

#endif // NEIGHBORLOCATOR_H
