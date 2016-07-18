// -*- C++ -*-
//
// DummyNeighborLocator.h: Fake do-nothing neighbor list used by OpenKIMcalculator
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


#ifndef DUMMYNEIGHBORLOCATOR_H
#define DUMMYNEIGHBORLOCATOR_H

#include "NeighborLocator.h"
#include <set>
#include <vector>
using std::set;
using std::vector;

namespace ASAPSPACE {

class Vec;

#define NBLDUM {throw AsapError("Invalid DummyNeighborLocator method called: ") << __FILE__ << ":" << __LINE__;}

PyAsap_NeighborLocatorObject *PyAsap_NewDummyNeighborLocator(Atoms *a,
     double rCut, double driftfactor = 0.0);

// This dummy neighbor locator is used when the KIM Model uses CLUSTER
// neighbor locator, i.e. it does not need/want any neighborlist info.
class DummyNeighborLocator : public NeighborLocator
{
protected:
  DummyNeighborLocator(Atoms *a, double rCut, double driftfactor = 0.05);

public:
  friend PyAsap_NeighborLocatorObject *PyAsap_NewDummyNeighborLocator(Atoms *a,
      double rCut, double driftfactor);

  // Destructor should be protected, but friend template functions do
  // not work with GNU C++.
  virtual ~DummyNeighborLocator();

  string GetName() const {return "DummyNeighborLocator";}

  /// Check if the neighbor list can still be reused, update if not.
  ///
  /// If an implementation does not support multithreading, it must
  /// throw an exception if called with non-default arguments.
  virtual bool CheckAndUpdateNeighborList() NBLDUM

  /// Check if the neighbor list can still be reused, update if not.
  ///
  /// This version is used when called from Python
  virtual bool CheckAndUpdateNeighborList(PyObject *atoms) NBLDUM

  /// Check the neighbor list.
  ///
  /// Check if the neighbor list can still be reused, return true if
  /// it should be updated.
  virtual bool CheckNeighborList();

  /// Update neighbor list
  virtual void UpdateNeighborList();
  
  /// Get wrapped positions of all the atoms
  virtual const vector<Vec> &GetWrappedPositions() const NBLDUM;

  virtual void GetWrappedPositions(vector<Vec> &wp) const NBLDUM;
  
  /// Get scaled positions of all the atoms
  virtual const vector<Vec> &GetScaledPositions() const NBLDUM;

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
			   int& size, double r = -1.0) const NBLDUM;

  /// Get the neighbors of atom n (half neighbor list).
  ///
  /// This version of GetNeighbors only returns the numbers of the neighbors.
  /// It is intended for the Python interface.
  virtual void GetNeighbors(int n, vector<int> &neighbors) const NBLDUM

  /// Return the guaranteed maximal length of a single atom's NB list.

  /// Call this before using GetNeighbors() to make sure the arrays
  /// are big enough.  The value may change when the neighbor list is
  /// updated. 
  virtual int MaxNeighborListLength() const NBLDUM

  /// Return the cutoff distance of this neighbor locator.
  virtual double GetCutoffRadius() const NBLDUM

  /// Return the cutoff distance including twice the drift.
  virtual double GetCutoffRadiusWithDrift() const NBLDUM
    
  /// Get the number of atoms in the corresponding list of atoms.
  virtual int GetNumberOfAtoms() const NBLDUM  // Used by the Python interface

  /// Return the atoms access object.  Used by a few tool functions.
  virtual Atoms *GetAtoms() const NBLDUM

  // The following methods are only implemented in a "Full" neighbor locator.

  virtual int GetFullNeighbors(int n, int *neighbors, Vec *diffs,
			       double *diffs2, int& size,
			       double r = -1.0) const NBLDUM

  virtual void GetFullNeighbors(int n, vector<int> &neighbors) const NBLDUM

  /// Print internal info about an atom
  virtual void print_info(int n) NBLDUM

  /// Print memory usage
  virtual long PrintMemory() const NBLDUM
  
protected:
  Atoms *atoms;
  double drift;
  int nAtoms;
  int nAllAtoms;
  vector<Vec> referencePositions;
};

} // end namespace

#endif // DUMMYNEIGHBORLOCATOR_H
