// -*- C++ -*-
//
// KimNeighborNEIGHRVECH.h: Kim RVEC half neighbor list.
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

#ifndef KIMNEIGHBORNEIGHRVECF_H
#define KIMNEIGHBORNEIGHRVECF_H

#include "KimNeighborNEIGHRVECH.h"

namespace ASAPSPACE {

PyAsap_NeighborLocatorObject *PyAsap_NewKimNeighborNEIGHRVECF(intptr_t* pkim,
                                                              KimAtoms *atoms,
                                                              double rCut);

class KimNeighborNEIGHRVECF : public KimNeighborNEIGHRVECH
{
protected:
  KimNeighborNEIGHRVECF(intptr_t* pkim, KimAtoms *atoms, double rCut):
    KimNeighborNEIGHRVECH(pkim, atoms, rCut) {};

  friend PyAsap_NeighborLocatorObject *PyAsap_NewKimNeighborNEIGHRVECF(intptr_t *pkim,
                                                                       KimAtoms *atoms,
                                                                       double rCut);

  friend void PyAsap_Dealloc<PyAsap_NeighborLocatorObject>(PyObject *self);

public:
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
                           int& size, double r = -1.0) const
    {throw AsapError("Trying to use an OpenKIM full neighbor list as a half list (KimNeighborNEIGHRVECF)");}

  virtual int GetFullNeighbors(int n, int *neighbors, Vec *diffs, double *diffs2,
                           int& size, double r = -1.0) const;
};

} // end namespace

#endif // !KIMNEIGHBORNEIGHRVECH_H

