// -*- C++ -*-
//
// KimNeighborNEIGHPUREH.cpp: Kim pure half neighbor list.
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

#include "KimNeighborNEIGHPUREH.h"
#include "KimAsapPython.h"
#include "KIM_API_C.h"
#include "KIM_API_status.h"
#include <math.h>

namespace ASAPSPACE {

PyAsap_NeighborLocatorObject *PyAsap_NewKimNeighborNEIGHPUREH(intptr_t* pkim,
                                                              KimAtoms *atoms,
                                                              double rCut)
{
  PyAsap_NeighborLocatorObject *self;

  self = (PyAsap_NeighborLocatorObject*) malloc(sizeof(PyAsap_NeighborLocatorObject));
  if (self == NULL)
    throw AsapError("malloc failed.");

  self->ob_refcnt = 1;
  self->weakrefs = NULL;
  self->fulllist = false;
  self->cobj = new KimNeighborNEIGHPUREH(pkim, atoms, rCut);
  if (self->cobj == NULL)
    {
      CHECKREF(self);
      Py_DECREF(self);
      throw AsapError("Failed to create a new NeighborList object.");
    }
  return self;
}

} // end namespace

KimNeighborNEIGHPUREH::KimNeighborNEIGHPUREH(intptr_t* pkim, KimAtoms *atoms,
                                             double rCut)
  : KimNeighborLocator(pkim, atoms)
{
  this->rcut = rCut;
  rcut2 = rCut*rCut;
}

int KimNeighborNEIGHPUREH::GetNeighbors(int n, int *neighbors, Vec *diffs, double *diffs2,
                                        int& size, double r /* = -1.0 */ ) const
{
  int currentAtom;
  int number;
  int *rawneighbors;
  double *Rij;  // Not used in this KIM neighborlist mode.
  int ier = KIM_API_get_neigh((void*)pkim, nbmode, check_iterator(n),
      &currentAtom, &number, &rawneighbors, &Rij);
  if (KIM_STATUS_OK != ier)
    throw AsapError("KIM_API_get_neigh failed ") << __FILE__ << ":" << __LINE__;
  assert(currentAtom == n);
  // Now construct the list of distance vectors
  const Vec *pos = atoms->GetPositions();
  const Vec *thispos = pos + n;
  int numnb = 0;
  double rcut2 = this->rcut2;
  if (r > 0)
    rcut2 = r*r;
  for (int i = 0; i < number; i++)
    {
      const Vec *otherpos = pos + rawneighbors[i];
      diffs[numnb] = *otherpos - *thispos;
      diffs2[numnb] = diffs[numnb] * diffs[numnb];
      if (diffs2[numnb] <= rcut2)
        {
          neighbors[numnb] = rawneighbors[i];
          numnb++;
        }
    }
  assert(numnb <= size);
  size -= numnb;
  return numnb;
}
