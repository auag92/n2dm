// -*- C++ -*-
//
// KimNeighborNEIGHRVECH.cpp: Kim RVEC half neighbor list.
//
// Copyright (C) 2001-2012 Jakob Schiotz and Center for Individual
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

#include "KimNeighborNEIGHRVECH.h"
#include "KimAsapPython.h"
#include "KIM_API_C.h"
#include "KIM_API_status.h"
#include "Asap.h"
#include "Debug.h"
#include <math.h>

namespace ASAPSPACE {

PyAsap_NeighborLocatorObject *PyAsap_NewKimNeighborNEIGHRVECH(intptr_t* pkim,
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
  self->cobj = new KimNeighborNEIGHRVECH(pkim, atoms, rCut);
  if (self->cobj == NULL)
    {
      CHECKREF(self);
      Py_DECREF(self);
      throw AsapError("Failed to create a new NeighborList object.");
    }
  return self;
}

} // end namespace

KimNeighborNEIGHRVECH::KimNeighborNEIGHRVECH(intptr_t* pkim, KimAtoms *atoms,
                                             double rCut)
  : KimNeighborLocator(pkim, atoms)
{
  CONSTRUCTOR;
  this->rcut = rCut;
  rcut2 = rCut*rCut;
}

KimNeighborNEIGHRVECH::~KimNeighborNEIGHRVECH()
{
  DESTRUCTOR;
}

int KimNeighborNEIGHRVECH::GetNeighbors(int n, int *neighbors, Vec *diffs, double *diffs2,
                                        int& size, double r /* = -1.0 */ ) const
{
  int currentAtom;
  int number;
  int *rawneighbors;
  double *Rij;
  assert(nbmode == 1);
  int ier = KIM_API_get_neigh((void*)pkim, nbmode, check_iterator(n),
      &currentAtom, &number, &rawneighbors, &Rij);
  if (KIM_STATUS_OK != ier)
    throw AsapError("KIM_API_get_neigh failed ") << __FILE__ << ":" << __LINE__;
  assert(currentAtom == n);
  // Now copy the list of distance vectors
  int numnb = 0;
  double rcut2 = this->rcut2;
  if (r > 0)
    rcut2 = r*r;
  for (int i = 0; i < number; i++)
    {
      diffs[numnb][0] = Rij[3*i];
      diffs[numnb][1] = Rij[3*i + 1];
      diffs[numnb][2] = Rij[3*i + 2];
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
