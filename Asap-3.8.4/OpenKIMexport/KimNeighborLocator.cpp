// -*- C++ -*-
//
// KimNeighborLocator.cpp: Common base class for KIM interface neighbor locators.
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

#include "KimNeighborLocator.h"
#include "KIM_API_C.h"
#include "KIM_API_status.h"
#include "Asap.h"
#include "Debug.h"

// IMPORTANT:
//
// Although the code below has been written to support both iterator and locator mode
// of the OpenKIM neighbor list, the iterator mode should not be used as it relies on
// the atoms being returned in order, which is not guaranteed by the KIM API.  For that
// reason, only Neigh_LocaAccess should be specified in the .kim file.


KimNeighborLocator::KimNeighborLocator(intptr_t *pkim, KimAtoms *atoms)
{
  CONSTRUCTOR;
  this->pkim = pkim;
  this->atoms = atoms;
  AsapAtoms_INCREF(atoms);
  nAtoms = nGhosts = 0;
  int ier;
  int loc_int = KIM_API_get_neigh_mode(pkim, &ier);
  if (KIM_STATUS_OK > ier)
    {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
      throw AsapError("KIM API error");
    }
  nbmode = (loc_int == 2);
}

KimNeighborLocator::~KimNeighborLocator()
{
  DESTRUCTOR;
  AsapAtoms_DECREF(atoms);
}

bool KimNeighborLocator::CheckNeighborList()
{
  bool result = (nAtoms != atoms->GetNumberOfAtoms() || nGhosts != atoms->GetNumberOfGhostAtoms());
  UpdateNeighborList();
  nAtoms = atoms->GetNumberOfAtoms();
  nGhosts = atoms->GetNumberOfAtoms();
  return result;
}

int KimNeighborLocator::check_iterator(int n) const
{
  if (nbmode == 1)
    return n;  // Locator mode, do nothing

  // Iterator mode.  Initialize or check for consistency.
  int particle;
  int numnei;
  int *nb;
  double *rij;
  if (n == 0)
    {
      // Initialize iteration
      int x = KIM_API_get_neigh(pkim, 0, 0, &particle, &numnei, &nb, &rij);
    }
  return 1;
}
