// -*- C++ -*-
//
// asap_kim_api.h: Common interfaces classes for Asap potentials in OpenKIM.
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


#ifndef ASAP_KIM_API_H
#define ASAP_KIM_API_H

#include "KimAtoms.h"
#include "NeighborLocator.h"

namespace ASAPSPACE {

class Potential;

typedef enum {
  nbl_cluster, nbl_pure_h, nbl_pure_f, nbl_miopbc_h, nbl_miopbc_f, nbl_rvec_h, nbl_rvec_f
  } nblisttype_t;

class AsapKimPotential
{
public:
  AsapKimPotential(intptr_t *pkim, const char* paramfile_names,
                   int nmstrlen, int numparamfiles,
                   bool supportvirial);
  virtual ~AsapKimPotential();

  static int compute_static(void *km);
  int compute(void *km);

  PyAsap_NeighborLocatorObject *CreateNeighborList(KimAtoms *atoms,
                                                   double cutoff,
                                                   double drift);

  bool UsesKimFullList();

public:
  Potential *potential;    // The ASAP potential being used
  intptr_t *pkim;          // Pointer to the KIM object
  char* paramfile_names;   // The original parameter file name list.  Possibly used by reinit.
  int nmstrlen;            // The original length of the file names.
  int numparamfiles;       // The original number of file names.
  bool supportvirial;      // This potential supports virials.

private:
  //Data
  KimAtoms *atoms;
  nblisttype_t nblisttype;
  bool need_contrib;
  const char *NBCstr;   // Neighbor list string, kept for error messages.
};

} // end namespace

#endif // not ASAP_KIM_API_H
