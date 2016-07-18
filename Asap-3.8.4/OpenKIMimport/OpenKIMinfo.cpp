// OpenKIMInfo.cpp - a class to collect information about OpenKIM models.
//
// This class is part of the optional Asap module to support OpenKIM
// models.  The class OpenKIMInfo collects information about a given
// OpenKIM model, so it can be opened using the best neighbor list etc.

// Copyright (C) 2014 Jakob Schiotz and Center for Individual
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

#include "Asap.h"
#include "OpenKIMinfo.h"
#include "KIM_API.h"
#include "KIM_API_status.h"
#include "Debug.h"
#include <iostream>
#include <cstdlib>

#define KIMASAPERROR(m, x) AsapError("KIM Error '") << get_kimerror(m, x) \
  << "' at " << __FILE__ << ":" << __LINE__

#define KIMCHECK(m, x) if (x < KIM_STATUS_OK) throw KIMASAPERROR(m,x)

#define KIMASAPERROR2(m, x, aux) AsapError("KIM Error '") << get_kimerror(m, x) \
  << "' at " << __FILE__ << ":" << __LINE__ << " (aux. info: " << aux << ")"

#define KIMCHECK2(m, x, aux) if (x < KIM_STATUS_OK) throw KIMASAPERROR2(m,x,aux)

static const char *get_kimerror(KIM_API_model *m, int x)
{
  const char *message;
  int err = m->get_status_msg(x, &message);
  if (err < KIM_STATUS_OK)
    throw AsapError("KIM double-fault: Error while getting error message.  Original error code = ") << x;
  return message;
}

static const char *nbc_names[] = {
    "CLUSTER",
    "NEIGH_PURE_H", "NEIGH_PURE_F",
    "NEIGH_RVEC_H", "NEIGH_RVEC_F",
    "MI_OPBC_H", "MI_OPBC_F",
    NULL
};

static const char *nb_access_names[] = {
    "Neigh_IterAccess", "Neigh_LocaAccess", "Neigh_BothAccess",
    NULL
};

OpenKIMinfo::OpenKIMinfo(const char* name) : particletypes(NULL)
{
  DEBUGPRINT;
  this->name = name;
  model = new KIM_API_model();
  int x = model->model_info(name);
  KIMCHECK(model, x);
  DEBUGPRINT;
}

OpenKIMinfo::~OpenKIMinfo()
{
  DEBUGPRINT;
  if (particletypes != NULL)
    std::free(particletypes);
  if (model != NULL)
    {
      int err;
      model->free(&err);  // Part of destructor?
      delete model;
    }
  DEBUGPRINT;
}

void OpenKIMinfo::GetSupportedSymbols(std::vector<const char *> &symbols)
{
  DEBUGPRINT;
  int numspecies;
  int maxlength;  // Not used
  int err;
  err = model->get_num_model_species(&numspecies, &maxlength);
  KIMCHECK(model, err);
  symbols.clear();
  for (int i = 0; i < numspecies; i++)
    {
      const char *sp;
      err = model->get_model_species(i, &sp);
      symbols.push_back(sp);
    }
  DEBUGPRINT;
}


const std::vector<const char *> &OpenKIMinfo::GetSupportedNBC()
{
  DEBUGPRINT;
  supported_nbc.clear();
  const char **candidate = nbc_names;
  while (*candidate != NULL)
    {
      int err;
      int idx = model->get_index(*candidate, &err);
      if ((err >= KIM_STATUS_OK) && (idx >= 0))
        supported_nbc.push_back(*candidate);
      candidate++;
    }
  DEBUGPRINT;
  return supported_nbc;
}

const std::vector<const char *> &OpenKIMinfo::GetSupportedAccess()
{
  DEBUGPRINT;
  supported_access.clear();
  const char **candidate = nb_access_names;
  while (*candidate != NULL)
    {
      int err;
      int idx = model->get_index(*candidate, &err);
      if ((err >= KIM_STATUS_OK) && (idx >= 0))
        supported_access.push_back(*candidate);
      candidate++;
    }
  DEBUGPRINT;
  return supported_access;
}

int OpenKIMinfo::GetAPIindex(const char *name)
{
  DEBUGPRINT;
  int err;
  int idx = model->get_index(name, &err);
  if (err < KIM_STATUS_OK)
    return -1;
  else
    return idx;
}
