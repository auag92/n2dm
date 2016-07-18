// -*- C++ -*-
//
// asap_emt_driver.cpp: OpenKIM Model Driver interface for EMT.
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

#include "KIM_API_C.h"
#include "KIM_API_status.h"
#include "Asap.h"
#include "Debug.h"
#include "asap_kim_api.h"
#include "asap_emt_driver.h"
#include "KimParameterProvider.h"

#include <stdlib.h>

static int asap_emt_driver_initmodel(AsapKimPotential *model);

KimEMT::KimEMT(AsapKimPotential *owner, EMTParameterProvider *provider) : EMT(NULL, NULL)
{
  CONSTRUCTOR;
  this->owner = owner;
  nblist = NULL;
  nblist_obj = NULL;
  provider_obj = NULL;  // Bypass EMT's Python-based memory management.
  this->provider = provider;
}

KimEMT::~KimEMT()
{
  DESTRUCTOR;
  assert(provider_obj == NULL);
  delete provider;
  if (nblist != NULL)
    delete nblist;  // Not deleted in the real potential.
}

void KimEMT::CreateNeighborList()
{
  PyAsap_NeighborLocatorObject *nbl = owner->CreateNeighborList(atoms, rNbCut, driftfactor);
  always_fullnblist = owner->UsesKimFullList();
  nblist = nbl->cobj;
  nblist_obj = (PyObject *) nbl;
  nblist->UpdateNeighborList();
}

/* Destruction function */

static int asap_emt_destroy(void *km)
{
  intptr_t* pkim = *((intptr_t**) km);
  int ier;
  AsapKimPotential *model = (AsapKimPotential *) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier)
    {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
      return ier;
    }
  delete model;

  // Remove the model pointer from the model buffer (if left, the API will
  // release it with free instead of delete).
  KIM_API_set_model_buffer(pkim, NULL, &ier);
  if (KIM_STATUS_OK > ier)
    {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_set_model_buffer", ier);
      return ier;
    }
  return KIM_STATUS_OK;
}

/* Reinit function */

static int asap_emt_reinit(void *km)
{
  intptr_t* pkim = *((intptr_t**) km);
  int ier;
  // Remove the old model
  AsapKimPotential *model = (AsapKimPotential *) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier)
    {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
      return ier;
    }

  // Add the new model
  AsapKimPotential *newmodel = new AsapKimPotential(pkim, model->paramfile_names, model->nmstrlen,
      model->numparamfiles, true);
  delete model;
  return asap_emt_driver_initmodel(newmodel);
}

/* Initialization function */
extern "C" int model_driver_init(void *km, char* paramfile_names, int* nmstrlen, int* numparamfiles)
{
  intptr_t* pkim = *((intptr_t**) km);
  int ier;

  if(*numparamfiles !=1)
  {
     ier = KIM_STATUS_FAIL;
     KIM_API_report_error(__LINE__, __FILE__, "Wrong number of parameter files", ier);
     return ier;
  }

  /* store pointer to compute function in KIM object */
  ier = KIM_API_set_method(pkim, "compute", 1, (func_ptr) &AsapKimPotential::compute_static);
  if (KIM_STATUS_OK > ier)
    {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_set_method", ier);
      return ier;
    }

  /* store pointer to destroy function in KIM object */
  ier = KIM_API_set_method(pkim, "destroy", 1, (func_ptr) &asap_emt_destroy);
  if (KIM_STATUS_OK > ier)
    {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_set_method", ier);
      return ier;
    }

  /* store pointer to reinit function in KIM object */
  ier = KIM_API_set_method(pkim, "reinit", 1, (func_ptr) &asap_emt_reinit);
  if (KIM_STATUS_OK > ier)
    {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_set_method", ier);
      return ier;
    }

  // Create the model
  AsapKimPotential *model = new AsapKimPotential(pkim, paramfile_names,
      *nmstrlen, *numparamfiles, true);
  return asap_emt_driver_initmodel(model);
}

static int asap_emt_driver_initmodel(AsapKimPotential *model)
{
  intptr_t *pkim = model->pkim;
  const char *paramfile_names = model->paramfile_names;
  double *model_cutoff;
  int ier;

  KimParameterProvider *provider = NULL;
  try
  {
      provider = new KimParameterProvider(paramfile_names, pkim);
  }
  catch (AsapError &e)
  {
      std::cerr << e.GetMessage() << std::endl;
      return KIM_STATUS_FAIL;
  }
  model->potential = new KimEMT(model, provider);

  // Store the model in the model buffer
  KIM_API_set_model_buffer(pkim, (void*) model, &ier);
  if (KIM_STATUS_OK > ier)
    {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_set_model_buffer", ier);
      return ier;
    }

  /* store maximal model cutoff in KIM object */
  model_cutoff = (double*) KIM_API_get_data(pkim, "cutoff", &ier);
  if (KIM_STATUS_OK > ier){
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_data", ier);
      return ier;
  }
  provider->CalcGammaEtc();
  *model_cutoff = provider->GetListCutoffDistance();
  //std::cerr << "Setting cutoff to " << *model_cutoff << std::endl;

  return KIM_STATUS_OK;
}

