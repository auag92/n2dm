// -*- C++ -*-
//
// asap_kim_api.cpp: Common interfaces classes for Asap potentials in OpenKIM.
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


#include "KimAsapPython.h"
#include "asap_kim_api.h"
#include "Potential.h"
#include "NeighborCellLocator.h"
#include "KimNeighborMIOPBCH.h"
#include "KimNeighborMIOPBCF.h"
#include "KimNeighborNEIGHPUREH.h"
#include "KimNeighborNEIGHPUREF.h"
#include "KimNeighborNEIGHRVECH.h"
#include "KimNeighborNEIGHRVECF.h"
#include "KIM_API_C.h"
#include "KIM_API_status.h"
#include "SymTensor.h"
#include "Debug.h"
#include <stdlib.h>
#include <string.h>

namespace ASAPSPACE {
  int verbose = 0;
}

AsapKimPotential::AsapKimPotential(intptr_t *pkim, const char* paramfile_names,
                                   int nmstrlen, int numparamfiles,
                                   bool supportvirial)
{
  CONSTRUCTOR;
  int ier;

  potential = NULL;
  atoms = NULL;
  this->pkim = pkim;
  // Store the parameter file name(s) for later use by reinit.
  this->paramfile_names = new char[nmstrlen * numparamfiles];
  memcpy(this->paramfile_names, paramfile_names, nmstrlen * numparamfiles);
  this->nmstrlen = nmstrlen;
  this->numparamfiles = numparamfiles;
  this->supportvirial = supportvirial;
  // Check for the neighbor list type.
  ier = KIM_API_get_NBC_method(pkim, &NBCstr);
  if (KIM_STATUS_OK > ier)
    {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_NBC_method", ier);
      exit(1);
    }
  if (!strcmp("CLUSTER",NBCstr))
    {
      nblisttype = nbl_cluster;
      need_contrib = false;
    }
  else if (!strcmp("MI_OPBC_H",NBCstr))
    {
      nblisttype = nbl_miopbc_h;
      need_contrib = true;
    }
  else if (!strcmp("MI_OPBC_F",NBCstr))
    {
      nblisttype = nbl_miopbc_f;
      need_contrib = false;
    }
  else if (!strcmp("NEIGH_PURE_H",NBCstr))
    {
      nblisttype = nbl_pure_h;
      need_contrib = true;
    }
  else if (!strcmp("NEIGH_PURE_F",NBCstr))
    {
      nblisttype = nbl_pure_f;
      need_contrib = false;
    }
  else if (!strcmp("NEIGH_RVEC_H",NBCstr))
    {
      nblisttype = nbl_rvec_h;
      need_contrib = true;
    }
  else if (!strcmp("NEIGH_RVEC_F",NBCstr))
    {
      nblisttype = nbl_rvec_f;
      need_contrib = false;
    }
  else
    {
      KIM_API_report_error(__LINE__, __FILE__, "Unknown NBC method", KIM_STATUS_FAIL);
      exit(1);
    }
}

AsapKimPotential::~AsapKimPotential()
{
  DESTRUCTOR;
  if (potential != NULL)
    delete potential;
  if (atoms != NULL)
    AsapAtoms_DECREF(atoms);
  delete paramfile_names;
}

int AsapKimPotential::compute_static(void* km)
{
  int ier;
  assert(km != NULL);
  intptr_t* pkim = *((intptr_t**) km);
  assert(pkim != NULL);
  AsapKimPotential *model = (AsapKimPotential *) KIM_API_get_model_buffer(pkim, &ier);
  if (KIM_STATUS_OK > ier)
    {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_get_model_buffer", ier);
      return ier;
    }
  assert(model != NULL);
  ier = model->compute(km);
  return ier;
}

int AsapKimPotential::compute(void* km)
{
  int ier;

  assert(potential != NULL);
  assert(pkim = *((intptr_t**) km));   // Sanity check
  double *cutoff = NULL;
  int *nAtoms = NULL;
  int *nTotalAtoms = NULL;
  int* particleSpecies = NULL;

  // Flags indicating what we need to compute.
  int comp_energy;
  int comp_force;
  int comp_particleEnergy;
  int comp_virial = 0;
  int comp_particleVirial = 0;

  // Quantities to be computed
  double *coords = NULL;
  double *energy = NULL;
  double *forces = NULL;
  double *particleEnergy = NULL;
  double *virial = NULL;
  double *particleVirial = NULL;

  /* check to see if we have been asked to compute the forces and particleEnergy */
  /* If we support virials, also check if we should calculate them */
  KIM_API_getm_compute(pkim, &ier, 5*3,
                       "energy",         &comp_energy,         1,
                       "forces",         &comp_force,          1,
                       "particleEnergy", &comp_particleEnergy, 1,
                       "virial",         &comp_virial,         supportvirial,
                       "particleVirial", &comp_particleVirial, supportvirial
  );
  if (KIM_STATUS_OK > ier)
    {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_compute", ier);
      return ier;
    }

  KIM_API_getm_data(pkim, &ier, 10*3,
                    "cutoff",                      &cutoff,          1,
                    "numberOfParticles",           &nTotalAtoms,     1,
                    "numberContributingParticles", &nAtoms,          need_contrib,
                    "particleSpecies",             &particleSpecies, 1,
                    "coordinates",                 &coords,          1,
                    "energy",                      &energy,          comp_energy,
                    "forces",                      &forces,          comp_force,
                    "particleEnergy",              &particleEnergy,  comp_particleEnergy,
                    "virial",                      &virial,          comp_virial,
                    "particleVirial",              &particleVirial,  comp_particleVirial
                    );
  if (KIM_STATUS_OK > ier)
    {
      KIM_API_report_error(__LINE__, __FILE__, "KIM_API_getm_data", ier);
      return ier;
    }
  if (!need_contrib)
    nAtoms = nTotalAtoms;

  if (atoms == NULL)
    {
      // First call, create the Atoms interface object
      atoms = new KimAtoms(pkim);
      assert(atoms != NULL);
      atoms->ReInit(*nAtoms, *nTotalAtoms - *nAtoms, coords, particleSpecies);
      potential->SetAtoms(NULL, atoms);
    }
  else
    {
      atoms->ReInit(*nAtoms, *nTotalAtoms - *nAtoms, coords, particleSpecies);
    }

  // Now do the actual computation
  try
  {
      if (comp_particleEnergy)
        {
          const vector<double> &energies_v = potential->GetPotentialEnergies(NULL);
          assert(energies_v.size() == *nAtoms);
          for (int i = 0; i < *nAtoms; i++)
            particleEnergy[i] = energies_v[i];
        }
      if (comp_energy)
        *energy = potential->GetPotentialEnergy(NULL);
      if (comp_particleVirial)
        {
          const vector<SymTensor> &virials = potential->GetVirials(NULL);
          assert(virials.size() == *nTotalAtoms);
          const double *virials_ptr = (double *) &virials[0];
          for (int i = 0; i < 6*(*nTotalAtoms); i++)
            particleVirial[i] = virials_ptr[i];
        }
      if (comp_virial)
        {
          SymTensor v = potential->GetVirial(NULL);
          for (int i = 0; i < 6; i++)
            virial[i] = v[i];
        }
      if (comp_force)
        {
          const vector<Vec> &forces_v = potential->GetForces(NULL);
          assert(forces_v.size() == *nTotalAtoms);
          const double *forces_v_ptr = (double *) &forces_v[0];
          for (int i = 0; i < 3*(*nTotalAtoms); i++)
            forces[i] = forces_v_ptr[i];
        }
  }
  catch (AsapError &e)
  {
      ier = KIM_STATUS_FAIL;
      string msg = e.GetMessage();
      // Will the following line store a pointer to something inside msg? Hopefully not!
      KIM_API_report_error(__LINE__, __FILE__, (char *) msg.c_str(), ier);
      return ier;
  }
  return KIM_STATUS_OK;
}


PyAsap_NeighborLocatorObject *AsapKimPotential::CreateNeighborList(KimAtoms *atoms,
                                                                   double cutoff,
                                                                   double drift)
{
  int ier;
  PyAsap_NeighborLocatorObject *nblist;
  if (nblisttype == nbl_cluster)
    {
      atoms->SetPBC(false, false, false);
      nblist = PyAsap_NewNeighborCellLocator(atoms, cutoff, drift);
      assert(nblist != NULL);
    }
  else if (nblisttype == nbl_miopbc_h)
    {
      atoms->SetPBC(true, true, true);
      nblist = PyAsap_NewKimNeighborMIOPBCH(pkim, atoms, cutoff);
    }
  else if (nblisttype == nbl_miopbc_f)
    {
      atoms->SetPBC(true, true, true);
      nblist = PyAsap_NewKimNeighborMIOPBCF(pkim, atoms, cutoff);
    }
  else if (nblisttype == nbl_pure_h)
    {
      atoms->SetPBC(true, true, true);
      nblist = PyAsap_NewKimNeighborNEIGHPUREH(pkim, atoms, cutoff);
    }
  else if (nblisttype == nbl_pure_f)
    {
      atoms->SetPBC(true, true, true);
      nblist = PyAsap_NewKimNeighborNEIGHPUREF(pkim, atoms, cutoff);
    }
  else if (nblisttype == nbl_rvec_h)
    {
      atoms->SetPBC(true, true, true);
      nblist = PyAsap_NewKimNeighborNEIGHRVECH(pkim, atoms, cutoff);
    }
  else if (nblisttype == nbl_rvec_f)
    {
      atoms->SetPBC(true, true, true);
      nblist = PyAsap_NewKimNeighborNEIGHRVECF(pkim, atoms, cutoff);
    }
  else
    {
      throw AsapError("Unknown NBC method: ") << NBCstr;
    }
  return nblist;
}

bool AsapKimPotential::UsesKimFullList()
{
  return ( (nblisttype == nbl_miopbc_f) ||
           (nblisttype == nbl_pure_f) ||
           (nblisttype == nbl_rvec_f) );
}
