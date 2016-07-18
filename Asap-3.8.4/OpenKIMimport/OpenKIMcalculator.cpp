// OpenKIMcalculator.cpp - interface to OpenKIM models.
//
// This class is part of the optional Asap module to support OpenKIM
// models.  The class OpenKIMcalculator does the actual interfacing
// to the model.

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
#include "OpenKIMcalculator.h"
#include "KIM_API.h"
#include "KIM_API_status.h"
#include "Atoms.h"
#include "NeighborList.h"
#include "DummyNeighborLocator.h"
#include "Timing.h"
//#define ASAPDEBUG
#include "Debug.h"
#include <iostream>
#include <cstdlib>

using std::cerr;
using std::endl;
using std::flush;

#if 1
#define VERB(x) if (verbose == 1) cerr << x
#else
#define VERB(x)
#endif

#define KIMASAPERROR(m, x) AsapError("KIM Error '") << get_kimerror(m, x) \
  << "' at " << __FILE__ << ":" << __LINE__

#define KIMCHECK(m, x) if (x < KIM_STATUS_OK) throw KIMASAPERROR(m,x)

static const char *get_kimerror(KIM_API_model *m, int x)
{
  const char *message;
  int err = m->get_status_msg(x, &message);
  if (err < KIM_STATUS_OK)
    throw AsapError("KIM double-fault: Error while getting error message.  Original error code = ") << x;
  return message;
}

extern "C" {
  static int get_neigh_wrapper(void **kimmdl, int *mode, int* request,
      int *particle, int *numnei, int **nei1particle,
      double **rij)
  {
    KIM_API_model *model = (KIM_API_model *) *kimmdl;
    if (model == NULL)
      return KIM_STATUS_API_OBJECT_INVALID;
    int x;
    OpenKIMcalculator *calc = (OpenKIMcalculator *) model->get_data("neighObject", &x);
    if (x != KIM_STATUS_OK)
      return x;  // return KIM_STATUS_API_OBJECT_INVALID ???
    return calc->get_neigh(mode, request, particle, numnei, nei1particle, rij);
  }
}

OpenKIMcalculator::OpenKIMcalculator(PyObject *self) : Potential(self)
{
  CONSTRUCTOR;
  DEBUGPRINT;
  model = NULL;
  model_initialized = false;
  counters.nblist = counters.compute_energy = counters.compute_force = -1;
  driftfactor = 0.05;  // Drift factor for the neighbor list.
  pls_alloc_n = 0;
  nAtoms = nSize = 0;
  nAtomsAlloc = nSizeAlloc = 0;
  nSpecies = 0;
  cutoff = 0.0;
  nblist_iterator = -10000;
  nblist = NULL;
  nblist_obj = NULL;
  this->self = self;
  // NB: No INCREF!  This is not a reference, but a pointer internal
  // to the Python object.
  DEBUGPRINT;
}

OpenKIMcalculator::~OpenKIMcalculator()
{
  DESTRUCTOR;
  DEBUGPRINT;
  if (model != NULL)
    {
      if (model_initialized)
        {
          model->model_destroy();
          int err;
          model->free(&err);  // Part of destructor?
        }
      delete model;
    }
  if (atoms != NULL)
    AsapAtoms_DECREF(atoms);
  Py_XDECREF(nblist_obj);
}

void OpenKIMcalculator::Initialize(const char *descr, const char *name)
{
  DEBUGPRINT;
  model = new KIM_API_model();
  if (model == NULL)
    throw AsapError("Failed to allocate OpenKIM model.  This should never happen!");
  int x = model->string_init(descr, name);
  if (x < KIM_STATUS_OK)
    throw AsapError("Failed to initialize OpenKIM model, see kim.log for details.");
  // Connect the most essential variables.
  model->setm_data(&x, 3*4,
      "numberOfParticles",      1,      &nSize,         1,
      "numberOfSpecies",        1,      &nSpecies,        1,
      "cutoff",                 1,      &cutoff,        1);
  KIMCHECK(model, x);

  // We need to initialize the model already now, so the neighbor list method
  // and the cutoff is available to the Python interface *before* it calls the
  // C++ version of SetAtoms2.
  model->model_init();
  KIMCHECK(model, x);
  model_initialized = true;
  // Identify the neighbor list type.
  DEBUGPRINT;
  const char *NBCstr;
  x = model->get_NBC_method(&NBCstr);
  KIMCHECK(model, x);
  DEBUGPRINT;
  if (!strcmp("CLUSTER",NBCstr))
     nblisttype = nbl_cluster;
  else if (!strcmp("MI_OPBC_H",NBCstr))
     nblisttype = nbl_miopbc_h;
  else if (!strcmp("MI_OPBC_F",NBCstr))
     nblisttype = nbl_miopbc_f;
  else if (!strcmp("NEIGH_PURE_H",NBCstr))
     nblisttype = nbl_pure_h;
  else if (!strcmp("NEIGH_PURE_F",NBCstr))
     nblisttype = nbl_pure_f;
  else if (!strcmp("NEIGH_RVEC_H",NBCstr))
     nblisttype = nbl_rvec_h;
  else if (!strcmp("NEIGH_RVEC_F",NBCstr))
     nblisttype = nbl_rvec_f;
  else
    throw AsapError("Unknown KIM neighbor list type.");
  DEBUGPRINT;
  nblist_half = model->is_half_neighbors(&x);
  KIMCHECK(model, x);
  if (nblist_half && (nblisttype != nbl_cluster))
    {
      x = model->set_data("numberContributingParticles", 1, &nAtoms);
      KIMCHECK(model, x);
    }
  DEBUGPRINT;
}

void OpenKIMcalculator::PleaseAllocate(string quantity, bool alloc)
{
  DEBUGPRINT;
  if (quantity == "energy")
    pls_alloc.energy = alloc;
  else if (quantity == "particleEnergy")
    pls_alloc.particleEnergy = alloc;
  else if (quantity == "forces")
    pls_alloc.forces = alloc;
  else if (quantity == "virial")
    pls_alloc.virial = alloc;
  else if (quantity == "particleVirial")
    pls_alloc.particleVirial = alloc;
  else
    throw AsapError("Unknown argument to OpenKIMcalculator::PleaseAllocate: ") << quantity;
  pls_alloc_n++;
  DEBUGPRINT;
}

void OpenKIMcalculator::SetAtoms(PyObject *pyatoms, Atoms* accessobj)
{
  DEBUGPRINT;
  if (atoms != NULL)
    {
      // SetAtoms should only do anything the first time it is called.
      // Subsequent calls should just check for accessobj being NULL.
      if (accessobj != NULL)
        throw AsapError("OpenKIMcalculator::SetAtoms called multiple times with accessobj != NULL");
      // SetAtoms should not do anything if called more than once!
      return;
    }

  // The first time SetAtoms is being called some initialization is done.
  if (accessobj != NULL)
    {
      atoms = accessobj;
      AsapAtoms_INCREF(atoms);
    }
  else
    atoms = new NormalAtoms();
  assert(atoms != NULL);
  atoms->Begin(pyatoms);
  nAtoms = atoms->GetNumberOfAtoms();
  nSize = nAtoms;   // nSize will be changed when the neighborlist is built.
  // Find the number of elements
  set<int> elements;
  atoms->GetListOfElements(elements);
  DEBUGPRINT;
  nSpecies = elements.size();
  atoms->End();
  DEBUGPRINT;
}

const vector<Vec> &OpenKIMcalculator::GetForces(PyObject *pyatoms){
  USETIMER("OpenKIMcalculator::GetForces")
  DEBUGPRINT;
  if (!pls_alloc.forces)
    throw AsapError("OpenKIMcalculator not prepared to calculate forces.");
  VERB(" Forces[");
  Calculate(pyatoms, false, true);
  assert(forces.size() == nSize);
  return forces;
}

const vector<double> &OpenKIMcalculator::GetPotentialEnergies(PyObject *pyatoms)
{
  USETIMER("OpenKIMcalculator::GetPotentialEnergies");
  DEBUGPRINT;
  if (!pls_alloc.particleEnergy)
    throw AsapError("OpenKIMcalculator not prepared to calculate particle energies.");
  VERB(" Energies[");
  Calculate(pyatoms, true, false);
  assert(particleEnergy.size() == nAtoms);
  return particleEnergy;
}

double OpenKIMcalculator::GetPotentialEnergy(PyObject *pyatoms)
{
  USETIMER("OpenKIMcalculator::GetPotentialEnergy");
  DEBUGPRINT;
  VERB(" Energy[");
  Calculate(pyatoms, true, false);
  return energy;
}


const vector<SymTensor> &OpenKIMcalculator::GetVirials(PyObject *pyatoms)
{
  USETIMER("OpenKIMcalculator::GetVirials");
  DEBUGPRINT;
  if (!pls_alloc.particleVirial)
    throw AsapError("OpenKIMcalculator not prepared to calculate particle virials.");
  VERB(" Virials[");
  Calculate(pyatoms, false, true);
  assert(particleVirial.size() == nSize);
  return particleVirial;
}

SymTensor OpenKIMcalculator::GetVirial(PyObject *pyatoms)
{
  USETIMER("OpenKIMcalculator::GetVirial");
  DEBUGPRINT;
  if (!pls_alloc.virial)
    throw AsapError("OpenKIMcalculator not prepared to calculate virial.");
  VERB(" Virial[");
  Calculate(pyatoms, false, true);
  return virial;
}

bool OpenKIMcalculator::CalcReq_Energy(PyObject *pyatoms)
{
  DEBUGPRINT;
  atoms->Begin(pyatoms);
  bool required = (!pls_alloc.energy || counters.compute_energy != atoms->GetPositionsCounter());
  atoms->End();
  return required;
}

bool OpenKIMcalculator::CalcReq_Forces(PyObject *pyatoms)
{
  DEBUGPRINT;
  atoms->Begin(pyatoms);
  bool required = (!pls_alloc.forces || counters.compute_force != atoms->GetPositionsCounter());
  atoms->End();
  return required;
}

bool OpenKIMcalculator::CalcReq_Virials(PyObject *pyatoms)
{
  DEBUGPRINT;
  atoms->Begin(pyatoms);
  bool required = (!pls_alloc.virial || counters.compute_force != atoms->GetPositionsCounter());
  atoms->End();
  return required;
}

bool OpenKIMcalculator::CheckNeighborList()
{
  RETURNIFASAPERROR2(false);
  USETIMER("OpenKIMcalculator::CheckNeighborList");
  DEBUGPRINT;
  assert(atoms != NULL);
  bool update = (nblist == NULL) || nblist->IsInvalid();  // Update if invalid
  if (!update && (counters.nblist != atoms->GetPositionsCounter()))
    {
      VERB("n");
      update = nblist->CheckNeighborList();
    }
  // May communicate
  update = atoms->UpdateBeforeCalculation(update,
                                          cutoff * (1 + driftfactor));
  counters.nblist = atoms->GetPositionsCounter();
  if (nblisttype == nbl_miopbc_f || nblisttype == nbl_miopbc_h)
    SetOrthoCell();
  DEBUGPRINT;
  return update;
}

void OpenKIMcalculator::UpdateNeighborList()
{
  RETURNIFASAPERROR;
  USETIMER("OpenKIMcalculator::UpdateNeighborList");
  DEBUGPRINT;
  VERB("N");
  DEBUGPRINT;
  RETURNIFASAPERROR;
  if (nblist)
    {
      DEBUGPRINT;
      nblist->UpdateNeighborList();
      RETURNIFASAPERROR;
      {
        if ((nAtoms != atoms->GetNumberOfAtoms())
            || (nSize - nAtoms != atoms->GetNumberOfGhostAtoms()))
          {
            nAtoms = atoms->GetNumberOfAtoms();
            nSize = nAtoms + atoms->GetNumberOfGhostAtoms();
            ghostatoms = atoms->HasGhostAtoms();
            Allocate();
          }
      }
    }
  else
    {
      // First call, create the neighbor list.
      DEBUGPRINT;
      CreateNeighborList();
      RETURNIFASAPERROR;
      {
        nAtoms = atoms->GetNumberOfAtoms();
        nSize = nAtoms + atoms->GetNumberOfGhostAtoms();
        ghostatoms = atoms->HasGhostAtoms();
        Allocate();
      }
    }
  DEBUGPRINT;
}

void OpenKIMcalculator::CreateNeighborList()
{
  DEBUGPRINT;
  RETURNIFASAPERROR;
  MEMORY;
  USETIMER("OpenKIMcalculator::CreateNeighborList");
  PyAsap_NeighborLocatorObject *nbl;
  int x;
  if (nblisttype == nbl_cluster)
    nbl = PyAsap_NewDummyNeighborLocator(atoms, cutoff, driftfactor);
  else
    nbl = PyAsap_NewNeighborList(atoms, cutoff, driftfactor);
  nblist = nbl->cobj;
  nblist_obj = (PyObject *) nbl;
  if (!nblist_half)
    {
      NeighborList *nblist2 = dynamic_cast<NeighborList *>(nblist);
      assert(nblist2 != NULL);
      nblist2->EnableFullNeighborLists();
    }
  nblist->UpdateNeighborList();
  switch (nblisttype)
  {
    case nbl_rvec_h:
    case nbl_rvec_f:
    case nbl_pure_h:
    case nbl_pure_f:
      model->setm_data(&x, 2*4,
           "get_neigh",     1,    &get_neigh_wrapper,   1,
           "neighObject",   1,    this,                 1);
      KIMCHECK(model, x);
      break;
    case nbl_miopbc_h:
    case nbl_miopbc_f:
      SetOrthoCell();
      model->setm_data(&x, 3*4,
           "get_neigh",      1,   &get_neigh_wrapper,   1,
           "neighObject",    1,   this,                 1,
           "boxSideLengths", 3,   orthocell,            1);
      KIMCHECK(model, x);
      break;
    case nbl_cluster:
      // No neighbor locator for this one, just a sanity check
      if (nAtoms != nSize)
        throw AsapError("KIM models with neighbor list type CLUSTER do not support ghost atoms.");
      break;
    default:
      throw AsapError("Unsupported OpenKIM neighborlist method: ") << nblisttype;
  }
#if 0
  nb_accessmode = model->get_neigh_mode(&x);
  KIMCHECK(model, x);
  if (nb_accessmode != 2)
    throw AsapError("Currently Asap only supports neighbor lists in Locator mode.");
#endif
  MEMORY;
}

void OpenKIMcalculator::Allocate()
{
  RETURNIFASAPERROR;
  assert(pls_alloc_n == 5);
  USETIMER("OpenKIMcalculator::Allocate");
  DEBUGPRINT;
  VERB(" Allocate[" << nAtoms << "," << nSize << "]" << flush);
  // WARNING: Resizing the vector may allocate way too much memory.  It
  // appears that calling reserve solves this problem.  For efficiency,
  // reserve is called with 5% extra space.  This is only necessary if the
  // atoms have ghosts, otherwise no reallocation will happen.

  // First, check if reallocation is necessary.
  if (nSize != nSizeAlloc || nAtoms != nAtomsAlloc)
    {
      DEBUGPRINT;
      // Do the reserve trick if the atoms have ghosts.
      if (ghostatoms)
        {
          // Atoms have ghosts.  Reserve a bit extra memory to minimize
          // reallocation due to migration.
          // If full neighbor lists are used, reserve extra space for
          // particle energies, as they will get resized when energy is calculated.
          int nSizeRes = nSize + nSize/20;
          if (pls_alloc.forces && forces.capacity() < nSize)
            forces.reserve(nSizeRes);
          if (pls_alloc.particleEnergy && particleEnergy.capacity() < nSize)
            particleEnergy.reserve(nSizeRes);
          if (pls_alloc.particleVirial && particleVirial.capacity() < nSize)
            particleVirial.reserve(nSizeRes);
        }
      DEBUGPRINT;
      // Resize the arrays.  Arrays that are not needed will be given the
      // size 1, so we can still take the address of the data in model->setm_data.
      species.resize(nSize);
      if (pls_alloc.forces)
        forces.resize(nSize);
      else
        forces.resize(1);
      if (pls_alloc.particleEnergy)
        particleEnergy.resize(nSize);
      else
        particleEnergy.resize(1);
      if (pls_alloc.particleVirial)
        particleVirial.resize(nSize);
      else
        particleVirial.resize(1);
      int x;
      model->setm_data(&x, 6*4,
          "particleSpecies", nSize,    &species[0],        1,
          "energy",          1,        &energy,            pls_alloc.energy,
          "particleEnergy",  nSize,    &particleEnergy[0], pls_alloc.particleEnergy,
          "forces",          3*nSize,  &forces[0],         pls_alloc.forces,
          "virial",          6,        &virial[0],         pls_alloc.virial,
          "particleVirial",  6*nSize,  &particleVirial[0], pls_alloc.particleVirial);
      KIMCHECK(model,x);
    }
  DEBUGPRINT;
}

void OpenKIMcalculator::Calculate(PyObject *pyatoms, bool calc_energy, bool calc_force)
{
  assert(atoms != NULL);
  atoms->Begin(pyatoms);
  recalc.nblist = CheckNeighborList();
  recalc.compute_energy = calc_energy && (counters.compute_energy != atoms->GetPositionsCounter());
  recalc.compute_force = calc_force && (counters.compute_force != atoms->GetPositionsCounter());
  if (recalc.compute_energy || recalc.compute_force)
    {
      DEBUGPRINT;
      DoCalculate();
      VERB("]" << flush);
    }
  else
    {
      assert(recalc.nblist == false);
      VERB("-]");
    }
  DEBUGPRINT;
  atoms->End();
}

void OpenKIMcalculator::DoCalculate()
{
  DEBUGPRINT;
  int x;
  if (recalc.nblist)
    UpdateNeighborList();
  if (pls_alloc.particleEnergy)
    {
      particleEnergy.resize(nSize);
      x = model->set_data("particleEnergy", nSize, &particleEnergy[0]);
      KIMCHECK(model, x);
    }
  // We need to map atomic numbers to the model type codes.
  const asap_z_int *z = atoms->GetAtomicNumbers();
  assert(species.size() == nSize);
  for (int i = 0; i < nSize; i++)
    species[i] = z_to_typecode[z[i]];
  x = model->set_data("coordinates", 3*nSize, (void *) atoms->GetPositions());
  KIMCHECK(model, x);
  model->setm_compute(&x, 5*3,
      "energy", (int) recalc.compute_energy, pls_alloc.energy,
      "particleEnergy", (int) recalc.compute_energy, pls_alloc.particleEnergy,
      "forces", (int) recalc.compute_force, pls_alloc.forces,
      "virial", (int) recalc.compute_force, pls_alloc.virial,
      "particleVirial", recalc.compute_force, pls_alloc.particleVirial);
  KIMCHECK(model, x);
  x = model->model_compute();
  KIMCHECK(model, x);
  if (recalc.compute_energy)
    counters.compute_energy = atoms->GetPositionsCounter();
  if (recalc.compute_force)
    counters.compute_force = atoms->GetPositionsCounter();
  if (pls_alloc.particleEnergy)
    particleEnergy.resize(nAtoms);
  DEBUGPRINT;
}

int OpenKIMcalculator::get_neigh(int *mode, int *request,
        int *particle, int *numnei, int **nei1particle,
        double **rij)
{
  //std::cerr << "get_neigh mode=" << *mode << " rq=" << *request << std::endl;
  int rq;  // Requested particle
  if (*mode == 1)
    {
      // Locator Mode
      rq = *request;
    }
  else if (*mode == 0)
    {
      if (*request == 0)
        {
          nblist_iterator = 0;
          *numnei = 0;
          return KIM_STATUS_NEIGH_ITER_INIT_OK;
        }
      else if (*request == 1)
        {
          rq = nblist_iterator++;
          if (rq >= nSize)
            return KIM_STATUS_NEIGH_ITER_PAST_END;
        }
      else
        return KIM_STATUS_NEIGH_INVALID_REQUEST;
    }
  else
    return KIM_STATUS_NEIGH_INVALID_MODE;
  if (rq < 0 || rq >= nSize)
    return KIM_STATUS_NEIGH_INVALID_REQUEST;
  int maxnb = nblist->MaxNeighborListLength();
  switch (nblisttype)
  {
    case nbl_rvec_h:
    case nbl_rvec_f:
      nb_buffer_n.resize(maxnb);
      nb_buffer_rij.resize(maxnb);
      nb_buffer_dist.resize(maxnb);
      if (nblist_half)
        *numnei = nblist->GetNeighbors(rq, &nb_buffer_n[0], &nb_buffer_rij[0],
            &nb_buffer_dist[0], maxnb);
      else
        *numnei = nblist->GetFullNeighbors(rq, &nb_buffer_n[0], &nb_buffer_rij[0],
            &nb_buffer_dist[0], maxnb);
      *particle = rq;
      *nei1particle = &nb_buffer_n[0];
      *rij = (double *) &nb_buffer_rij[0];
      break;
    case nbl_miopbc_h:
    case nbl_miopbc_f:
    case nbl_pure_h:
    case nbl_pure_f:
      if (nblist_half)
        nblist->GetNeighbors(rq, nb_buffer_n);
      else
        nblist->GetFullNeighbors(rq, nb_buffer_n);
      *numnei = nb_buffer_n.size();
      *particle = rq;
      *nei1particle = &nb_buffer_n[0];
      break;
    default:
      std::cerr << "Asap internal error: Unsupported KIM NBL method in OpenKIMcalculator::get_neigh()"
          << std::endl;
      return KIM_STATUS_NEIGH_INVALID_MODE;
  }
  return KIM_STATUS_OK;
}

const char *OpenKIMcalculator::GetNBCmethod()
{
  DEBUGPRINT;
  const char *result;
  int err = model->get_NBC_method(&result);

  KIMCHECK(model, err);
  DEBUGPRINT;
  return result;
}

void OpenKIMcalculator::SetOrthoCell()
{
  DEBUGPRINT;
  const Vec *cell = atoms->GetCell();
  for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
        if ((i != j ) && fabs(cell[i][j]) > 1e-13)
            throw AsapError("Attemting to use MI_OPBC boundary conditions with non-orthogonal cell.");
      orthocell[i] = cell[i][i];
    }
  DEBUGPRINT;
}


void OpenKIMcalculator::GetSupportedSymbols(std::vector<const char *> &symbols)
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


int OpenKIMcalculator::GetParticleTypeCode(const char *symbol)
{
  DEBUGPRINT;
  int err;
  int code = model->get_species_code(symbol, &err);
  KIMCHECK(model, err);
  DEBUGPRINT;
  return code;
}

