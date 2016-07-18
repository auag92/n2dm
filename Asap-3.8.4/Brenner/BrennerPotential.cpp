// BrennerPotential.cpp - Asap implementation of the Brenner Potential.
//
// Copyright (C) 2010-2011 Jakob Schiotz and Center for Individual
// Nanoparticle Functionality, Department of Physics, Technical
// University of Denmark.  Email: schiotz@fysik.dtu.dk
//
// This file is part of Asap version 3.
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

#include "BrennerPotential.h"
#include "AtomPairInfoState.h"
#include "asapbrenner.h"
#include "expand.h"
#include <iostream>

// SMARTCUTOFF allows for a smaller cutoff if only some elements are present.
#define SMARTCUTOFF

// The following two constants are taken from the original defs.h
#define PSEED              0.6955
//#define RLL                2.5

#if 1
#define VERB(x) if (verbose == 1) std::cerr << x
#else
#define VERB(x)
#endif

int BrennerPotential::ktype_to_z[5];
int BrennerPotential::z_to_ktype[MAXATNO+1];
double BrennerPotential::rb2[4][4];     
double BrennerPotential::rmax[4][4];    
double BrennerPotential::rmax_nosq;     


BrennerPotential::BrennerPotential(PyObject *self) : Potential(self)
{
  apis = new AtomPairInfoState;
  atoms = NULL;
  nblist = NULL;
  nblist_obj = NULL;
  counter = 0;
  counter_z = 0;
}

BrennerPotential::~BrennerPotential()
{
  Py_XDECREF(nblist_obj);
  if (atoms != NULL)
    AsapAtoms_DECREF(atoms);
  delete apis;
}

void BrennerPotential::Initialize()
{
  DEBUGPRINT;
  // Initialize ktype_to_z and z_to_ktype
  ktype_to_z[0] = 0;
  ktype_to_z[1] = 6;  // Carbon
  ktype_to_z[2] = 1;  // Hydrogen
  ktype_to_z[3] = 14; // Silicon
  ktype_to_z[4] = 32; // Germanium
  for (int i = 0; i < MAXATNO+1; i++)
    z_to_ktype[i] = 0;
  for (int i = 0; i < 5; i++)
    z_to_ktype[ktype_to_z[i]] = i;

  for(int i = 0; i < 4; ++i)
    for(int j = 0; j < 4; ++j)
      rb2[i][j] = 0;   // Looks like otherwise may be uninitialized!

  init_c();
  init_xh();
  init_in2();
  init_in3();
  
  // Initialize rmax and RLIST (not used?).  Uses rb2 initialized in init_c().
  //double rll = RLL;
  // RLIST is apparently the square of the neighborlist cutoff.  XXX Check!
  rmax_nosq = 0;
  for(int i = 0; i < 4; ++i)
    {
      for(int j = 0; j < 4; ++j)
        {
          Float r2 = rb2[i][j];
          //RLIST[i+1][j+1] = (r2+rll)*(r2+rll);
          rmax[i][j] = r2*r2;
	  if (rmax_nosq < r2)
	    rmax_nosq = r2;
        }
    }
  
  si_ge_init();
}


void BrennerPotential::SetAtoms(PyObject *a, Atoms* accessobj)
{
  VERB(" SetAtoms");
  if (accessobj != NULL)
    throw AsapError("BrennerPotential::SetAtoms called with accessobj != NULL");
  if (atoms == NULL)
    atoms = new NormalAtoms();
  assert(atoms != NULL);
}

double BrennerPotential::GetPotentialEnergy(PyObject *a)
{
  VERB(" Energy[");
  Calculate(a);
  VERB("]");
  return Epot;
}

const vector<Vec> &BrennerPotential::GetForces(PyObject *a)
{
  VERB(" Force[");
  Calculate(a);
  VERB("]");
  return force;
}

void BrennerPotential::Calculate(PyObject *a)
{
  assert(atoms != NULL);
  atoms->Begin(a);
  z = atoms->GetAtomicNumbers();
  positions = atoms->GetPositions();
  nAtoms = atoms->GetNumberOfAtoms();
  if (counter != atoms->GetPositionsCounter())
    {
      // The atoms have been modified.  Do a calculation.
      Epot = 0.0;
      force.resize(nAtoms);
      for (vector<Vec>::iterator i = force.begin();
	   i != force.end(); ++i)
	(*i)[0] = (*i)[1] = (*i)[2] = 0.0;
      
      if (counter_z != atoms->GetNumbersCounter())
	{
	  CountAtoms();
	  counter_z = atoms->GetNumbersCounter();
	}
      CheckAndUpdateNeighborList();
      VERB("c");
      Epot = caguts();  // Do the actual calculation.
      counter = atoms->GetPositionsCounter();
    }
  atoms->End();
}

void BrennerPotential::CheckAndUpdateNeighborList()
{
  VERB("n");
  if (nblist == NULL)
    {
      // Make a new neighbor list.
      VERB("N");
#ifdef SMARTCUTOFF
      rmax_nosq = 0.0;
      set<int> elements;
      atoms->GetListOfElements(elements);
      for (set<int>::const_iterator z1 = elements.begin();
	   z1 != elements.end(); ++z1)
	for (set<int>::const_iterator z2 = elements.begin();
	     z2 != elements.end(); ++z2)
	  {
	    int k1 = z_to_ktype[*z1] - 1;
	    int k2 = z_to_ktype[*z2] - 1;
	    double r = sqrt(rmax[k1][k2]);
	    if (r > rmax_nosq)
	      rmax_nosq = r;
	  }
#endif
      double driftfactor = 0.1;
      PyAsap_NeighborLocatorObject *nbl =
	PyAsap_NewNeighborList(atoms, rmax_nosq, driftfactor);
      nblist = dynamic_cast<NeighborList *>(nbl->cobj);
      assert(nblist != NULL);
      nblist_obj = (PyObject *) nbl;
      nblist->EnableFullNeighborLists();
      nblist->CheckAndUpdateNeighborList();
      return;
    }
  bool update = (nblist->IsInvalid() || nblist->CheckNeighborList());
  if (update)
    {
      VERB("U");
      nblist->UpdateNeighborList();
    }
}

void BrennerPotential::CountAtoms()
{
  VERB("+");
  for (int i = 0; i < 5; i++)
    noa[i] = 0;
  for (int i = 0; i < nAtoms; i++)
    {
      int zz = z[i];
      if (zz < 1 || zz > MAXATNO)
	throw AsapError("Invalid atomic number: z[") << i << "]=" << zz;
      noa[z_to_ktype[zz]]++;
    }
  if (noa[0])
    throw AsapError("BrennerPotential only supports Hydrogen, Carbon, Silicon and Germanium.");
}

