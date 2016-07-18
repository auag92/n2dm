// Copyright (C) 2008-2011 Jakob Schiotz and Center for Individual
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

#include "ImagePotential.h"
#include "Vec.h"
#include "ImageAtoms.h"
#include "Timing.h"
//#define ASAPDEBUG
#include "Debug.h"
#include <set>
using std::cerr;
using std::endl;
using std::flush;
using std::set;

ImagePotential::ImagePotential(PyObject *self, Potential *p) : Potential(self)
{
  DEBUGPRINT;
  CONSTRUCTOR;
  img_atoms = NULL;
  potential = p;
  assert(potential != NULL);
  virials_collect_cnt = force_collect_cnt = 0;
}

ImagePotential::~ImagePotential()
{
  DEBUGPRINT;
  DESTRUCTOR;
  if (img_atoms != NULL)
    AsapAtoms_DECREF(img_atoms);
}

void ImagePotential::SetAtoms(PyObject *pyatoms,
		              Atoms* accessobj /* = NULL */)
{
  DEBUGPRINT;
  if (atoms != NULL)
    {
      // SetAtoms should only do anything the first time it is called.
      // Subsequent calls should just check for accessobj being NULL.
      if (accessobj != NULL && accessobj != atoms)
        throw AsapError("EMT::SetAtoms called multiple times with accessobj != NULL");
      // SetAtoms should not do anything if called more than once!
      return;
    }
  if (accessobj == NULL)
    accessobj = new NormalAtoms();
  else
    AsapAtoms_INCREF(accessobj);
  img_atoms = new ImageAtoms(accessobj);
  AsapAtoms_DECREF(accessobj); // img_atoms now owns it.
  atoms = img_atoms;
  // Since we are wrapping the potential within the Python implementation
  // of set_atoms/SetAtoms in the potential being wrapped, we have to bypass
  // the Python layer of SetAtoms, and call the C++ layer directly.
  potential->SetAtoms(pyatoms, img_atoms);
  DEBUGPRINT;
}

const vector<Vec> &ImagePotential::GetForces(PyObject *a)
{
  USETIMER("ImagePotential::GetForces");
  DEBUGPRINT;
  img_atoms->Begin(a, true);  // Allow reopen.
  int cnt = img_atoms->GetPositionsCounter();
  if (force_collect_cnt != cnt)
    {
      DEBUGPRINT;
      const vector<Vec> &orig_forces = potential->GetForces(a);
      forces = orig_forces;
      CollectFromImages<Vec>(forces);
      force_collect_cnt = cnt;
    }
  img_atoms->End();
  return forces;
}

const vector<SymTensor> &ImagePotential::GetVirials(PyObject *a)
{
  USETIMER("ImagePotential::GetStresses");
  DEBUGPRINT;
  img_atoms->Begin(a, true);  // Allow reopen.
  int cnt = img_atoms->GetPositionsCounter();
  if (virials_collect_cnt != cnt)
    {
      DEBUGPRINT;
      const vector<SymTensor> &orig_virials = potential->GetVirials(a);
      virials = orig_virials;
      CollectFromImages<SymTensor>(virials);
      virials_collect_cnt = cnt;
    }
  img_atoms->End();
  DEBUGPRINT;
  return virials;
}

template <class T>
void ImagePotential::CollectFromImages(vector<T> &data)
{
  DEBUGPRINT;
  int nOrig = img_atoms->GetOriginalNumberOfAtoms();
  const vector<int>  orig_map = img_atoms->GetOriginalAtomsMap();

  for (int i = 0; i < orig_map.size(); i++)
    {
      data[orig_map[i]] += data[i + nOrig];
    }
  data.resize(nOrig);
  DEBUGPRINT;
}

