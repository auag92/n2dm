// -*- C++ -*-
//
// DummyNeighborLocator.h: Fake do-nothing neighbor list used by OpenKIMcalculator
//
// Copyright (C) 2001-2011 Jakob Schiotz and Center for Individual
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

#include "AsapPython.h"
#include "DummyNeighborLocator.h"
#include "NeighborLocatorInterface.h"

namespace ASAPSPACE {


PyAsap_NeighborLocatorObject *PyAsap_NewDummyNeighborLocator(Atoms *atoms,
    double rCut, double driftfactor)
{
  PyAsap_NeighborLocatorObject *self;

  self = PyObject_NEW(PyAsap_NeighborLocatorObject,
                      &PyAsap_NeighborLocatorType);
  if (self == NULL)
    throw AsapError("OOPS XXXX");

  self->weakrefs = NULL;
  self->fulllist = false;
  self->cobj = new DummyNeighborLocator(atoms, rCut, driftfactor);
  if (self->cobj == NULL)
    {
      CHECKREF(self);
      Py_DECREF(self);
      throw AsapError("Failed to create a new NeighborCellLocator object.");
    }
  return self;
}


DummyNeighborLocator::DummyNeighborLocator(Atoms *a, double rCut, double driftfactor) : NeighborLocator()
{
  atoms = a;
  AsapAtoms_INCREF(atoms);
  drift = rCut * (1 + driftfactor);
}

DummyNeighborLocator::~DummyNeighborLocator()
{
  AsapAtoms_DECREF(atoms);
}

bool DummyNeighborLocator::CheckNeighborList()
{
  if (invalid)
    return true;

  int n_at = atoms->GetNumberOfAtoms();
  if (nAtoms != n_at || nAllAtoms != n_at + atoms->GetNumberOfGhostAtoms())
    return true;

  bool updateRequired = false;
  const Vec *positions = atoms->GetPositions();
  double max2 = drift * drift;
  for (int n = 0; n < nAtoms; n++)
    if (Length2(positions[n] - referencePositions[n]) > max2)
      {
        updateRequired = true;
        break;
      }
  return updateRequired;
}

void DummyNeighborLocator::UpdateNeighborList()
{
  invalid = false;
  nAtoms = atoms->GetNumberOfAtoms();
  nAllAtoms = nAtoms + atoms->GetNumberOfGhostAtoms();
  const Vec *positions = atoms->GetPositions();
  referencePositions.resize(nAllAtoms);
  for (int i = 0; i < nAllAtoms; i++)
    referencePositions[i] = positions[i];
}

} // namespace
