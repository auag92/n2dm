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


#include "AsapPython.h"
#include "Asap.h"
#include "Potential.h"
#include "ParallelAtoms.h"
#include "ParallelAtomsInterface.h"
#include "ExceptionInterface.h"

namespace ASAPSPACE {

PyObject *PyAsap_DistributeAtoms(PyObject *noself, PyObject *atoms)
{
  ParallelAtoms *patoms;
  try {
    patoms = new ParallelAtoms(atoms);
  }
  CATCHEXCEPTION;
  assert(patoms != NULL);
  try {
    patoms->Begin(atoms);
    CHECKNOASAPERROR;
    patoms->Distribute();
    PROPAGATEASAPERROR;
    patoms->End();
  }
  CATCHEXCEPTION_CMD(if (patoms) AsapAtoms_DECREF(patoms););
  AsapAtoms_DECREF(patoms);
  Py_RETURN_NONE;
};

} // end namespace
