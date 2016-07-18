// PotentialInterface.h: Python interface to the calculator objects.
//
// Copyright (C) 2008 Jakob Schiotz and Center for Individual
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

#ifndef _POTENTIAL_INTERFACE_H
#define _POTENTIAL_INTERFACE_H

#include "AsapPython.h"
#include "Asap.h"
#include "Potential.h"

namespace ASAPSPACE {

extern PyTypeObject PyAsap_PotentialType;

int PyAsap_InitPotentialInterface(PyObject *module);

void InitPotentialType(PyTypeObject &type, bool mark=true);

PyObject *PyAsap_PotentialUseImageAtoms(PyAsap_PotentialObject *self,
                                        PyObject *noargs);

extern char PyAsap_PotentialUseImageAtoms_Docstring[];

} // end namespace

#endif // _POTENTIAL_INTERFACE_H
