// NeighborLocatorInterface.h: Python interface to the EMT object.
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


#ifndef _NEIGHBORLOCATOR_INTERFACE_H
#define _NEIGHBORLOCATOR_INTERFACE_H

#include "AsapPython.h"
#include "Asap.h"

namespace ASAPSPACE {

PyObject *PyAsap_NewNeighborList_Py(PyObject *noself, PyObject *args,
				    PyObject *kwargs);

PyObject *PyAsap_NewNeighborCellLocator_Py(PyObject *noself, PyObject *args,
                                           PyObject *kwargs);

PyObject *PyAsap_NewFullNeighborList(PyObject *noself, PyObject *args,
				     PyObject *kwargs);

extern char PyAsap_NewNeighborList_Docstring[];
extern char PyAsap_NewNeighborCellLocator_Docstring[];
extern char PyAsap_NewFullNeighborList_Docstring[];

extern PyTypeObject PyAsap_NeighborLocatorType;

int PyAsap_InitNeighborLocatorInterface(PyObject *module);

} // end namespace

#endif // _NEIGHBORLOCATOR_INTERFACE_
