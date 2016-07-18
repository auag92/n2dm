// OpenKIMinterface.cpp - Python interface to OpenKIM models.
//
// This file is part of the optional Asap module to support OpenKIM
// models.  Defines the Python interface to the modules in the
// other files in OpenKIMimport.

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


#ifndef _OPENKIMINFO_INTERFACE_H
#define _OPENKIMINFO_INTERFACE_H

#include "AsapPython.h"
#include "Asap.h"

namespace ASAPSPACE {

extern char OpenKIMinfo_Docstring[];
extern char OpenKIMcalculator_Docstring[];

PyObject *PyAsap_NewOpenKIMinfo(PyObject *noself, PyObject *args,
                                PyObject *kwargs);

PyObject *PyAsap_NewOpenKIMcalculator(PyObject *noself, PyObject *args,
                                      PyObject *kwargs);

int PyAsap_InitOpenKIMInterface(PyObject *module);

} // end namespace

#endif // _OPENKIMINFO_INTERFACE_
