// -*- C++ -*-
// ExceptionInterface.h: Python interface to Asap exeptions.
//
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


#ifndef _INTERFACE_H
#define _INTERFACE_H

#include "AsapPython.h"
#include "Exception.h"
#ifdef CATCHFLOAT
#define _GNU_SOURCE 1
#include <fenv.h>
#ifndef ALLOCATEEXCEPTION
extern
#endif // not ALLOCATEEXPRESSION
fenv_t PyAsap_savedfloatenv; 
#endif // CATCHFLOAT

namespace ASAPSPACE {

#ifndef ALLOCATEEXCEPTION
extern
#endif
PyObject *PyAsap_ErrorObject;

int PyAsap_InitExceptionInterface(PyObject *module);

#define CATCHEXCEPTION_CMD(cmd) catch (AsapError &e)	\
  {							\
    cmd                                                 \
    string msg = e.GetMessage();			\
    PyErr_SetString(PyAsap_ErrorObject, msg.c_str());	\
    return NULL;					\
  }							\
  catch (AsapPythonError &e)				\
  {							\
    cmd                                                 \
    return NULL;					\
  }							\
  catch (AssertionFailed &e)				\
  {							\
    cmd                                                 \
    string msg = e.GetMessage();			\
    PyErr_SetString(PyExc_AssertionError, msg.c_str());	\
    return NULL;					\
  }

#define CATCHEXCEPTION CATCHEXCEPTION_CMD()

// For potentials: Release memory if an exception occured during a calculation.
#define POTCATCHEXCEPTION CATCHEXCEPTION_CMD(self->cobj->RecoverAfterException();)

#ifdef CATCHFLOAT
#define FP_EXCEPT_ON {fegetenv(&PyAsap_savedfloatenv); \
    feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);}

#define FP_EXCEPT_OFF fesetenv(&PyAsap_savedfloatenv)
#else
#define FP_EXCEPT_ON
#define FP_EXCEPT_OFF
#endif  // CATCHFLOAT

} // end namespace

#endif // _INTERFACE_H


