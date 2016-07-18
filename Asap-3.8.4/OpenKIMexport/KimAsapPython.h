// -*- C++ -*-
//
// KimAsapPython.h: Fake AsapPython.h file defining fake Python objects.
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



// Fake AsapPython.h file defining fake Python objects!

#ifndef KIMASAPPYTHON_H
#define KIMASAPPYTHON_H

#include "KIM_API_C.h"
#include "KIM_API_status.h"

#define PyObject_HEAD int ob_refcnt;

typedef struct {
  PyObject_HEAD
} PyObject;


#define Py_INCREF(op) ((PyObject*)(op))->ob_refcnt++

#define Py_DECREF(op)                                   \
    do {                                                \
        if (--((PyObject*)(op))->ob_refcnt == 0)        \
          free(op);                                     \
    } while (0)

#define Py_XINCREF(op) do { if ((op) == NULL) ; else Py_INCREF(op); } while (0)
#define Py_XDECREF(op) do { if ((op) == NULL) ; else Py_DECREF(op); } while (0)

#define CHECKREF(x)

#define Py_None NULL

#define PyErr_SetString(a, b) KIM_API_report_error(__LINE__, __FILE__, (char *) b, KIM_STATUS_FAIL)

// The following looks wrong on a 64-bit machine, but is actually OK
// as all stuff depending on this type being 32-bit is gone in the KIM
// wrapper case.
typedef int npy_int32;

#endif // not KIMASAPPYTHON_H
