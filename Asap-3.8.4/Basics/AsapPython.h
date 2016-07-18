/* -*- C -*-                                                             */
/* AsapPython.h: Includes Python.h and numarrays correctly               */
/*                                                                       */
/* Copyright (C) 2008 Jakob Schiotz and Center for Individual            */
/* Nanoparticle Functionality, Department of Physics, Technical          */
/* University of Denmark.  Email: schiotz@fysik.dtu.dk                   */
/*                                                                       */
/* This file is part of Asap version 3.                                  */
/*                                                                       */
/* This program is free software: you can redistribute it and/or         */
/* modify it under the terms of the GNU Lesser General Public License    */
/* version 3 as published by the Free Software Foundation.               */
/* Permission to use other versions of the GNU Lesser General Public     */
/* License may granted by Jakob Schiotz or the head of department of the */
/* Department of Physics, Technical University of Denmark, as described  */
/* in section 14 of the GNU General Public License.                      */
/*                                                                       */
/* This program is distributed in the hope that it will be useful,       */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/* GNU General Public License for more details.                          */
/*                                                                       */
/* You should have received a copy of the GNU General Public License     */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */


/* This file should be included by any file accessing Python objects,    */
/* except the AsapModule.cpp which should import the NumPy module.       */
/*                                                                       */
/* This file MUST remain valid C as well as valid C++ !                  */
/*                                                                       */
/* This file must be included BEFORE any other header file !             */


#ifndef _ASAPPYTHON_H
#define _ASAPPYTHON_H

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL Asap_Array_API
#define NO_IMPORT_ARRAY
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#ifdef __cplusplus  // Allow inclusion from mpi.c
#include "Exception.h"

// Take a PyObject and check that it is a reasonable PyArrayObject
inline PyArrayObject *AsPyArray(PyObject *p) {
  if ((p != NULL) && !PyArray_Check(p))
    throw AsapError("Expected a NumPy array, got something else!");
  return (PyArrayObject *) p;
}
#endif // __cplusplus

#define CHECKREF(x) assert((x)->ob_refcnt >= 1 && (x)->ob_refcnt <= 100);
#define XCHECKREF(x) assert((x) == NULL || ((x)->ob_refcnt >= 1 && (x)->ob_refcnt <= 100));
#define PRINTREF(x) cerr << __FILE__ << ":" << __LINE__ << " Refcount=" << (x)->ob_refcnt << endl;

/* Python 2.3 support */
#ifndef Py_RETURN_NONE
#define Py_RETURN_NONE return Py_INCREF(Py_None), Py_None
#define Py_RETURN_TRUE return Py_INCREF(Py_True), Py_True
#define Py_RETURN_FALSE return Py_INCREF(Py_False), Py_False
#define Py_CLEAR(op)				\
        do {                            	\
                if (op) {			\
                        PyObject *tmp = (PyObject *)(op);	\
                        (op) = NULL;		\
                        Py_DECREF(tmp);		\
                }				\
        } while (0)
#endif

/* Python 2.3 and 2.4 support */
#if (PY_VERSION_HEX < 0x02050000)
typedef Py_ssize_t (*lenfunc)(PyObject *);
typedef PyObject *(*ssizeargfunc)(PyObject *, Py_ssize_t);
typedef PyObject *(*ssizessizeargfunc)(PyObject *, Py_ssize_t, Py_ssize_t);
typedef int(*ssizeobjargproc)(PyObject *, Py_ssize_t, PyObject *);
typedef int(*ssizessizeobjargproc)(PyObject *, Py_ssize_t, Py_ssize_t, PyObject *);
#endif
#if 0
typedef PyObject * (*unaryfunc)(PyObject *);
typedef PyObject * (*binaryfunc)(PyObject *, PyObject *);
typedef PyObject * (*ternaryfunc)(PyObject *, PyObject *, PyObject *);
typedef int (*inquiry)(PyObject *);
typedef int (*coercion)(PyObject **, PyObject **);
typedef int(*intobjargproc)(PyObject *, int, PyObject *);
typedef int(*intintobjargproc)(PyObject *, int, int, PyObject *);
typedef int(*objobjargproc)(PyObject *, PyObject *, PyObject *);
#endif

#endif /* _ASAPPYTHON_H */
