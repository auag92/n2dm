// -*- C++ -*-
// EMTParameterProviderInterface.cpp: Python interface to the
// EMTParameterProvider objects.
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

#include "AsapPython.h"

#include "EMTParameterProviderInterface.h"
#include "ExceptionInterface.h"
#include "PythonConversions.h"
#include "Templates.h"
#include "EMTDefaultParameterProvider.h"
#include "EMTRasmussenParameterProvider.h"
#include "EMTPythonParameterProvider.h"
//#define ASAPDEBUG
#include "Debug.h"

// All EMT parameter provider are of the type PyAsap_EMTParamProvType,
// but the cobj pointer may point to C++ objects of different type,
// depending on how the PyAsap_EMTParamProvObject was created.
// Different factory functions provide PyAsap_EMTParamProvObject
// objects with different contents.  Objects created the standard way
// use a EMTPythonParameterProvider (XXX NOT IMPLEMENTED YET) which
// should be subclassed in Python.

// The PyAsap_EMTParamProvObject is defined in Basics/EMTParameterProvider.h

namespace ASAPSPACE {

PyTypeObject PyAsap_EMTParamProvType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "asapserial3.EMTParameterProvider",
  sizeof(PyAsap_EMTParamProvObject),
  // The rest are initialized by name for reliability.
};

static char EMTParamProv_Docstring[] =
  "Generic EMT parameter provider.\n";
  
char EMTDefaultParamProv_Docstring[] =
  "Provides default EMT parameters.\n";

char EMTRasmussenParamProv_Docstring[] =
  "Provides the alternative EMT parameters optimized by Torben Rasmussen.\n";


// The default allocator, it created an EMTPythonParameterProvider.
// The following factory function creates other kinds of parameter
// providers.
static PyObject *PyAsap_NewParameterProvider(PyTypeObject *type,
					     PyObject *args, PyObject *kwargs)
{
  DEBUGPRINT;
  // Ignore the arguments.
  PyAsap_EMTParamProvObject *self =
    (PyAsap_EMTParamProvObject *)type->tp_alloc(type, 0);
  if (self == NULL)
    return NULL;
  self->weakrefs = NULL;
  self->cobj = new EMTPythonParameterProvider((PyObject *) self);
  if (self->cobj == NULL)
    {
      CHECKREF(self);
      Py_DECREF(self);
      return NULL;
    }
  return (PyObject *) self;
}

static int PyAsap_InitParameterProvider(PyObject *self,
					PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "", kwlist))
    return -1;
  return 0;
}

// Factory function for creating a default parameter provider
PyObject *PyAsap_EMTDefaultParamProvNew(PyObject *noself, PyObject *noargs)
{
  PyAsap_EMTParamProvObject *self;
  DEBUGPRINT;
  
  self = PyObject_New(PyAsap_EMTParamProvObject, &PyAsap_EMTParamProvType);
  if (self == NULL)
    return NULL;
  self->weakrefs = NULL;
  self->cobj = new EMTDefaultParameterProvider();
  if (self->cobj == NULL)
    {
      CHECKREF(self);
      Py_DECREF(self);
      return NULL;
    }
  return (PyObject *) self;
}

// Factory function for creating a Rasmussen parameter provider
PyObject *PyAsap_EMTRasmussenParamProvNew(PyObject *noself, PyObject *noargs)
{
  PyAsap_EMTParamProvObject *self;
  DEBUGPRINT;
  
  self = PyObject_New(PyAsap_EMTParamProvObject, &PyAsap_EMTParamProvType);
  if (self == NULL)
    return NULL;
  self->weakrefs = NULL;
  self->cobj = new EMTRasmussenParameterProvider();
  if (self->cobj == NULL)
    {
      CHECKREF(self);
      Py_DECREF(self);
      return NULL;
    }
  return (PyObject *) self;
}

static PyObject *PyAsap_ParamProvMaxCutoff(PyAsap_EMTParamProvObject *self,
                                           PyObject *noargs)
{
  double cutoff = self->cobj->GetMaxListCutoffDistance();
  return PyFloat_FromDouble(cutoff);
}


static PyMethodDef PyAsap_ParamProvMethods[] = {
  {"get_max_cutoff_beforeinit", (PyCFunction)PyAsap_ParamProvMaxCutoff,
   METH_NOARGS, "The maximal value of the neighbor list cutoff that may later be set."},
   {NULL}  // Sentinel
};


int PyAsap_InitEMTParameterProviderInterface(PyObject *module)
{
  PyAsap_EMTParamProvType.tp_dealloc =
    PyAsap_Dealloc<PyAsap_EMTParamProvObject>;
  PyAsap_EMTParamProvType.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE;
  PyAsap_EMTParamProvType.tp_repr =
    PyAsap_Representation<PyAsap_EMTParamProvObject>;
  PyAsap_EMTParamProvType.tp_new = PyAsap_NewParameterProvider; 
  PyAsap_EMTParamProvType.tp_init = PyAsap_InitParameterProvider; 
  PyAsap_EMTParamProvType.tp_doc = EMTParamProv_Docstring;
  PyAsap_EMTParamProvType.tp_methods = PyAsap_ParamProvMethods;

  if (PyType_Ready(&PyAsap_EMTParamProvType) < 0)
    return -1;
  PyModule_AddObject(module, "EMTParameters",
		     (PyObject *) &PyAsap_EMTParamProvType);
  return 0;
}

} // end namespace
