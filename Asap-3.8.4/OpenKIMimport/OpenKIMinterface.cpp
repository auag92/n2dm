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

#include "OpenKIMinterface.h"
#include "OpenKIMinfo.h"
#include "OpenKIMcalculator.h"
#include "Templates.h"
#include "ExceptionInterface.h"
#include "PotentialInterface.h"
//#define ASAPDEBUG
#include "Debug.h"
#include <cstdlib>


namespace ASAPSPACE {

//////////////////////////////
//
//  OpenKIMinfo object
//
//////////////////////////////

static PyTypeObject PyAsap_OpenKIMinfoType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "asapserial3.OpenKIMinfo",
  sizeof(PyAsap_OpenKIMinfoObject),
  // The rest are initialized by name for reliability.
};

char OpenKIMinfo_Docstring[] =
    "Informational about an OpenKIM model\n\n\
  Parameters:\n\
    name: The name of the model.\n\n";

PyObject *PyAsap_NewOpenKIMinfo(PyObject *noself, PyObject *args,
                                PyObject *kwargs)
{
  static char *kwlist[] = {"name", NULL};
  DEBUGPRINT;
  const char *name = NULL;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,  "s:NewOpenKIMinfo",
      kwlist, &name))
    return NULL;
  try {
      DEBUGPRINT;
      PyAsap_OpenKIMinfoObject *self = PyAsap_NewOpenKIMinfoObject(name);
      DEBUGPRINT;
      return (PyObject *) self;
  }
  CATCHEXCEPTION;
}

PyAsap_OpenKIMinfoObject *PyAsap_NewOpenKIMinfoObject(const char *name)
{
  DEBUGPRINT;
  PyAsap_OpenKIMinfoObject *self;

  self = PyObject_NEW(PyAsap_OpenKIMinfoObject,
                      &PyAsap_OpenKIMinfoType);
  if (self == NULL)
    throw AsapError("Failed to create OpenKIMinfo object ");

  self->weakrefs = NULL;
  self->cobj = new OpenKIMinfo(name);
  DEBUGPRINT;
  if (self->cobj == NULL)
    {
      CHECKREF(self);
      Py_DECREF(self);
      throw AsapError("Failed to create a new OpenKIMinfo object.");
    }
  DEBUGPRINT;
  return self;
}


static PyObject *PyAsap_OpenKIMinfoGetSupportedTypes(PyAsap_OpenKIMinfoObject
                                                     *self, PyObject *noargs)
{
  DEBUGPRINT;
  std::vector<const char *> symbols;
  self->cobj->GetSupportedSymbols(symbols);
  int num = symbols.size();
  PyObject *result = PyTuple_New(num);
  if (result == NULL)
    return NULL;
  for (int i = 0; i < num; i++)
    {
      PyTuple_SET_ITEM(result, i, PyString_FromString(symbols[i]));
    }
  DEBUGPRINT;
  return result;
}

static PyObject *PyAsap_OpenKIMinfoGetSupportedNbList(PyAsap_OpenKIMinfoObject
                                                   *self, PyObject *noargs)
{
  DEBUGPRINT;
  PyObject *result;
  try {
      const std::vector<const char *> &pbc_types = self->cobj->GetSupportedNBC();

      int n = pbc_types.size();
      result = PyTuple_New(n);
      if (result == NULL)
        return NULL;
      for (int i = 0; i < n; i++)
        {
          PyTuple_SET_ITEM(result, i, PyString_FromString(pbc_types[i]));
        }
  }
  CATCHEXCEPTION;
  DEBUGPRINT;
  return result;
}

static PyObject *PyAsap_OpenKIMinfoGetSupportedAccess(PyAsap_OpenKIMinfoObject
                                                     *self, PyObject *noargs)
{
  DEBUGPRINT;
  PyObject *result;
  try {
      const std::vector<const char *> &pbc_types = self->cobj->GetSupportedAccess();

      int n = pbc_types.size();
      result = PyTuple_New(n);
      if (result == NULL)
        return NULL;
      for (int i = 0; i < n; i++)
        {
          PyTuple_SET_ITEM(result, i, PyString_FromString(pbc_types[i]));
        }
  }
  CATCHEXCEPTION;
  DEBUGPRINT;
  return result;
}

static PyObject *PyAsap_OpenKIMinfoGetAPIindex(PyAsap_OpenKIMinfoObject
                                               *self, PyObject *args,
                                               PyObject *kwargs)
{
  DEBUGPRINT;
  static char *kwlist[] = {"name", NULL};

  const char *name = NULL;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,  "s:get_API_index",
      kwlist, &name))
    return NULL;
  try {
      int result = self->cobj->GetAPIindex(name);
      return Py_BuildValue("i", result);
  }
  CATCHEXCEPTION;
  DEBUGPRINT;
}

static PyMethodDef PyAsap_OpenKIMinfoMethods[] = {
  {"get_supported_types", (PyCFunction) PyAsap_OpenKIMinfoGetSupportedTypes,
   METH_NOARGS, "Get the supported elements as a tuple of chemical symbols"},
  {"get_supported_nblist", (PyCFunction) PyAsap_OpenKIMinfoGetSupportedNbList,
   METH_NOARGS, "Get the supported boundary condition / neighbor list methods"},
  {"get_supported_access", (PyCFunction) PyAsap_OpenKIMinfoGetSupportedAccess,
    METH_NOARGS, "Get the supported neighbor list access methods"},
  {"get_API_index", (PyCFunction) PyAsap_OpenKIMinfoGetAPIindex,
   METH_VARARGS | METH_KEYWORDS, "Get API index of property - or -1 if not supported"},
  {NULL}
};



//////////////////////////////
//
//  OpenKIMcalculator object
//
//////////////////////////////

static PyTypeObject PyAsap_OpenKIMcalculatorType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "asapserial3.OpenKIMcalculator",
  sizeof(PyAsap_PotentialObject),
  // The rest are initialized by name for reliability.
};

char OpenKIMcalculator_Docstring[] =
    "Internal interface to an OpenKIM model\n\n\
  Parameters:\n\
    descr: OpenKIM descriptor string.\n\
    name: The name of the model.\n\n";


static int PyAsap_OpenKIMcalculatorInit(PyAsap_PotentialObject *self, PyObject *args,
                                        PyObject *kwargs)
{
  if (PyAsap_PotentialType.tp_init((PyObject *)self, args, kwargs) < 0)
    return -1;
  try
    {
      self->cobj = new OpenKIMcalculator((PyObject *) self);
      self->orig_cobj = self->cobj;
    }
  catch (AsapError &e)
    {
      string msg = e.GetMessage();
      PyErr_SetString(PyAsap_ErrorObject, msg.c_str());
      return -1;
    }
  catch (AsapPythonError &e)
    {
      return -1;
    }
  if (self->cobj == NULL)
    return -1;
  return 0;
}

static PyObject *PyAsap_OpenKIMcalcInitialize(PyAsap_PotentialObject *self, PyObject *args,
                                              PyObject *kwargs)
{
  static char *kwlist[] = {"descr", "name", NULL};
  const char *descr = NULL;
  const char *name = NULL;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,  "ss:NewOpenKIMcalculator",
                                   kwlist, &descr, &name))
    return NULL;
  OpenKIMcalculator *cobj = dynamic_cast<OpenKIMcalculator *>(self->orig_cobj);
  assert(cobj != NULL);
  try {
    cobj->Initialize(descr, name);
  }
  catch (AsapError &e)
    {
      string msg = e.GetMessage();
      PyErr_SetString(PyAsap_ErrorObject, msg.c_str());
      return NULL;
    }
  catch (AsapPythonError &e)
    {
      return NULL;
    }
  Py_RETURN_NONE;
}


static PyObject *PyAsap_OpenKIMcalcGetNBCmethod(PyAsap_PotentialObject
                                               *self, PyObject *noargs)
{
  OpenKIMcalculator *cobj = dynamic_cast<OpenKIMcalculator *>(self->orig_cobj);
  assert(cobj != NULL);
  const char *pbc;
  try {
      pbc = cobj->GetNBCmethod();
  }
  catch (AsapError &e)
    {
      string msg = e.GetMessage();
      PyErr_SetString(PyAsap_ErrorObject, msg.c_str());
      return NULL;
    }
  catch (AsapPythonError &e)
    {
      return NULL;
    }
  PyObject *result = Py_BuildValue("s", pbc);
  return result;
}

static PyObject *PyAsap_OpenKIMcalcPleaseAlloc(PyAsap_PotentialObject *self,
                                               PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"quantity", "alloc", NULL};
  const char *quantity = NULL;
  int alloc = 0;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,  "si:please_allocate", kwlist,
                                   &quantity, &alloc))
    return NULL;
  OpenKIMcalculator *cobj = dynamic_cast<OpenKIMcalculator *>(self->orig_cobj);
  assert(cobj != NULL);
  try
  {
    cobj->PleaseAllocate(quantity, alloc);
  }
  catch (AsapError &e)
    {
      string msg = e.GetMessage();
      PyErr_SetString(PyAsap_ErrorObject, msg.c_str());
      return NULL;
    }
  catch (AsapPythonError &e)
    {
      return NULL;
    }
  Py_RETURN_NONE;
}

static PyObject *PyAsap_OpenKIMcalcSetTranslation(PyAsap_PotentialObject *self,
                                                 PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"translation", NULL};
  PyObject *translation = NULL;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,  "O!:set_translation", kwlist,
                                   &PyDict_Type, &translation))
    return NULL;
  OpenKIMcalculator *cobj = dynamic_cast<OpenKIMcalculator *>(self->orig_cobj);
  assert(cobj != NULL);
  cobj->ClearTranslation();
  PyObject *key, *value;
  Py_ssize_t i = 0;
  while (PyDict_Next(translation, &i, &key, &value))
    {
      int z = PyInt_AsLong(key);
      int code = PyInt_AsLong(value);
      if (z == -1 || code == -1)
        return PyErr_Format(PyExc_ValueError,
            "Illegal translation %i -> %i (or non-integer type)", z, code);
      cobj->AddTranslation(z, code);
    }
  Py_RETURN_NONE;
}


static PyObject *PyAsap_OpenKIMcalcGetSupportedTypes(PyAsap_PotentialObject
                                                     *self, PyObject *noargs)
{
  DEBUGPRINT;
  std::vector<const char *> symbols;
  OpenKIMcalculator *cobj = dynamic_cast<OpenKIMcalculator *>(self->orig_cobj);
  assert(cobj != NULL);
  cobj->GetSupportedSymbols(symbols);
  int num = symbols.size();
  PyObject *result = PyTuple_New(num);
  if (result == NULL)
    return NULL;
  for (int i = 0; i < num; i++)
    {
      PyTuple_SET_ITEM(result, i, PyString_FromString(symbols[i]));
    }
  DEBUGPRINT;
  return result;
}

static PyObject *PyAsap_OpenKIMcalcGetTypeCode(PyAsap_PotentialObject
                                               *self, PyObject *args,
                                               PyObject *kwargs)
{
  DEBUGPRINT;
  static char *kwlist[] = {"symbol", NULL};
  const char *symbol = NULL;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,  "s:get_type_code",
      kwlist, &symbol))
    return NULL;
  OpenKIMcalculator *cobj = dynamic_cast<OpenKIMcalculator *>(self->orig_cobj);
  assert(cobj != NULL);
  try {
      int result = cobj->GetParticleTypeCode(symbol);
      return Py_BuildValue("i", result);
  }
  CATCHEXCEPTION;
  DEBUGPRINT;
}

static PyMethodDef PyAsap_OpenKIMcalculatorMethods[] = {
  {"_initialize", (PyCFunction) PyAsap_OpenKIMcalcInitialize,
    METH_VARARGS | METH_KEYWORDS, "Initialize the OpenKIM model by name and descriptor."},
  {"get_NBC_method", (PyCFunction) PyAsap_OpenKIMcalcGetNBCmethod,
   METH_NOARGS, "Get name of the neighborlist method"},
  {"please_allocate", (PyCFunction) PyAsap_OpenKIMcalcPleaseAlloc,
   METH_VARARGS | METH_KEYWORDS, "Enable specific property"},
  {"set_translation", (PyCFunction) PyAsap_OpenKIMcalcSetTranslation,
   METH_VARARGS | METH_KEYWORDS, "Set Z->typecode translation table."},
  {"_use_imageatoms", (PyCFunction) PyAsap_PotentialUseImageAtoms,
   METH_NOARGS, PyAsap_PotentialUseImageAtoms_Docstring},
  {"get_supported_types", (PyCFunction) PyAsap_OpenKIMcalcGetSupportedTypes,
   METH_NOARGS, "Get the supported elements as a tuple of chemical symbols"},
  {"get_type_code", (PyCFunction) PyAsap_OpenKIMcalcGetTypeCode,
   METH_VARARGS | METH_KEYWORDS, "Get type code of an element"},
  {NULL}
};



//////////////////////////////
//
//  Module initialization
//
//////////////////////////////


int PyAsap_InitOpenKIMInterface(PyObject *module)
{

  InitPotentialType(PyAsap_OpenKIMcalculatorType);
  PyAsap_OpenKIMcalculatorType.tp_init = (initproc) PyAsap_OpenKIMcalculatorInit;
  PyAsap_OpenKIMcalculatorType.tp_doc = OpenKIMcalculator_Docstring;
  PyAsap_OpenKIMcalculatorType.tp_methods = PyAsap_OpenKIMcalculatorMethods;
  if (PyType_Ready(&PyAsap_OpenKIMcalculatorType) < 0)
    return -1;
  Py_INCREF(&PyAsap_OpenKIMcalculatorType);
  PyModule_AddObject(module, "OpenKIMcalculator", (PyObject *) &PyAsap_OpenKIMcalculatorType);

  PyAsap_OpenKIMinfoType.tp_new = NULL;  // Use factory functions
  PyAsap_OpenKIMinfoType.tp_dealloc =
    PyAsap_Dealloc<PyAsap_OpenKIMinfoObject>;
  PyAsap_OpenKIMinfoType.tp_flags = Py_TPFLAGS_DEFAULT;
  PyAsap_OpenKIMinfoType.tp_methods = PyAsap_OpenKIMinfoMethods;
  PyAsap_OpenKIMinfoType.tp_repr =
    PyAsap_Representation<PyAsap_OpenKIMinfoObject>;
  PyAsap_OpenKIMinfoType.tp_doc = OpenKIMinfo_Docstring;
  if (PyType_Ready(&PyAsap_OpenKIMinfoType) < 0)
    return -1;

  return 0;
}

} // namespace
