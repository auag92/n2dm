// -*- C++ -*-
// ToolsInterface.cpp: Python interface to simple Tools functions.
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

#include "AsapPython.h"
#include "ToolsInterface.h"
#include "CoordinationNumbers.h"
#include "ExceptionInterface.h"
#include "PythonConversions.h"
#include "Atoms.h"
#include "CNA.h"
#include "FullCNA.h"
#include "Templates.h"

namespace ASAPSPACE {

PyObject *PyAsap_CoordinationNumbers(PyObject *noself, PyObject *args)
{
  // Arguments are atoms and rCut.
  PyObject *py_atoms;
  double rCut;
  if (!PyArg_ParseTuple(args, "Od:CoordinationNumber", &py_atoms, &rCut))
    return NULL;
  PyObject *py_result = NULL;
  try {
    vector<int> result;
    CHECKNOASAPERROR;
    CoordinationNumbers(py_atoms, rCut, result);
    PROPAGATEASAPERROR;
    py_result = PyAsap_ArrayFromVectorInt(result);
  }
  CATCHEXCEPTION;
  return py_result;
}

PyObject *PyAsap_RestrictedCNA(PyObject *noself, PyObject *args)
{
  // Arguments are atoms and rCut.
  PyObject *py_atoms;
  double rCut;
  if (!PyArg_ParseTuple(args, "Od:RestrictedCNA", &py_atoms, &rCut))
    return NULL;
  PyObject *py_result;
  try {
    vector<char> cna;
    CHECKNOASAPERROR;
    CNA(py_atoms, rCut, cna);
    PROPAGATEASAPERROR;
    py_result = PyAsap_ArrayFromVectorChar(cna);
  }
  CATCHEXCEPTION;
  return py_result;
}


// The FullCNA object

static PyTypeObject PyAsap_FullCNAType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "asapserial3.FullCNA",
  sizeof(PyAsap_FullCNAObject),
  // The rest are initialized by name for reliability.
};

static char FullCNA_Docstring[] = "FullCNA object (internal use only).\n";

static int PyAsap_FullCNAInit(PyAsap_FullCNAObject *self, PyObject *args,
                              PyObject *kwargs)
{
  static char *kwlist[] = {"atoms", "cutoff", NULL};

  PyObject *atoms;
  double cutoff;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,  "Od:FullCNA",
      kwlist, &atoms, &cutoff))
    return -1;
  if (cutoff <= 0.0)
    {
      PyErr_SetString(PyExc_ValueError,
          "FullCNA: Cutoff must be greater than zero.");
      return -1;
    }
  assert(self->cobj == NULL);
  try
  {
      CHECKNOASAPERROR;
      self->cobj = new FullCNA(atoms, cutoff);
      PROPAGATEASAPERROR;
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

static PyObject *PyAsap_FullCNA_GetRawCNA(PyAsap_FullCNAObject *self, PyObject *noargs)
{
  assert(self->cobj != NULL);
  return self->cobj->GetRawCNA();
}

static PyObject *PyAsap_FullCNA_GetPerAtomCNA(PyAsap_FullCNAObject *self, PyObject *noargs)
{
  assert(self->cobj != NULL);
  return self->cobj->GetPerAtomCNA();
}

static PyObject *PyAsap_FullCNA_GetTotalCNA(PyAsap_FullCNAObject *self, PyObject *noargs)
{
  assert(self->cobj != NULL);
  return self->cobj->GetTotalCNA();
}

static PyMethodDef PyAsap_FullCNAMethods[] = {
    {"get_per_atom_cna", (PyCFunction)PyAsap_FullCNA_GetPerAtomCNA,
        METH_NOARGS,  "Return the per-atom ('normal') CNA data."},
    {"get_total_cna", (PyCFunction)PyAsap_FullCNA_GetTotalCNA,
        METH_NOARGS,  "Return the total CNA data."},
    {"get_raw_cna", (PyCFunction)PyAsap_FullCNA_GetRawCNA,
        METH_NOARGS,  "Return the raw CNA data."},
    {NULL}  // Sentinel
};

int PyAsap_InitToolsInterface(PyObject *module)
{
  PyAsap_FullCNAType.tp_new = PyType_GenericNew;
  PyAsap_FullCNAType.tp_dealloc = PyAsap_Dealloc<PyAsap_FullCNAObject>;
  PyAsap_FullCNAType.tp_flags = Py_TPFLAGS_DEFAULT;
  PyAsap_FullCNAType.tp_methods = PyAsap_FullCNAMethods;
  PyAsap_FullCNAType.tp_repr = PyAsap_Representation<PyAsap_FullCNAObject>;
  PyAsap_FullCNAType.tp_init = (initproc) PyAsap_FullCNAInit;
  PyAsap_FullCNAType.tp_doc = FullCNA_Docstring;
  if (PyType_Ready(&PyAsap_FullCNAType) < 0)
    return -1;
  Py_INCREF(&PyAsap_FullCNAType);
  PyModule_AddObject(module, "FullCNA", (PyObject *) &PyAsap_FullCNAType);
  return 0;
}

} // end namespace
