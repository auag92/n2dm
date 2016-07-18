// -*- C++ -*-
// EMTPythonParameterProvider.cpp:  Get EMT parameters from a Python class.
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

#include "EMTPythonParameterProvider.h"
#include "Asap.h"
#include "Exception.h"
// #define ASAPDEBUG
#include "Debug.h"
#include <iostream>

EMTPythonParameterProvider::EMTPythonParameterProvider(PyObject *self)
{
  DEBUGPRINT;
  this->self = self;
  // NB: No INCREF!  This is not a reference, but a pointer internal
  // to the Python object.
}

static double extract(PyObject *dict, const char *name)
{
  DEBUGPRINT;
  PyObject *x = PyDict_GetItemString(dict, name);  // Borrowed ref!
  if (x == NULL || !PyFloat_Check(x))
    {
      Py_DECREF(dict);  // We also abort GetNewParameters()
      throw AsapError("EMT parameter dictionary had no (or non-float) element ")
	<< name;
    }
  return PyFloat_AS_DOUBLE(x);
}

static int extract_int(PyObject *dict, const char *name)
{
  DEBUGPRINT;
  PyObject *x = PyDict_GetItemString(dict, name);
  if (x == NULL || !PyInt_Check(x))
    {
      Py_DECREF(dict);  // We also abort GetNewParameters()
      throw AsapError("EMT parameter dictionary had no (or non-int) element ")
	<< name;
    }
  return (int) PyInt_AS_LONG(x);
}

static char *extract_str(PyObject *dict, const char *name)
{
  DEBUGPRINT;
  PyObject *x = PyDict_GetItemString(dict, name);
  if (x == NULL || !PyString_Check(x))
    {
      Py_DECREF(dict);  // We also abort GetNewParameters()
      throw
	AsapError("EMT parameter dictionary had no (or non-string) element ")
	<< name;
    }
  // "Intern" the string (the element name).  Interning the string
  // assures that the Python string object never goes away, but does
  // not introduce a memory leak as multiple identical strings are mapped
  // to the same string (no leak as long as the number of distinct element
  // names remains limited!).  The alternative, copying the string,
  // introduce a tiny leak as the pointer cannot be deallocated (usually
  // it refers to static memory).
  PyString_InternInPlace(&x);
  char *str = PyString_AsString(x);
  // Note that by failing to DECREF the string here we make the interned
  // string immortal, in contrast to what the docs say this does not
  // happen automatically when interning.
  // It might be better to make it immortal through the C API (undocumented!)
  // or to keep it alive in a dictionary for as long as this parameter provider
  // is alive.
  return str;
}
      
emt_parameters *EMTPythonParameterProvider::GetNewParameters(int element)
{
  DEBUGPRINT;
  PyObject *newparam = PyObject_CallMethod(self, "get_parameters", "(i)",
					   element);
  if (newparam == NULL)
    throw AsapPythonError();
  if (!PyDict_Check(newparam))
    throw AsapError("get_parameters did not return a dictionary");
  
  DEBUGPRINT;
  emt_parameters *p = new emt_parameters;
  p->e0 = extract(newparam, "E0");
  p->seq = extract(newparam, "S0");
  p->neq = extract(newparam, "n0");
  p->V0 = extract(newparam, "V0");
  p->eta2 = extract(newparam, "eta2");
  p->kappa = extract(newparam, "kappa");
  p->lambda = extract(newparam, "lambda");
  p->mass = extract(newparam, "mass");
  p->Z = extract_int(newparam, "Z");
  p->name = extract_str(newparam, "name");
  p->invmass = 1.0 / p->mass;
  p->gamma1 = p->gamma2 = 0.0;
  p->lengthscale = 0.0;  // Remove
  assert(element == p->Z);

  DEBUGPRINT;
  Py_DECREF(newparam);
  return p;
}


void EMTPythonParameterProvider::CalcGammaEtc()
{
  int n = params.size();
  PyObject *data = PyObject_CallMethod(self, "get_gammas_etc", "");
  if (data == NULL)
    throw AsapPythonError();
  if (!PyTuple_Check(data))
    throw AsapError("get_gammas_etc did not return a tuple");
  PyObject *py_gammas;
  PyArrayObject *py_chi;
  if (!PyArg_Parse(data, "((ddd)OO!)", &cutoff, &cutslope, &listcutofffactor,
		   &py_gammas, &PyArray_Type, (PyObject *) &py_chi))
    throw AsapPythonError();
  if (!PyList_Check(py_gammas) || PyList_GET_SIZE(py_gammas) != n)
    {
      Py_DECREF(data);
      throw AsapError("get_gammas_etc returned improper gammas.");
    }
  if (PyArray_NDIM(py_chi) != 2
      || PyArray_DIM(py_chi, 0) != n
      || PyArray_DIM(py_chi, 1) != n
      || PyArray_TYPE(py_chi) != NPY_DOUBLE
      || !PyArray_ISCARRAY_RO(py_chi))
    {
      Py_DECREF(data);
      throw AsapError("get_gammas_etc returned improper chi.");
    }
  for (int i = 0; i < n; i++)
    {
      PyObject *gam = PyList_GET_ITEM(py_gammas, i);
      if (gam == NULL || !PyArg_Parse(gam, "(dd)", &(params[i]->gamma1),
				      &(params[i]->gamma2)))
	{
	  Py_DECREF(data);
	  throw AsapError("Failed to parse gammas - item ") << i;
	}
    }
  chi = new TinyDoubleMatrix(n,n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      {
	(*chi)[i][j] = *((double*) PyArray_GETPTR2(py_chi, i, j));
      }
  Py_DECREF(data);
#if 0
  std::cerr << "cutoff " << cutoff << std::endl;
  std::cerr << "Beta " << Beta << std::endl;
  std::cerr << "cutslope " << cutslope << std::endl;
  std::cerr << "shell0 " << shell0 << std::endl;
  std::cerr << "shell1 " << shell1 << std::endl;
#endif
}

double EMTPythonParameterProvider::GetMaxListCutoffDistance()  // Max value, useful before initialization.
{
  PyObject *data = PyObject_CallMethod(self, "get_maximal_cutoff", "");
  if (data == NULL)
    throw AsapPythonError();
  if (!PyFloat_Check(data))
    throw AsapError("get_maximal_cutoff did not return a float");
  double result = PyFloat_AsDouble(data);
  Py_DECREF(data);
  return result;
}

