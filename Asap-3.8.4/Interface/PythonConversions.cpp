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

#include "PythonConversions.h"

namespace ASAPSPACE {

int PyAsap_VectorIntFromArray(vector<int> &to, PyObject *from)
{
  // It is impossible to predict if we get an array of 32 or 64 bit integers.  Cast to large
  // type and convert down (this function is never critical for performance).
  PyArrayObject *array = (PyArrayObject *) PyArray_ContiguousFromObject(from, NPY_LONG, 1, 1);
  if (array == NULL)
    {
      PyErr_SetString(PyExc_TypeError,
		      "Not compatible with 1D array of integers."); 
      return -1;
    };
  int n = PyArray_DIM(array, 0);
  to.resize(n);
  npy_long *p = (npy_long *) PyArray_DATA(array);
  for (int i = 0; i < n; i++)
    to[i] = (int) p[i];
  CHECKREF(array);
  Py_DECREF(array);
  return 0;
}

int PyAsap_VectorDoubleFromArray(vector<double> &to, PyObject *from)
{
  PyArrayObject *array = (PyArrayObject *) PyArray_ContiguousFromObject(from, NPY_DOUBLE, 1, 1);
  if (array == NULL)
    {
      PyErr_SetString(PyExc_TypeError,
		      "Not compatible with 1D array of double."); 
      return -1;
    };
  to.resize(PyArray_DIM(array, 0));
  memcpy(&to[0], PyArray_DATA(array), PyArray_DIM(array, 0)*sizeof(double));
  CHECKREF(array);
  Py_DECREF(array);
  return 0;
}

int PyAsap_TinyMatrixDoubleFromArray(TinyMatrix<double> &to, PyObject *from)
{
  PyArrayObject *array = (PyArrayObject *) PyArray_ContiguousFromObject(from, NPY_DOUBLE, 2, 2);
  if (array == NULL)
    {
      PyErr_SetString(PyExc_TypeError,
		      "Not compatible with 2D array of double."); 
      return -1;
    };
  to.Allocate(PyArray_DIM(array, 0), PyArray_DIM(array, 1));
  int n = PyArray_DIM(array, 0) * PyArray_DIM(array, 1);
  memcpy(to[0], PyArray_DATA(array), n*sizeof(double));
  Py_DECREF(array);
  return 0;
}

int PyAsap_SetIntFromArray(set<int> &to, PyObject *from)
{
  to.clear();
  PyArrayObject *array = (PyArrayObject *) PyArray_ContiguousFromObject(from, NPY_INT, 1, 1);
  if (array == NULL)
    {
      PyErr_SetString(PyExc_TypeError,
		      "Not compatible with 1D array of integers."); 
      return -1;
    };
  for (int i = 0; i < PyArray_DIM(array, 0); i++)
    to.insert(*((int *)PyArray_GETPTR1(array, i)));
  CHECKREF(array);
  Py_DECREF(array);
  return 0;
}

} // end namespace
