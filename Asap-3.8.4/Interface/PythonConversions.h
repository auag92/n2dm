// -*- C++ -*-
// PythonConversions.h: Helper functions creating Python objects.
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

#ifndef _PYTHONCONVERSIONS_H
#define _PYTHONCONVERSIONS_H

#include "AsapPython.h"
#include "Asap.h"
#include "Vec.h"
#include "SymTensor.h"
#include "TinyMatrix.h"
#include <vector>
#include <set>
using std::vector;
using std::set;

namespace ASAPSPACE {

inline PyObject *PyAsap_ArrayFromVectorDouble(const vector<double> &data)
{
  npy_intp size = data.size();
  PyObject *res = PyArray_SimpleNew(1, &size, NPY_DOUBLE);
  if (res == NULL)
    return NULL;
  assert(PyArray_NBYTES((PyArrayObject *) res) ==  size*sizeof(double));
  memcpy(PyArray_DATA((PyArrayObject *) res), &data[0], size*sizeof(double));
  return res;
}

inline PyObject *PyAsap_ArrayFromVectorInt(const vector<int> &data)
{
  npy_intp size = data.size();
  PyObject *res = PyArray_SimpleNew(1, &size, NPY_INT);
  if (res == NULL)
    return NULL;
  assert(PyArray_NBYTES((PyArrayObject *) res) ==  size*sizeof(int));
  if (size > 0)
    memcpy(PyArray_DATA((PyArrayObject *) res), &data[0], size*sizeof(int));
  return res;
}

inline PyObject *PyAsap_ArrayFromVectorChar(const vector<char> &data)
{
  npy_intp size = data.size();
  PyObject *res = PyArray_SimpleNew(1, &size, NPY_BYTE);
  if (res == NULL)
    return NULL;
  assert(PyArray_NBYTES((PyArrayObject *) res) ==  size*sizeof(char));
  memcpy(PyArray_DATA((PyArrayObject *) res), &data[0], size*sizeof(char));
  return res;
}

inline PyObject *PyAsap_ArrayFromVectorVec(const vector<Vec> &data)
{
  npy_intp size[2];
  size[0] = data.size();
  size[1] = 3;
  PyObject *res = PyArray_SimpleNew(2, size, NPY_DOUBLE);
  if (res == NULL)
    return NULL;
  assert(PyArray_NBYTES((PyArrayObject *) res) ==  size[0]*sizeof(Vec));
  memcpy(PyArray_DATA((PyArrayObject *) res), &data[0], size[0]*sizeof(Vec));
  return res;
}

inline PyObject *PyAsap_ArrayFromVectorSymTensor(const vector<SymTensor> &data)
{
  npy_intp size[2];
  size[0] = data.size();
  size[1] = 6;
  PyObject *res = PyArray_SimpleNew(2, size, NPY_DOUBLE);
  if (res == NULL)
    return NULL;
  assert(PyArray_NBYTES((PyArrayObject *) res) ==  size[0]*sizeof(SymTensor));
  memcpy(PyArray_DATA((PyArrayObject *) res), &data[0], size[0]*sizeof(SymTensor));
  return res;
}

int PyAsap_VectorIntFromArray(vector<int> &to, PyObject *from);

int PyAsap_VectorDoubleFromArray(vector<double> &to, PyObject *from);

int PyAsap_TinyMatrixDoubleFromArray(TinyMatrix<double> &to, PyObject *from);

int PyAsap_SetIntFromArray(set<int> &to, PyObject *from);

} // end namespace

#endif // _PYTHONCONVERSIONS_H
