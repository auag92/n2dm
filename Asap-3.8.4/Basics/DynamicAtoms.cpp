// DynamicAtoms.cpp  --  Access the atoms from a dynamics object
// -*- c++ -*-
//
// Copyright (C) 2001-2011 Jakob Schiotz and Center for Individual
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

#include "DynamicAtoms.h"
#include "Exception.h"
#include <vector>
using std::vector;

DynamicAtoms::DynamicAtoms(PyObject *py_atoms)
{
    atoms = py_atoms;
    arrays = PyObject_GetAttrString(atoms, "arrays");
    if (arrays == NULL)
        throw AsapError("Atoms object has no 'arrays' attribute");
    if (!PyDict_Check(arrays))
    {
        Py_DECREF(arrays);
        throw AsapError("Atoms.arrays is not a dictionary!");
    }
    Py_INCREF(atoms);
    // Now, we need to get the masses of all atoms.
    // They are extracted from ase.data
    PyObject *datamodule = PyImport_ImportModule("ase.data");
    if (datamodule == NULL)
        throw AsapPythonError();
    PyArrayObject *atomic_masses = AsPyArray(PyObject_GetAttrString(datamodule,
        "atomic_masses"));
    Py_DECREF(datamodule);
    if (atomic_masses == NULL)
        throw AsapPythonError();
    if (PyArray_NDIM(atomic_masses) != 1           // one-dimensional
	  || PyArray_TYPE(atomic_masses) != NPY_DOUBLE  // array of double
	  || !PyArray_ISCARRAY_RO(atomic_masses))       // Contiguous etc.
    {
        Py_DECREF(atomic_masses);
        throw AsapError("ase.data.atomic_masses has unexpected type");
    }
    int n = PyArray_DIM(atomic_masses, 0);
    masses.resize(n);
    invmasses.resize(n);
    double *m = (double *) PyArray_DATA(atomic_masses);
    for (int i = 0; i < n; i++)
    {
        masses[i] = m[i];
        invmasses[i] = 1.0 / m[i];
    }
    Py_DECREF(atomic_masses);
}

DynamicAtoms::~DynamicAtoms()
{
    Py_DECREF(arrays);
    Py_DECREF(atoms);
}

Vec *DynamicAtoms::GetVecData(const char *name)
{
    // Get a borrowed reference to the array.  We do not keep 
    // a reference ourselves as the array cannot go away before
    // we are done using it.
    PyArrayObject *py_data = AsPyArray(PyDict_GetItemString(arrays, name));
    if (py_data == NULL)
        throw AsapError("Atoms.arrays has no ") << name;
    if (PyArray_NDIM(py_data) != 2           // two-dimensional
	  || PyArray_DIM(py_data, 1) != 3         // Second dim is 3.
	  || PyArray_TYPE(py_data) != NPY_DOUBLE  // array of double
	  || !PyArray_ISCARRAY_RO(py_data))       // Contiguous etc.
        throw AsapError("Atoms data has unexpected type: ") << name;
    return (Vec *) PyArray_DATA(py_data);
}

int DynamicAtoms::GetNAtoms()
{
  PyArrayObject *py_data = AsPyArray(PyDict_GetItemString(arrays, "positions"));
  if (py_data == NULL)
    throw AsapError("DynamicsAtoms::GetNAtoms: Atoms.arrays has no positions");
  if (PyArray_NDIM(py_data) != 2           // two-dimensional
        || PyArray_DIM(py_data, 1) != 3         // Second dim is 3.
        || PyArray_TYPE(py_data) != NPY_DOUBLE  // array of double
        || !PyArray_ISCARRAY_RO(py_data))       // Contiguous etc.
      throw AsapError("Atoms positions have unexpected type: ");
  return PyArray_DIM(py_data, 0);
}

double *DynamicAtoms::GetDoubleData(PyObject *name)
{
  PyArrayObject *py_data = AsPyArray(PyDict_GetItem(arrays, name));
  if (py_data == NULL)
      throw AsapError("Atoms.arrays has no ") << PyString_AsString(name);
  if (PyArray_Check(py_data) != 1           // A NumPy array
        || PyArray_TYPE(py_data) != NPY_DOUBLE  // array of double
        || !PyArray_ISCARRAY_RO(py_data))       // Contiguous etc.
      throw AsapError("Atoms data has unexpected type: ") << name;
  return (double *) PyArray_DATA(py_data);
}

double *DynamicAtoms::GetDoubleDataMaybe(PyObject *name)
{
  if (PyDict_Contains(arrays, name))
    return GetDoubleData(name);
  else
    return NULL;
}

// Helper function for DynamicAtoms::GetAtomicNumbers
template<class T>
static void copynum(vector<asap_z_int> &num, 
                    PyArrayObject *py_numbers)
{
    T *from = (T *) PyArray_DATA(py_numbers);
    vector<asap_z_int>::iterator j = num.begin();
    for (int i = 0; i < PyArray_DIM(py_numbers, 0); i++)
        *j++ = (asap_z_int) from[i];
}

const asap_z_int *DynamicAtoms::GetAtomicNumbers()
{
    // Get a borrowed reference to the array.  We do not keep 
    // a reference ourselves as the array cannot go away before
    // we are done using it.
    PyArrayObject *py_data = AsPyArray(PyDict_GetItemString(arrays, "numbers"));
    if (py_data == NULL)
        throw AsapError("Atoms.arrays has no numbers");
    if (PyArray_NDIM(py_data) != 1           // one-dimensional
  	  || !PyArray_ISCARRAY_RO(py_data))       // Contiguous etc.
        throw AsapError("Atoms data 'numbers' has unexpected shape");
    
    // Now, we need to handle the different integer types that may be used
    // for the atomic numbers, depending on 32-bit vs 64 bit architectures 
    // etc.
    int tn = PyArray_TYPE(py_data);
    if (PyArray_EquivTypenums(tn, ASAP_Z_ARRAYTYPE))
        return (asap_z_int *) PyArray_DATA(py_data);  // Lucky, no conversion.
    // The data needs to be converted
    conv_numbers.resize(PyArray_DIM(py_data, 0));
    if (PyArray_EquivTypenums(tn, NPY_INT32))
        copynum<npy_int32>(conv_numbers, py_data);
    else if (PyArray_EquivTypenums(tn, NPY_INT64))
	    copynum<npy_int64>(conv_numbers, py_data);
    else if (PyArray_EquivTypenums(tn, NPY_INT8))
	    copynum<npy_int8>(conv_numbers, py_data);
    else if (PyArray_EquivTypenums(tn, NPY_INT16))
	    copynum<npy_int16>(conv_numbers, py_data);
    else
	    throw AsapError("Atomic numbers are an unsupported integer type.");
	return &conv_numbers[0];
}

