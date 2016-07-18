// MolecularDynamics.cpp  --  The Velocity Verlet molecular dynamics algorithm.
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

#include "MolecularDynamics.h"
#include "Potential.h"
#include "DynamicAtoms.h"
//#define ASAPDEBUG
#include "Debug.h"
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

// Keep a reference while this object is in scope.
class keep_pyreference
{
public:
  keep_pyreference(PyObject *p)
  {
    ref = p;
    Py_INCREF(ref);
  }
  ~keep_pyreference()
  {
    Py_DECREF(ref);
  }
private:
  PyObject *ref;
};

class threadchecker
{
public:
  threadchecker(PyObject *atoms)
  {
#ifdef _OPENMP
    n = omp_get_max_threads();
    active = false;
    if (n > 1)
      {
        Py_ssize_t natoms = PySequence_Length(atoms);
        if (natoms == -1)
          throw AsapError("Atoms do not appear to be a sequence.");
        if (natoms < 1000)
          {
            active = true;
            std::cerr << "\nASAP Warning: disabling OpenMP parallelization, too few atoms ("
                << natoms << ")" << std::endl;
            omp_set_num_threads(1);
          }
      }
#endif
  }
#ifdef _OPENMP
  ~threadchecker()
  {
    if (active)
      omp_set_num_threads(n);
  }
private:
  int n;
  bool active;
#endif
};

MolecularDynamics::MolecularDynamics(PyObject *py_atoms, Potential *calc,
    double timestep)
{
  DEBUGPRINT;
  Py_INCREF(py_atoms);
  this->py_atoms = py_atoms;
  atoms = new DynamicAtoms(py_atoms);
  calculator = calc;
  //std::cout << "MolecularDynamics using the calculator " << calc << std::endl;
  this->timestep = timestep;
  getforcestring = PyString_FromString("get_forces");
  steps = 0;
  DEBUGPRINT;
}

MolecularDynamics::~MolecularDynamics()
{
  Py_DECREF(py_atoms);
  delete atoms;
  DEBUGPRINT;
}
void MolecularDynamics::Run(int nsteps, PyObject *observers, PyObject *self)
{
  threadchecker th = threadchecker(py_atoms);
  keep_pyreference keep_self = keep_pyreference(self);
  Run2(nsteps, observers, self);
}

const Vec *MolecularDynamics::GetForces()
{
  // Get forces, either directly from the Asap Potential or through Python.
  if (calculator != NULL)
  {
    const vector<Vec> &forces = calculator->GetForces(py_atoms);
    nAtoms = forces.size();
    assert(nAtoms == atoms->GetNAtoms());
    return &forces[0];
  }
  else
  {
    const vector<Vec> &forces = GetForcesThroughPython();
    nAtoms = forces.size();
    assert(nAtoms == atoms->GetNAtoms());
    return &forces[0];
  }
}

const vector<Vec> &MolecularDynamics::GetForcesThroughPython()
{
  PyArrayObject *py_forces = AsPyArray(PyObject_CallMethodObjArgs(py_atoms,
      getforcestring, NULL));
  if (py_forces == NULL)
    throw AsapPythonError();
  if (PyArray_NDIM(py_forces) != 2           // two-dimensional
      || PyArray_DIM(py_forces, 1) != 3         // Second dim is 3
      || PyArray_TYPE(py_forces) != NPY_DOUBLE  // array of double
      || !PyArray_ISCARRAY_RO(py_forces))       // Contiguous etc.
    throw AsapError("The forces array has a wrong type or shape.");
  int n = PyArray_DIM(py_forces, 0);
  const Vec *f = (Vec *) PyArray_DATA(py_forces);
  forces_from_python.resize(n);
  for (int i = 0; i < n; i++)
    forces_from_python[i] = f[i];
  Py_DECREF(py_forces);
  return forces_from_python;
}

void MolecularDynamics::ParseObservers(PyObject *observers)
{
  if (obs_callable.size() != 0)
    CleanupObservers();  // Probably leftovers due to an exception in an observer.
  if (observers == Py_None)
    return;
  PyObject *obs = PySequence_Fast(observers, "List of observers is not a sequence.");
  if (obs == NULL)
    throw AsapPythonError();
  Py_ssize_t n_obs = PySequence_Fast_GET_SIZE(obs);
  if (n_obs == 0)
    {
      Py_DECREF(obs);
      return;
    }
  PyObject **items = PySequence_Fast_ITEMS(obs);
  for (int i = 0; i < n_obs; i++)
    {
      PyObject *py_line = PySequence_Fast(items[i], "Observer tuple is not a sequence.");
      if (py_line == NULL)
        {
          Py_DECREF(obs);
          CleanupObservers();
          throw AsapPythonError();
        }
      if (PySequence_Fast_GET_SIZE(py_line) != 4)
        {
          Py_DECREF(obs);
          Py_DECREF(py_line);
          CleanupObservers();
          throw AsapError("Each observer should be described by a 4-tuple.");
        }
      PyObject **line = PySequence_Fast_ITEMS(py_line);
      int interval = (int) PyInt_AsLong(line[1]);
      if (interval < 0)
        {
          Py_DECREF(obs);
          Py_DECREF(py_line);
          CleanupObservers();
          throw AsapError("Observer has negative or non-integer interval.");
        }
      // Now we are ready to add the info about the observer
      obs_callable.push_back(line[0]);
      Py_INCREF(line[0]);
      obs_interval.push_back(interval);
      obs_args.push_back(line[2]);
      Py_INCREF(line[2]);
      obs_kwargs.push_back(line[3]);
      Py_INCREF(line[3]);
      Py_DECREF(py_line);
    }
  Py_DECREF(obs);
}

void MolecularDynamics::CleanupObservers()
{
  for (std::vector<PyObject *>::iterator i = obs_callable.begin();
      i != obs_callable.end(); ++i)
    {
      Py_DECREF(*i);
    }
  obs_callable.clear();
  obs_interval.clear();
  for (std::vector<PyObject *>::iterator i = obs_args.begin();
      i != obs_args.end(); ++i)
    {
      Py_DECREF(*i);
    }
  obs_args.clear();
  for (std::vector<PyObject *>::iterator i = obs_kwargs.begin();
      i != obs_kwargs.end(); ++i)
    {
      Py_DECREF(*i);
    }
  obs_kwargs.clear();
}

void MolecularDynamics::CallObservers(PyObject *self)
{
  bool stepsupdated = false;
  int n = obs_callable.size();
  assert(n == obs_kwargs.size());
  for (int i = 0; i < n; i++)
    if (steps % obs_interval[i] == 0)
      {
        if (!stepsupdated)
          {
            UpdateStepsInPython(self);
            stepsupdated = true;
          }
        PyObject *retval = PyObject_Call(obs_callable[i], obs_args[i], obs_kwargs[i]);
        if (retval == NULL)
          throw AsapPythonError();
      }
}

void MolecularDynamics::UpdateStepsInPython(PyObject *self)
{
  PyObject *n = PyInt_FromLong(steps);
  int x = PyObject_SetAttrString(self, "nsteps", n);
  Py_DECREF(n);
  if (x == -1)
    throw AsapPythonError();
}

