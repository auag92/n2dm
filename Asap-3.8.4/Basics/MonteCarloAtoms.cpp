// -*- C++ -*-
// MonteCarloAtoms.cpp:  Adds Monte Carlo optimization to the atoms interface.
//
// Copyright (C) 2009-2011 Jakob Schiotz and Center for Individual
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
#include "MonteCarloAtoms.h"
#include "Asap.h"
#include "Exception.h"
//#define ASAPDEBUG
#include "Debug.h"


bool MonteCarloAtoms::check_numbers(PyArrayObject *py_numbers, PyArrayObject *py_gh_num,
				    bool step_count_atoms)
{
  assert(py_gh_num == NULL);  // Must be a serial simulation.
  PyArrayObject *mc_optim = AsPyArray(PyObject_GetAttrString(py_atoms, "mc_optim"));
  if (mc_optim == NULL)
    throw AsapError("Not a Monte Carlo enabled Atoms object!");
  if (PyArray_NDIM(mc_optim) != 1           // two-dimensional
      || PyArray_TYPE(mc_optim) != NPY_INT     // array of integers
      || !PyArray_ISCARRAY_RO(mc_optim))       // Contiguous etc.
    {
      DEBUGPRINT;
      throw AsapError("The mc_optim array has a wrong type or shape.");
    }
  int *mc_optim_data = (int *) PyArray_DATA(mc_optim);
  int nmodif = mc_optim_data[0];
  if (nmodif > 100 || step_count_atoms)
    {
      // MC optimization irrelevant
      mc_optim_relevant = false;
      return NormalAtoms::check_numbers(py_numbers, py_gh_num, step_count_atoms);
    }
  else if (nmodif > 0)
    {
      // A small number of atoms have changed: MC optimizations relevant.
      mc_optim_relevant = true;
      modified_atoms.clear();
      long *from = (long *) PyArray_DATA(py_numbers);
      bool step_count_numbers = false;
      for (int i = 0; i < nmodif; i++)
	{
	  int n = *++mc_optim_data;
	  modified_atoms.insert(n);
	  int newz = (int) from[n];
	  if (numbers[n] != newz)
	    {
	      step_count_numbers = true;
	      numbers[n] = newz;
	    }
	}
      return step_count_numbers;
    }
  else
    {
      // Atoms are not modified!  Signal to check_positions that this is the case.
      mc_optim_relevant = true;
      modified_atoms.clear();
      return false;
    }
}
bool MonteCarloAtoms::check_positions(PyArrayObject *py_positions,
				      PyArrayObject *py_gh_pos,
				      bool step_count_atoms_or_cell)
{
  assert(py_gh_pos == NULL);  // Must be a serial simulation.
  if (mc_optim_relevant)
    {
      if (modified_atoms.size() == 0)
	{
	  // Atoms are not modified! 
	  mc_optim_relevant = false;
	  return false;
	}
      Vec *from = (Vec *) PyArray_DATA(py_positions);
      for (set<int>::const_iterator i = modified_atoms.begin();
	   i != modified_atoms.end(); ++i)
	{
	  positions[*i] = from[*i];
	}
      return true;
    }
  else
    {
      return NormalAtoms::check_positions(py_positions, py_gh_pos,
				    step_count_atoms_or_cell);
    }
}
