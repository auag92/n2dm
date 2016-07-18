// -*- C++ -*-
// Atoms.cpp:  The interface to the ase Atoms object.
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

#include "NormalAtoms.h"
#include "Asap.h"
#include "Exception.h"
//#define ASAPDEBUG
#include "Debug.h"
#include <cstring>
#include <iostream>
#include <math.h>
#include <cstdio>

using std::cerr;
using std::endl;
using std::flush;

#define PARANOID

NormalAtoms::NormalAtoms()
{
  CONSTRUCTOR;
  DEBUGPRINT;
  nAtoms = 0;
  nGhosts = 0;
  active = false;
  firsttime = true;
  count_atoms = count_cell = count_positions = count_momenta
    = count_numbers = count_inverse_cell = 1;
  py_arrays = py_atoms = NULL;
  momenta_ready = momenta_nonzero = false;
  getmasses_pyname = PyString_FromString("get_masses");
  assert(getmasses_pyname != NULL);
  DEBUGPRINT;
}

NormalAtoms::~NormalAtoms()
{
  DESTRUCTOR;
  DEBUGPRINT;
  if (active > 1)
    std::cerr << "ASAP INTERNAL ERROR: Atoms in nested open when deallocated.!!" << std::endl;
  while (active)
    End();
  Py_XDECREF(getmasses_pyname);
  DEBUGPRINT;
}

void NormalAtoms::Begin(PyObject *pyatoms, bool allow_reopen /* = false */)
{
  DEBUGPRINT;
  // Grab a reference to the atoms while we access it, as we may
  // delete object in End() if it is called.
  assert(pyatoms != NULL);
  py_atoms = pyatoms;
  Py_INCREF(py_atoms); // DECREFed by End() - or DoBegin if trouble.
  if (active)
    {
      assert(active > 0);  // If negative something is REALLY wrong.
      // Atoms are already open - this may be unexpected!
      if (expect_reopen < active)
        {
          // Yes, unexpected reopen - may happen after an exception has
          // been caught.  Print a warning, then close atoms and continue
          // with lowest-level opening.
          cerr << endl
               << "Asap warning: NormalAtoms::Begin called twice without a call to End."
               << endl
               << "This may signify a programming error, but is OK if the previous"
               << endl
               << "calculation ended with an exception."
               << endl << endl;
          active=1;  // Makes sure End() really ends.
          End();
        }
      else
        {
          // Expected reopen
          active++;  // Count it.
          if (allow_reopen)
            expect_reopen = active;
          Py_DECREF(py_atoms);
          return;  // Nothing more to do!
        }
    }
  // OK, atoms are not open, open them.
  assert(active == 0);
  DoBegin(pyatoms);   // also sets active=1 if no exception occurs.
  if (allow_reopen)
    expect_reopen = active;
  else
    expect_reopen = 0;
}

void NormalAtoms::DoBegin(PyObject *pyatoms)
{
  DEBUGPRINT;
  PyArrayObject *py_positions = NULL;
  PyArrayObject *py_pbc = NULL;
  PyArrayObject *py_cell = NULL;
  PyArrayObject *py_numbers = NULL;
  PyObject *py_ghosts = NULL;
  PyArrayObject *py_gh_pos = NULL;
  PyArrayObject *py_gh_num = NULL;
  py_masses = NULL;

  bool step_count_atoms = false;
  bool step_count_cell = false;
  bool step_count_positions = false;
  bool step_count_numbers = false;
  bool changed_boundary_conditions = false;

  // We need an exception handler, so we can DECREF stuff if
  // something goes wrong.
  try
    {
      py_arrays = PyObject_GetAttrString(pyatoms, "arrays");
      py_positions = AsPyArray(PyMapping_GetItemString(py_arrays, "positions"));
      py_numbers = AsPyArray(PyMapping_GetItemString(py_arrays, "numbers"));
      py_pbc = AsPyArray(PyObject_GetAttrString(pyatoms, "pbc"));
      py_cell = AsPyArray(PyObject_GetAttrString(pyatoms, "cell"));
      if ((py_arrays == NULL) || (py_positions == NULL) || (py_numbers == NULL)
	  || (py_pbc == NULL) || (py_cell == NULL))
	throw AsapError("Failed to extract data from the atoms object.");
      DEBUGPRINT;

      // Check that the positions object is sensible.
#ifdef PARANOID
      if (PyArray_NDIM(py_positions) != 2           // two-dimensional
	  || PyArray_DIM(py_positions, 1) != 3         // Second dim is 3
	  || PyArray_TYPE(py_positions) != NPY_DOUBLE  // array of double
	  || !PyArray_ISCARRAY_RO(py_positions))       // Contiguous etc.
	{
	  DEBUGPRINT;
	  throw AsapError("The positions array has a wrong type or shape.");
	}
#endif
      DEBUGPRINT;

      // Check the number of atoms
      int newNatoms = PyArray_DIM(py_positions, 0);
      DEBUGPRINT;

#ifdef PARANOID
      // Check that the atomic numbers are sensible.
      if (PyArray_NDIM(py_numbers) != 1           // One-dimensional
	  || PyArray_DIM(py_numbers, 0) != newNatoms // One per atom
	  || !PyArray_ISINTEGER(py_numbers)          // array of integers
	  || !PyArray_ISCARRAY_RO(py_numbers))      // Contiguous etc.
	throw AsapError("The atomic numbers array has a wrong type or shape.");
#endif
      
      DEBUGPRINT;
      // Check for ghost atoms.
      py_ghosts = PyObject_GetAttrString(pyatoms, "ghosts");
      if (py_ghosts == NULL)
	{
	  PyErr_Clear();
	  hasGhosts = false;
	}
      else
	{
          hasGhosts = true;
	  py_gh_pos = AsPyArray(PyMapping_GetItemString(py_ghosts, "positions"));
	  py_gh_num = AsPyArray(PyMapping_GetItemString(py_ghosts, "numbers"));
	  if ((py_gh_pos == NULL) || (py_gh_num == NULL))
	    {
	      throw AsapError("Failed to extract ghost atom information.");
	    }
#ifdef PARANOID
	  if (PyArray_NDIM(py_gh_pos) != 2           // two-dimensional
	      || PyArray_DIM(py_gh_pos, 1) != 3         // Second dim is 3
	      || PyArray_TYPE(py_gh_pos) != NPY_DOUBLE  // array of double
	      || !PyArray_ISCARRAY_RO(py_gh_pos))       // Contiguous etc.
	    {
	      DEBUGPRINT;
	      throw AsapError("The ghost positions array has a wrong type or shape.");
	    }
#endif
	  nGhosts = PyArray_DIM(py_gh_pos, 0);
#ifdef PARANOID
	  // Check that the atomic numbers are sensible.
	  if (PyArray_NDIM(py_gh_num) != 1            // One-dimensional
	      || PyArray_DIM(py_gh_num, 0) != nGhosts // Correct number
	      || PyArray_TYPE(py_gh_num) != PyArray_TYPE(py_numbers)     // array of longs
	      || !PyArray_ISCARRAY_RO(py_gh_num))       // Contiguous etc.
	    throw AsapError("The ghost atomic numbers array has a wrong type or shape.");
#endif
	}  

      // A change in the number of ghost atoms is not a change in the object.
      if (newNatoms != nAtoms)
	{
	  step_count_atoms = true;
	  nAtoms = newNatoms;
	  // We do not check if nGhosts has changed, since nGhosts
	  // will already have been updated by the ParallelAtoms class
	  // when it changed the sizes of these arrays.  Also,
	  // changing the number of ghosts is not a change in the
	  // atoms, as that would cause extra NB list updates (?)
	  // Instead, we just resize the arrays if needed.
	}
	
      // Check and copy the boundary conditions
      check_boundary_conditions(py_pbc, step_count_atoms,
				changed_boundary_conditions);
      all_periodic_boundaries = (pbc[0] && pbc[1] && pbc[2]);
      all_free_boundaries = !(pbc[0] || pbc[1] || pbc[2]);
      
      // Check and copy the unit cell
      step_count_cell = check_unit_cell(py_cell);

      DEBUGPRINT;
      // Check and copy the atomic numbers
      step_count_numbers = check_numbers(py_numbers, py_gh_num,
					 step_count_atoms);
      
      // Check and copy the positions
      step_count_positions =
	check_positions(py_positions, py_gh_pos,
			step_count_atoms || step_count_cell);
      DEBUGPRINT;
    }
  catch(AsapError)
    {
      DEBUGPRINT;
      CHECKREF(py_atoms);
      Py_CLEAR(py_atoms);
      XCHECKREF(py_pbc);
      Py_XDECREF(py_pbc);
      XCHECKREF(py_cell);
      Py_XDECREF(py_cell);
      XCHECKREF(py_positions);
      Py_XDECREF(py_positions);
      XCHECKREF(py_arrays);
      Py_CLEAR(py_arrays);
      XCHECKREF(py_numbers);
      Py_XDECREF(py_numbers);
      XCHECKREF(py_ghosts);
      Py_XDECREF(py_ghosts);
      XCHECKREF(py_gh_pos);
      Py_XDECREF(py_gh_pos);
      XCHECKREF(py_gh_num);
      Py_XDECREF(py_gh_num);
      throw;
    }

  // Now step the counters
  if (step_count_atoms)
    {
      count_atoms++;
      count_cell++;
      count_positions++;
      count_numbers++;
    }
  else if (step_count_cell)
    {
      count_cell++;
      count_positions++;
    }
  else
    {
      if (step_count_positions || step_count_numbers)
	count_positions++;
      if (step_count_numbers)
	count_numbers++;
    }
  CHECKREF(py_pbc);
  Py_DECREF(py_pbc);
  CHECKREF(py_cell);
  Py_DECREF(py_cell);
  CHECKREF(py_positions);
  Py_DECREF(py_positions);
  XCHECKREF(py_numbers);
  Py_XDECREF(py_numbers);
  XCHECKREF(py_ghosts);
  Py_XDECREF(py_ghosts);
  XCHECKREF(py_gh_pos);
  Py_XDECREF(py_gh_pos);
  XCHECKREF(py_gh_num);
  Py_XDECREF(py_gh_num);

  // We have not checked the momenta.
  momenta_ready = false;
  // Mark that the atoms are now open
  active = 1;

  // ParallelAtoms must do something special if boundary conditions
  // have changed.
  if (changed_boundary_conditions && !firsttime)
    NewBoundaryConditions(pbc);
  firsttime = false;
  DEBUGPRINT;
}

void NormalAtoms::check_boundary_conditions(PyArrayObject *py_pbc, bool &step_count_atoms,
				      bool &changed_boundary_conditions)
{
#ifdef PARANOID
  if (PyArray_NDIM(py_pbc) != 1         // one-dimensional
      || PyArray_DIM(py_pbc, 0) != 3       // Three data point 3
      || PyArray_TYPE(py_pbc) != NPY_BOOL  // array of booleans
      || !PyArray_ISCARRAY_RO(py_pbc))     // Contiguous etc.
    throw AsapError("The boundary conditions array has a wrong type or shape.");
#endif
  DEBUGPRINT;
  npy_bool *new_pbc = (npy_bool*) PyArray_DATA(py_pbc);
  for (int i = 0; i < 3; ++i)
    if (pbc[i] != (bool) new_pbc[i])
      {
	step_count_atoms = true;
	changed_boundary_conditions = true;
	pbc[i] = (bool) new_pbc[i];
      }
}

bool NormalAtoms::check_unit_cell(PyArrayObject *py_cell)
{
  bool step_count_cell = false;
#ifdef PARANOID
  if (PyArray_NDIM(py_cell) != 2           // two-dimensional
      || PyArray_DIM(py_cell, 0) != 3         // First dim is 3
      || PyArray_DIM(py_cell, 1) != 3         // Second dim is 3
      || PyArray_TYPE(py_cell) != NPY_DOUBLE  // array of double
      || !PyArray_ISCARRAY_RO(py_cell))       // Contiguous etc.
    throw AsapError("The unit cell has a wrong type or shape.");
  assert(PyArray_NBYTES(py_cell) == 3 * sizeof(Vec));
#endif
  if (memcmp(cell, PyArray_DATA(py_cell), 3*sizeof(Vec)) != 0)
    {
      step_count_cell = true;
      memcpy(cell, PyArray_DATA(py_cell), 3*sizeof(Vec));
    }
  return step_count_cell;
}

// Helper function for NormalAtoms::check_numbers
template<class T>
static void copynum(vector<asap_z_int> &num, 
    PyArrayObject *py_numbers, PyArrayObject *py_gh_num)
{
  T *from = (T *) PyArray_DATA(py_numbers);
  vector<asap_z_int>::iterator j = num.begin();
  for (int i = 0; i < PyArray_DIM(py_numbers, 0); i++)
    *j++ = (asap_z_int) from[i];
  if (py_gh_num != NULL)
    {
      from = (T *) PyArray_DATA(py_gh_num);
      for (int i = 0; i < PyArray_DIM(py_gh_num,0); i++)
	*j++ = (asap_z_int) from[i];
    }
  assert(j == num.end()); 
}

template<class T>
static bool chknum(vector<asap_z_int> &num, 
    PyArrayObject *py_numbers, PyArrayObject *py_gh_num)
{
  bool step_count_numbers = false;
  T *from = (T *) PyArray_DATA(py_numbers);
  vector<asap_z_int>::iterator j = num.begin();
  for (int i = 0; i < PyArray_DIM(py_numbers,0); i++)
    {
      asap_z_int z = (asap_z_int) from[i];
      step_count_numbers |= (*j != z);
      *j++ = z;
    }
  if (py_gh_num)
    {
      from = (T *) PyArray_DATA(py_gh_num);
      for (int i = 0; i < PyArray_DIM(py_gh_num,0); i++)
	*j++ = (asap_z_int) from[i];
    }
  assert(j == num.end());
  return step_count_numbers;
}

bool NormalAtoms::check_numbers(PyArrayObject *py_numbers, PyArrayObject *py_gh_num,
			  bool step_count_atoms)
{
  bool step_count_numbers = false;
  if (numbers.size() != nAtoms+nGhosts)
    numbers.resize(nAtoms+nGhosts);
  if (step_count_atoms)
    {
      // Copy the atomic numbers
      step_count_numbers = true;
      int tn = PyArray_TYPE(py_numbers);
      if (PyArray_EquivTypenums(tn, ASAP_Z_ARRAYTYPE))
	copynum<asap_z_int>(numbers, py_numbers, py_gh_num);
      else if (PyArray_EquivTypenums(tn, NPY_INT32))
	copynum<npy_int32>(numbers, py_numbers, py_gh_num);
      else if (PyArray_EquivTypenums(tn, NPY_INT64))
	copynum<npy_int64>(numbers, py_numbers, py_gh_num);
      else if (PyArray_EquivTypenums(tn, NPY_INT8))
	copynum<npy_int8>(numbers, py_numbers, py_gh_num);
      else if (PyArray_EquivTypenums(tn, NPY_INT16))
	copynum<npy_int16>(numbers, py_numbers, py_gh_num);
      else
	throw AsapError("Atomic numbers are an unsupported integer type.");
    }
  else
    {
      // Check if the atomic numbers have changed
      int tn = PyArray_TYPE(py_numbers);
      if (PyArray_EquivTypenums(tn, ASAP_Z_ARRAYTYPE))
	step_count_numbers = chknum<asap_z_int>(numbers, py_numbers, py_gh_num);
      else if (PyArray_EquivTypenums(tn, NPY_INT32))
	step_count_numbers = chknum<npy_int32>(numbers, py_numbers, py_gh_num);
      else if (PyArray_EquivTypenums(tn, NPY_INT64))
	step_count_numbers = chknum<npy_int64>(numbers, py_numbers, py_gh_num);
      else if (PyArray_EquivTypenums(tn, NPY_INT8))
	step_count_numbers = chknum<npy_int8>(numbers, py_numbers, py_gh_num);
      else if (PyArray_EquivTypenums(tn, NPY_INT16))
	step_count_numbers = chknum<npy_int16>(numbers, py_numbers, py_gh_num);
      else
	throw AsapError("Atomic numbers are an unsupported integer type.");
    }
  return step_count_numbers;
}

bool NormalAtoms::check_positions(PyArrayObject *py_positions, PyArrayObject *py_gh_pos,
			    bool step_count_atoms_or_cell)
{
  DEBUGPRINT;
  bool step_count_positions;
  bool pos_resized = (positions.size() != nAtoms+nGhosts);
  if (pos_resized)
    positions.resize(nAtoms+nGhosts);
  step_count_positions = (step_count_atoms_or_cell ||
			  memcmp(&positions[0], PyArray_DATA(py_positions),
				 nAtoms*sizeof(Vec)) != 0);
  if (step_count_positions || pos_resized)
    {
      memcpy(&positions[0], PyArray_DATA(py_positions), nAtoms*sizeof(Vec));
      if (py_gh_pos && nGhosts > 0)
	memcpy(&positions[nAtoms], PyArray_DATA(py_gh_pos),
	       nGhosts*sizeof(Vec));
#if 0
      for (int aa=0; aa < nAtoms; aa++)
	{
	  assert(!isnan(positions[aa][0]));
	  assert(!isnan(positions[aa][1]));
	  assert(!isnan(positions[aa][2]));
	}
      for (int aa=nAtoms; aa < nAtoms+nGhosts; aa++)
	{
	  assert(!isnan(positions[aa][0]));
	  assert(!isnan(positions[aa][1]));
	  assert(!isnan(positions[aa][2]));
	}
#endif
    }
  return step_count_positions;
}

void NormalAtoms::End()
{
  DEBUGPRINT;
  if (active <= 0)
    throw AsapError("NormalAtoms::End() called without a previous call to Begin()");
  active--;
  if (expect_reopen > active)
    expect_reopen = active;
  if (active == 0)
    {
      // Really close the atoms.
      Py_XDECREF(py_masses);
      py_masses = NULL;
      CHECKREF(py_atoms);
      Py_CLEAR(py_atoms);
      CHECKREF(py_arrays)
      Py_CLEAR(py_arrays);
      DEBUGPRINT;
    }
}

void NormalAtoms::GetPositions(vector<Vec> &pos, bool ghosts /* = false */) const
{
  DEBUGPRINT;
  assert(active);
  pos.clear();
  if (ghosts || nGhosts == 0)
    {
      assert(positions.size() == nAtoms + nGhosts);
      if (pos.capacity() < nAtoms + nGhosts)
        pos.reserve(nAtoms + nGhosts + (nAtoms + nGhosts)/25);  // 4% extra
      pos.insert(pos.begin(), positions.begin(), positions.end());
    }
  else
    {
      if (pos.capacity() < nAtoms)
        pos.reserve(nAtoms + nAtoms/25);  // 4% extra
      pos.insert(pos.begin(), positions.begin(), positions.begin() + nAtoms);
      assert(pos.size() == nAtoms);
    }
  DEBUGPRINT;
}

void NormalAtoms::GetScaledPositions(vector<Vec> &pos, bool ghosts /* = false */)
{
  DEBUGPRINT;
  int n = nAtoms;
  if (ghosts)
    n += nGhosts;
  assert(positions.size() >= n);
  const Vec *inv = GetInverseCell();
  if (pos.capacity() < n)
    pos.reserve(n + n/25);  // Reserve 4% extra.
  pos.resize(n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < 3; j++)
      pos[i][j] = positions[i][0] * inv[0][j]   
                + positions[i][1] * inv[1][j]
                + positions[i][2] * inv[2][j];
  DEBUGPRINT;
}

void NormalAtoms::GetScaledPositions(vector<Vec> &scaledpos, const set<int> &which)
{
  DEBUGPRINT;
  assert(scaledpos.size() == which.size());
  const Vec *inv = GetInverseCell();
  vector<Vec>::iterator spi = scaledpos.begin();
  for (set<int>::const_iterator i = which.begin(); i != which.end(); ++i,++spi)
    for (int j = 0; j < 3; j++)
      (*spi)[j] = positions[*i][0] * inv[0][j]
                + positions[*i][1] * inv[1][j]
                + positions[*i][2] * inv[2][j];
  DEBUGPRINT;
}

void NormalAtoms::GetListOfElements(set<int> &elements) const
{
  DEBUGPRINT;
  const asap_z_int *atomicnumbers = GetAtomicNumbers();

  elements.clear();
  for (int i = 0; i < nAtoms; i++)
    {
      int z = atomicnumbers[i];
      if (elements.find(z) == elements.end())
        elements.insert(z);
    }
  DEBUGPRINT;
}

const Vec *NormalAtoms::GetMomenta()
{
  DEBUGPRINT;
  (void) GetMomentaCounter();  // Make sure momenta are updated
  if (momenta_ready)
    return &momenta[0];
  else
    return NULL;
}

int NormalAtoms::GetMomentaCounter()
{
  DEBUGPRINT;
  assert(active);
  if (!momenta_ready)
    {
      PyArrayObject *py_momenta = AsPyArray(PyDict_GetItemString(py_arrays, "momenta"));
      if (py_momenta != NULL)
        {
#ifdef PARANOID
          // Check that the positions object is sensible.
          if (PyArray_NDIM(py_momenta) != 2           // two-dimensional
              || PyArray_DIM(py_momenta, 0) != nAtoms    // Fisrt dim is nAtoms
              || PyArray_DIM(py_momenta, 1) != 3         // Second dim is 3
              || PyArray_TYPE(py_momenta) != NPY_DOUBLE  // array of double
              || !PyArray_ISCARRAY_RO(py_momenta))       // Contiguous etc.
            throw AsapError("The momenta array has a wrong type or shape.");
#endif
          bool step_mom_ctr = false;
          if (momenta.size() != nAtoms)
            {
              step_mom_ctr = true;
              momenta.resize(nAtoms);
            }
          else
            step_mom_ctr = (memcmp(&momenta[0], PyArray_DATA(py_momenta),
                nAtoms*sizeof(Vec)) != 0);
          if (step_mom_ctr)
            {
              count_momenta++;
              memcpy(&momenta[0], PyArray_DATA(py_momenta), nAtoms*sizeof(Vec));
            }
          momenta_nonzero = true;
          momenta_ready = true;
        }
      else
        {
          // No momenta, make them zero
          PyErr_Clear();  // Clear the KeyError from PyDict_GetItemString
          if ((momenta.size() != 0) && momenta_nonzero)
            count_momenta++;  // OK, this is a change
          momenta.resize(nAtoms);
          memset(&momenta[0], 0, nAtoms*sizeof(Vec));
          momenta_nonzero = false;
          momenta_ready = true;
        }
    }
  DEBUGPRINT;
  return count_momenta;
}
  
const double *NormalAtoms::GetMasses()
{
  if (py_masses == NULL)
    {
      assert(active);
      py_masses = AsPyArray(PyObject_CallMethodObjArgs(py_atoms, getmasses_pyname, NULL));
      if (py_masses == NULL)
        throw AsapPythonError();
#ifdef PARANOID
      if (PyArray_NDIM(py_masses) != 1              // one-dimensional
          || PyArray_DIM(py_masses, 0) < nAtoms       // Correct size
          || PyArray_TYPE(py_masses) != NPY_DOUBLE  // array of double
          || !PyArray_ISCARRAY_RO(py_masses)
          )       // Contiguous etc.
        {
          DEBUGPRINT;
          std::cerr << PyString_AsString(PyObject_Repr((PyObject *) py_masses)) << std::endl;
          throw AsapError("The masses array has a wrong type or shape.");
        }
#endif
    }
  return (const double *) PyArray_DATA(py_masses);
}

double NormalAtoms::GetVolume() const
{
  DEBUGPRINT;
  double det;
  assert(active);
  det = -cell[0][2]*cell[1][1]*cell[2][0] +
    cell[0][1]*cell[1][2]*cell[2][0] +
    cell[0][2]*cell[1][0]*cell[2][1] -
    cell[0][0]*cell[1][2]*cell[2][1] -
    cell[0][1]*cell[1][0]*cell[2][2] +
    cell[0][0]*cell[1][1]*cell[2][2];
  DEBUGPRINT;
  return fabs(det);
}

const double *NormalAtoms::GetCellHeights()
{
  DEBUGPRINT;
  if (count_inverse_cell < count_cell)
    invert_cell();
  return heights;
}

const Vec *NormalAtoms::GetInverseCell()
{
  DEBUGPRINT;
  if (count_inverse_cell < count_cell)
    invert_cell();
  return inverse;
}

void NormalAtoms::invert_cell()
{
  DEBUGPRINT;
  assert(active);
  count_inverse_cell = count_cell;
  double determinant = Cross(cell[0], cell[1]) * cell[2];
  // Find heights
  for (int i = 0; i < 3; i++)
    {
      Vec inv = Cross(cell[(i + 1) % 3], cell[(i + 2) % 3]);
      heights[i] = fabs(determinant) / sqrt(Length2(inv));
    }
  // Invert matrix.  I_i,j = { C_j-1,i-1 C_j+1,i+1 - C_j+1,i-1 C_j-1,i+1 } / det
  for (int i = 0; i < 3; i++)
    {
      int ip = (i + 1) % 3;
      int im = (i + 2) % 3;
      for (int j = 0; j < 3; j++)
        {
          int jp = (j + 1) % 3;
          int jm = (j + 2) % 3;
          inverse[i][j] = (cell[jm][im]*cell[jp][ip] -
              cell[jp][im]*cell[jm][ip]) / determinant;
        }
    }
  DEBUGPRINT;
}

// Update flag and counters across processors.
//
// Called by a Potential with a flag indicating if the neighborlist
// should be updated.  In a serial simulation, just return the flag.
//
// In a parallel simulation, communicate across processors so the
// counters of this Atoms object and the flag passed as an argument
// are updated if just one processor thinks an update is necessary.
//
// This version is overloaded in parallel simulations.
bool NormalAtoms::UpdateBeforeCalculation(bool flag, double range /* NOT USED */)
{
  DEBUGPRINT;
  return flag;
}

// Set a python array on the atoms.
void NormalAtoms::SetData(const char *name, PyObject *data)
{
  assert(py_arrays != NULL);
  if (PyDict_SetItemString(py_arrays, name, data) == -1)
    throw AsapPythonError();
}

// Get a python array from the atoms
PyObject *NormalAtoms::GetData(const char *name)
{
  assert(py_arrays != NULL);
  PyObject *res = PyDict_GetItemString(py_arrays, name);
  if (res == NULL)
    throw AsapError("Failed to get array from atoms: ") << name;
  Py_INCREF(res);
  return res;
}

// Remove a python array from the atoms
void NormalAtoms::DeleteData(const char *name)
{
  assert(py_arrays != NULL);
  if (PyDict_DelItemString(py_arrays, name) == -1)
    throw AsapError("Failed to delete array from atoms: ") << name;
}
    
long NormalAtoms::PrintMemory() const
{
  long mem = 0;  // Count the big stuff.
  mem += positions.capacity() * sizeof(Vec);
  mem += momenta.capacity() * sizeof(Vec);
  mem += numbers.capacity() * sizeof(int);
  mem = (mem + 512*1024)/(1024*1024);
  char buffer[500];
  snprintf(buffer, 500, "*MEM* Atoms/C++  %ld MB.", mem);
  cerr << buffer << endl;
  return mem;
}
