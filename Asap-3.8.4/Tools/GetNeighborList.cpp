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
#include "GetNeighborList.h"
#include "Atoms.h"
#include "Potential.h"
#include "NeighborCellLocator.h"

namespace ASAPSPACE {
// Return a neighborlist object and an Atoms access object (either new
// or reused from the potential).  Begin() will have been called on the atoms.

void GetNbList_FromAtoms(PyObject *pyatoms, double rCut,
			      PyObject **pynblist, Atoms **atoms)
{
  PyObject *py_pot = NULL;
  Potential *pot = NULL;
  NeighborLocator *nblist = NULL;
  py_pot = PyObject_CallMethod(pyatoms, "get_calculator", "");
  if (py_pot && PyAsap_PotentialCheck(py_pot))
    {
      // Got an ASAP potential.
      pot = ((PyAsap_PotentialObject *)py_pot)->cobj;
      *pynblist = pot->GetNeighborList();
      if (*pynblist)
	{
	  nblist = ((PyAsap_NeighborLocatorObject *)*pynblist)->cobj;
	  if (nblist->GetCutoffRadius() >= rCut)
	    {
	      // Use this neighbor list
	      Py_INCREF(*pynblist);
	      *atoms = nblist->GetAtoms();
	      AsapAtoms_INCREF(*atoms);
	      (*atoms)->Begin(pyatoms);
	      nblist->CheckAndUpdateNeighborList();
	    }
	  else
	    {
	      // Discard neighbor list
	      nblist = NULL;
	      *pynblist = NULL;
	    }
	}
    }
  Py_XDECREF(py_pot);
  if (nblist == NULL)
    {
      // Create interface object and neighborlist object
      *atoms = new NormalAtoms();
      (*atoms)->Begin(pyatoms);
      *pynblist = (PyObject *) PyAsap_NewNeighborCellLocator(*atoms, rCut, .0);
      nblist = ((PyAsap_NeighborLocatorObject *)*pynblist)->cobj;
      nblist->CheckAndUpdateNeighborList();
    }
  assert(*atoms != NULL);
  assert(*pynblist != NULL);
}
	  
} // end namespace
