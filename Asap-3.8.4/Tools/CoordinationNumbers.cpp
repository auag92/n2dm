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
#include "CoordinationNumbers.h"
#include "Atoms.h"
#include "NeighborLocator.h"
#include "GetNeighborList.h"
#include "Vec.h"
#include "Timing.h"
#include "Exception.h"
#include <vector>
using std::vector;

namespace ASAPSPACE {

void CoordinationNumbers(PyObject *py_atoms, double rCut, vector<int> &coord)
{
  USETIMER("CoordinationNumbers");
  Atoms *atoms = NULL;
  PyObject *py_nblist = NULL;
  // Get an Atoms access object and a neighbor list with cutoff rCut or larger.
  GetNbList_FromAtoms(py_atoms, rCut, &py_nblist, &atoms);
  NeighborLocator *nl = ((PyAsap_NeighborLocatorObject*) py_nblist)->cobj;
  assert(nl != NULL);

  int nAtoms = atoms->GetNumberOfAtoms();
  int nTotal = nAtoms + atoms->GetNumberOfGhostAtoms();
  
  // The neighbor list is unlikely to support full lists, and it is
  // likely to have a too long cutoff (if we reuse the one from the
  // potential).

  coord.clear();
  coord.resize(nAtoms);
  int maxlistlen = nl->MaxNeighborListLength();
  vector<int> nb_buffer(maxlistlen);
  vector<Vec> diffs(maxlistlen);
  vector<double> diffs2(maxlistlen);
  double rC2 = rCut * rCut;
  for (int a1 = 0; a1 < nTotal; a1++)
    {
      int size = maxlistlen;
      int nnb =  nl->GetNeighbors(a1, &nb_buffer[0], &diffs[0], &diffs2[0],
				  size, rCut);
      for (int n = 0; n < nnb; n++)
	{
	  int a2 = nb_buffer[n];
	  if (diffs2[n] <= rC2)
	    {
	      if (a1 < nAtoms)
		coord[a1]++;
	      if (a2 < nAtoms)
		coord[a2]++;
	    }
	}
    }
  atoms->End();
  Py_DECREF(py_nblist);
  AsapAtoms_DECREF(atoms);
}

} // end namespace
