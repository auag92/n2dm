// -*- C++ -*-
// GetNeighborList.h: A utility function.
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

#ifndef _GETNEIGHBORLISTS_H
#define _GETNEIGHBORLISTS_H

#include "AsapPython.h"

namespace ASAPSPACE {

class Atoms;

/// Return a neighborlist object and an Atoms access object (either new
/// or reused from the potential).  Begin() will have been called on the atoms.
/// The returned neighborlist will have a cutoff of at least rCut, but
/// it may be larger.

void GetNbList_FromAtoms(PyObject *pyatoms, double rCut,
			 PyObject **nblist, Atoms **atoms);

} // end namespace

#endif // _GETNEIGHBORLISTS_H
