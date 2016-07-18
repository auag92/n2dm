// -*- C++ -*-
// Atoms.h:  The interface to the ase Atoms object.
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
//
// The Atoms interface object is the interface between the C++ code
// and the Python Atoms object.
//
// It has been split into an abstract base class (Atoms, in AtomsBasis.h)
// and an implementation (NormalAtoms, in NormalAtoms.h).  This makes it
// possible not only to extend the atoms with inheritance (as in
// ParallelAtoms), but also to use the Decorator pattern to extend the
// functionality (as in ImageAtoms).

#include "AtomsBasis.h"
#include "NormalAtoms.h"
