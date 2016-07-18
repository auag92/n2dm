// -*- C++ -*-
// AsapObject.h: Provides a name and __repr__ for all Python objects.
//
// AsapObject should be a base class for all C++ classes that will
// have a life as Python object, to support among others the __repr__
// method.
//
// Copyright (C) 2008-2011 Jakob Schiotz and Center for Individual
// Nanoparticle Functionality, Department of Physics, Technical
// University of Denmark.  Email: schiotz@fysik.dtu.dk
//
// This file is part of Asap version 3.
// Asap is released under the GNU Lesser Public License (LGPL) version 3.
// However, the parts of Asap distributed within the OpenKIM project
// (including this file) are also released under the Common Development
// and Distribution License (CDDL) version 1.0.
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

#ifndef _ASAPOBJECT_H
#define _ASAPOBJECT_H

#include "AsapNamespace.h"
#include <string>
using std::string;

namespace ASAPSPACE {

class AsapObject
{
public:
  virtual ~AsapObject() {};
  virtual string GetName() const = 0;
  virtual string GetRepresentation() const;
};

} // End namespace

#endif // _ASAPOBJECT_H
