// OpenKIMInfo.h - a class to collect information about OpenKIM models.
//
// This class is part of the optional Asap module to support OpenKIM
// models.  The class OpenKIMInfo collects information about a given
// OpenKIM model, so it can be opened using the best neighbor list etc.

// Copyright (C) 2014 Jakob Schiotz and Center for Individual
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

#ifndef _OPENKIMINFO_H
#define _OPENKIMINFO_H

#include "AsapPython.h"
#include "AsapNamespace.h"
#include "AsapObject.h"
#include <string>
#include <vector>

class KIM_API_model;

namespace ASAPSPACE {

class OpenKIMinfo;  // Defined later in this file.

/// The Python object corresponding to a OpenKIMinfo object.
typedef struct {
  PyObject_HEAD
  OpenKIMinfo *cobj;
  PyObject *weakrefs;
} PyAsap_OpenKIMinfoObject;

PyAsap_OpenKIMinfoObject *PyAsap_NewOpenKIMinfoObject(const char *name);

class OpenKIMinfo : public AsapObject
{
protected:
  OpenKIMinfo(const char *name);

  friend PyAsap_OpenKIMinfoObject *PyAsap_NewOpenKIMinfoObject(const char *name);

public:
  ~OpenKIMinfo();

  virtual std::string GetName() const {return "OpenKIMinfo";}

  /// Get information about supported elements.
  void GetSupportedSymbols(std::vector<const char *> &symbols);

  /// Get information about supported NBC methods.
  const std::vector<const char *> &GetSupportedNBC();

  /// Get information about supported neighbor access methods.
  const std::vector<const char *> &GetSupportedAccess();

  /// Check if a property exists.  Return index >= 0 if it exists, -1 if not.
  int GetAPIindex(const char *name);

private:
  const char* name;
  KIM_API_model *model;
  char *particletypes;
  std::vector<const char *> supported_nbc;
  std::vector<const char *> supported_access;
};

} // End namespace

#endif // _OPENKIMINFO_H
