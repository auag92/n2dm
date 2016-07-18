// -*- C++ -*-
// KimParameterProvider.h:  Get parameters from a string (for OpenKIM).
//
// Copyright (C) 2008-2012 Jakob Schiotz and Center for Individual
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

#ifndef _EMTPYTHONPARAMETERPROVIDER_H
#define _EMTPYTHONPARAMETERPROVIDER_H

#include "EMTDefaultParameterProvider.h"

namespace ASAPSPACE {

class KimParameterProvider : public EMTDefaultParameterProvider
{
public:
  KimParameterProvider(const char *parameter_filename, intptr_t *pkim);
  //virtual ~EMTPythonParameterProvider() {};  // Not changed.
  
  virtual string GetName() const {return "KimParameterProvider";}

protected:
  virtual emt_parameters *GetNewParameters(int element);
};

} // end namespace

#endif // _EMTPYTHONPARAMETERPROVIDER_H
