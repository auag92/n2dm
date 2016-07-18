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

#include "AsapPython.h"
#include "Exception.h"
#include <iostream>
#include "Debug.h"
#ifdef _OPENMP
#include <omp.h>
#endif

AsapError::AsapError(const char *m) : message() 
{
  message << m;
}

AsapError::AsapError(const AsapError& ex) : message() 
{
  message << ex.GetMessage();
}

AsapError::~AsapError() 
{
}
string AsapError::GetMessage() const 
{
  return message.str();
}

AsapNotImplementedError::AsapNotImplementedError(const char *m)
{
  PyErr_SetString(PyExc_NotImplementedError, m);
}

AssertionFailed::AssertionFailed(const char *expression, const char *file,
				 int line, const char *func)
  : message() 
{
  DEBUGPRINT;
  message << file << ":" << line << ": ";
  if (func)
    message << func << ": ";
  message << "Assertion '" << expression << "' failed.";
  std::cerr << message.str() << std::endl;
  DEBUGPRINT;
}

AssertionFailed::AssertionFailed(const AssertionFailed& ex)
  : message() 
{
  DEBUGPRINT;
  message << ex.GetMessage();
}

AssertionFailed::~AssertionFailed() 
{
}
string AssertionFailed::GetMessage() const 
{
  DEBUGPRINT;
  return message.str();
}

#ifdef _OPENMP

namespace ASAPSPACE {

void AsapThrowHelper(const AsapError &err)
{
#pragma omp critical
  {
    if(AsapGlobalException == NULL)
      {
        AsapGlobalException = new AsapError(err);
        std::cerr << "\n\nASAP EXCEPTION IN OPEN-MP CODE:" << std::endl;
        std::cerr << err.GetMessage() << std::endl;
        std::cerr << "ATTEMPTING TO PROPAGATE EXCEPTION (May fail due to OpenMP)" << std::endl;
      }
  }
}

void AsapThrowHelper(const AsapPythonError &err)
{
#pragma omp critical
  {
    if(AsapGlobalException == NULL)
      {
        AsapGlobalException = new AsapPythonError(err);
        std::cerr << "\n\nASAP EXCEPTION IN OPEN-MP CODE:" << std::endl;
        std::cerr << "  -  a Python error occurred! " << std::endl;
        std::cerr << "ATTEMPTING TO PROPAGATE EXCEPTION (May fail due to OpenMP)" << std::endl;
      }
  }
}

} // end namespace

#endif // _OPENMP
