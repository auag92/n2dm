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
#include "Debug.h"
#include "Exception.h"
#include <stdio.h>
#include <iostream>
using std::cerr;
using std::endl;
#include <stdlib.h>
#ifdef __DARWIN_UNIX03
/* Allows for special MaxOS magic */
#include <malloc/malloc.h>
#endif
#ifdef __linux__
/* stdlib.h does not define mallinfo (it should!) */
#include <malloc.h>
#endif

namespace ASAPSPACE {

#ifdef ASAPMEMORYDEBUG

void AsapPrintMemory(const char *file, int line)
{
  static char buf[500];
  snprintf(buf, 500, "%s:%d", file, line);
  buf[499] = '\0';  // Paranoia
  PyObject *modules = PyImport_GetModuleDict();
  if (modules == NULL)
    throw AsapPythonError();
  PyObject *asapmodule = PyMapping_GetItemString(modules, "asap3");
  if (asapmodule == NULL)
    throw AsapPythonError();
  PyObject *printmemory = PyObject_GetAttrString(asapmodule, "print_memory");
  Py_DECREF(asapmodule);
  if (printmemory == NULL)
    throw AsapPythonError();
  PyObject *retval = PyObject_CallFunction(printmemory, "s", buf);
  Py_DECREF(printmemory);
  if (retval == NULL)
    throw AsapPythonError();
  Py_DECREF(retval);
}

#endif //ASAPMEMORYDEBUG

/* get heap memory using mallinfo.
   There is a UNIX version and a Mac OS X version is not well tested
   but seems to give credible values in simple tests.*/
double heap_mallinfo()
{
  double heap;
#ifdef __linux__
  unsigned int mmap, arena, small;
  struct mallinfo mi; /* structure in bytes */

  mi = mallinfo();
  mmap = mi.hblkhd;
  arena = mi.uordblks;
  small = mi.usmblks;
  heap = ((double)(mmap + arena + small))/1024.0; /* convert to KB */
#elif defined(__DARWIN_UNIX03)
  /* Mac OS X specific hack */
  struct malloc_statistics_t mi; /* structure in bytes */

  malloc_zone_statistics(NULL, &mi);
  heap = ((double)(mi.size_in_use))/1024.0; /* convert to KB */
#else
  heap = -1;
#endif
  return heap;
}

} // end namespace
