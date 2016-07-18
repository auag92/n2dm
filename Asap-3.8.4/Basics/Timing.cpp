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
#include "Timing.h"
#include "TimingResults.h"
#include <map>
#include <string>

#ifdef TIMING
static PyObject *maketimingresults(std::map<string, Timing_timer *> &pool,
                                   int level = 0)
{
  if (level == 3) {
    PyErr_SetString(PyExc_RuntimeError,
                    "Too deep recursion in maketimingresults");
    return NULL;
  }
  int nelem = pool.size();
  if (nelem == 0) {
    Py_INCREF(Py_None);
    return Py_None;
  }
  PyObject *result = PyDict_New();
  if (!result)
    return NULL;
  std::map<string, Timing_timer *>::const_iterator p;
  for (p = pool.begin(); p != pool.end(); ++p)
    {
      std::string name;
      double wall, cpu, wall_c, cpu_c;
      int count;
      p->second->Report(name, wall, cpu, wall_c, cpu_c, count);
      PyObject *key = Py_BuildValue("s", name.c_str());
      PyObject *value = Py_BuildValue("(ddddiN)", wall, cpu, wall_c, cpu_c,
                                      count,
                                      maketimingresults(p->second->masters,
                                                        level+1));
      if (!key || !value) {
        Py_DECREF(result); Py_XDECREF(key); Py_XDECREF(value);
        return NULL;
      }
      int x = PyDict_SetItem(result, key, value);
      Py_DECREF(key);
      Py_DECREF(value);
      if (x == -1) {
        Py_DECREF(result);
        return NULL;
      }
    }
  return result;
}  
#endif //TIMING

namespace ASAPSPACE {

PyObject *PyAsap_TimingResults(PyObject *noself, PyObject *noargs) {
#ifdef TIMING
  return maketimingresults(Timing_timerpool);
#else // TIMING
  Py_INCREF(Py_None);
  return Py_None;
#endif //TIMING
}

} // end namespace
