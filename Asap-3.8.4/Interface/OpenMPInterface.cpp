// OpenMPInterface.cpp  --  Iterface controlling number of OMP threads
// -*- c++ -*-
//
// Copyright (C) 2001-2011 Jakob Schiotz and Center for Individual
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

#include "OpenMPInterface.h"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace ASAPSPACE {

char PyAsap_support_openmp_doc[] =
    "Check if ASAP is compiled with OpenMP support.";

PyObject *PyAsap_support_openmp(PyObject *noself, PyObject *noargs)
{
#ifdef _OPENMP
  Py_RETURN_TRUE;
#else
  Py_RETURN_FALSE;
#endif
}

char PyAsap_get_num_procs_doc[] =
    "Return the number of processors available to OpenMP.";

PyObject *PyAsap_get_num_procs(PyObject *noself, PyObject *noargs)
{
  int n = 1;
#ifdef _OPENMP
  n = omp_get_num_procs();
#endif // _OPENMP
  return Py_BuildValue("i", n);
}

#ifdef _OPENMP
#ifdef ASAP_FIX_AFFINITY
char PyAsap_set_num_threads_doc[] =
    "Set number of OpenMP threads.\n\n"
    "The first parameter (n) sets the number of threads.\n"
    "If the optional second parameter (noaffinity) is true, processor\n"
    "affinity is turned off.\n";
#else // ASAP_FIX_AFFINITY
char PyAsap_set_num_threads_doc[] =
    "Set number of OpenMP threads.\n\n"
    "The first parameter (n) sets the number of threads.\n"
    "The second parameter (noaffinity) is ignored, as ASAP is not compiled\n"
    "with a compiler that per default binds the process to CPU 0.\n";
#endif // ASAP_FIX_AFFINITY
#else // _OPENMP
char PyAsap_set_num_threads_doc[] =
    "Set number of OpenMP threads.\n\n"
    "As ASAP is not compiled with OpenMP support, the first parameter (n)\n"
    "must be zero, and an optional second parameter (noaffinity) is ignored.\n";
#endif // _OPENMP

PyObject *PyAsap_set_num_threads(PyObject *noself, PyObject *args,
                                 PyObject *kwargs)
{
  static char *kwlist[] = {"n", "noaffinity", NULL};
  int n;
  int noaffinity = 0;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,  "i|i:set_num_threads", kwlist,
      &n, &noaffinity))
    return NULL;
#ifdef _OPENMP
  omp_set_num_threads(n);
#ifdef ASAP_FIX_AFFINITY
  if (noaffinity)
    {
      // Some compilers (Open64) may have bound the process to CPU number 0.
      // This will totally break MPI performance and must be turned off.
      cpu_set_t cpumask;
      CPU_ZERO(&cpumask);
      int nmax = omp_get_num_procs();
      for (int i = 0; i < nmax; i++)
        CPU_SET(i, &cpumask);
      sched_setaffinity(0, sizeof(cpu_set_t), &cpumask);
    }
#endif // ASAP_FIX_AFFINITY
#else // _OPENMP
  if (n != 1)
    {
      PyErr_SetString(PyExc_ValueError,
          "No OpenMP support: Cannot set number of threads different from 1.");
      return NULL;
    }
#endif // _OPENMP
  Py_RETURN_NONE;
}

} // end namespace

