/// \file
// -*- C++ -*-
// Asap.h  -- Global Asap header file
//
// Asap.h contains global configuration macros and a few widely used
// macros.
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



#ifndef _ASAP_H
#define _ASAP_H

///
/// Global Asap header file, containing a few defines etc

// A C++ namespace is defined for all of Asap to reduce the likelihood
// of symbol clashes.  It is defined in a separate header file so it can
// be changed when Asap source files are used in e.g. OpenKIM
#include "AsapNamespace.h"

/// The batch-buffer size in EMT.cpp.
#define BUFLEN 10000

/// Maximal number of elements in a single simulation (EMT.cpp)
#define NMAXELEMENTS 8

/// Is timing enabled?
//#define TIMING

/// Do we want to catch floating point exceptions
///
/// Recommended on architectures where Asap compiles with this option
/// enabled, as it will cause a crash immediately rather than letting
/// a simulation continue forever while more and more data becomes
/// Not-A-Number.  There is no measurable performance penalty when
/// enabling this option.  Does not work on Mac OS X
///
/// CURRENTLY DISABLED: This causes problems with OpenKIM.

#ifndef _WIN32
#ifndef __APPLE__
//#define CATCHFLOAT
#endif
#endif

/// Do we want to use our own assert macro?
#ifndef NDEBUG        // Never if NDEBUG is defined
#define ASAPASSERT
#endif

/// Do we want asserts that cause a (small) performance penalty
///
/// The performance penalty is approximately 2 percent.
#define SLOWASSERT

// Some openMP compilers (Open64 for example) turns on processor
// affinity per default.  This will break MPI calculations where all
// tasks will be forced to share the same CPU.  Defining ASAP_FIX_AFFINITY
// explicitly turns off affinity when asap loads, unless the user has set
// environment variables controlling affinity (handled in Python).
// Defining ASAP_FIX_AFFINITY probably only works on Linux, and will be
// ignored if OpenMP is not used.
#ifdef __OPEN64__
#ifdef __linux__
#define ASAP_FIX_AFFINITY
#endif // __linux__
#endif // __OPEN64__

/// Is debugging turned on globally?
//#define ASAPDEBUG

/// Do we want debugging info printed when objects are created or destroyed?
//#define DEBUGOBJECTS

/// Do we want to debug memory usage
//#define ASAPMEMORYDEBUG

// The following symbol definitions should not be changed
#ifdef SLOWASSERT
#define ASSERT(x) assert(x);
#else
#define ASSERT(x)
#endif

#define CHECKNAN(x) assert(!isnan(x) && !isinf(x))

// Define the size of the integer used to hold atomic numbers
#define asap_z_int npy_int32
#define ASAP_Z_ARRAYTYPE NPY_INT32
#define ASAP_Z_SIZE NPY_SIZEOF_INT32


// A single global variable
namespace ASAPSPACE {
  extern int verbose;
}
using ASAPSPACE::verbose;

// Macros for controling visibility of symbols, for use with -fvisibility=hidden
#if defined _WIN32 || defined __CYGWIN__
  // Prepare for the unlikely case that someone will try to compile for Windoze.
  #define ASAP_PUBLIC __declspec(dllexport)
  #define ASAP_LOCAL
#else  // _WIN32
  // Test for GNU C++ version 4.0 or higher, or compatibel Intel compiler.
  #if __GNUC__ >= 4
    #define ASAP_PUBLIC __attribute__ ((visibility("default")))
    #define ASAP_LOCAL  __attribute__ ((visibility("hidden")))
  #else
    #define ASAP_PUBLIC
    #define ASAP_LOCAL
  #endif //GCC 4
#endif  // Win or unix

// Not all platforms has defined pi using this macro.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace ASAPSPACE;

#endif // _ASAP_H
