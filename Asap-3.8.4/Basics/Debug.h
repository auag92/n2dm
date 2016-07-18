// -*- C++ -*-

// Debug.h - helps debugging.

// If ASAPDEBUG is set in Asap.h, all functions print debugging info.
//
// If ASAPDEBUG is not set in Asap.h, but is set in a .cpp file just
// before including this header, methods in that file print debugging
// info.
//
// If DEBUGOBJECTS is set in Asap.h, debugging messages are printed
// when central objects are created or destroyed.
//
// If ASAPMEMORYDEBUG is set in Asap.h or in a .cpp file just before
// including this header, memory usage is printed when MEMORY is
// called.
//
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


#ifndef _DEBUG_H
#define _DEBUG_H

namespace ASAPSPACE {

#ifdef ASAPDEBUG

#ifdef _OPENMP
#import <omp.h>
#define OPENMPTHR << omp_get_thread_num() << " "
#else
#define OPENMPTHR
#endif

#ifdef __GNUC__
#define DEBUGPRINT std::cerr OPENMPTHR << "Now in " << __PRETTY_FUNCTION__ << " (" << __FILE__ << ":" << __LINE__ << ")" << std::endl << std::flush;
#else // __GNUC__
#define DEBUGPRINT std::cerr OPENMPTHR << "Now in " << __FILE__ << ":" << __LINE__ << std::endl << std::flush;
#endif // __GNUC_

#else  //ASAPDEBUG

#define DEBUGPRINT

#endif //ASAPDEBUG





#ifdef DEBUGOBJECTS

#ifdef __GNUC__
#define CONSTRUCTOR std::cerr << "Constructing object: " << __PRETTY_FUNCTION__  << " " << __FILE__ << " (" << this << ")" << std::endl << std::flush;
#define DESTRUCTOR std::cerr << "Destroying object: " << __PRETTY_FUNCTION__ << " (" << this << ")" << std::endl << std::flush;
#else // __GNUC__
#define CONSTRUCTOR std::cerr << "Constructing object: " << __FILE__ << ":" << __LINE__ << " (" << this << ")" << std::endl << std::flush;
#define DESTRUCTOR std::cerr << "Destroying object: " << __FILE__ << ":" << __LINE__ << " (" << this << ")" << std::endl << std::flush;
#endif // __GNUC_

#else // DEBUGOBJECTS

#define CONSTRUCTOR
#define DESTRUCTOR

#endif // DEBUGOBJECTS


#ifdef ASAPMEMORYDEBUG

void AsapPrintMemory(const char *file, int line);
#define MEMORY AsapPrintMemory(__FILE__, __LINE__)

#else // ASAPMEMORYDEBUG

#define MEMORY

#endif // ASAPMEMORYDEBUG


/* get heap memory using mallinfo.
   There is a UNIX version and a Mac OS X version is not well tested
   but seems to give credible values in simple tests.*/

double heap_mallinfo();

} // end namespace

#else // _DEBUG_H

#warning "Multiple inclusion of Debug.h"

#endif // _DEBUG_H
