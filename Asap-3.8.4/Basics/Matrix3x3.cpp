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

// -*- C++ -*-
// Asap:  Copyright (C) 2008 CINF/CAMD and Jakob Schiotz

// Helper functions for 3x3 matrices defined as Vec[3].

#include "Matrix3x3.h"

namespace ASAPSPACE {

void matrixMultiply3x3(Vec product[3], const Vec a[3], const Vec b[3])
{
  product[0] = a[0][0] * b[0] + a[0][1] * b[1] + a[0][2] * b[2];
  product[1] = a[1][0] * b[0] + a[1][1] * b[1] + a[1][2] * b[2];
  product[2] = a[2][0] * b[0] + a[2][1] * b[1] + a[2][2] * b[2];
}

}
