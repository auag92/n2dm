// -*- C++ -*-
// SymTensor.h: Symmetric tensor objects (containing six doubles).
//
// Copyright (C) 2008 Jakob Schiotz and Center for Individual
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

#ifndef _SYMTENSOR_H
#define _SYMTENSOR_H

#include <iostream>
using std::istream;
using std::ostream;

namespace ASAPSPACE {

class SymTensor
{
public:
  SymTensor() {};
  /// Set to zero
  void clear();
  /// const indexing
  double operator[](int n) const;
  /// Non-const indexing
  double& operator[](int n);
  /// Add a SymTensor to this one.
  SymTensor& operator+=(const SymTensor& v);
  /// Multiply a SymTensor with a scalar.
  SymTensor& operator*=(double x);
  /// Print a SymTensor
  friend ostream& operator<<(ostream& out, const SymTensor& v);

private:
  double x[6];  ///< The actual data.
};

inline void SymTensor::clear()
{
  x[0] = x[1] = x[2] = x[3] = x[4] = x[5] = 0.0;
}

inline double SymTensor::operator[](int n) const
{
  return x[n];
}

inline double& SymTensor::operator[](int n)
{
  return x[n];
}

inline SymTensor& SymTensor::operator+=(const SymTensor& v)
{
  x[0] += v.x[0]; x[1] += v.x[1]; x[2] += v.x[2];
  x[3] += v.x[3]; x[4] += v.x[4]; x[5] += v.x[5];
  return *this;
}

inline SymTensor& SymTensor::operator*=(double v)
{
  x[0] *= v; x[1] *= v; x[2] *= v;
  x[3] *= v; x[4] *= v; x[5] *= v;
  return *this;
}

inline ostream& operator<<(ostream& out, const SymTensor& v)
{
  out << "(" << v[0] << ", " << v[1] << ", " << v[2] << ", "
      << v[3] << ", " << v[4] << ", " << v[5] << ")";
  return out;
}



} // end namespace

#endif // _SYMTENSOR_H
