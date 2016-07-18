// -*- C++ -*-
// IVec.h: Integer version of the 3-vector class in Vec.h
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


#ifndef __IVEC_H__
#define __IVEC_H__

#include <iostream>
using std::istream;
using std::ostream;

namespace ASAPSPACE {

/// An integer 3-vector.

/// The only data is the three positions (and there are no virtual
/// functions), so the memory layout of an array of IVecs will be x0,
/// y0, z0, x1, y1, z1, x2, ...
///
/// Almost all operations are inline for speed.

class IVec
{
public:
  /// Dummy constructor needed by STL containers.
  IVec() {};
  /// Construct a 3-vector from three ints.
  IVec(int x0, int x1, int x2);
  /// Dot product
  int operator*(const IVec& v) const;
  /// Multiplication with scalar
  IVec operator*(const int& s) const;
  /// Add two IVecs
  IVec operator+(const IVec& v) const;
  /// Subtract two IVecs
  IVec operator-(const IVec& v) const;
  /// Unary minus
  IVec operator-() const;
  /// Add a IVec to this one.
  IVec& operator+=(const IVec& v);
  /// Subtract a vec from this one.
  IVec& operator-=(const IVec& v);
  /// Multiply this vec with a scalar
  IVec& operator*=(int s);
  /// Divide this IVec with a scalar.
  IVec& operator/=(int s);
  /// IVec equality
  bool operator==(const IVec &v) const;
  /// const indexing
  int operator[](int n) const;
  /// Non-const indexing
  int& operator[](int n);
  /// Cross product of two IVecs.
  friend IVec Cross(const IVec& v1, const IVec& v2);
  /// The length of a IVec
  friend int Length2(const IVec& v);
  /// Increment y with a times x.
  friend void Vaxpy(int a, const IVec& x, IVec& y);
  /// Print a Vec
  friend ostream& operator<<(ostream& out, const Vec& v);
public:
  int x[3];  ///< The actual data.
};

inline IVec::IVec(int x0, int x1, int x2)
{
  x[0] = x0;
  x[1] = x1;
  x[2] = x2;
}

inline int IVec::operator*(const IVec& v) const
{
  return x[0] * v.x[0] + x[1] * v.x[1] + x[2] * v.x[2];
}

inline IVec IVec::operator*(const int& s) const
{
  return IVec(s * x[0], s * x[1], s * x[2]);
}

inline IVec IVec::operator+(const IVec& v) const
{
  return IVec(x[0] + v.x[0], x[1] + v.x[1], x[2] + v.x[2]);
}

inline IVec IVec::operator-(const IVec& v) const
{
  return IVec(x[0] - v.x[0], x[1] - v.x[1], x[2] - v.x[2]);
}

inline IVec IVec::operator-() const
{
  return IVec(-x[0], -x[1], -x[2]);
}

inline IVec& IVec::operator+=(const IVec& v)
{
  x[0] += v.x[0]; x[1] += v.x[1]; x[2] += v.x[2];
  return *this;
}

inline IVec& IVec::operator-=(const IVec& v)
{
  x[0] -= v.x[0]; x[1] -= v.x[1]; x[2] -= v.x[2];
  return *this;
}

inline IVec& IVec::operator*=(int s)
{
  x[0] *= s; x[1] *= s; x[2] *= s;
  return *this;
}

inline IVec& IVec::operator/=(int s)
{
  x[0] /= s; x[1] /= s; x[2] /= s;
  return *this;
}

inline bool IVec::operator==(const IVec &v) const
{
  return (x[0] == v.x[0]) && (x[1] == v.x[1]) && (x[2] == v.x[2]);
}

inline int IVec::operator[](int n) const
{
  return x[n];
}

inline int& IVec::operator[](int n)
{
  return x[n];
}

inline IVec Cross(const IVec& v1, const IVec& v2)
{
  return IVec(v1.x[1] * v2.x[2] - v1.x[2] * v2.x[1],
             v1.x[2] * v2.x[0] - v1.x[0] * v2.x[2],
             v1.x[0] * v2.x[1] - v1.x[1] * v2.x[0]);
}

inline void Vaxpy(int a, const IVec& x, IVec& y)
{
  y.x[0] += a * x.x[0];
  y.x[1] += a * x.x[1];
  y.x[2] += a * x.x[2];
}

inline int Length2(const IVec& v)
{
  return v.x[0] * v.x[0] + v.x[1] * v.x[1] + v.x[2] * v.x[2];
}

inline ostream& operator<<(ostream& out, const IVec& v)
{
  out << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")";
  return out;
}

} // end namespace

#endif // __IVEC_H__

