// -*- C++ -*-
// TinyMatrix.h  ---  Implements a rudimentary matrix class
//
// Copyright (C) 2001-2011 Jakob Schiotz and Center for Individual
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


#ifndef _TINYMATRIX_H
#define _TINYMATRIX_H

namespace ASAPSPACE {

template<class Type>
// A TinyMatrix is a 2D array where size is only known at runtime.
class TinyMatrix
{
public:
  /// Default constructor: Create an empty TinyMatrix
  TinyMatrix() {data = 0;}

  /// Construct a 2D array of known size.
  TinyMatrix(int rows, int columns) {Allocate(rows, columns);}

  /// Copy constructor
  TinyMatrix(const TinyMatrix<Type> &original) {CopyFrom(original);}

  void CopyFrom(const TinyMatrix<Type> &original)
  {
    Allocate(original.r, original.c);
    for (int i = 0; i < r*c; i++)
      data[i] = original.data[i];
  }

  void Allocate(int rows, int columns) {r = rows; c = columns; data = new Type[r*c];}
  ~TinyMatrix() {delete[] data;}
  Type *operator[](int row) {return data + c*row;}
  const Type *operator[](int row) const {return data + c*row;}

  void operator=(const TinyMatrix<Type> &from) {CopyFrom(from);}

protected:
  int r, c;
  Type *data;
};

/// A TinyMatrixOfVector is a TinyMatrix where each element is a vector<X>
template<class T>
class TinyMatrixOfVector : public TinyMatrix<T>
{
public:
  TinyMatrixOfVector(int rows, int columns, int size) {
    TinyMatrix<T>::Allocate(rows, columns);
    for (int i = 0; i < rows * columns; i++)
      TinyMatrix<T>::data[i].resize(size);
  }
};

typedef TinyMatrix<double> TinyDoubleMatrix;

} // end namespace

#endif // _TINYMATRIX_H
