// -*- C++ -*-
// KimAtoms.cpp:  Interface to KIM pretending to be the atoms.
//
// Copyright (C) 2012-2013 Jakob Schiotz and the Department of Physics,
// Technical University of Denmark.  Email: schiotz@fysik.dtu.dk
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

// This is a reimplementation of the ASAP Atoms object.  It does not
// inherit from the usual Atoms object, as Atoms.h include Python.h
// which is not available in OpenKIM.

#include "KimAtoms.h"
#include "Exception.h"
#include "Debug.h"
#include <math.h>

KimAtoms::KimAtoms(intptr_t* pkim)
{
  CONSTRUCTOR;
  assert(pkim != NULL);
  refcount = 1;
  this->pkim = pkim;
  counter = 2;
  count_inverse_cell = 1;
}

KimAtoms::~KimAtoms()
{
  DESTRUCTOR;
}

void KimAtoms::ReInit(int nAtoms, int nGhosts, double *pos, int *z)
{
  this->nAtoms = nAtoms;
  this->nGhosts = nGhosts;
  nAllAtoms = nAtoms + nGhosts;
  positions.resize(nAllAtoms);
  numbers.resize(nAllAtoms);
  for (int i = 0; i < nAllAtoms; i++)
    {
      positions[i] = ((Vec *)pos)[i];
      numbers[i] = z[i];
    }
  counter++;
  for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
        cell[i][j] = (i == j) * 50.0;
    }
}

void KimAtoms::GetListOfElements(set<int> &elements) const
{
  DEBUGPRINT;
  const asap_z_int *atomicnumbers = GetAtomicNumbers();

  elements.clear();
  for (int i = 0; i < nAllAtoms; i++)
    {
      int z = atomicnumbers[i];
      if (elements.find(z) == elements.end())
        elements.insert(z);
    }
  DEBUGPRINT;
}

/// Get the cartesian momenta
const Vec *KimAtoms::GetMomenta()
{
  throw AsapError("GetMomenta called in KIM mode.");
  return NULL;
}


/// Get the cartesian momenta
const double *KimAtoms::GetMasses()
{
  throw AsapError("GetMasses called in KIM mode.");
  return NULL;
}

double KimAtoms::GetVolume() const
{
  DEBUGPRINT;
  double det;
  det = -cell[0][2]*cell[1][1]*cell[2][0] +
    cell[0][1]*cell[1][2]*cell[2][0] +
    cell[0][2]*cell[1][0]*cell[2][1] -
    cell[0][0]*cell[1][2]*cell[2][1] -
    cell[0][1]*cell[1][0]*cell[2][2] +
    cell[0][0]*cell[1][1]*cell[2][2];
  DEBUGPRINT;
  return fabs(det);
}

void KimAtoms::GetPositions(vector<Vec> &pos, bool ghosts /* = false */) const
{
  DEBUGPRINT;
  pos.clear();
  int n = nAtoms;
  if (ghosts)
    n = nAllAtoms;

  pos.reserve(n + n/25);  // 4% extra
  for (int i = 0; i < n; i++)
    pos[i] = positions[i];

  DEBUGPRINT;
}

void KimAtoms::GetScaledPositions(vector<Vec> &pos, bool ghosts /* = false */)
{
  DEBUGPRINT;
  int n = nAtoms;
  if (ghosts)
    n += nGhosts;
  const Vec *inv = GetInverseCell();
  if (pos.capacity() < n)
    pos.reserve(n + n/25);  // Reserve 4% extra.
  pos.resize(n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < 3; j++)
      pos[i][j] = positions[i][0] * inv[0][j]
                + positions[i][1] * inv[1][j]
                + positions[i][2] * inv[2][j];
  DEBUGPRINT;
}

void KimAtoms::GetScaledPositions(vector<Vec> &scaledpos, const set<int> &which)
{
  DEBUGPRINT;
  assert(scaledpos.size() == which.size());
  const Vec *inv = GetInverseCell();
  vector<Vec>::iterator spi = scaledpos.begin();
  for (set<int>::const_iterator i = which.begin(); i != which.end(); ++i,++spi)
    for (int j = 0; j < 3; j++)
      (*spi)[j] = positions[*i][0] * inv[0][j]
                + positions[*i][1] * inv[1][j]
                + positions[*i][2] * inv[2][j];
  DEBUGPRINT;
}

const double *KimAtoms::GetCellHeights()
{
  DEBUGPRINT;
  if (count_inverse_cell < counter)
    invert_cell();
  return heights;
}

const Vec *KimAtoms::GetInverseCell()
{
  DEBUGPRINT;
  if (count_inverse_cell < counter)
    invert_cell();
  return inverse;
}

void KimAtoms::invert_cell()
{
  DEBUGPRINT;
  count_inverse_cell = counter;
  double determinant = Cross(cell[0], cell[1]) * cell[2];
  // Find heights
  for (int i = 0; i < 3; i++)
    {
      Vec inv = Cross(cell[(i + 1) % 3], cell[(i + 2) % 3]);
      heights[i] = fabs(determinant) / sqrt(Length2(inv));
    }
  // Invert matrix.  I_i,j = { C_j-1,i-1 C_j+1,i+1 - C_j+1,i-1 C_j-1,i+1 } / det
  for (int i = 0; i < 3; i++)
    {
      int ip = (i + 1) % 3;
      int im = (i + 2) % 3;
      for (int j = 0; j < 3; j++)
        {
          int jp = (j + 1) % 3;
          int jm = (j + 2) % 3;
          inverse[i][j] = (cell[jm][im]*cell[jp][ip] -
              cell[jp][im]*cell[jm][ip]) / determinant;
        }
    }
  DEBUGPRINT;
}

void KimAtoms::SetDiagonalCell(double d[3])
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      if (i == j)
        cell[i][j] = d[i];
      else
        cell[i][j] = 0.0;
  count_inverse_cell = 0;  // Must be recalculated
}

