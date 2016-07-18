// -*- C++ -*-
// ImageAtoms.cpp:  Beyond the Minimum Image Convention
//
// Copyright (C) 2014 Jakob Schiotz and Center for Individual
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

#include "ImageAtoms.h"
//#define ASAPDEBUG
#include "Debug.h"
#include <cmath>

using std::sqrt;

const bool ImageAtoms::bc_are_free[3] = {false, false, false};

ImageAtoms::ImageAtoms(Atoms *atoms)
{
  CONSTRUCTOR;
  nAtoms = nGhosts = nImages = nSize = 0;
  //bc_are_free[0] = bc_are_free[1] = bc_are_free[2] = false;
  initialized = false;
  realatoms = atoms;
  last_range = 0.0;
  AsapAtoms_INCREF(realatoms);
}

ImageAtoms::~ImageAtoms()
{
  DESTRUCTOR;
  AsapAtoms_DECREF(realatoms);
  assert(refcount == 0);
}

bool ImageAtoms::UpdateBeforeCalculation(bool flag, double range)
{
  DEBUGPRINT;
  bool update = realatoms->UpdateBeforeCalculation(flag, range);
  if (update)
    {
      make_images(range);
      last_range = range;
      update_images();
    }
  DEBUGPRINT;
  return update;
}

void ImageAtoms::Begin(PyObject *pyatoms, bool expect_reopen /*=false*/)
{
  DEBUGPRINT;
  realatoms->Begin(pyatoms, expect_reopen);
  if ((nAtoms != realatoms->GetNumberOfAtoms()) || (nGhosts != realatoms->GetNumberOfGhostAtoms()))
    make_images(last_range);
  update_images();
  DEBUGPRINT;
}

void ImageAtoms::GetPositions(vector<Vec> &pos, bool ghosts /* = false */) const
{
  DEBUGPRINT;
  pos.clear();
  int nTot = ghosts ? (nAtoms + nGhosts + nImages) : nAtoms;
  if (pos.capacity() < nTot)
    pos.reserve(nTot + nTot/25);  // 4% extra
  assert(allpositions.size() >= nTot);
  pos.insert(pos.begin(), allpositions.begin(), allpositions.begin() + nTot);
  assert(pos.size() == nTot);
  DEBUGPRINT;
}


void ImageAtoms::GetScaledPositions(vector<Vec> &pos, bool ghosts /* = false */)
{
  DEBUGPRINT;
  int n = nAtoms;
  if (ghosts)
    n += nGhosts + nImages;
  assert(allpositions.size() >= n);
  const Vec *inv = GetInverseCell();
  if (pos.capacity() < n)
    pos.reserve(n + n/25);  // Reserve 4% extra.
  pos.resize(n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < 3; j++)
      pos[i][j] = allpositions[i][0] * inv[0][j]
                + allpositions[i][1] * inv[1][j]
                + allpositions[i][2] * inv[2][j];
  DEBUGPRINT;
}

void ImageAtoms::GetScaledPositions(vector<Vec> &scaledpos, const set<int> &which)
{
  DEBUGPRINT;
  assert(scaledpos.size() == which.size());
  const Vec *inv = GetInverseCell();
  vector<Vec>::iterator spi = scaledpos.begin();
  for (set<int>::const_iterator i = which.begin(); i != which.end(); ++i,++spi)
    for (int j = 0; j < 3; j++)
      (*spi)[j] = allpositions[*i][0] * inv[0][j]
                + allpositions[*i][1] * inv[1][j]
                + allpositions[*i][2] * inv[2][j];
  DEBUGPRINT;
}

void ImageAtoms::make_images(double range)
{
  DEBUGPRINT;
  initialized = true;
  translations.clear();
  first_atom.clear();
  last_atom.clear();
  original_atom.clear();
  nImages = 0;

  // Find how many replications we need in each direction.
  const bool *pbc = realatoms->GetBoundaryConditions();
  const double *heights = realatoms->GetCellHeights();
  int nbox[3];
  for (int i = 0; i < 3; i++)
    {
      if (pbc[i])
        nbox[i] = (int) ceil(range / heights[i]);
      else
        nbox[i] = 0;
    }


  // Now create the images.  We loop over each replica of the system
  // that may contain images.
  const Vec *positions = realatoms->GetPositions();
  nAtoms = realatoms->GetNumberOfAtoms();
  nGhosts = realatoms->GetNumberOfGhostAtoms();
  const Vec *cell = realatoms->GetCell();
  Vec toppoint(cell[0] + cell[1] + cell[2]);

  // Normalize and make sure the sign is right (we cannot assume
  // a right-handed unit cell.
  Vec directions[3];
  for (int i = 0; i < 3; i++)
    {
      int j = (i + 1) % 3;
      int k = (i + 2) % 3;
      directions[i] = Cross(cell[j], cell[k]);
      directions[i] *= 1.0/sqrt(Length2(directions[i]));
      if (directions[i] * cell[i] < 0.0)
        directions[i] *= -1.0;
    }

  const Vec *r = realatoms->GetPositions();
  for (int i = -nbox[0]; i <= nbox[0]; i++)
    for (int j = -nbox[1]; j <= nbox[1]; j++)
      for (int k = -nbox[2]; k <= nbox[2]; k++)
        {
          if ((i == 0) && (j == 0) && (k == 0))
            continue;   // Central box.  Not images
          IVec tranl(i, j, k);
          translations.push_back(tranl);
          first_atom.push_back(original_atom.size());
          for (int n = 0; n < nAtoms + nGhosts; n++)
            {
              Vec rr = r[n] + tranl[0] * cell[0] + tranl[1] * cell[1] + tranl[2] * cell[2];
              Vec r2 = rr - toppoint;
              // Check the six limits, skip atom if any test fails
              // Check a lower limit, direction 0
              if ( (i < 0) && (rr * directions[0] < - range) )
                continue;
              // Check upper limit, direction 0
              if ( (i > 0 ) && (r2 * directions[0] > range) )
                continue;
              // Check a lower limit, direction 1
              if ( (j < 0) && (rr * directions[1] < - range) )
                continue;
              // Check upper limit, direction 1
              if ( (j > 0 ) && (r2 * directions[1] > range) )
                continue;
              // Check a lower limit, direction 2
              if ( (k < 0) && (rr * directions[2] < - range) )
                continue;
              // Check upper limit, direction 2
              if ( (k > 0 ) && (r2 * directions[2] > range) )
                continue;
              // Passed all tests, include this atom.
              original_atom.push_back(n);
            }
          last_atom.push_back(original_atom.size());
        }
  nImages = original_atom.size();
  // Reserve the necessary space in allpositions and allnumbers
  nSize = nAtoms + nGhosts + nImages;
  if (allpositions.capacity() < nSize)
    allpositions.reserve(nSize + nSize/25);
  if (allnumbers.capacity() < nSize)
    allnumbers.reserve(nSize + nSize/25);
  DEBUGPRINT;
}

void ImageAtoms::update_images()
{
  DEBUGPRINT;
  if (!initialized)
    {
      DEBUGPRINT;
      return;  // This happens when Begin is called from the initial SetAtoms call.
    }
  // Get the positions and atomic numbers of the real atoms and the ghosts.
  realatoms->GetPositions(allpositions, true);
  assert(allpositions.size() <= nSize);
  allpositions.resize(nSize);
  allnumbers.resize(nSize);
  const asap_z_int *z = realatoms->GetAtomicNumbers();
  for (int i = 0; i < nAtoms + nGhosts; i++)
    allnumbers[i] = z[i];
  // Now update the images.  Both positions and atomic numbers may have changed.
  int target = nAtoms + nGhosts;
  int source = 0;
  int ngroups = translations.size();
  const Vec *r = realatoms->GetPositions();
  const Vec *cell = realatoms->GetCell();

  for (int i = 0; i < ngroups; i++)
    {
      const IVec &tr = translations[i];
      for (int n = first_atom[i]; n < last_atom[i]; n++)
        {
          allpositions[target] = r[original_atom[source]] + tr[0] * cell[0] + tr[1] * cell[1] + tr[2] * cell[2];
          allnumbers[target] = z[original_atom[source]];
          source++;
          target++;
        }
    }
  assert(target == nSize);
  DEBUGPRINT;
}
