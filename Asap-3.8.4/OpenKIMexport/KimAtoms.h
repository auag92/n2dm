// -*- C++ -*-
// KimAtoms.h:  Interface to KIM pretending to be the atoms.
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


#ifndef KIMATOMS_H
#define KIMATOMS_H

#include "KimAsapPython.h"
#include "Asap.h"
#include "Vec.h"
#include "IVec.h"
#include <set>
#include <vector>

using std::set;
using std::vector;

namespace ASAPSPACE {

class NeighborLocator;

// Supported neighbor list types
typedef enum {asapkim_cluster, asapkim_opbc_h} asapkim_nbtype;

class KimAtoms  // Will be renamed Atoms in a typedef !
{
protected:
  /// Delete the interface object.
  ///
  /// Protected: may only be called from AsapAtoms_DECREF() friend function.
  virtual ~KimAtoms();

public:
  /// Create the interface object.
  KimAtoms(intptr_t* pkim = NULL);

  /// Pass new KIM pointers to this object.
  void ReInit(int nAtoms, int nGhost, double *pos, int *z);

  // Reference counting
  friend void AsapAtoms_INCREF(KimAtoms *atoms);
  friend void AsapAtoms_DECREF(KimAtoms *atoms);

  /// Begin and End do nothing
  virtual void Begin(void *pyatoms, bool expect_reopen=false) {};
  virtual void End() {};

  /// Atoms are always active
  inline bool IsActive() const {return true;}

  /// Get the number of atoms
  inline int GetNumberOfAtoms() const {return nAtoms;}

  inline int GetNumberOfGhostAtoms() const {return nGhosts;}

  /// Get total number of atoms in parallel simulation.
  ///
  /// In a serial simulation, same as GetNumberOfAtoms, in a parallel simulation
  /// it is the sum over these over all processors.
  virtual int GetTotalNumberOfAtoms() const {return GetNumberOfAtoms();}

  // Note: HasGhostAtoms may have to be improved when parallel simulations
  // are supported, see ticket #56, which is fixed for "normal" Asap use,
  // but not necessarily for OpenKIM use (will depend on how parallelization
  // is finally handled).
  inline bool HasGhostAtoms() const {return nGhosts != 0;}

  /// Get the cartesian positions.  Ghost positions, if any, are at the end.
  inline const Vec *GetPositions() const {return &positions[0];}

  void GetPositions(vector<Vec> &pos, bool ghosts=false) const;

  void GetScaledPositions(vector<Vec> &scaledpos, bool ghosts=false);

  /// Get a copy of some positions, converted to scaled space.
  void GetScaledPositions(vector<Vec> &scaledpos, const set<int> &which);

  /// Get the atomic numbers
  const asap_z_int *GetAtomicNumbers() const {return &numbers[0];}

  /// Print memory usage
  virtual long PrintMemory() const {return 0;}

  /// Get a set of all elements present in the simulations.
  virtual void GetListOfElements(set<int> &elements) const;

  int GetPositionsCounter() const {return counter;}
  int GetMomentaCounter() const {return counter;}
  int GetNumbersCounter() const {return counter;}
  int GetCellCounter() const {return counter;}

  /// Get the cartesian momenta
  const Vec *GetMomenta();

  const double *GetMasses();

  bool UpdateBeforeCalculation(bool flag, double range) {return flag;}

  virtual void CommunicateData(double *address, int n = 1) {}

  /// Return three integers specifying the CPU layout.
  virtual IVec GetNumberOfCells() const {return IVec(1,1,1);}

  /// Get the supercell
  const Vec *GetCell() const {return cell;}

  /// Get the supercell.   No sanity check!
  const Vec *GET_CELL() const {return cell;}

  /// Get the volume of the supercell
  double GetVolume() const;

  /// Get the height of the supercell
  const double *GetCellHeights();

  /// Get the inverse supercell
  const Vec *GetInverseCell();

  const bool *GetBoundaryConditions() const {return pbc;}

  void SetPBC(bool x, bool y, bool z) {pbc[0] = x; pbc[1] = y; pbc[2] = z;}

  void SetDiagonalCell(double d[3]);

private:
  void invert_cell();

public:
  // Public data (naughty!)
  int refcount;           ///< Number of references to this object.
  NeighborLocator *kim_interface_nblist;

private:
  // Data
  intptr_t* pkim;   // KIM API object.
  int nAtoms;  // The number of atoms
  int nGhosts;
  int nAllAtoms;
  vector<Vec> positions;  ///< A copy of the positions of the atoms.
  vector<asap_z_int> numbers;    ///< A copy of the atomic numbers.
  int counter, count_inverse_cell;
  Vec cell[3];            ///< The unit cell
  Vec inverse[3];         ///< The inverse unit cell
  double heights[3];      ///< Heights of the unit cell
  bool pbc[3];
  asapkim_nbtype nbtype;
};


/// Increase the reference count of an Atoms object
inline void AsapAtoms_INCREF(KimAtoms *atoms)
{
  atoms->refcount++;
}

/// Decrease the reference count of Atoms, deallocate if it reaches zero.
inline void AsapAtoms_DECREF(KimAtoms *atoms)
{
  atoms->refcount--;
  if (atoms->refcount == 0)
    delete atoms;
}

} // end namespace

// Fake the real Atoms object
#define Atoms KimAtoms
#define NormalAtoms KimAtoms

#endif // KIM_ATOMS_H
