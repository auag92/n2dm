// -*- C++ -*-
// Atoms.h:  The interface to the ase Atoms object.
//
// Copyright (C) 2008-2011 Jakob Schiotz and Center for Individual
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


#ifndef ATOMS_H
#define ATOMS_H

#include "AsapPython.h"
#include "Asap.h"
#include <vector>
#include <set>
#include "Vec.h"
#include "IVec.h"

using std::vector;
using std::set;

namespace ASAPSPACE {

/// The Atoms object interfaces to the ase-3 Python Atoms object.
///
/// Before accessing the atoms, open them by calling
/// Atoms::Begin(pyatoms), this extracts information from the Python
/// object.  The remaining methods will give access to this
/// information until Atoms::End() is called.
///
/// It is an error (but not a serious one) to forget calling End()
/// before the next time Begin() is called.  This should only happen
/// if a previous access was aborted by an exception, but is allowed
/// in that case.  Calling End() twice in a row is a serious error, as
/// it is a symptom of two different objects accessing (possibly
/// different) atoms though this object simultaneously.
///
/// It gives the best performance if there is a one-to-one mapping
/// between Atoms objects and Python Atoms objects, but several Atoms
/// objects may interface a python object simultaneously, and the same
/// Atoms object may be used to access different python objects
/// sequentially.
///
/// This object is refcounted, see AsapAtoms_INCREF() and AsapAtoms_DECREF().
class Atoms
{
protected:
  /// Delete the interface object.
  ///
  /// Protected: may only be called from AsapAtoms_DECREF() friend function.
  virtual ~Atoms() {assert(refcount == 0);}
  
public:
  /// Create the interface object.
  Atoms() {refcount = 1;}

  // Reference counting
  friend void AsapAtoms_INCREF(Atoms *atoms);
  friend void AsapAtoms_DECREF(Atoms *atoms);
  
  /// Start accessing the Python atoms.
  ///
  /// Until End() is called, references are kept to the Python atoms
  /// and some of their attributes.
  /// If expect_reopen=true, further calls to Begin will be allowed as
  /// long as they are matched by corresponding calls to End().
  virtual void Begin(PyObject *pyatoms, bool expect_reopen=false) = 0;

  /// Finish accessing the Python atoms, freeing resources.
  virtual void End() = 0;
  
  /// Are the atoms active (opened by Begin()) ?
  virtual bool IsActive() const = 0;
  
  /// Get the number of atoms
  virtual int GetNumberOfAtoms() const = 0;

  /// Get the number of ghost atoms
  virtual int GetNumberOfGhostAtoms() const = 0;
  
  /// Get total number of atoms in parallel simulation.
  ///
  /// In a serial simulation, same as GetNumberOfAtoms, in a parallel simulation
  /// it is the sum over these over all processors.
  virtual int GetTotalNumberOfAtoms() const = 0;

  /// Do we have ghost atoms? True if arrays are present, even if empty.
  virtual bool HasGhostAtoms() const = 0;

  /// Return three integers specifying the CPU layout.
  virtual IVec GetNumberOfCells() const = 0;

  /// Get the cartesian positions.  Ghost positions, if any, are at the end.
  virtual const Vec *GetPositions() const = 0;

  /// Get a copy of the cartesian positions.
  ///
  /// Include ghosts if possible and requested.
  virtual void GetPositions(vector<Vec> &pos, bool ghosts=false) const = 0;
  
  /// Get a copy of the positions converted to scaled space.
  ///
  /// Include ghosts if possible and requested.
  virtual void GetScaledPositions(vector<Vec> &scaledpos, bool ghosts=false) = 0;

  /// Get a copy of some positions, converted to scaled space.
  virtual void GetScaledPositions(vector<Vec> &scaledpos, const set<int> &which) = 0;

  /// Get the atomic numbers
  virtual const asap_z_int *GetAtomicNumbers() const = 0;
  
  /// Get the cartesian momenta
  virtual const Vec *GetMomenta() = 0;

  /// Get the boundary conditions.
  virtual const bool *GetBoundaryConditions() const = 0;

  virtual bool AllFreeBoundaries() const = 0;

  virtual bool AllPeriodicBoundaries() const = 0;

  /// Get the supercell
  virtual const Vec *GetCell() const = 0;

  /// Get the supercell.   No sanity check!
  virtual const Vec *GET_CELL() const = 0;

  /// Get the volume of the supercell
  virtual double GetVolume() const = 0;
  
  /// Get the height of the supercell
  virtual const double *GetCellHeights() = 0;
  
  /// Get the inverse supercell
  virtual const Vec *GetInverseCell() = 0;
  
  /// Get a set of all elements present in the simulations.
  virtual void GetListOfElements(set<int> &elements) const = 0;

  /// Get the masses of the atoms
  virtual const double *GetMasses() = 0;

  /// Set a python array on the atoms.
  virtual void SetData(const char *name, PyObject *data) = 0;

  /// Get a python array on the atoms.  Returns a new reference.
  virtual PyObject *GetData(const char *name) = 0;

  /// Remove a python array from the atoms
  virtual void DeleteData(const char *name) = 0;
  
  /// Update flag and counters across processors.
  ///
  /// Called by a Potential with a flag indicating if the neighborlist
  /// should be updated.  In a serial simulation, just return the flag.
  ///
  /// In a parallel simulation, communicate across processors so the
  /// counters of this Atoms object and the flag passed as an argument
  /// are updated if just one processor thinks an update is necessary.
  /// If the flag is true, a migration will also be triggered in a
  /// parallel simulation.
  ///
  /// The second argument is the range of the neighbor locator
  /// including the drift.  It is used when attaching ghost atoms in
  /// parallel simulations.
  virtual bool UpdateBeforeCalculation(bool flag, double range) = 0;
  
  /// Communicate some kind of data to the ghost atoms.
  ///
  /// This is used if a potential needs to communicate some kind of
  /// date during the energy/force calculations.  If the data should
  /// be communicated at the beginning of the calculation, it should
  /// instead be placed in a ghost array.
  ///
  /// NB: In a serial simulation, do nothing except checking that no
  /// ghosts are present.
  virtual void CommunicateData(double *address, int n = 1) = 0;
  
  /// Counter showing changes to the number of atoms or boundary conditions
  virtual int GetAtomsCounter() const = 0;

  /// Counter showing changes to the unit cell
  ///
  /// Also updated if AtomsCounter has changed.
  virtual int GetCellCounter() const = 0;

  /// Counter showing changes to the atomic numbers.
  ///
  /// Also updated if AtomsCounter has changed.
  virtual int GetNumbersCounter() const = 0;
  
  /// Counter showing changes to the positions.
  ///
  /// Also updated if CellCounter, NumbersCounter or AtomsCounter is updated.
  virtual int GetPositionsCounter() const = 0;

  /// Counter showing changes to the momenta.
  ///
  /// Reading this counter will cause the momenta to be extracted from Python.
  virtual int GetMomentaCounter() = 0;

  /// Print memory usage
  virtual long PrintMemory() const = 0;
  
protected:
  int refcount;           ///< Number of references to this object.
};

/// Increase the reference count of an Atoms object
inline void AsapAtoms_INCREF(Atoms *atoms)
{
  atoms->refcount++;
}

/// Decrease the reference count of Atoms, deallocate if it reaches zero.
inline void AsapAtoms_DECREF(Atoms *atoms)
{
  atoms->refcount--;
  if (atoms->refcount == 0)
    delete atoms;
}

} // end namespace

#endif // ATOMS_H
