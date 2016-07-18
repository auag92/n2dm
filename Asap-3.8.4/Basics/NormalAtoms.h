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


#ifndef NORMALATOMS_H
#define NORMALATOMS_H

#include "AsapPython.h"
#include "Asap.h"
#include "AtomsBasis.h"
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
class NormalAtoms : public Atoms
{
protected:
  /// Delete the interface object.
  ///
  /// Protected: may only be called from AsapAtoms_DECREF() friend function.
  virtual ~NormalAtoms();
  
public:
  /// Create the interface object.
  NormalAtoms();

  // Reference counting
  friend void AsapAtoms_INCREF(Atoms *atoms);
  friend void AsapAtoms_DECREF(Atoms *atoms);
  
  /// Start accessing the Python atoms.
  ///
  /// Until End() is called, references are kept to the Python atoms
  /// and some of their attributes.
  /// If expect_reopen=true, further calls to Begin will be allowed as
  /// long as they are matched by corresponding calls to End().
  virtual void Begin(PyObject *pyatoms, bool expect_reopen=false);

  /// Finish accessing the Python atoms, freeing resources.
  virtual void End();
  
  /// Are the atoms active (opened by Begin()) ?
  virtual bool IsActive() const {return active;}
  
  /// Get the number of atoms
  virtual int GetNumberOfAtoms() const {ASSERT(active) return nAtoms;}

  /// Get the number of ghost atoms
  virtual int GetNumberOfGhostAtoms() const {ASSERT(active) return nGhosts;}
  
  /// Get total number of atoms in parallel simulation.
  ///
  /// In a serial simulation, same as GetNumberOfAtoms, in a parallel simulation
  /// it is the sum over these over all processors.
  virtual int GetTotalNumberOfAtoms() const {return GetNumberOfAtoms();}

  /// Do we have ghost atoms? True if arrays are present, even if empty.
  virtual bool HasGhostAtoms() const {return hasGhosts;}

  /// Return three integers specifying the CPU layout.
  virtual IVec GetNumberOfCells() const {return IVec(1,1,1);}

  /// Get the cartesian positions.  Ghost positions, if any, are at the end.
  virtual const Vec *GetPositions() const {ASSERT(active) return &positions[0];}

  /// Get a copy of the cartesian positions.
  ///
  /// Include ghosts if possible and requested.
  virtual void GetPositions(vector<Vec> &pos, bool ghosts=false) const;
  
  /// Get a copy of the positions converted to scaled space.
  ///
  /// Include ghosts if possible and requested.
  virtual void GetScaledPositions(vector<Vec> &scaledpos, bool ghosts=false);

  /// Get a copy of some positions, converted to scaled space.
  virtual void GetScaledPositions(vector<Vec> &scaledpos, const set<int> &which);

  /// Get the atomic numbers
  virtual const asap_z_int *GetAtomicNumbers() const
  {ASSERT(active) assert(numbers.size() > 0); return &numbers[0];}
  
  /// Get the cartesian momenta
  virtual const Vec *GetMomenta();

  /// Get the boundary conditions.
  virtual const bool *GetBoundaryConditions() const
  {ASSERT(active) return pbc;}

  virtual bool AllFreeBoundaries() const
  {ASSERT(active) return all_free_boundaries;}

  virtual bool AllPeriodicBoundaries() const
  {ASSERT(active) return all_periodic_boundaries;}

  /// Get the supercell
  virtual const Vec *GetCell() const {ASSERT(active) return cell;}

  /// Get the supercell.   No sanity check!
  virtual const Vec *GET_CELL() const {return cell;}

  /// Get the volume of the supercell
  virtual double GetVolume() const;
  
  /// Get the height of the supercell
  virtual const double *GetCellHeights();
  
  /// Get the inverse supercell
  virtual const Vec *GetInverseCell();
  
  /// Get a set of all elements present in the simulations.
  virtual void GetListOfElements(set<int> &elements) const;

  /// Get the masses of the atoms
  virtual const double *GetMasses();

  /// Set a python array on the atoms.
  virtual void SetData(const char *name, PyObject *data);

  /// Get a python array on the atoms.  Returns a new reference.
  virtual PyObject *GetData(const char *name);

  /// Remove a python array from the atoms
  virtual void DeleteData(const char *name);
  
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
  virtual bool UpdateBeforeCalculation(bool flag, double range);
  
  /// Communicate some kind of data to the ghost atoms.
  ///
  /// This is used if a potential needs to communicate some kind of
  /// date during the energy/force calculations.  If the data should
  /// be communicated at the beginning of the calculation, it should
  /// instead be placed in a ghost array.
  ///
  /// NB: In a serial simulation, do nothing except checking that no
  /// ghosts are present.
  virtual void CommunicateData(double *address, int n = 1)
  {assert(nGhosts == 0);}
  
  /// Counter showing changes to the number of atoms or boundary conditions
  virtual int GetAtomsCounter() const {ASSERT(active) return count_atoms;}

  /// Counter showing changes to the unit cell
  ///
  /// Also updated if AtomsCounter has changed.
  virtual int GetCellCounter() const {ASSERT(active) return count_cell;}

  /// Counter showing changes to the atomic numbers.
  ///
  /// Also updated if AtomsCounter has changed.
  virtual int GetNumbersCounter() const {ASSERT(active) return count_numbers;}
  
  /// Counter showing changes to the positions.
  ///
  /// Also updated if CellCounter, NumbersCounter or AtomsCounter is updated.
  virtual int GetPositionsCounter() const {ASSERT(active) return count_positions;}

  /// Counter showing changes to the momenta.
  ///
  /// Reading this counter will cause the momenta to be extracted from Python.
  virtual int GetMomentaCounter();

  /// Print memory usage
  virtual long PrintMemory() const;
  
protected:
  /// Start accessing the Python atoms.
  ///
  /// Until End() is called, references are kept to the Python atoms
  /// and some of their attributes.
  /// The actual work of opening the atoms is done here.
  virtual void DoBegin(PyObject *pyatoms);

  /// Calculate the inverse supercell.
  virtual void invert_cell();

  /// Inform the object about new boundary conditions
  virtual void NewBoundaryConditions(const bool newperiodic[3]) {};

  virtual void check_boundary_conditions(PyArrayObject *py_pbc,
					 bool &step_count_atoms,
					 bool &changed_boundary_conditions);

  virtual bool check_unit_cell(PyArrayObject *py_cell);

  virtual bool check_numbers(PyArrayObject *py_numbers, PyArrayObject *py_gh_num,
			     bool step_count_atoms);

  virtual bool check_positions(PyArrayObject *py_positions, PyArrayObject *py_gh_pos,
			       bool step_count_atoms_or_cell);

protected:
  bool firsttime;         ///< True the first time Begin is called.
  bool hasGhosts;         ///< We have ghost atoms.
  
  int active;             ///< 1 (or larger) if currently processing atoms.
  int expect_reopen;      ///< Allow new calls to Begin while active==expect_reopen.
  int nAtoms;             ///< The number of atoms.
  int nGhosts;            ///< The number of ghost atoms.
  vector<Vec> positions;  ///< A copy of the positions of the atoms.
  vector<Vec> momenta;    ///< A copy of the momenta.
  vector<asap_z_int> numbers;    ///< A copy of the atomic numbers.
  
  Vec cell[3];            ///< The unit cell
  bool pbc[3];            ///< Boundary conditions
  bool all_free_boundaries;      ///< All pbc[] are false
  bool all_periodic_boundaries;  ///< All pbc[] are true
  Vec inverse[3];         ///< Inverse unit cell
  double heights[3];      ///< Heights of the unit cell

  PyObject *py_atoms;     ///< Stored reference to the atoms object.
  PyObject *py_arrays;    ///< Stored reference to arrays dictionary.
  PyArrayObject *py_masses;    ///< NumPy array containing the masses
  PyObject *getmasses_pyname;  ///< Python string "get_masses" (optimization).

  bool momenta_ready;   ///< The momenta have been read and checked.
  bool momenta_nonzero; ///< The momenta may be non-zero.

  // The counters.  Remember to update ParallelAtoms::UpdateBeforeCalculation
  // when adding new counters.
  int count_atoms;     ///< Counts rare, significant changes to atoms.
  int count_cell;      ///< Counts changes to the supercell.
  int count_positions; ///< Counts changes to the positions
  int count_numbers;   ///< Counts changes to the atomic numbers.
  int count_momenta;   ///< Counts changes to the momenta.
  int count_inverse_cell; ///< When was the inverse cell last updated?
};

} // end namespace

#endif // NORMALATOMS_H
