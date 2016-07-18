// -*- C++ -*-
// ImageAtoms.h:  Beyond the Minimum Image Convention
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


#ifndef IMAGEATOMS_H
#define IMAGEATOMS_H

#include "AsapPython.h"
#include "Asap.h"
#include "Atoms.h"

namespace ASAPSPACE {

/// ImageAtoms: Handle boundary conditions beyond the minimum-image convention.
///
/// The ImageAtoms object is a decorator for the Atoms and ParallelAtoms
/// objects, in the Design Pattern sense.  It presents the atoms to
/// other objects by hiding any periodic boundary conditions and instead
/// adding "image atoms", extra atoms representing the images of the real
/// atoms through the periodic boundary conditions.  This is not limited
/// to the minimum-image convention, but can handle arbitrarily small
/// unit cells.
///
/// ImageAtoms has a sister class, ImagePotential, that reintegrates the
/// forces on the image atoms into the real atoms.
///
/// The vast majority of the methods just forward to the real Atoms or
/// ParallelAtoms object.  This is implemented as a Decorator pattern instead
/// of a subclass because we only know at run time if we will be extending
/// an Atoms or a ParallelAtoms object.  Since this object is intercalated
/// between the real Potential and the ParallelAtoms potential, the
/// communication between ParallelPotential and its ParallelAtoms object
/// does not pass through ImageAtoms.
///
/// This object is refcounted, see AsapAtoms_INCREF() and AsapAtoms_DECREF().
class ImageAtoms : public Atoms
{
protected:
  /// Delete the interface object.
  ///
  /// Protected: may only be called from AsapAtoms_DECREF() friend function.
  virtual ~ImageAtoms();
  
public:
  /// Create the interface object.
  ImageAtoms(Atoms *atoms);

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

  /// Get the number of ghost atoms
  virtual int GetNumberOfGhostAtoms() const {return nGhosts + nImages;}
  
  /// Do we have ghost atoms? True if arrays are present, even if empty.
  virtual bool HasGhostAtoms() const {return (nGhosts + nImages > 0);}

  /// Get the cartesian positions.  Ghost positions, if any, are at the end.
  virtual const Vec *GetPositions() const {return &allpositions[0];}

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
   {assert(allnumbers.size() > 0); return &allnumbers[0];}

  /// Get the boundary conditions.
  ///
  /// These will be faked to free boundary conditions.
  virtual const bool *GetBoundaryConditions() const {return bc_are_free;}

  virtual bool AllFreeBoundaries() const {return true;}

  virtual bool AllPeriodicBoundaries() const {return false;}


  /// ALL THE FOLLOWING METHODS JUST FORWARD TO THE REAL ATOMS OBJECT.

  /// Finish accessing the Python atoms, freeing resources.
  virtual void End() {realatoms->End();}

  /// Are the atoms active (opened by Begin()) ?
  virtual bool IsActive() const {return realatoms->IsActive();}
  
  /// Get the number of atoms
  virtual int GetNumberOfAtoms() const {return realatoms->GetNumberOfAtoms();}

  /// Get total number of atoms in parallel simulation.
  ///
  /// In a serial simulation, same as GetNumberOfAtoms, in a parallel simulation
  /// it is the sum over these over all processors.
  virtual int GetTotalNumberOfAtoms() const {return realatoms->GetTotalNumberOfAtoms();}

  /// Return three integers specifying the CPU layout.
  virtual IVec GetNumberOfCells() const {return realatoms->GetNumberOfCells();}

  /// Get the cartesian momenta
  virtual const Vec *GetMomenta() {return realatoms->GetMomenta();}

  /// Get the supercell
  virtual const Vec *GetCell() const {return realatoms->GetCell();}

  /// Get the supercell.   No sanity check!
  virtual const Vec *GET_CELL() const {return realatoms->GET_CELL();}

  /// Get the volume of the supercell
  virtual double GetVolume() const {return realatoms->GetVolume();}
  
  /// Get the height of the supercell
  virtual const double *GetCellHeights() {return realatoms->GetCellHeights();}
  
  /// Get the inverse supercell
  virtual const Vec *GetInverseCell() {return realatoms->GetInverseCell();}
  
  /// Get a set of all elements present in the simulations.
  virtual void GetListOfElements(set<int> &elements) const {realatoms->GetListOfElements(elements);}

  /// Get the masses of the atoms
  virtual const double *GetMasses() {return realatoms->GetMasses();}

  /// Set a python array on the atoms.
  virtual void SetData(const char *name, PyObject *data) {realatoms->SetData(name, data);}

  /// Get a python array on the atoms.  Returns a new reference.
  virtual PyObject *GetData(const char *name) {return realatoms->GetData(name);}

  /// Remove a python array from the atoms
  virtual void DeleteData(const char *name) {realatoms->DeleteData(name);}
  

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
  {realatoms->CommunicateData(address, n);}
  
  /// Counter showing changes to the number of atoms or boundary conditions
  virtual int GetAtomsCounter() const {return realatoms->GetAtomsCounter();}

  /// Counter showing changes to the unit cell
  ///
  /// Also updated if AtomsCounter has changed.
  virtual int GetCellCounter() const {return realatoms->GetCellCounter();}

  /// Counter showing changes to the atomic numbers.
  ///
  /// Also updated if AtomsCounter has changed.
  virtual int GetNumbersCounter() const {return realatoms->GetNumbersCounter();}
  
  /// Counter showing changes to the positions.
  ///
  /// Also updated if CellCounter, NumbersCounter or AtomsCounter is updated.
  virtual int GetPositionsCounter() const {return realatoms->GetPositionsCounter();}

  /// Counter showing changes to the momenta.
  ///
  /// Reading this counter will cause the momenta to be extracted from Python.
  virtual int GetMomentaCounter() {return realatoms->GetMomentaCounter();}

  /// Print memory usage
  virtual long PrintMemory() const {return realatoms->PrintMemory();}
  
  // Return the map from image atoms to original atoms
  virtual const vector<int> &GetOriginalAtomsMap() const {return original_atom;}

  virtual int GetOriginalNumberOfAtoms() const {return nAtoms + nGhosts;}

private:
  void make_images(double range);

  void update_images();

protected:
  Atoms *realatoms;

  int nAtoms;
  int nGhosts;
  int nImages;
  int nSize;

  double last_range;

  vector<Vec> allpositions;
  vector<asap_z_int> allnumbers;

  static const bool bc_are_free[3];

  // A translation and two indices for each block of images that share the same
  // translation vector.
  vector<IVec> translations;
  vector<int> first_atom;
  vector<int> last_atom;

  // A list of atoms (ghost or real) on which the images are based.
  vector<int> original_atom;

  bool initialized;
};

} // end namespace

#endif // IMAGEATOMS_H
