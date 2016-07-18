// -*- C++ -*-
// Potential.h: Abstract base class for all potentials.
//
// Copyright (C) 2001-2012 Jakob Schiotz and Center for Individual
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


#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "AsapPython.h"
#include <vector>

#include "Exception.h"
#include "AsapObject.h"
#include "SymTensor.h"

using std::vector;

namespace ASAPSPACE {

class Atoms;
class Vec;
class SuperCell;
class NeighborLocator;

class Potential; // Defined later in this file.

/// The Python object corresponding to a Potential object.
typedef struct {
  PyObject_HEAD
  Potential *cobj;      // Access through this object, either same as orig_cobj or a wrapper.
  Potential *orig_cobj; // The actual potential object
  PyObject *weakrefs;
  bool setatoms_called; // True if SetAtoms has been called.
} PyAsap_PotentialObject;

bool PyAsap_PotentialCheck(PyObject *obj);  // Implemented in PotentialInterface.c

//typedef double symTensor[6];

/// Abstract base class of all potentials.

/// A Potential calculates forces, energies and stresses for a list of
/// atoms.  A given instance of a Potential is associated with a
/// specific Atoms object on a one-to-one bases, this is established
/// when Atoms.SetCalculator() is called.
///
/// Four types of classes are derived from Potential.
///   - Concrete implementations of potentials, such as EMT
///     and MoPotential, LJPotential, BrennerPotential.
///   - Potentials wrapping concrete implementations, but providing
///     special functionalily, such as ParallelPotential (implementing
///     parallel simulations) or QCPotential (implementing the
///     QuasiContinuum method).
class Potential : public AsapObject
{
public:
  Potential(PyObject *self) {
    stresses_cnt = stresses_mom_cnt = -1;
    atoms=NULL;
    this->self = self;
  }

  virtual ~Potential() {}
  
  /// Set the atoms belonging to this potential.

  /// This is called automatically by Atoms.SetCalculator() and should
  /// not be called elsewhere.
  ///
  /// Any wrapping Potential (such as ParallelPotential) that needs to
  /// chain to another Potential's SetAtoms method MUST instead call
  /// Potential::SetAtoms_ThroughPython.
  virtual void SetAtoms(PyObject *a, Atoms* accessobj = NULL) = 0;

  virtual void SetAtoms_ThroughPython(PyObject *a, Atoms* accessobj = NULL);

  /// Calculate the total energy of the system.
  virtual double GetPotentialEnergy(PyObject *a) = 0;

  /// Calculate the forces on all atoms and return the result.
  virtual const vector<Vec> &GetForces(PyObject *a) = 0;

  /// Calculate the "virials" of the atoms.
  ///
  /// The virial divided by the atomic volume is the stress of an atom.
  virtual const vector<SymTensor> &GetVirials(PyObject *a) = 0;

  /// Calculate the total "virial" of the system.
  ///
  /// The virial gives the stress of the sytem when divided by the system volume.
  virtual SymTensor GetVirial(PyObject *a);

  /// Calculate the stress on all atoms.
  virtual const vector<SymTensor> &GetStresses(PyObject *a);

  /// Calculate the total stress of the system.
  virtual SymTensor GetStress(PyObject *a);

  /// Calculate the energy of all atoms.
  virtual const vector<double> &GetPotentialEnergies(PyObject *a) = 0;

  // The following three functions are used to check if recalculations
  // are needed.  If a potential does not contain any logic to prevent
  // recalculations, these functions should return True.

  /// Is work required to calculate the energy?
  virtual bool CalcReq_Energy(PyObject *pyatoms) {return true;}

  /// Is work required to calculate the forces?
  virtual bool CalcReq_Forces(PyObject *pyatoms) {return true;}

  /// Is work required to calculate the stress?
  virtual bool CalcReq_Virials(PyObject *pyatoms) {return true;}

  /// Is work required to calculate the stress?
  virtual bool CalcReq_Stress(PyObject *pyatoms);

  // Check if Neighbor lists are up to date.
  //virtual void CheckNeighborLists() = 0;

  /// Return the neighbor list.

  /// Return a BORROWED reference to the Python object containing
  /// neighbor list for the
  /// potential, if this type of potential supports it, and if it is
  /// defined now.  Otherwise, return NULL without setting a Python
  /// error.
  virtual PyObject *GetNeighborList() const {return NULL;}

  /// Return the cutoff radius used in the potential.
  virtual double GetCutoffRadius() const = 0;

  /// Return the lattice constant of the material, if well-defined.

  /// If a lattice constant of the material can be defined, return it
  /// in Angstrom, otherwise throw an exception.
  virtual double GetLatticeConstant() const = 0;

  /// Can this potential be used in parallel simulations?

  ///
  /// ParallelPotential checks this and protests if false.  Note that
  /// ParallelPotential itself returns false if this method is called,
  /// as it cannot be passed to another instance of ParallelPotential.
  virtual bool Parallelizable() const {return false;}

  /// Print memory usage
  virtual long PrintMemory() const = 0;

  /// Clean up after an exception.  Called by the interface code.
  void RecoverAfterException();

  /// Return the atomic volumes that should be used to calculate the stresses
  ///
  /// The atoms should already be open, so no atoms parameter is passed.
  /// If an empty volumes vector is returned, GetStress will divide the unit
  /// cell evenly amongst the atoms.
  virtual void GetAtomicVolumes(vector<double> &volumes) {volumes.clear();}

protected:
  Atoms *atoms;
  PyObject *self;              ///< This objects Python counterpart

private:  // These *really* need to be private, not protected!
  int stresses_cnt;
  int stresses_mom_cnt;
  vector<SymTensor> stresses;
  SymTensor stress;
};

} // end namespace

#endif // POTENTIAL_H
