// -*- C++ -*-
// ImagePotential.h: Support for going beyond the minimum image convention.
//
// Copyright (C) 2001-2011 Jakob Schiotz and Center for Individual
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


#include "AsapPython.h"
#include "Potential.h"
#include <assert.h>
#include <vector>
#include <utility>
using std::vector;
using std::pair;

namespace ASAPSPACE {

class Atoms;
class ParallelAtoms;
class Communicator;
class ImageAtoms;

class ImagePotential : public Potential
{
public:
  /// Create a parallel wrapper for the serial potential p
  ImagePotential(PyObject *self, Potential *p);

  virtual ~ImagePotential();

  virtual string GetName() const {return "ImagePotential";}

  /// Set the atoms belonging to this potential.
  ///
  /// This is called automatically by Atoms.SetCalculator() and should
  /// not be called elsewhere.
  virtual void SetAtoms(PyObject *a, Atoms* accessobj = NULL);

  /// Calculate the forces on all atoms and return the result.
  virtual const vector<Vec> &GetForces(PyObject *a);

  /// Calculate the stress on all atoms.
  virtual const vector<SymTensor> &GetVirials(PyObject *a);


  /// THE FOLLOWING METHODS JUST FORWARDS TO THE REAL POTENTIAL

  /// Print memory usage
  virtual long PrintMemory() const {return potential->PrintMemory();}

  virtual double GetPotentialEnergy(PyObject *a) {return potential->GetPotentialEnergy(a);}

  /// Calculate the energy of all atoms.
  virtual const vector<double> &GetPotentialEnergies(PyObject *a) {return potential->GetPotentialEnergies(a);}

  virtual bool CalcReq_Energy(PyObject *pyatoms) {return potential->CalcReq_Energy(pyatoms);}

  virtual bool CalcReq_Forces(PyObject *pyatoms) {return potential->CalcReq_Forces(pyatoms);}


  /// Return the neighbor list.

  /// Return the Python object containing neighbor list for the
  /// potential, if this type of potential supports it, and if it is
  /// defined now.  Otherwise, return NULL without setting a Python
  /// error.
  virtual PyObject *GetNeighborList() const {return potential->GetNeighborList();}

  /// Return the cutoff radius used in the potential.
  virtual double GetCutoffRadius() const {return potential->GetCutoffRadius();}

  /// Return the lattice constant of the material, if well-defined.

  /// If a lattice constant of the material can be defined, return it
  /// in Angstrom, otherwise throw an exception.
  virtual double GetLatticeConstant() const {return potential->GetLatticeConstant();}

  virtual void GetAtomicVolumes(vector<double> &volumes) {potential->GetAtomicVolumes(volumes);}

protected:
  template <class T>
  void CollectFromImages(vector<T> &data);

private:
  Potential *potential;      ///< The wrapped potential (C++ object)
  ImageAtoms *img_atoms;     ///< Access to atoms with images.
  vector<Vec> forces;        ///< Forces after gathering contributions from ghosts.
  vector<SymTensor> virials; ///< Stresses after gathering contributions from ghosts.
  int virials_collect_cnt;   ///< Counter when virials were collected.
  int force_collect_cnt;     ///< Counter when forces were collected.
};

} // end namespace


