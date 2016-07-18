// -*- C++ -*-
// ParalelPotential.h: Support for MPI communication in potentials.
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


class ParallelPotential : public Potential
{
public:
  /// Create a parallel wrapper for the serial potential p
  ParallelPotential(PyObject *self, PyObject *p);

  virtual ~ParallelPotential();

  virtual string GetName() const {return "ParallelPotential";}

  /// Print memory usage
  virtual long PrintMemory() const;

  /// Set the atoms belonging to this potential.

  /// This is called automatically by Atoms.SetCalculator() and should
  /// not be called elsewhere.
  virtual void SetAtoms(PyObject *a, Atoms* accessobj = NULL);

  virtual double GetPotentialEnergy(PyObject *a);

  /// Calculate the forces on all atoms and return the result.
  virtual const vector<Vec> &GetForces(PyObject *a);

  /// Calculate the stress on all atoms.
  virtual const vector<SymTensor> &GetVirials(PyObject *a);

  /// Calculate the total stress of the system.
  virtual SymTensor GetStress(PyObject *a);
  virtual SymTensor GetVirial(PyObject *a);

  /// Calculate the energy of all atoms.
  virtual const vector<double> &GetPotentialEnergies(PyObject *a);

  virtual bool CalcReq_Energy(PyObject *pyatoms);
  virtual bool CalcReq_Forces(PyObject *pyatoms);
  virtual bool CalcReq_Stress(PyObject *pyatoms);

  /// Return the neighbor list.

  /// Return the Python object containing neighbor list for the
  /// potential, if this type of potential supports it, and if it is
  /// defined now.  Otherwise, return NULL without setting a Python
  /// error.
  virtual PyObject *GetNeighborList() const;

  /// Return the cutoff radius used in the potential.
  virtual double GetCutoffRadius() const;

  /// Return the lattice constant of the material, if well-defined.

  /// If a lattice constant of the material can be defined, return it
  /// in Angstrom, otherwise throw an exception.
  virtual double GetLatticeConstant() const;

  virtual void GetAtomicVolumes(vector<double> &volumes);

private:
  PyObject *py_potential;   ///< The wrapped potential (Python object)
  Potential *potential;     ///< The wrapped potential (C++ object)
  ParallelAtoms *par_atoms; ///< Parallel access to atoms.
  Communicator *mpi;        ///< The communicator object.
  vector<Vec> forces;       ///< Forces after gathering contributions from ghosts.
  vector<SymTensor> stresses; ///< Stresses after gathering contributions from ghosts.
  int stress_collect_cnt;   ///< Counter when stresses were collected.
  int force_collect_cnt;    ///< Counter when forces were collected.
};

} // end namespace


