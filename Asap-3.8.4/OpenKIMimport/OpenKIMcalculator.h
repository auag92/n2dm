// OpenKIMcalculator.h - interface to OpenKIM models.
//
// This class is part of the optional Asap module to support OpenKIM
// models.  The class OpenKIMcalculator does the actual interfacing
// to the model.

// Copyright (C) 2014 Jakob Schiotz and Center for Individual
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

#ifndef _OPENKIMCALCULATOR_H
#define _OPENKIMCALCULATOR_H

#include "AsapPython.h"
#include "Asap.h"
#include "Potential.h"
#include <map>
#include <vector>

class KIM_API_model;

namespace ASAPSPACE {

class OpenKIMcalculator : public Potential
{
public:
  OpenKIMcalculator(PyObject *self);
  ~OpenKIMcalculator();

  void Initialize(const char *descr, const char *name);
  /// Indicate if a quantity should be allocated in the KIM API.
  void PleaseAllocate(string quantity, bool alloc);

  /// This potential can be used in parallel simulations
  virtual bool Parallelizable() const {return true;}

  /// Check if neighbor lists need an update.
  ///
  /// Return true if the neigbor list needs an update.  In parallel
  /// simulations, the atoms will cause a communication so it returns
  /// true if an update is needed on any processor.  In parallel
  /// simulations, it also instructs the atoms to communicate their
  /// update status.

  virtual bool CheckNeighborList();
  /// Update the neighbor lists
  virtual void UpdateNeighborList();

  /// Get information about supported elements.
  void GetSupportedSymbols(std::vector<const char *> &symbols);

  void ClearTranslation() {z_to_typecode.clear();}
  void AddTranslation(int z, int typecode) {z_to_typecode[z] = typecode;}

  /// Get the particle type code (number) from the particle symbol
  int GetParticleTypeCode(const char *symbol);

  virtual void SetAtoms(PyObject *pyatoms, Atoms* accessobj = NULL);

  virtual const vector<Vec> &GetForces(PyObject *pyatoms);
  virtual const vector<double> &GetPotentialEnergies(PyObject *pyatoms);
  virtual double GetPotentialEnergy(PyObject *pyatoms);
  virtual const vector<SymTensor> &GetVirials(PyObject *a);
  virtual SymTensor GetVirial(PyObject *a);

  /// Is work required to calculate the energy?
  virtual bool CalcReq_Energy(PyObject *pyatoms);

  /// Is work required to calculate the forces?
  virtual bool CalcReq_Forces(PyObject *pyatoms);

  /// Is work required to calculate the stress?
  virtual bool CalcReq_Virials(PyObject *pyatoms);

  virtual std::string GetName() const {return "OpenKIMcalculator";}

  const char *GetNBCmethod();

  /// Return the cutoff radius used in the potential.
  virtual double GetCutoffRadius() const {return cutoff;}

  /// Return the lattice constant of the material, if well-defined.

  /// If a lattice constant of the material can be defined, return it
  /// in Angstrom, otherwise throw an exception.
  virtual double GetLatticeConstant() const
  {throw AsapError("OpenKIMcalculator::GetLatticeConstant not supported.");}

  /// Print memory usage
  virtual long PrintMemory() const
  {throw AsapError("OpenKIMcalculator::PrintMemory not supported.");}

  // The OpenKIM neighbor list interface function passed to the model.
  int get_neigh(int *mode, int* request,
        int *particle, int *numnei, int **nei1particle,
        double **rij);

protected:
  typedef enum {
    nbl_cluster, nbl_pure_h, nbl_pure_f, nbl_miopbc_h, nbl_miopbc_f, nbl_rvec_h, nbl_rvec_f
    } nblisttype_t;

  /// Allocate the neighbor list.
  virtual void CreateNeighborList();

  /// (Re)allocate storage for forces, energies and intermediate results.
  virtual void Allocate();

  /// Calculate stuff
  virtual void Calculate(PyObject *pyatoms, bool calc_energy, bool calc_force);
  virtual void DoCalculate();

  /// Set orthocell from atoms cell, and check orthogonality of the latter.
  void SetOrthoCell();

private:
  KIM_API_model *model;
  bool model_initialized;
  NeighborLocator *nblist;     ///< The neighborlist object.
  PyObject *nblist_obj;
  double driftfactor;             ///< Drift factor for the neighbor list.

  int nAtoms;      ///< Number of particles without ghost atoms
  int nSize;       ///< Number of particles with ghost atoms
  int nAtomsAlloc, nSizeAlloc;  ///< Values of nAtoms and nSize at last allocation.
  bool ghostatoms;          ///< True if atoms have ghosts.
  int nSpecies;
  double cutoff;      // The cutoff, will be set by OpenKIM model.
  vector<int> species;
  vector<Vec> forces;
  vector<double> particleEnergy;
  vector<SymTensor> particleVirial;
  double energy;
  SymTensor virial;

  nblisttype_t nblisttype;
  bool nblist_half;
#if 0
  int nb_accessmode;
#endif
  vector<int> nb_buffer_n;
  vector<Vec> nb_buffer_rij;
  vector<double> nb_buffer_dist;
  double orthocell[3];       // Cell size for MI_OPBC boundary conditions.
  int nblist_iterator;       // Next neighbor (used in iterator mode).

  std::map<asap_z_int,int> z_to_typecode;  // Translation of atomic number to KIM typecode.

  struct {
    bool energy;
    bool particleEnergy;
    bool forces;
    bool virial;
    bool particleVirial;
  } pls_alloc;
  int pls_alloc_n;  // Counts that they have all been set.

  /// A structure of counters to check if recalculations are necessary.
  struct {
    int nblist;
    int compute_energy;
    int compute_force;
  } counters;

  /// A structure of bools, to remember if recalculations are necessary.
  struct {
    int nblist;
    int compute_energy;
    int compute_force;
  } recalc;

};

} // namespace



#endif // _OPENKIMCALCULATOR_H
