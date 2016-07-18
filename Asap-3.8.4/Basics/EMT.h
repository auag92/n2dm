// EMT.h  --  Implements the EMT potential.  
// -*- c++ -*-
//
// Copyright (C) 2001-2011 Jakob Schiotz and Center for Individual
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

#ifndef _EMT_H
#define _EMT_H

#include "AsapPython.h"
#include "Asap.h"
#include "Atoms.h"
#include "Vec.h"
#include "SymTensor.h"
#include "Potential.h"
#include "EMTParameterProvider.h"
#include "TinyMatrix.h"
#include <vector>
using std::vector;

namespace ASAPSPACE {

class NeighborLocator;
struct emt_parameters;

/// The Effective Medium Theory (EMT) potential.

/// The Effective Medium Theory (EMT) as originally implemented in Per
/// Stoltze's ARTwork program, and documented in
///   - K. W. Jacobsen, P. Stoltze and J. K. Norskov,
///     Surf. Sci. vol. 366, p. 394-402 (1996).
class EMT : public Potential
{
public:
  /// Create an EMT potential optionally using a parameter provider.
  EMT(PyObject *self, PyObject *prov);

  /// Delete the EMT potential.
  virtual ~EMT();

  virtual string GetName() const {return "EMT";}
  virtual string GetRepresentation() const;

  /// Whether we subtract E0 from atomic energies.
  ///
  /// Defines the zero of potential energy; if we don't subtract, zero
  /// corresponds to infinite seperation. If we do, it is the bulk fcc
  /// crystal at the equilibrium lattice constant.
  void SetSubtractE0(bool subtractE0) { this->subtractE0 = subtractE0; }

  virtual void SetAtoms(PyObject *pyatoms, Atoms* accessobj = NULL);
  virtual const vector<Vec> &GetForces(PyObject *pyatoms);
  virtual const vector<double> &GetPotentialEnergies(PyObject *pyatoms);
  virtual double GetPotentialEnergy(PyObject *pyatoms);
  virtual const vector<SymTensor> &GetVirials(PyObject *pyatoms);
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

  virtual bool CalcReq_Energy(PyObject *pyatoms);
  virtual bool CalcReq_Forces(PyObject *pyatoms);
  virtual bool CalcReq_Virials(PyObject *pyatoms);

  virtual double GetCutoffRadius() const {return rNbCut;}
  virtual double GetLatticeConstant() const;
  virtual int GetNumberOfAtoms() const {return nAtoms;}

  /// Return a pointer to the EMT "density" sigma1.
  const vector<vector<double> > &GetSigma1() {return sigma1;}

  /// Return a pointer to the EMT "density" sigma2.
  const vector<vector<double> > &GetSigma2() {return sigma2;}

  /// Print the EMT parameters.
  void PrintParameters();

  /// Return the neighbor list.

  /// Return a BORROWED reference to the Python object containing
  /// neighbor list for the
  /// potential, if this type of potential supports it, and if it is
  /// defined now.  Otherwise, return NULL without setting a Python
  /// error.
  PyObject *GetNeighborList() const {return nblist_obj;}

  /// This potential can be used in parallel simulations
  virtual bool Parallelizable() const {return true;}

  /// Mark that the atoms have changed the boundary conditions.
  virtual void BoundaryConditionsChanged();

  /// Print memory usage
  virtual long PrintMemory() const;
  
  /// Return the atomic volumes that should be used to calculate the stresses
  virtual void GetAtomicVolumes(vector<double> &volumes);

protected:
  /// Initialization of the EMT parameters.
  virtual void InitParameters();
  /// (Re)allocate storage for forces, energies and intermediate results.
  virtual void Allocate();
  /// (Re)allocate storage for stresses
  virtual void AllocateStress();
  /// Calculate type numbers from the atomic numbers.
  virtual void CalculateIDs();
  
  /// Calculate sigma1 and perhaps sigma2.
  virtual void CalculateSigmas(bool calculatesigma2);

  /// Calculate energies
  virtual void CalculateEnergies();

  /// Calculate energies once sigma1 is known

  /// If Epot is NULL, energies are not calculated, only derivatives
  /// needed for the force.
  virtual void CalculateEnergiesAfterSigmas(bool calcEpot);

  /// Calculate the forces
  virtual void CalculateForces();

  /// Calculate forces in a multicomponent system.
  virtual void CalculateForcesAfterEnergies();
  
  /// Calculate forces in a system with only one element.
  virtual void CalculateForcesAfterEnergiesSingle();
  
  /// Calculate stresses.
  virtual void CalculateVirials();

  /// Allocate the neighbor list.
  virtual void CreateNeighborList();
  
  /// sigma_batch does the hard work in CalculateSigmas().
  virtual void sigma_batch(int *self, int *other, Vec rnb[],
                           double *sq_dist, int zs, int zo, int n,
                           bool calculatesigma2, bool partialupdate = false);

  /// force_batch does the hard work in CalculateForcesAfterEnergy().
  virtual void force_batch(const int *self, const int *other, const Vec rnb[],
                           const double sq_dist[], const double dEdss[],
                           const double dEdso[], int zs, int zo, int n);
  void distribute_force(const int *self, const int *other,
      const double *df, const Vec *rnb, int n);
  
protected:  // Data
  // Atoms *atoms;             ///< The atoms we are working on
  bool ghostatoms;          ///< True if atoms have ghosts.
  int nAtoms;               ///< The number of (real) atoms.
  int nSize;  ///< Number of atoms including ghost atoms.
  NeighborLocator *nblist;     ///< The neighborlist object.
  PyObject *nblist_obj;
  double driftfactor;             ///< Drift factor for the neighbor list.
  EMTParameterProvider *provider; ///< The source of the EMT parameters
  PyObject *provider_obj;         ///< The Python object behind provider.
  bool subtractE0;          ///< Whether we subtract E0 from atomic energies.
  const emt_parameters *singleelement; ///< Element if only one element, 0 otherwise.
  std::vector<const emt_parameters *> parameters; ///< The EMT parameters
  const TinyDoubleMatrix *chi;  ///< The Chi matrix of EMT.
  int nelements;    ///< The number of different elements in the simulation
  /// Cutoff parameters (from EMTParameterProvider).
  double rFermi, rNbCut;  
  /// Cutoff slope
  double cutoffslope;
  bool sigma2isvalid;              ///< For consistency checks.
  bool initialized;                ///< SetAtoms has been called.
  bool always_fullnblist;          ///< Allways use full neighbor lists (beyond minimum-image).

  /// Temporary data for the atoms
  
  /// Each atom has two sigmas for each possible neighboring element.
  /// So these are arrays (nelements long).
  vector<vector<double> > sigma1;
  vector<vector<double> > sigma2;

  //@{
  /// Each atom has a single Ec, Eas, radius etc.
  vector<double> Ec;
  vector<double> Eas;
  vector<double> Epot;
  vector<double> radius;
  vector<double> dEds;
  vector<Vec> force;
  /// The stresses.
  vector<SymTensor> virials;
  //@}

  /// Temporary data
  vector<double> tmp_double;
  vector<double> ex2;  // Data stored between calls to CalculateEnergiesAfterSigma (with and without actual energy calculation)
  
  int nAtomsRes, nSizeRes;
  ///< If there are ghostatoms, some extra space is reserved for the arrays

  /// The atomic numbers are translated into IDs, integers in [0, nelements-1]
  vector<int> id;
    
  /* The buffers for batch processing */
  // static const int BUFLEN;  // The batch--buffer size.  Avoid powers of 2
  // static const int NMAXELEMENTS;
  
  /// A structure of counters to check if recalculations are necessary.
  struct {
    int ids;
    int nblist;
    int sigma1;
    int sigma2;
    int beforeforces;
    int energies;
    int forces;
    int virials;
  } counters;

  /// A structure of bools, to remember if recalculations are necessary.
  struct {
    bool ids;
    bool nblist;
    bool sigma1;
    bool sigma2;
    bool beforeforces;
    bool energies;
    bool forces;
    bool virials;
  } recalc;

  /// A flag indicating that Atoms has been opened in a child class.
  bool skip_begin;
}; 

} // end namespace

#endif // ! _EMT_H
