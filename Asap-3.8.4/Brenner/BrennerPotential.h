// -*- C++ -*-
// BrennerPotential.h - Asap implementation of the Brenner Potential.
//
// Copyright (C) 2010-2011 Jakob Schiotz and Center for Individual
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

#ifndef BRENNER_POTENTIAL_H
#define BRENNER_POTENTIAL_H

#include "Potential.h"
#include "Atoms.h"
#include "NeighborList.h"

#define BRENNERNTAB 10000

namespace ASAPSPACE {

struct AtomPairInfoState;

class BrennerPotential : public Potential
{
public:
  BrennerPotential(PyObject *self);
  virtual ~BrennerPotential();

  // Initialize global data before using this kind of object!
  static void Initialize();  

  virtual string GetName() const {return "BrennerPotential";}

  virtual void SetAtoms(PyObject *a, Atoms* accessobj = NULL);

  /// Calculate the total energy of the system.
  virtual double GetPotentialEnergy(PyObject *a);

  /// Calculate the forces on all atoms and return the result.
  virtual const vector<Vec> &GetForces(PyObject *a);

  /// Calculate the stress on all atoms.
  virtual const vector<SymTensor> &GetVirials(PyObject *a)
  {
    throw AsapNotImplementedError("BrennerPotential does not support stress calculations.");
  }

  /// Calculate the energy of all atoms.
  virtual const vector<double> &GetPotentialEnergies(PyObject *a)
  {
    throw AsapNotImplementedError("BrennerPotential does not support per-atom energies.");
  }

  // Check if Neighbor lists are up to date.
  //virtual void CheckNeighborLists();

  /// Return the neighbor list.

  /// Return the Python object containing neighbor list for the
  /// potential, if this type of potential supports it, and if it is
  /// defined now.  Otherwise, return NULL without setting a Python
  /// error.
  virtual PyObject *GetNeighborList() const {return nblist_obj;}

  /// Return the cutoff radius used in the potential.
  virtual double GetCutoffRadius() const {return rmax_nosq;}

  /// Return the lattice constant of the material, if well-defined.

  /// If a lattice constant of the material can be defined, return it
  /// in Angstrom, otherwise throw an exception.
  virtual double GetLatticeConstant() const
  { throw AsapNotImplementedError("BrennerPotential::GetLatticeConstant not implemented.");}

  /// Can this potential be used in parallel simulations?
  virtual bool Parallelizable() const {return false;}

  /// Print memory usage
  virtual long PrintMemory() const {return 0;}  // XXXX

private:
  /// Calculate energies and forces.
  void Calculate(PyObject *a);

  void CheckAndUpdateNeighborList();

  void CountAtoms();

  double caguts();

  double pibond();

  void calc1side(const int k_this,
		 const int kother,
		 const int index1,
		 const int index2,
		 const int j,
		 const int c_sign,
		 double *xni,
		 const struct AtomPairInfo *const atom_pairj,
		 const struct AtomPairInfo *const atom_pairk,
		 const Vec cj,
		 double *xsij, double *xhc[2], double *cfuni, double *dcfuni,
		 double* conk, double* dctjk, double* dctij, double* dctik,
		 double* ssumk, double* sdalik, double* xsjk,
		 double* xsik,
		 const double sij,
		 const double rsqij,
		 Vec *xk /* xl */, double *cosk, double *sink,
		 double *exnij,
		 double *dexni);

  void calc_dihedral_terms(int i, int j, int jn,
			   const struct AtomPairInfoState *apis,
			   double xnt1,
			   double xnt2,
			   double conjug, double sij,
			   double *sink, double *sinl,
			   double *dctik, double *dctjk,
			   double *dctil, double *dctij,
			   double *dctjl, double *dctji,
			   double *cosk, double *cosl, Vec cj,
			   double vatt,
			   const Vec *xk, const Vec *xl,
			   double *btot,
			   double *dradi,
			   double *dradj,
			   double *drdc);

  void calc_many_body(const struct AtomPairInfoState *apis,
		      int i_side, int index1, int index2,
		      int indexj,
		      double vdbdi,
		      double *xsik,
		      double *dexni,
		      double vdrdi, double vdrdc, double *cfuni, double sdalik,
		      double *xsjk, Vec *xk, double *dcfuni, double conk);
  
  double BCUINT(int KJ, double XX1, double XX2,
		double *ansy1_ptr, double *ansy2_ptr);

  double TOR(double XNT1, double XNT2, double CONJUG,
	     double *drdl_ptr, double *drdm_ptr, double *drdn_ptr);

  double RADIC(int KI, int KJ, double XNT1, double XNT2, double CONJUG,
	       double *drdl_ptr, double *drdm_ptr, double *drdn_ptr);

  double sili_germ();
  
  inline int getKtype(int n) {return z_to_ktype[z[n]];}

  inline void transForce(int i, int j, const Vec &f)
  {
    force[i] += f;
    force[j] -= f;
  }   
  
  // Get number of atoms of a given ktype.
  int getNoa(int ktype) const {return noa[ktype];}

  // Initialize data used in caguts.cpp
  static void init_c();

  // Initialize data used in pibond.cpp
  static void init_xh();

  static void init_in2();

  static void init_in3();

  // Initialize data used in sili_germ.cpp
  static void si_ge_init();

  // Generate lookup tables.
  static void mtable(double drtable[4][4][BRENNERNTAB], double ddtab[4][4],
		     double rtable[4][4][BRENNERNTAB],
		     double atable[4][4][BRENNERNTAB],
		     double datable[4][4][BRENNERNTAB], 
		     double tabfc[4][4][BRENNERNTAB],
		     double tabdfc[4][4][BRENNERNTAB]);


private:
  Atoms *atoms;
  AtomPairInfoState *apis;
  NeighborList *nblist;
  PyObject *nblist_obj;

  // Info about the atoms
  const Vec *positions;
  const asap_z_int *z;  // Atomic numbers
  int nAtoms;
  vector<Vec> force;
  double Epot;
  int noa[5];  // Number of atoms

  // Counters to detect changes in atoms
  int counter;
  int counter_z;
  
  // Static data initialized by Initialize()
  static const int MAXATNO = 93;     // Plutonium not supported :-)
  static int ktype_to_z[5];
  static int z_to_ktype[MAXATNO+1];
  static double rb2[4][4];     // Initialized in mtable.
  static double rmax[4][4];    // Initialized in Initialize() using rb2.
  static double rmax_nosq;     // Max value of sqrt(rmax). CUTOFF!
};

} // end namespace

#endif // BRENNER_POTENTIAL_H
