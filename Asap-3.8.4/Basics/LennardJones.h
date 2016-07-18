// -*- C++ -*-
// LennardJones.h: The Lennard-Jones potential.
//
// Copyright (C) 2007-2011 David Landis, Jakob Schiotz and Center for
// Individual Nanoparticle Functionality, Department of Physics,
// Technical University of Denmark.  Email: schiotz@fysik.dtu.dk
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

#ifndef LENNARDJONES_H
#define LENNARDJONES_H

#include "AsapPython.h"
#include "Asap.h"
#include "Potential.h"
#include "Vec.h"
#include "SymTensor.h"
#include <vector>

using std::vector;

namespace ASAPSPACE {

class Atoms;
class NeighborList;

/**
 * Implements the Lennard Jones potential 
 *
 * Optimizations that remain to do:
 *   - currently the user specified whether a sparse matrix
 *     or a map for atomic numbers and sigma/epsilon tables
 *     should be used when calling the constructor. The argument
 *     that defines this behaviour is <code>map</code>.<br>
 *     In case <code>map=true</code> turned out to be much better
 *     on every architecture and setup map should be transformed 
 *     into a compiler directive (#define) to gain more performance.
 *
 * <b>Warning:</b> It is not expected in the threaded version that
 *          that a user calculates forces, stresses or 
 *          energies concurrently. I think it might be
 *          possible though.
 *
 */

class LennardJones : public Potential
{
 public:
  /** \brief Constructs a new Lennard Jones Potential instance.
    
     \param numElements the size of the *elements* array (mustn't be zero) 
     \param elements specifies all atom numbers you're using.
     \param epsilon lower triangular matrix giving the epsilons in the order of *elements*.
     \param sigma lower triangular matrix giving the epsilons in the order of *elements*.
     \param masses an array of the masses of the used atoms in the order of *elements*.
     \param modified specifies whether V0 should be added (true: added, false, not added)

     <b>Example:</b>
           
     Warning!!! The numbers in theses examples are random. They don't
     correspond to the real values
     <pre>
      elements = [29, 37, 62]

                   29    37    62
                   ¦     ¦     ¦
                   v     v     v
      epsilon  = [0.15, 0.00, 0.00,  <-- 29
                  0.20, 0.19, 0.00,  <-- 37
                  0.70, 0.25, 0.05]  <-- 62

      explanation of the epsilon matrix interpretation
        epsilon(29, 37)=epsilon(37, 29) is 0.19
        epsilon(37, 62)=epsilon(62, 37) is 0.25
      The values in the upper triangle of the matrice are ignored!!


                   29    37    62
                   ¦     ¦     ¦
                   v     v     v
      sigma    = [2.70, 0.00, 0.00, <-- 29
                  1.50, 3.70, 0.00, <-- 37
                  1.70, 1.25, 4.70] <-- 62

                   29    37    62
                   ¦     ¦     ¦
                   v     v     v
      masses   = [63,    78,   33]
    </pre>
  */
  LennardJones(PyObject *self,
               int numElements, const std::vector<int> &elements,
	       const std::vector<double> &epsilon,
	       const std::vector<double> &sigma,
	       const std::vector<double> &masses, double rCut=-1.0,
	       bool modified=true);
  ~LennardJones();

  virtual string GetName() const {return "LennardJones";}

  /// Print memory usage
  virtual long PrintMemory() const;

  void SetAtoms(PyObject *atoms, Atoms* accessobj = NULL); 
  double GetPotentialEnergy(PyObject *atoms);
  const vector<Vec> &GetForces(PyObject *atoms);
  const vector<SymTensor> &GetVirials(PyObject *atoms);
  const vector<double> &GetPotentialEnergies(PyObject *atoms);
  void CheckNeighborLists();
  double GetCutoffRadius() const {return rCut;}
  double GetLatticeConstant() const {return latticeConstant;}
  /// Return the neighbor list.

  /// Return the Python object containing neighbor list for the
  /// potential, if this type of potential supports it, and if it is
  /// defined now.  Otherwise, return NULL without setting a Python
  /// error.
  PyObject *GetNeighborList() const {return neighborList_obj;}

  /// This potential can be used in parallel simulations
  virtual bool Parallelizable() const {return true;}

private:
  /** \brief sorts the values in items and adds them up*/
  double Add(vector<double> &items);
  /** \brief return abs(d1)<abs(d2) */
  static bool LessFabs(const double d1, const double d2);
  /** \brief Allocates the memory needed by the data structures */
  void Allocate();
  /** \brief calculates the total energy and individual energies of the atoms */
  double CalculateEnergyAndEnergies();
  /** \brief helper function for CalculateEnergyAndEnergies*/
  void CalculateEnergyAndEnergies(vector<double>& atomicEnergies);

  /** \brief helper function for GetCartesianForces*/
  /** \param start the start atom id
      \param end   the end atom id (exclusive)*/
  void GetCartesianForces(vector<Vec>& forces);
  /** \brief  private member function to calculate the stress for a range for atoms with id <code>start</code> to <code>end</code>*/
  /** \param start the start atom id
      \param end   the end atom id (exclusive)*/
  void GetVirials(const Vec *momenta, SymTensor *stresses);
  void GetAtomicVolumes(vector<double> &volumes);

  /** \brief calculates a r-cut from the sigmas */
  double CalculateRCut(int numElements, const vector<double> &sigma); 
  /** \brief converts the parameters to the internal form */
  /** Fills the local epsilon and sigma matrices (stored as an array)

      Depending on the paramter *map* either sparse matrices are 
      created or a lookup array to map atom numbers to ids.
      <pre>
      The given arrays are interpreted as follows:
 
                           n (j as counter)
                          <->
      epsilon  =  0.15 | 0.00 | 0.00 |  29   ^
                  0.20 | 0.19 | 0.00 |  37   |n   (i as counter)
                  0.70 | 0.25 | 0.05 |  62   v

                  So the position of an element at [i][j] is
                  n * i + j
      </pre>
  */
  void Internalize(int numElements, const std::vector<int> &elements,
      const std::vector<double> &epsilon, const std::vector<double> &sigma,
      const std::vector<double> &masses);
  /** \brief prints internal information */
  void Print();
  // Atoms  *atoms;  // Defined in Potential.h
  /** \brief reference to the neighborlist */
  NeighborList *neighborList;
  PyObject *neighborList_obj;
  /** \brief the number of real atoms */
  int    nAtoms;   
  int nSize;
  /** \brief V0 table */
  vector<double> v0;
  /** \brief epsilons  */
  vector<double> epsilon;
  // Atomic volumes
  vector<double> atomicvolume;
  //sigmas (by the 6th power)  */
  vector<double> sigma6;  
  //sigmas (by the 12th power)  */
  vector<double> sigma12;  
  /** \brief cut of radius  */
  double rCut;
  double driftfactor;
  /** \brief shoud V0 be added? */
  bool   modified;
  /** \brief the element order in epsilon and sigma6, if map==true */
  vector<int> elements; 
  /** \brief the masses of the elements; masses[29] returns the mass of Cu if it was specified in the constructor */
  vector<double> masses;

  /** \brief the # of distinct elements used in this simulation */ 
  int    numElements; 
  /** \brief the lattice constant */ 
  double latticeConstant;
  /** \brief atomic energies for atoms*/
  vector<double> atomicEnergies;
  /** \brief the function GetStresses returns pointer into this one.*/
  vector<SymTensor> virials;
  /** \brief GetCartesianForces returns a pointer into this one.*/
  vector<Vec> forces;

  /** \brief A structure of counters to check if recalculations are necessary.*/
  struct {
    int ids;
    int nblist;
    int energies;
    int forces;
    int stresses;
  } counters;
};

} // end namespace


#endif //LENNARDJONES
