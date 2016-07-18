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

/*
 * Calculates the Morse potential
 */


// This implementation supports parallel simulations. Due to the
// simplicity of Morse, everything is handled by the ParallelAtoms
// and ParallelPotential classes, this potential only has to be
// prepared for parallel simulations in two ways:
//
// 1. After a neighborlist update, the number of atoms may have
//    changed, and arrays may have to be reallocated.
//
// 2. The neighborlist may return atoms with number higher than
//    nAtoms, these are ghosts of atoms residing on other processor,
//    and only have positions and atomic numbers. It must be avoided
//    to update energies, forces and stresses for these, as that would
//    overwrite the end of the arrays.

#include "AsapPython.h"
#include "Asap.h"
#include "Morse.h"
#include "NeighborList.h"
#include "Atoms.h"
#include "Vec.h"
#include "Timing.h"
//#define ASAPDEBUG
#include "Debug.h"
#include <math.h>
#include <stdio.h>
#include <algorithm>

#ifndef _WIN32
#include <sys/wait.h>
#endif

using std::cerr;
using std::endl;
using std::flush;
using std::less;

extern int verbose;

#define SPARSE_MATRIX_SIZE 92  // No Pu, please!  

// Use same volume for all atoms in stress calculation.
#define SHAREVOLUME   

// Default driftfactor
#define DRIFTFACTOR 0.05


#if 0
#define VERB(x) if (verbose == 1) cerr << x
#else
#define VERB(x)
#endif


// SunOS compiler: sqrt does not work with a literal number (cannot decide if
// it is a double or a long double).
#ifdef SUN_SQRT_LITERAL_HACK
#define sqrt(x) sqrt((double) x)
#endif

// Standard mapping of the six independent parts of the stress tensor to
// vector notation
const static int stresscomp[3][3] = {{0, 5, 4}, {5, 1, 3}, {4, 3, 2}};



//******************************************************************************
//                                    Morse
//******************************************************************************
Morse::Morse(PyObject *self,
             const std::vector<int> &elements,
             const std::vector<double> &epsilon,
             const std::vector<double> &alpha,
             const std::vector<double> &rmin,
             double rCut, bool modified) : Potential(self)
{
  DEBUGPRINT;

  atoms = NULL;
  this->nAtoms = 0;
  neighborList = 0;
  neighborList_obj = 0;
  driftfactor = DRIFTFACTOR;
  this->modified = modified;
  numElements = elements.size();
  this->latticeConstant = 0.0;

  memset(&counters, 0x0, sizeof(counters));

  if(rCut<0)
    this->rCut = CalculateRCut(numElements*numElements, alpha, rmin);
  else
    this->rCut = rCut;

  Internalize(elements, epsilon, alpha, rmin);
  //Print();
  DEBUGPRINT;
}

Morse::~Morse()
{
  DEBUGPRINT;
  Py_XDECREF(neighborList_obj);
  if (atoms != NULL)
    AsapAtoms_DECREF(atoms);
}


//******************************************************************************
//                                   Print
//******************************************************************************
void Morse::Print() {

  std::cout<<"   ***************************************"<<endl;
  std::cout<<"   atoms is        "<<atoms<<endl;
  std::cout<<"   neighborList is "<<neighborList<<endl;
  std::cout<<"   nAtoms is       "<<nAtoms<<endl;
  if(atoms&&neighborList) {
    std::cout<<"   atoms->GetNumberOfAtoms() is       "<<atoms->GetNumberOfAtoms()<<endl;
    std::cout<<"   neighborList->MaxNeighborListLength() "<<neighborList->MaxNeighborListLength()<<endl;
  }
  std::cout<<"   modified is     "<<modified<<endl;
  std::cout<<"   numElements is  "<<numElements<<endl;
  std::cout<<"   latticeConstant is    "<<latticeConstant<<endl;
  std::cout<<"   SPARSE_MATRIX_SIZE is "<<SPARSE_MATRIX_SIZE<<endl;
  std::cout<<endl;
  std::cout<<"   ***************************************"<<endl<<endl;
  //std::cout<<"Epsilon table:"<<endl;
  //for(int i=0;i<SPARSE_MATRIX_SIZE;i++)
  //  for(int j=0;j<SPARSE_MATRIX_SIZE;j++)
  //	std::cout<<"epsilon["<<i<<","<<j<<"]:"<<epsilon[i*SPARSE_MATRIX_SIZE+j]<<endl;
}


//******************************************************************************
//                                Internalize
//******************************************************************************
void Morse::Internalize(const std::vector<int> &p_elements,
                        const std::vector<double> &p_epsilon,
                        const std::vector<double> &p_alpha,
                        const std::vector<double> &p_rmin)
{
  DEBUGPRINT;
  //Print();

  // 1. Copy the masses.  REMOVED: May be added back for stress calculations.

  // 2. Setting the inverse volume.  
  inversevolume.resize(SPARSE_MATRIX_SIZE);
  memset(&inversevolume[0], -1, sizeof(double) * SPARSE_MATRIX_SIZE);
  //double volumefactor = 0.94;  // 1.1^3 / sqrt(2) ## Skal det slettes? ##
  for(int i = 0; i < numElements; i++) {
    //double sig = p_sigma[i*numElements+i]; ## Skal det slettes? ##
    inversevolume[p_elements[i]] = 1.0;
  }

  // 3. Allocate the matrices
  epsilon.resize(SPARSE_MATRIX_SIZE * SPARSE_MATRIX_SIZE);
  alpha.resize(SPARSE_MATRIX_SIZE * SPARSE_MATRIX_SIZE);
  rmin.resize(SPARSE_MATRIX_SIZE * SPARSE_MATRIX_SIZE);
  v0.resize(SPARSE_MATRIX_SIZE * SPARSE_MATRIX_SIZE);

  // 4. Fill the matrix with the values we know.
  //   a) We know that A[i][j] = A[j][i]
  //   b) The position in the array is calculated as follows:
  //      (element number of atom i) * SPARSE_MATRIX_SIZE + (element number of atom j)
  for(int i = 0; i < numElements; i++) {
    for(int j = 0; j <= i; j++) {
      int position1 = p_elements[i] + p_elements[j] * SPARSE_MATRIX_SIZE;
      int position2 = p_elements[j] + p_elements[i] * SPARSE_MATRIX_SIZE;

      epsilon[position1] = p_epsilon[i * numElements + j];
      epsilon[position2] = p_epsilon[i * numElements + j];
      
      alpha[position1] = p_alpha[i * numElements + j];
      alpha[position2] = p_alpha[i * numElements + j];
      
      rmin[position1] = p_rmin[i * numElements + j];
      rmin[position2] = p_rmin[i * numElements + j];

      if (modified) {
	v0[position1] = epsilon[position1] * (exp(-2*alpha[position1] * (rCut - rmin[position1])) \
                                              - 2*exp(-alpha[position1] * (rCut - rmin[position1])));
	v0[position2] = epsilon[position2] * (exp(-2*alpha[position2] * (rCut - rmin[position2])) \
                                              - 2*exp(-alpha[position2] * (rCut - rmin[position2])));
      } else {
	v0[position1] = 0.0;
	v0[position2] = 0.0;
      }
    }       
  }

  DEBUGPRINT;
}


void Morse::SetAtoms(PyObject *pyatoms, Atoms* accessobj /* = NULL */)
{
  // First time SetAtoms is called, it should create an Atoms access
  // object, or use the one provided.  On subsequent calls, it should
  // test that no access object is provided.
  if (atoms == NULL) 
    {
      // The first time SetAtoms is being called some initialization is done.
      if (accessobj != NULL) {
	atoms = accessobj;
	AsapAtoms_INCREF(atoms);
      } else {
	atoms = new NormalAtoms();
      }
    }
  else
    {
      // This is at least the second time SetAtoms is called.
      if (accessobj != NULL)
	throw AsapError("Morse::SetAtoms called multiple times with accessobj != NULL");
    }
  assert(atoms != NULL);
  
  // The Morse potential is often used for thin gasses or the like.
  // In those cases, the neighbor list's "drift factor" should be
  // increased for more efficient reuse of the list.  Experiments show
  // that a broad optimum occurs when there is on average around 0.1
  // atom per cell in the NeighborCellLocator object used by
  // NeighborList.  This may pt. only be done for serial simulations

  atoms->Begin(pyatoms);
  IVec CPUlayout = atoms->GetNumberOfCells();
  if (CPUlayout[0] * CPUlayout[1] * CPUlayout[2] == 1)
    {
      //Serial simulation.  Test density.
      int n = atoms->GetNumberOfAtoms();
      double vol = atoms->GetVolume();
      double atomspercell = n / vol * pow(rCut,3);
      if (atomspercell < 0.05)
	{
	  // Low density, increase driftfactor
	  driftfactor = (pow(0.1 / atomspercell, 1.0/3) - 1.0) / 2;
	  assert(driftfactor >= DRIFTFACTOR);
	}
      else
	driftfactor = DRIFTFACTOR;  // High or medium density.
    }
  else
    driftfactor = DRIFTFACTOR;
  atoms->End();
}

//******************************************************************************
//                            Calculate cutoff radius
//******************************************************************************
/*
     Returns the lowest cutoff radius, where all potentials have droped below
     0.5 % of the minimum value.
*/
double Morse::CalculateRCut(int numElements2,
                            const vector<double> &alpha,
                            const vector<double> &rmin) 
{
  double min_alpha = alpha[0];
  double max_rmin = rmin[0];

  for (int i = 0; i < numElements2; i++) {
    if (alpha[i] < min_alpha) min_alpha = alpha[i];
    if (rmin[i] > max_rmin) max_rmin = rmin[i];
  }
  return (max_rmin + 6 / min_alpha);
}

//******************************************************************************
//                                 Allocate
//******************************************************************************
void Morse::Allocate()
{
  DEBUGPRINT;
  if (verbose)
    cerr << "Allocate(" << nAtoms << ") " << endl;
  assert(nAtoms != 0);
  atomicEnergies.resize(nAtoms);
  forces.resize(nSize);
  //stresses.resize(nAtoms);
  DEBUGPRINT;
}

//******************************************************************************
//                             CheckNeighborLists
//******************************************************************************
void Morse::CheckNeighborLists()
{
  DEBUGPRINT;
  if (counters.nblist == atoms->GetPositionsCounter() && neighborList != NULL &&
      !neighborList->IsInvalid())
    return;

  if (neighborList) {
      DEBUGPRINT;
      bool update = neighborList->CheckNeighborList();
      update = atoms->UpdateBeforeCalculation(update, rCut * (1 + driftfactor));
      if (update)
        neighborList->UpdateNeighborList();
      if ((nAtoms != atoms->GetNumberOfAtoms()) ||
          (nSize != nAtoms + atoms->GetNumberOfGhostAtoms())) {
          DEBUGPRINT;
          assert(update);
          nAtoms = atoms->GetNumberOfAtoms();
          nSize = nAtoms + atoms->GetNumberOfGhostAtoms();
          Allocate();
      }
  } else {
      DEBUGPRINT;
      atoms->UpdateBeforeCalculation(true, rCut * (1 + driftfactor));
      PyAsap_NeighborLocatorObject *nbl = PyAsap_NewNeighborList(atoms, rCut, driftfactor);
      neighborList_obj = (PyObject *)nbl;
      neighborList = nbl->cobj;
      assert(neighborList != NULL);
      neighborList->CheckAndUpdateNeighborList();
      nAtoms = atoms->GetNumberOfAtoms();
      nSize = nAtoms + atoms->GetNumberOfGhostAtoms();
      Allocate();
  }
  counters.nblist = atoms->GetPositionsCounter();
  DEBUGPRINT;
}

//******************************************************************************
//                                GetStresses
//******************************************************************************
const vector<SymTensor> &Morse::GetVirials(PyObject *pyatoms)
{
  throw AsapNotImplementedError("Morse: Stresses not implemented");
}


//******************************************************************************
//                                 GetForces 
//******************************************************************************
const vector<Vec> &Morse::GetForces(PyObject *pyatoms)
{
  assert(atoms != NULL);
  atoms->Begin(pyatoms);
  CheckNeighborLists();
  memset(&(forces[0]), 0, nSize * sizeof(Vec)); //zero the forces
  GetCartesianForces(forces);
  atoms->End();
  return forces;
}


//******************************************************************************
//                           GetCartesianForces
//******************************************************************************
//
// Calculates
//  
//                                  dV(r_{i,j})   dr_{i,j}
//   Fn_{k} = - Sum_{i} Sum_{j<>i} (----------- * --------)
//                                   dr_{i,j}     drn_{k}  
//
//                          dV(r_{k,i})   (rn_{k} - rn_{i})
//          = - Sum_{i<>k} (----------- * -----------------)
//                           dr_{k,i}          r_{k,i}
//
// where
//  
//   dV(r_{k,i})                     -alpha*(r_{k,i} - rmin)    -2*alpha*(r_{k,i} - rmin)
//   ----------- = 2*alpha*epsilon*(e                        - e                         )
//    dr_{k,i}
//
// and 
//
//   dr_{i,j}   (rn_{i} - rn_{j})
//   -------- = ----------------- * (d_{i,k} - d_{j,k})
//   drn_{k}         r_{i,j}
//
//   n = {x,y,z} (the dimension)
//
void Morse::GetCartesianForces(vector<Vec> &forces)
{
  DEBUGPRINT
  const asap_z_int *z = atoms->GetAtomicNumbers();
  int maxNeighbors = neighborList->MaxNeighborListLength();
  vector<int> neighbors(maxNeighbors);
  vector<Vec> diffs(maxNeighbors);
  vector<double> diffs2(maxNeighbors);

  DEBUGPRINT;

  // Iterate over all atoms
  for(int k=0; k<nAtoms; k++) {
    int size = maxNeighbors;
    int numNeighbors = neighborList->GetNeighbors(k, &neighbors[0], &diffs[0],
                                                  &diffs2[0], size);

    // Iterate over all neighbors inside rcut (in priciple this should be all other atoms)
    for(int a=0; a<numNeighbors; a++) {
      int i = neighbors[a];
      int pos = z[k] * SPARSE_MATRIX_SIZE + z[i]; //position in the sigma/epsilon array

      double r = sqrt(diffs2[a]);
      double dx = alpha[pos]*(r - rmin[pos]);
      double expdx = exp(-dx);
      double dF = alpha[pos]*epsilon[pos]*(expdx - expdx*expdx) / r;
      if (i < nAtoms)
        dF *= 2.0;
      forces[k] += dF * diffs[a];
      forces[i] -= dF * diffs[a];
    }
  }
  DEBUGPRINT
}

//******************************************************************************
//                             GetPotentialEnergy
//******************************************************************************
double Morse::GetPotentialEnergy(PyObject *pyatoms)
{
  DEBUGPRINT;
  assert(atoms != NULL);
  atoms->Begin(pyatoms);
  CheckNeighborLists();
  DEBUGPRINT;
  double e = CalculateEnergyAndEnergies();
  atoms->End();
  return e;
}

const vector<double> &Morse::GetPotentialEnergies(PyObject *pyatoms)
{
  DEBUGPRINT;
  assert(atoms != NULL);
  atoms->Begin(pyatoms);
  CheckNeighborLists();
  CalculateEnergyAndEnergies();
  atoms->End();
  DEBUGPRINT;
  return atomicEnergies;
}

//******************************************************************************
//                           CalculateEnergyAndEnergies
//******************************************************************************
double Morse::CalculateEnergyAndEnergies() {
  if (counters.energies != atoms->GetPositionsCounter()) {
    memset(&atomicEnergies[0], 0, nAtoms * sizeof(double));
    CalculateEnergyAndEnergies(atomicEnergies);
    counters.energies = atoms->GetPositionsCounter();
  }

  double energy = 0.0;
  assert(atomicEnergies.size() == nAtoms);
  for (int a = 0; a < nAtoms; a++) {
    energy +=  atomicEnergies[a];
  }
  return energy;
}

//******************************************************************************
//                           CalculateEnergyAndEnergies
//******************************************************************************
//  
//  V(r_ij) = epsilon *[exp(-2*alpha*(r_ij - r_min)) - 2*exp(-alpha*(r_ij - r_min))]
//  
void Morse::CalculateEnergyAndEnergies(vector<double>& atomicEnergies)
{
  DEBUGPRINT;
  const asap_z_int *z = atoms->GetAtomicNumbers();
  int maxNeighbors = neighborList->MaxNeighborListLength();
  vector<int> neighbors(maxNeighbors);
  vector<double> diffs2(maxNeighbors);
  vector<Vec> diffs(maxNeighbors);

  // Iterate over the real number of atoms
  for(int i = 0; i < nAtoms; i++) {
    int size = maxNeighbors;
    int numNeighbors = neighborList->GetNeighbors(i, &neighbors[0], &diffs[0],
                                                  &diffs2[0], size);
    for(int j = 0; j < numNeighbors; j++) {
      int n = neighbors[j];
      int position = z[i] * SPARSE_MATRIX_SIZE + z[n];

      double r = sqrt(diffs2[j]);
      double dx = alpha[position]*(r - rmin[position]);
      double expdx = exp(-dx);
      double dv = epsilon[position]*(expdx*expdx - 2*expdx) - v0[position];

      atomicEnergies[i] += 0.5 * dv;
      if (n < nAtoms)
        atomicEnergies[n] += 0.5 * dv;
	  }
  } 
  
  DEBUGPRINT;
}

long Morse::PrintMemory() const
{
  cerr << "*MEM*  Morse: Memory estimate not supported." << endl;
  return 0;
}

//******************************************************************************
//                                   Add
//******************************************************************************
//
// adds up the items in the vector, sorting them first by their absolute values
// (the result is pretty exact, since small numbers (1e-20 etc.) do not just
// drop out)
//
double Morse::Add(vector<double> &items) {
  int numItems = items.size();
  double result = 0.0;

  sort(items.begin(), items.end(), Morse::LessFabs);
  for(int a=0; a<numItems; a++) {
    result += items[a];
  }
  return result;
}

//******************************************************************************
//                                   LessFabs
//******************************************************************************
//compares the absolute values of a number
bool Morse::LessFabs(const double d1, const double d2) {
  return fabs(d1) < fabs(d2);
}

