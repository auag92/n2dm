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
 * Calculates the Lennard Jones potential
 */


// This implementation supports parallel simulations.  Due to the
// simplicity of Lennard-Jones, everything is handled by the
// ParallelAtoms and ParallelPotential classes, this potential only
// has to be prepared for parallel simulations in two ways:
//
// 1. After a neighborlist update, the number of atoms may have
//    changed, and arrays may have to be reallocated.
//
// 2. The neighborlist may return atoms with number higher than
//    nAtoms, these are ghosts of atoms residing on other processor,
//    and only have positions and atomic numbers.  It must be avoided
//    to update energies, forces and stresses for these, as that would
//    overwrite the end of the arrays.

#include "AsapPython.h"
#include "Asap.h"
#include "LennardJones.h"
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

#if 0
#define VERB(x) if (verbose == 1) cerr << x
#else
#define VERB(x)
#endif


#ifdef SUN_SQRT_LITERAL_HACK
// SunOS compiler: sqrt does not work with a literal number (cannot decide if
// it is a double or a long double).
#define sqrt(x) sqrt((double) x)
#endif

// Standard mapping of the six independent parts of the stress tensor to
// vector notation
const static int stresscomp[3][3] = {{0, 5, 4}, {5, 1, 3}, {4, 3, 2}};



//******************************************************************************
//                             LennardJones
//******************************************************************************
LennardJones::LennardJones(PyObject *self,
                           int numElements, const std::vector<int> &elements,
			   const std::vector<double> &epsilon,
			   const std::vector<double> &sigma,
			   const std::vector<double> &masses, double rCut,
			   bool modified) : Potential(self)
{
  DEBUGPRINT;
  atoms = NULL;
  this->nAtoms      = 0;
  neighborList= 0;
  neighborList_obj= 0;
  driftfactor = 0.05;
  this->modified    = modified;
  this->numElements = numElements;
  this->latticeConstant = 0.0;
  
  memset(&counters, 0x0, sizeof(counters));

  if(rCut<0)
    this->rCut = CalculateRCut(numElements*numElements, sigma);
  else
    this->rCut = rCut;

  Internalize(numElements, elements, epsilon, sigma, masses);
  //Print();
  DEBUGPRINT;
}

LennardJones::~LennardJones()
{
  DEBUGPRINT;
  Py_XDECREF(neighborList_obj);
  if (atoms != NULL)
    AsapAtoms_DECREF(atoms);
}


//******************************************************************************
//                             Print
//******************************************************************************
void LennardJones::Print() {

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
//                             Internalize
//******************************************************************************
void LennardJones::Internalize(int p_numElements,
			       const std::vector<int> &p_elements,
			       const std::vector<double> &p_epsilon,
			       const std::vector<double> &p_sigma,
			       const std::vector<double> &p_masses) 
{
  DEBUGPRINT;
  //Print();
  double rCut6 = pow(rCut, 6.0);   //cCut^6
  double rCut12 = pow(rCut, 12.0); //cCut^12

  //0. complete upper triangular matrices
  std::vector<double> full_epsilon; full_epsilon.resize(p_numElements*p_numElements);
  std::vector<double> full_sigma;   full_sigma.resize(p_numElements*p_numElements);
  for(int i=0;i<p_numElements;i++) {
    for(int j=0;j<=i;j++) {
        full_epsilon[j*p_numElements+i]=p_epsilon[i*p_numElements+j];
        full_sigma[j*p_numElements+i]=p_sigma[i*p_numElements+j];
        full_epsilon[i*p_numElements+j]=p_epsilon[i*p_numElements+j];
        full_sigma[i*p_numElements+j]=p_sigma[i*p_numElements+j];
    }
  }
  //copy the masses
  masses.resize(SPARSE_MATRIX_SIZE);
  memset(&masses[0], -1, sizeof(double)*SPARSE_MATRIX_SIZE);
  for(int i=0; i<p_numElements; i++) 
    masses[p_elements[i]]   = p_masses[i];
  atomicvolume.resize(SPARSE_MATRIX_SIZE);
  memset(&atomicvolume[0], -1, sizeof(double)*SPARSE_MATRIX_SIZE);
  double volumefactor = 0.94;  // 1.1^3 / sqrt(2)
  for(int i=0; i<p_numElements; i++)
    {
      double sig = p_sigma[i*p_numElements+i];
      atomicvolume[p_elements[i]]   = (volumefactor * sig * sig * sig);
    }
  //1. allocate the matrices
  epsilon.resize(SPARSE_MATRIX_SIZE*SPARSE_MATRIX_SIZE);
  sigma6.resize(SPARSE_MATRIX_SIZE*SPARSE_MATRIX_SIZE);
  sigma12.resize(SPARSE_MATRIX_SIZE*SPARSE_MATRIX_SIZE);
  v0.resize(SPARSE_MATRIX_SIZE*SPARSE_MATRIX_SIZE);

  //2. Fill the matrix with the values we know.
  //   a) We know that sigma[i][j] = sigma[j][i]
  //   b) The position in the array is calculated as follows:
  //      (element number of atom i)*SPARSE_MATRIX_SIZE + (element number of atom j)
  for(int i=0; i<p_numElements; i++) {
     for(int j=0; j<=i; j++) {
       double newEpsilon = full_epsilon[i*p_numElements+j];
       double newSigma1 = full_sigma[i*p_numElements+j];
       double newSigma3 = newSigma1 * newSigma1 * newSigma1;
       double newSigma6   = newSigma3 * newSigma3 ;
       int position1 = p_elements[i]+p_elements[j]*SPARSE_MATRIX_SIZE;
       int position2 = p_elements[j]+p_elements[i]*SPARSE_MATRIX_SIZE;
       epsilon[position1] = newEpsilon;
       epsilon[position2] = newEpsilon;
       //sigma^6
       sigma6  [position1] = newSigma6;
       sigma6  [position2] = newSigma6;
       //sigma^12
       sigma12 [position1] = newSigma6*newSigma6;
       sigma12 [position2] = newSigma6*newSigma6;
       //v0
       if (modified)
	 {
	   v0[position1] = 2*epsilon[position1]*(sigma12[position1]/rCut12-sigma6[position1]/rCut6);
	   v0[position2] = 2*epsilon[position2]*(sigma12[position2]/rCut12-sigma6[position2]/rCut6);
	 }
       else
	 {
	   v0[position1] = 0.0;
	   v0[position2] = 0.0;
	 }
     }       
  }
  DEBUGPRINT;
}


void LennardJones::SetAtoms(PyObject *pyatoms, Atoms* accessobj /* = NULL */)
{
  if (atoms != NULL)
    {
      // SetAtoms should only do anything the first time it is called.
      // Subsequent calls should just check for accessobj being NULL.
      if (accessobj != NULL)
	throw AsapError("LennardJones::SetAtoms called multiple times with accessobj != NULL");
      // SetAtoms should not do anything if called more than once!
      return;
    }

  // The first time SetAtoms is being called some initialization is done.
  if (accessobj != NULL)
    {
      atoms = accessobj;
      AsapAtoms_INCREF(atoms);
    }
  else
    atoms = new NormalAtoms();
  assert(atoms != NULL);
}

//******************************************************************************
//                             CalculateIDs
//******************************************************************************
/*
   calculates or defines the default rCut
*/
double LennardJones::CalculateRCut(int numElements,
				   const vector<double> &sigma) 
{
   //we return 3*min(sigma) 
   double result = sigma[0];
   int i = 0;
   while(++i<numElements)
     if(sigma[i]>result)
       result = sigma[i];
    return 3*result;
}

//******************************************************************************
//                             Allocate
//******************************************************************************
void LennardJones::Allocate()
{
  DEBUGPRINT;
  if (verbose)
    cerr << "Allocate(" << nAtoms << ") " << endl;
  assert(nAtoms != 0);
  atomicEnergies.resize(nAtoms);
  forces.resize(nSize);
  virials.resize(nSize);
  DEBUGPRINT;
}

//******************************************************************************
//                             CheckNeighborLists
//******************************************************************************
void LennardJones::CheckNeighborLists()
{
  DEBUGPRINT;
  if (counters.nblist == atoms->GetPositionsCounter() && neighborList != NULL
      && !neighborList->IsInvalid())
    return;
  if (neighborList)
    {
      DEBUGPRINT;
      bool update = neighborList->CheckNeighborList();
      update = atoms->UpdateBeforeCalculation(update,
					  rCut * (1 + driftfactor));
      if (update)
	neighborList->UpdateNeighborList();
      if ((nAtoms != atoms->GetNumberOfAtoms()) ||
          (nSize != nAtoms + atoms->GetNumberOfGhostAtoms()))
	{
	  DEBUGPRINT;
	  assert(update);
	  nAtoms = atoms->GetNumberOfAtoms();
	  nSize = nAtoms + atoms->GetNumberOfGhostAtoms();
	  Allocate();
	}
    }
  else
    {
      DEBUGPRINT;
      atoms->UpdateBeforeCalculation(true, rCut * (1 + driftfactor));
      PyAsap_NeighborLocatorObject *nbl = PyAsap_NewNeighborList(atoms, rCut,
								 driftfactor);
      neighborList_obj = (PyObject *)nbl;
      neighborList = dynamic_cast<NeighborList*>(nbl->cobj);
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
//                             GetPotentialEnergy
//******************************************************************************
double LennardJones::GetPotentialEnergy(PyObject *pyatoms)
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

const vector<double> &LennardJones::GetPotentialEnergies(PyObject *pyatoms)
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
//                               GetStresses 
//******************************************************************************
const vector<SymTensor> &LennardJones::GetVirials(PyObject *pyatoms)
{
  assert(atoms != NULL);
  atoms->Begin(pyatoms);
  CheckNeighborLists();
  if (nAtoms != virials.size())
    virials.resize(nSize);
  memset(&virials[0], 0, nSize * sizeof(SymTensor));
  GetVirials(atoms->GetMomenta(), &virials[0]);
  atoms->End();
  return virials;
}

//******************************************************************************
//                               Stress 
//******************************************************************************
/*
                 p_{i,alpha} * p_{i, beta}              dV(r_{ij}   r_{ij,alpha)*r_{ij,beta)
      Calculates ------------------------- - Sum( 0.5 * --------- * ------------------------ )
                          m_{i}                          dr_{ij}            r_{ij)


      where

  
       dV(r_ij)                       sigma^{12}      sigma^{6}
       -------- =   4*epsilon*[(-12)*(----------) +6*(---------) ]
       dr_{ij}                        r_ij^{13}        r_ij^{7}


      resulting in

                 p_{i,alpha} * p_{i, beta}                         sigma^{12}    sigma^{6}    
                 ------------------------- - 12*Sum(epsilon*[(-2)*(----------) + (---------) ] * r_{ij,alpha)*r_{ij,beta)
                          m_{i}                                     r_ij^{14}     r_ij^{8}             
*/


//******************************************************************************
//                               GetStresses
//******************************************************************************
void LennardJones::GetVirials(const Vec *momenta, SymTensor *stresses)
{
  DEBUGPRINT;

  int maxNeighbors=neighborList->MaxNeighborListLength();//the max number of neighbors in neighborlist
  int numNeighbors; //the current number of neighbors in the iteration

  vector<Vec> diffs(maxNeighbors);
  vector<int> neighbors(maxNeighbors);
  vector<double> diffs2(maxNeighbors);

  const asap_z_int *atomicNumberPtr;
  int  matrixLength;

  atomicNumberPtr = atoms->GetAtomicNumbers();
  matrixLength    = SPARSE_MATRIX_SIZE;
  
  for(int i=0; i<nAtoms; i++) {
     for(int alpha = 0; alpha < 3; alpha++)
        for(int beta = alpha; beta < 3; beta++)
           {
	     int j = stresscomp[alpha][beta];
	     double result;
	     int MaxNB = maxNeighbors;
	     numNeighbors = neighborList->GetNeighbors(i, &neighbors[0],
						       &diffs[0], &diffs2[0],
						       MaxNB);
	     
	     for(int a=0; a<numNeighbors; a++) {//iterate over the neighbors
	       int neighbor = neighbors[a];
	       
	       double dist4 =diffs2[a]*diffs2[a];  //=dist^4
	       double dist8 =dist4*dist4;          //=dist^8
	       double dist14=dist8*dist4*diffs2[a];//=dist^14
	       
	       int arrayPos; //position in the sigma/epsilon array
	       
	       arrayPos = atomicNumberPtr[i]*matrixLength+atomicNumberPtr[neighbor];
	       
	       result = epsilon[arrayPos]*(sigma6[arrayPos]/dist8 -
					      2*sigma12[arrayPos]/dist14) *
		 (diffs[a][alpha]) * (diffs[a][beta]);
	       if (neighbor < nAtoms)
	         result *= 12.0;
	       else
	         result *= 6.0;
	       stresses[i][j] += result;
	       stresses[neighbor][j] += result;
	     }
           }
  }
  DEBUGPRINT;
}

void LennardJones::GetAtomicVolumes(vector<double> &volumes)
{
#ifdef SHAREVOLUME
  volumes.clear();
#else
  volumes.resize(nAtoms);
  for (int i = 0; i < nAtoms; i++)
    volumes[i] = atomicvolume[z[i]];
#endif
}

//******************************************************************************
//                               GetCartesianForces 
//******************************************************************************
const vector<Vec> &LennardJones::GetForces(PyObject *pyatoms)
{
  assert(atoms != NULL);
  atoms->Begin(pyatoms);
  CheckNeighborLists();
  assert(nSize >= nAtoms);
  assert(forces.size() == nSize);
  memset(&(forces[0]), 0, nSize * sizeof(Vec)); //zero the forces
  GetCartesianForces(forces);
  atoms->End();
  return forces;
}

//******************************************************************************
//                           GetCartesianForces
//******************************************************************************
//
//                                       dV(r_{in})         dr_{in}
// Calculates F_{alpha}(n) = Sum_{i<>n} (---------- * ( -------------- ))
//                                        dr_{in)         dr_{n,alpha}  
//
//
// where
//
//  
//  dV(r_in)                       sigma^{12}      sigma^6
//  -------- =   4*epsilon*[(-12)*(----------) +6*(-------) ]
//  dr_{in}                        r_in^{13}        r_in^7
//
// and 
//
//       dr_{in}                                  - (r_{i0} - r_{n,alpha})
//  ( -------------- )) = ---------------------------------------------------------------------------
//     dr_{n,alpha}        ((r_{n0} - r_{i0})^2 + (r_{n0} - r_{i0})^2 + (r_{n0} - r_{i0})^2 ) ^ (1/2)
//                          
//                      =  - (r_{i0} - r_{n,alpha})  
//                         ------------------------
//                               | r_{n,i} |
//
//   alpha = {0,1,2} (the dimension)
//
//
// Simplified: 
// ===========
//
//                                             sigma^{12}      sigma^6
// F_{alpha}(n) = -24*Sum_{i<>n} epsilon*[ -2*(----------) + *(-------) ]*(r_{i0} - r_{n,alpha}) ]
//                                             r_in^{14}        r_in^8
//
void LennardJones::GetCartesianForces(vector<Vec>& forces)
{
  DEBUGPRINT
  const asap_z_int *z = atoms->GetAtomicNumbers();
  int maxNeighbors=neighborList->MaxNeighborListLength();
  vector<int> neighbors(maxNeighbors);
  vector<Vec> diffs(maxNeighbors);
  vector<double> diffs2(maxNeighbors);

  int size;

  DEBUGPRINT;
  const asap_z_int *zz;
  int msize;
  zz = z;
  msize = SPARSE_MATRIX_SIZE;
  for(int n=0; n<nAtoms; n++) { //iterate over all atoms
    size = maxNeighbors;
    int numNeighbors = neighborList->GetNeighbors(n, &neighbors[0], &diffs[0],
						  &diffs2[0], size);
    for(int a=0; a<numNeighbors; a++) {//iterate over all neighbors
      int i = neighbors[a];
      double dist4 =diffs2[a]*diffs2[a]; //=dist^4
      double dist8 =dist4*dist4;         //=dist^12
      double dist14=dist8*dist4*diffs2[a];
      int arrayPos; //position in the sigma/epsilon array
      double dV;

      arrayPos = zz[n]*msize + zz[i];
      
      dV = epsilon[arrayPos]*(sigma6[arrayPos]/dist8 - 2*sigma12[arrayPos]/dist14);///sqrt(diffs2[a]);
      if (i < nAtoms)
        dV *= -24.0;
      else
        dV *= -12.0;
      forces[n] -= dV * diffs[a];
      forces[i] += dV * diffs[a];
    }
  }
  DEBUGPRINT
}

//******************************************************************************
//                                   Add
//******************************************************************************
//
// adds up the items in the vector, sorting them first by their absolute values
// (the result is pretty exact, since small numbers (1e-20 etc.) do not just
// drop out)
//
double LennardJones::Add(vector<double> &items) {
  int numItems = items.size();
  double result = 0.0;

  sort(items.begin(), items.end(), LennardJones::LessFabs);
  for(int a=0;a<numItems;a++) {
    result += items[a];
  }
  return result;
}

//******************************************************************************
//                                   LessFabs
//******************************************************************************
//compares the absolute values of a number
bool LennardJones::LessFabs(const double d1, const double d2) {
  return fabs(d1)<fabs(d2);
}


//******************************************************************************
//                           CalculateEnergyAndEnergies
//******************************************************************************
double LennardJones::CalculateEnergyAndEnergies() {
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
///  
///                         sigma         sigma
///  V(r_ij) = 4*epsilon*[(-------)^12 - (------)^6 ]
///                         r_{ij}         r_{ij}
///  
void LennardJones::CalculateEnergyAndEnergies(vector<double>& atomicEnergies)
{
  DEBUGPRINT;
      const asap_z_int *z; //lookup for atomic numbers
      int size;
      int cols;     //the number of columns in sigma6, epsilon matrices
      int maxNeighbors=neighborList->MaxNeighborListLength();
      vector<int> neighbors(maxNeighbors);
      vector<double> diffs2(maxNeighbors);
      vector<Vec> diffs(maxNeighbors);

      z = atoms->GetAtomicNumbers();
      cols = SPARSE_MATRIX_SIZE;

      for(int i = 0;i<nAtoms;i++) {//iterate over the real number of atoms
	  size = maxNeighbors;
	  int n = neighborList->GetNeighbors(i, &neighbors[0], &diffs[0],
					     &diffs2[0], size);
	    for(int j = 0;j<n;j++) {
  	         int nn = neighbors[j];
		 int position = z[i]*cols+z[nn];
		 double s6_div_r6 = sigma6[position]/diffs2[j]/diffs2[j]/diffs2[j];
		 double s12_div_r12  = s6_div_r6 * s6_div_r6 ;
		 double result = 2.0*epsilon[position]*(s12_div_r12-s6_div_r6)-v0[position];
		 atomicEnergies[i] += result;
		 if (nn < nAtoms)
		   atomicEnergies[nn] += result;
	    }
      } 
      DEBUGPRINT;
}

long LennardJones::PrintMemory() const
{
  cerr << "*MEM*  LennardJones: Memory estimate not supported." << endl;
  return 0;
}
