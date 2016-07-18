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

#include "FullCNA.h"
#include "Atoms.h"
#include "Exception.h"
#include "NeighborList.h"
#include "GetNeighborList.h"
#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
#include <assert.h>

#define CNABUF 1000

// TO DO: Change from NeighborList to NeighborList2, and then use the
// full neigbor list support instead of building a full list
// ourselves.  Milestone: Asap 3.



using namespace std;

/*  does a full Common Neighbor Analysis, following
    Clarke and Jonsson in PRE, vol 47, p 3975 (1993)
 */


FullCNA::FullCNA(PyObject *py_atoms, double cnaCutoff)
{
  this->py_atoms = py_atoms;
  Py_INCREF(py_atoms);
  this->cnaCutoff = cnaCutoff;
  this->usingCutoffMap = false;
  atoms = NULL;
  py_nblist = NULL;
  ready = false;
}

FullCNA::~FullCNA()
{
  if (atoms != NULL)
    AsapAtoms_DECREF(atoms);
  Py_XDECREF(py_nblist);  // Already cleared except in case of errors
  Py_DECREF(py_atoms);
  // Delete the map of Python-representations of CNA tuples.
  for (pyCNA_map::const_iterator i = pythonCNAindices.begin();
      i != pythonCNAindices.end(); ++i)
    {
      Py_DECREF(i->second);
    }
}

void FullCNA::SetMultipleCutoffs(std::map< intPair, double> &cutoffMap)
{
  this->cutoffMap.clear();
  map< intPair, double>::iterator mpIter;
  double maxCutoff = 0.;

  for( mpIter = cutoffMap.begin(); mpIter != cutoffMap.end(); mpIter++ )
    {
      intPair key = (*mpIter).first;
      double value = (*mpIter).second;
      this->cutoffMap[key] = value;
      if(key.first != key.second)
        {
          intPair reverseKey = make_pair(key.second, key.first);
          this->cutoffMap[reverseKey] = value;
        }
      if ((*mpIter).second > maxCutoff)
        maxCutoff = (*mpIter).second;
    }
  cnaCutoff = maxCutoff;
  usingCutoffMap = true;
}

void FullCNA::MakeCNA()
{
  // Don't call this function again
  ready = true;

  // Get Atoms access object (open) and neighbor list
  GetNbList_FromAtoms(py_atoms, cnaCutoff, &py_nblist, &atoms);

  nAtoms = atoms->GetNumberOfAtoms();

  if (atoms->GetNumberOfGhostAtoms() > 0)
    {
      // Here we need to make a copy of the atoms object in order
      // to get access to ghost neighbors of ghost atoms.  This
      // will be handled later.
      throw AsapNotImplementedError("Full CNA is not yet supported in parallel simulations.");
    }
  int nAtomsIncGhosts = nAtoms;

  // clear maps from previous calls
  cnaIndices.clear();

  NeighborLocator *nblistBonds = ((PyAsap_NeighborLocatorObject *) py_nblist)->cobj;

  // We make our own full neighbor list by folding the original list for two
  // reasons.  First, it makes it easier to reuse the original list, which is
  // likely to be half, while allowing us to have per-element-pair cutoffs.
  // Second, this is the way it was already coded.
  vector<vector<int> > fullNeighbors(nAtomsIncGhosts);
  int other[CNABUF] ;
  Vec rnb[CNABUF] ;
  double sqdist[CNABUF] ;
  const asap_z_int *atomicNumbers = atoms->GetAtomicNumbers();

  for(int atIdx = 0; atIdx < nAtomsIncGhosts; ++atIdx)
    {
      int listSpaces = CNABUF;
      int numNbrs = nblistBonds->GetNeighbors(atIdx, other, rnb, sqdist, listSpaces,
          cnaCutoff);
      int thisType = atomicNumbers[atIdx];

      if(!usingCutoffMap)
        for(int nbrIdx = 0; nbrIdx < numNbrs; ++nbrIdx)
          {
            int nbrAt = other[nbrIdx];
            fullNeighbors[atIdx].push_back(nbrAt);
            fullNeighbors[nbrAt].push_back(atIdx);
          }
      else
        {
          for(int nbrIdx = 0; nbrIdx < numNbrs; ++nbrIdx)
            {
              int nbrAt = other[nbrIdx];
              int nbrType = atomicNumbers[nbrAt];
              if( sqrt(sqdist[nbrIdx]) < cutoffMap[make_pair(thisType, nbrType)])
                {
                  fullNeighbors[atIdx].push_back(nbrAt);
                  fullNeighbors[nbrAt].push_back(atIdx);
                }
            }
        }
    }

  // now we do CNA analysis on all pairs separated by less than cnaCutoff
  for(int atIdx = 0; atIdx < nAtoms; ++atIdx)
    {
      for(vector<int>::const_iterator nbrIdx = fullNeighbors[atIdx].begin();
          nbrIdx != fullNeighbors[atIdx].end(); ++nbrIdx)
        {
          int nbrAt = *nbrIdx;
          if (nbrAt < atIdx)
            {
              intPair orderedPair = make_pair(nbrAt,atIdx);
              // Do the hard work
              cna_int cnaIndex = CNAonPair(orderedPair, fullNeighbors);
              cnaIndices.push_back(make_pair(orderedPair, cnaIndex));
            }
        }

    }

  // Release the neigbor list by releasing the Python object.
  Py_CLEAR(py_nblist);
  atoms->End();
}



/*
    -- finds the bonds in bondsToProcess containing atom
    -- removes them from bondsToProcess
    -- adds atom to atomsProcessed (already removed from atomsToProcess)
    -- adds all atoms in the found bonds not in atomsProcessed to atomList
    -- returns the bonds found
 */

vector<intPair> FullCNA::GetAdjacentBonds(int atom, vector<intPair> &bondsToProcess, vector<int> &atomsToProcess, vector<int> &atomsProcessed)
{

  vector<intPair> adjacentBonds;

  vector<intPair>::iterator bondIter;

  for( bondIter = bondsToProcess.begin(); bondIter != bondsToProcess.end(); 	bondIter++)
    {
      intPair bond = *bondIter;
      if(atom == bond.first || atom == bond.second)
        {
          adjacentBonds.push_back(bond);
        }
    }

  // erase adjacent bonds from bondsToProcess
  for( bondIter = adjacentBonds.begin(); bondIter != adjacentBonds.end(); bondIter++ )
    {
      std::vector<intPair>::iterator bond = find(bondsToProcess.begin(),bondsToProcess.end(),*bondIter);
      assert(bond != bondsToProcess.end());
      bondsToProcess.erase(bond);
      //bondsToProcess.remove(*bondIter);
    }

  atomsProcessed.push_back(atom);
  for( bondIter = adjacentBonds.begin(); bondIter != adjacentBonds.end(); bondIter++ )
    {
      intPair bond = *bondIter;
      int nextAtom = bond.first;
      if(find(atomsProcessed.begin(),atomsProcessed.end(),nextAtom) ==  atomsProcessed.end() && find(atomsToProcess.begin(),atomsToProcess.end(),nextAtom) == atomsToProcess.end())
        atomsToProcess.push_back(nextAtom);
      nextAtom = bond.second;
      if(find(atomsProcessed.begin(),atomsProcessed.end(),nextAtom) ==  atomsProcessed.end() && find(atomsToProcess.begin(),atomsToProcess.end(),nextAtom) == atomsToProcess.end())
        atomsToProcess.push_back(nextAtom);
    }
  return adjacentBonds;
}

bool FullCNA::Bonded(vector< vector<int> > &fullNeighbors, int at1, int at2)
{
  bool result = false;
  vector<int> &nbrs1 = fullNeighbors[at1];
  for( unsigned int idx = 0; idx < nbrs1.size(); ++idx )
    {
      if(nbrs1[idx] == at2)
        {
          result = true;
          break;
        }

    }
  return result;
}

cna_int FullCNA::CNAonPair(intPair pPair, vector< vector<int> > &fullNeighbors)
{
  int at1 = pPair.first;
  int at2 = pPair.second;

  // this will be the actual CNA index
  int index[3];


  vector<int> commonNeighbors;
  vector<int> &nbrs1 = fullNeighbors[at1];
  vector<int> &nbrs2 = fullNeighbors[at2];
  int nNbrs1 = nbrs1.size();
  int nNbrs2 = nbrs2.size();


  for(int nbrIdx1 = 0; nbrIdx1 < nNbrs1; ++nbrIdx1 )
    {
      int count = 0;
      for(int nbrIdx2 = 0; nbrIdx2 < nNbrs2; ++nbrIdx2 )
        {
          if(nbrs1[nbrIdx1] == nbrs2[nbrIdx2])
            {
              commonNeighbors.push_back(nbrs1[nbrIdx1]);
              count += 1;
            }
        }
      assert(count <= 1);
    }

  // the first index is the number of common neighbors
  int numCommonNbrs = commonNeighbors.size();
  index[0] = numCommonNbrs;

  vector<intPair> bondsToProcess;
  for(int cn1 = 0; cn1 <numCommonNbrs; ++cn1)
    {
      for(int cn2 = cn1; cn2 < numCommonNbrs; ++cn2 )
        {
          int nb1 = commonNeighbors[cn1];
          int nb2 = commonNeighbors[cn2];

          if(Bonded(fullNeighbors, nb1,nb2))
            {
              if(nb2 < nb1)
                {int tmp = nb1; nb1 = nb2; nb2 = tmp;}
              // nb1,nb2 is now an ordered pair
              intPair thisBond = make_pair(nb1,nb2);
              bondsToProcess.push_back(thisBond);
            }
        }
    }
  // the second index is the number of bonds among common neighbors
  int numCNBonds = bondsToProcess.size();
  index[1] = numCNBonds;

  // now group the common bonds into clusters
  int maxClusterSize = 0;
  int totalBondsInClusters = 0;
  int nBondsLeft = bondsToProcess.size();

  while(nBondsLeft > 0)
    {
      // make a new cluster starting with the first remaining bond to
      // be processed
      vector<intPair> thisCluster;
      intPair newBond = bondsToProcess.back();
      bondsToProcess.pop_back();
      thisCluster.push_back(newBond);
      vector<int> atomsToProcess, atomsProcessed;
      atomsToProcess.push_back(newBond.first);
      atomsToProcess.push_back(newBond.second);

      int nAtomsLeft = atomsToProcess.size();
      while(nAtomsLeft != 0)
        {
          int nextAtom = atomsToProcess.back();
          atomsToProcess.pop_back();
          vector<intPair> adjacentBonds = GetAdjacentBonds(nextAtom,bondsToProcess,atomsToProcess, atomsProcessed);
          thisCluster.insert(thisCluster.end(),adjacentBonds.begin(),adjacentBonds.end());
          nAtomsLeft = atomsToProcess.size();
        }

      int thisClusterSize = thisCluster.size();
      totalBondsInClusters += thisClusterSize;
      // is it bigger than the maximum?
      if(thisClusterSize > maxClusterSize)
        maxClusterSize = thisCluster.size();

      nBondsLeft = bondsToProcess.size();
    }
  // sanity check
  assert(totalBondsInClusters == numCNBonds);

  // the third index is the size of the largest cluster
  index[2] = maxClusterSize;

  // Now code the numbers into a single integer
  int codeindex = (index[0] * CNA_CODE_MULT + index[1]) * CNA_CODE_MULT
      + index[2];
  return codeindex;
}

// Return a 3xN NumPy array of int32 with all the CNA data. For each bond, three
// numbers are given.  The first two are the indices of the two atoms involved in
// the bond, the third is the CNA encoded as a single 32-bit integer with the three
// CNA digits as the three least significant bytes.
PyObject *FullCNA::GetRawCNA()
{
  if (!ready)
    MakeCNA();
  npy_intp dims[2];
  dims[0] = cnaIndices.size();
  dims[1] = 3;
  PyObject *result = PyArray_SimpleNew(2, dims, NPY_INT32);
  if (result == NULL)
    return NULL;
  npy_int32 *ptr = (npy_int32 *) PyArray_DATA((PyArrayObject *) result);
  for (int i = 0; i < dims[0]; ++i)
    {
      *(ptr++) = cnaIndices[i].first.first;
      *(ptr++) = cnaIndices[i].first.second;
      *(ptr++) = cnaIndices[i].second;
    }
  assert(3 * dims[0] == ptr - (npy_int32 *) PyArray_DATA((PyArrayObject *) result));
  return result;
}

// Return a Python list with CNA data for each atom.  Each element is
// a dictionary with the CNA tuples as keys and bond numbers as values.
PyObject *FullCNA::GetPerAtomCNA()
{
  typedef std::map<cna_int, int> atomdata;
  class Abort {};

  if (!ready)
    MakeCNA();

  // First, collect the data into C++ structures, translate to Python later.
  vector<atomdata> buffer(nAtoms);
  for (vector<CNA_data_item>::const_iterator i = cnaIndices.begin();
      i != cnaIndices.end(); ++i)
    {
      intPair atomnumbers = i->first;
      cna_int cna = i->second;
      buffer[atomnumbers.first][cna] += 1;
      buffer[atomnumbers.second][cna] += 1;
    }

  // Now, convert to Python.
  PyObject *result = PyList_New(nAtoms);
  if (result == NULL)
    return NULL;
  try
  {
      for (int i = 0; i < nAtoms; i++)
        {
          PyObject *dict = PyDict_New();
          if (dict == NULL) throw Abort();
          PyList_SET_ITEM(result, i, dict);
          const atomdata &ccna = buffer[i];
          for(atomdata::const_iterator j = ccna.begin();
              j != ccna.end(); ++j)
            {
              PyObject *pycna = PyCNAindex(j->first);
              PyObject *pycount = PyInt_FromLong((long) j->second);
              int x = PyDict_SetItem(dict, pycna, pycount);
              Py_DECREF(pycna);
              Py_DECREF(pycount);
              if (x) throw Abort();
            }
        }
  }
  catch (Abort &a)
  {
      Py_DECREF(result);
      result = NULL;
  }
  return result;
}

// Return a dictionary with the CNA for all the atoms.
PyObject *FullCNA::GetTotalCNA()
{
  typedef std::map<cna_int, int> atomdata;

  if (!ready)
    MakeCNA();

  // First, collect the data into C++ structures, translate to Python later.
  atomdata totalcna;
  for (vector<CNA_data_item>::const_iterator i = cnaIndices.begin();
      i != cnaIndices.end(); ++i)
    {
      totalcna[i->second] += 1;
    }

  // Now, convert to Python.
  PyObject *result = PyDict_New();
  if (result == NULL)
    return NULL;
  for(atomdata::const_iterator j = totalcna.begin();
      j != totalcna.end(); ++j)
    {
      PyObject *pycna = PyCNAindex(j->first);
      if (pycna == NULL)
        {
          Py_DECREF(result);
          return NULL;
        }
      PyObject *pycount = PyInt_FromLong((long) j->second);
      int x = PyDict_SetItem(result, pycna, pycount);
      Py_DECREF(pycna);
      Py_DECREF(pycount);
      if (x) {
          Py_DECREF(result);
          return NULL;
      }
    }
  return result;
}

#if 0
// Return a dictionary with pairs of atomic numbers as keys, and dictionaries
// with CNA for all such bonds as values.
PyObject *FullCNA::GetPerZTotalCNA()
{
  typedef std::pair<asap_z_int, asap_z_int> zpair;
  typedef std::map<cna_int, int> cnadata;
  typedef std::map<zpair, cnadata> bonddata;

  if (!ready)
    MakeCNA();

  // First, collect the data into C++ structures, translate to Python later.
  bonddata totalcna;
  const asap_z_int *atomicnumbers = atoms->GetAtomicNumbers();
  for (vector<CNA_data_item>::const_iterator i = cnaIndices.begin();
      i != cnaIndices.end(); ++i)
    {
      asap_z_int z1 = atomicnumbers[i->first.first];
      asap_z_int z2 = atomicnumbers[i->first.second];
      zpair z = zpair(z1, z2);
      if (z2 < z1)
        z = zpair(z2, z1);
      // z is now an ordered pair
      totalcna[z][i->second] += 1;
    }

  // Now, convert to Python
  PyObject *result = PyDict_New();
  if (result == NULL)
    return NULL;
  for (bonddata::const_iterator i = totalcna.begin();
      i != totalcna.end(); ++i)
    {
      PyObject *pyz = Py_BuildValue("ii", i->first.first, i->first.second);
      if (pyz == NULL)
        {
          Py_DECREF(result);
          return NULL;
        }
      PyObject *pydict = PyDict_New();
      if (pydict == NULL)
        {
          Py_DECREF(pyz);
          Py_DECREF(result);
          return NULL;
        }
      int x = PyDict_SetItem(result, pyz, pydict);
      Py_DECREF(pyz);
      if (x)
        {
          Py_DECREF(pydict);
          Py_DECREF(result);
          return NULL;
        }
      // Translate the CNA dictionary.
      for (cnadata::const_iterator j = i->second.begin();
          j != i->second.end(); ++j)
        {
          PyObject *pycna = PyCNAindex(j->first);
          if (pycna == NULL)
            {
              Py_DECREF(result);
              Py_DECREF(pydict);
              return NULL;
            }
          PyObject *pycount = PyInt_FromLong((long) j->second);
          int x = PyDict_SetItem(pydict, pycna, pycount);
          Py_DECREF(pycount);
          Py_DECREF(pycna);
          if (x)
            {
              Py_DECREF(pydict);
              Py_DECREF(result);
              return NULL;
            }
        }
      Py_DECREF(pydict);
    }
  return result;
}
#endif

// Convert a coded CNA to a Python 3-tuple, making sure to cache and reuse it.
PyObject *FullCNA::PyCNAindex(cna_int idx)
{
  pyCNA_map::iterator found = pythonCNAindices.find(idx);
  if (found != pythonCNAindices.end())
    {
      // Found the Python representation of the CNA triplet
      Py_INCREF(found->second);
      return found->second;
    }
  else
    {
      // Make the Python representation
      int index[3];
      index[2] = idx & CNA_CODE_MASK;
      index[1] = (idx >> CNA_CODE_BITS) & CNA_CODE_MASK;
      index[0] = (idx >> (2*CNA_CODE_BITS)) & CNA_CODE_MASK;
      PyObject *result = Py_BuildValue("iii", index[0], index[1], index[2]);
      if (result != NULL)  // Failure handled by calling function!
        {
          pythonCNAindices.insert(pair<cna_int,PyObject *>(idx, result));
          Py_INCREF(result); // We store it
        }
      return result;
    }
}
