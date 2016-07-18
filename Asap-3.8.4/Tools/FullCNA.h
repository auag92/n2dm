// -*- C++ -*-
// FullCNA.h: Common Neighbor Analysis.
//
// Copyright (C) 2003-2011 Nicholas Bailey, Jakob Schiotz and Center
// for Individual Nanoparticle Functionality, Department of Physics,
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

#ifndef _FULLCNA_H
#define _FULLCNA_H

#include "AsapPython.h"
#include "Asap.h"
#include "AsapObject.h"
#include <vector>
#include <list>
#include <map>

namespace ASAPSPACE {

class Atoms;

// The CNA is bitmapped into an integer of at least 32 bits
typedef int cna_int;
#define CNA_CODE_BITS 8
#define CNA_CODE_MULT 256
#define CNA_CODE_MASK 0xFF;

typedef std::pair<int, int> intPair;
typedef std::pair<intPair, cna_int> CNA_data_item;
typedef std::map<intPair, cna_int > pair_CNA_map ;  // REMOVE?
typedef std::map<intPair, double> pair_dist_map;
typedef std::map<cna_int, PyObject *> pyCNA_map;

// The FullCNA calculator performs a Common Neighbor Analysis and returns
// the result in various ways.  It can only be used for a single calculation,
// and must be destroyed before control is returned to the user's Python
// script, as it keeps the atoms object open.
class FullCNA : public AsapObject
{
public:
  // Initialize the CNA calculator.
  FullCNA(PyObject *py_atoms, double cnaCutoff);

  ~FullCNA();

  virtual string GetName() const {return "FullCNA_internal";}

  // Call before calling MakeCNA to use per-element cutoffs.
  void SetMultipleCutoffs(std::map< intPair, double> &cutoffMap);

  // Return the raw CNA data
  PyObject *GetRawCNA();

  // Get CNA per atom
  PyObject *GetPerAtomCNA();

  // Get CNA globally
  PyObject *GetTotalCNA();

#if 0
  // Get Per-z CNA globally
  PyObject *GetPerZTotalCNA();
#endif

private:
  // Do the CNA calculation.
  void MakeCNA();

  std::vector<intPair> GetAdjacentBonds(int atom, std::vector<intPair> &bondsToProcess, std::vector<int> &atomsToProcess, std::vector<int> &atomsProcessed);
	
  bool Bonded(std::vector< std::vector<int> > &fullNeighbors, int at1, int at2);
	
  cna_int CNAonPair(intPair pPair, std::vector< std::vector<int> > &fullNeighbors);

  PyObject *PyCNAindex(cna_int idx);

  std::vector<CNA_data_item> cnaIndices;  // The raw CNA data
  pair_dist_map cnaDistances;

  double cnaCutoff;
  bool usingCutoffMap;
  bool ready;
  std::map< intPair, double> cutoffMap;

  Atoms *atoms;
  int nAtoms;
  PyObject *py_atoms;
  PyObject *py_nblist;
  pyCNA_map pythonCNAindices;
};

} // end namespace

#endif //_FULLCNA_H

