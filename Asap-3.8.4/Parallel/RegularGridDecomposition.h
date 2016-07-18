// -*- C++ -*-
//
// RegularGridDecomposition.h: Domain decomposition based on a regular
// grid of processors.
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



#ifndef _REGULARGRIDDECOMPOSITION_H
#define _REGULARGRIDDECOMPOSITION_H

#include "DomainDecomposition.h"
#include "Vec.h"
#include <vector>
using std::vector;

namespace ASAPSPACE {

class ParallelAtoms;
class Communicator;


/// Domain decomposition based on a regular grid of processors.

///
/// \copydoc DomainDecomposition
///
/// RegularGridDecomposition assumes that the density of atoms is
/// approximately homogeneous.  Space is then divided in a regular
/// grid, and each processor becomes responsible for a cell on the
/// grid.
class RegularGridDecomposition : public DomainDecomposition
{
public:
  /// Create a decomposition from a SuperCell and a CPU layout (ncells).
  RegularGridDecomposition(const Vec superCell[3], const bool periodic[3],
			   const int nCells[3], Communicator *comm);
  virtual ~RegularGridDecomposition() {};
  void whichProcessor(ParallelAtoms *atoms,
                      vector< vector<int> > &processor,
                      vector<int> &migrateAway);
  void makeGhostExportLists(ParallelAtoms *atoms, double rCut,
                            vector< vector< pair<int, int> > > &ghostLists);
  const vector<int> *GetSendList() {return &sendlist;}
  const vector<int> *GetRecvList() {return &recvlist;}

  void GetTranslations(ParallelAtoms *atoms, vector<Vec> &translations) const;
  
private:
  void smallestBox(vector<Vec> &positions, const bool *periodic, int nAtoms, 
                   Vec &min, Vec &size);
  void InitializeNeighborsTable();
  void makeSendRecvLists();
  
  Vec minimum;  // Storing result from smallestBox between calls of 
  Vec size;     // whichProcessor and makeGhostExportList.
  Vec superCell[3];
  Communicator *comm;
  bool originalPeriodicity[3];
  int nCells[3];
  int nTotalCells[4];
  int nProcs;   // Total number of processors.
  int thisProcessor;  // The number of this processor.
  vector< pair<int, int> > neighbors[27];
  vector<int> sendlist;
  vector<int> recvlist;
};

} // end namespace

#endif //_REGULARGRIDDECOMPOSITION_H
