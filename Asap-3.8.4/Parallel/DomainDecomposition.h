// -*- C++ -*-
//
// DomainDecomposition.h: Abstract base class for domain decompositions.
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


#ifndef _DOMAINDECOMPOSITION_H
#define _DOMAINDECOMPOSITION_H

#include <Asap.h>
#include <Vec.h>
#include <vector>
#include <utility>
using std::vector;
using std::pair;

namespace ASAPSPACE {

class ParallelAtoms;
class Communicator;

/// Abstract base class for a domain decomposition.

/// A DomainDecomposition is responsible for deciding which atoms
/// belong on which processors.  It also decides if a given processor
/// needs to know about a given atom.  Only one implementation exists
/// so far: RegularGridDecomposition.
class  DomainDecomposition
{
public:
  virtual ~DomainDecomposition() {};
  /// Create a list telling which processor should receive which atoms.

  /// Given a ParallelAtoms object, the DomainDecomposition decides
  /// which atoms should be sent to which processors.  The returned
  /// list processor is a list with one element for each processor,
  /// the element is itself a list containing the indices of all atom
  /// which should migrate to that processor.  Atoms which should stay
  /// on this processor are not included, i.e. the list corresponding
  /// to this processor is empty.
  ///
  /// The last argument is a list of all atoms migrating away,
  /// i.e. the union of all the above-mentioned lists.
  virtual void whichProcessor(ParallelAtoms *atoms,
                              vector< vector<int> > &processor,
                              vector<int> &migrateAway) = 0;
  /// Make the ghost export list (see ParallelPotential::ghosts)
  virtual void makeGhostExportLists(ParallelAtoms *atoms, double rCut,
                                    vector< vector< pair<int, int> > > &ghostLists) = 0;
  /// The list of processors to which data should be sent (in that order!)
  virtual const vector<int> *GetSendList() = 0;
  /// The list of processors from which data should be received (in that order).
  virtual const vector<int> *GetRecvList() = 0;

  virtual void GetTranslations(ParallelAtoms *atoms,
			       vector<Vec> &translations) const = 0;

};

} // end namespace
#endif //_DOMAINDECOMPOSITION_H

