// -*- C++ -*-
// MonteCarloEMT.h  --  EMT with Monte Carlo optimizations.  
//
// Copyright (C) 2006-2011 Jakob Schiotz and Center for Individual
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


#ifndef _MONTECARLOEMT_H
#define _MONTECARLOEMT_H

#include "AsapPython.h"
#include "Asap.h"
#include "EMT.h"
#include <set>
#include <map>
using std::set;
using std::map;

namespace ASAPSPACE {

class MonteCarloAtoms;

class MonteCarloEMT : public EMT
{
public:
  /// Create an EMT potential optionally using a parameter provider.
  MonteCarloEMT(PyObject *self, PyObject *prov = 0);

  /// Delete the EMT potential.
  virtual ~MonteCarloEMT();

  /// Set the atoms, and enable Monte Carlo optimizations in the atoms.
  virtual void SetAtoms(PyObject *pyatoms, Atoms* accessobj = NULL);

  /// Calculate the energies
  ///
  /// In a Monte Carlo simulation, energies will (if possible) only be
  /// recalculated for atoms that have been moved since last energy
  /// calculation, or atoms where a neighbor has moved.
  virtual const vector<double> &GetPotentialEnergies(PyObject *pyatoms);

  /// Unlike EMT, this potential can NOT be used in parallel simulations
  virtual bool Parallelizable() const {return false;}

protected:
  virtual void CreateNeighborList();

  /// Handle alchemy in Monte Carlo simulations.
  void PartialCalculateIDs(const set<int> &changedatoms);

  /// Update sigmas in Monte Carlo simulations.
  void PartialCalculateSigmas(const set<int> &changedatoms);

  /// Calculate the energies of affected atoms in MC simulations.
  void PartialCalculateEnergiesAfterSigmas(const set<int> &changedatoms);

private:
  MonteCarloAtoms *mc_atoms;  // Secondary reference to the atoms.
  map<int, int> zmap;  ///< Map from atomic numbers to ids.
};

} // end namespace

#endif // _MONTECARLOEMT_H
