// -*- C++ -*-
// MonteCarloAtoms.h:  Adds Monte Carlo optimization to the atoms interface.
//
// Copyright (C) 2009-2011 Jakob Schiotz and Center for Individual
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

#ifndef MONTECARLOATOMS_H
#define MONTECARLOATOMS_H
#include "Atoms.h"
#include <set>

namespace ASAPSPACE {

class MonteCarloAtoms : public NormalAtoms
{
public:
  MonteCarloAtoms() : NormalAtoms() {};

  bool GetMonteCarloRelevant() const {return mc_optim_relevant;}

  const set<int> &GetModifiedAtoms() const {return modified_atoms;}

  void MonteCarloEnd() const {};  // Pt. do nothing.
  
protected:
  virtual bool check_numbers(PyArrayObject *py_numbers, PyArrayObject *py_gh_num,
			     bool step_count_atoms);

  virtual bool check_positions(PyArrayObject *py_positions, PyArrayObject *py_gh_pos,
			       bool step_count_atoms_or_cell);
  
private:
  bool mc_optim_relevant;
  set<int> modified_atoms;
};

} // end namespace

#endif // MONTECARLOATOMS_H
