// VelocityVerlet.h  --  The Velocity Verlet molecular dynamics algorithm.
// -*- c++ -*-
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

#ifndef _VELOCITYVERLET_H
#define _VELOCITYVERLET_H

#include "AsapPython.h"
#include "Asap.h"
#include "AsapObject.h"
#include "Vec.h"
#include "MolecularDynamics.h"
#include <vector>

namespace ASAPSPACE {

class Potential;
class DynamicAtoms;

class VelocityVerlet : public MolecularDynamics
{
public:
  VelocityVerlet(PyObject *py_atoms, Potential *calc, double timestep);
  ~VelocityVerlet();

  virtual string GetName() const {return "VelocityVerlet";}

protected:
  // Run the dynamics.  nsteps is the number of steps, observers is a Python
  // list of observers, each observer described by a tuple
  // (function, interval, args, kwargs) and function(*args, **kwargs) is
  // called for each interval time steps.
  //virtual void Run(int nsteps);
  virtual void Run2(int nsteps, PyObject *observers, PyObject *self);

private:
  PyObject *fixatoms_name;
};

} // end namespace

#endif // _VELOCITYVERLET_H
