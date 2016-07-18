// Langevin.h  --  The Velocity Verlet molecular dynamics algorithm.
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

#ifndef _LANGEVIN_H
#define _LANGEVIN_H

#include "AsapPython.h"
#include "Asap.h"
#include "AsapObject.h"
#include "Vec.h"
#include "MolecularDynamics.h"
#include <vector>

namespace ASAPSPACE {

class Potential;
class DynamicAtoms;
class AsapRandomThread;

class Langevin : public MolecularDynamics
{
public:
  Langevin(PyObject *py_atoms, Potential *calc, double timestep,
           PyObject *sdpos_name, PyObject *sdmom_name, PyObject *c1_name,
           PyObject *c2_name, bool fixcm, unsigned int seed);
  virtual ~Langevin();

  virtual string GetName() const {return "Langevin";}

  // Constants are scalars and given here.
  void SetScalarConstants(double act0, double c3, double c4, double pmcor,
                          double cnst);
  // Constants are vectors on the atoms.
  void SetVectorConstants(PyObject *act0_name, PyObject *c3_name,
                          PyObject *c4_name, PyObject *pmcor_name,
                          PyObject *cnst_name);

  // Get the random numbers
  void GetRandom(std::vector<Vec> &x1, std::vector<Vec> &x2, bool gaussian=true);

protected:
  // Run the dynamics.  nsteps is the number of steps, observers is a Python
  // list of observers, each observer described by a tuple
  // (function, interval, args, kwargs) and function(*args, **kwargs) is
  // called for each interval time steps.
  //virtual void Run(int nsteps);
  virtual void Run2(int nsteps, PyObject *observers, PyObject *self);
  void ClearPyNames();  // Release Python variable names.

private:
  bool fixcm;
  bool vectorconstants;  // Constants are vectors
  double act0;
  double c3;
  double c4;
  double pmcor;
  double cnst;
  PyObject *sdpos_name;
  PyObject *sdmom_name;
  PyObject *c1_name;
  PyObject *c2_name;
  PyObject *act0_name;
  PyObject *c3_name;
  PyObject *c4_name;
  PyObject *pmcor_name;
  PyObject *cnst_name;
  AsapRandomThread *random;
};

} // end namespace

#endif // _LANGEVIN_H
