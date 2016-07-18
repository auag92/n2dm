// VelocityVerlet.cpp  --  The Velocity Verlet molecular dynamics algorithm.
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

#include "VelocityVerlet.h"
#include "Potential.h"
#include "DynamicAtoms.h"
//#define ASAPDEBUG
#include "Debug.h"

#include <iostream>

VelocityVerlet::VelocityVerlet(PyObject *py_atoms, Potential *calc, double timestep) :
  MolecularDynamics(py_atoms, calc, timestep)
{
  fixatoms_name = PyString_FromString("FixAtoms_mult_double");
  assert(fixatoms_name != NULL);  // Probably cannot fail.
}

VelocityVerlet::~VelocityVerlet()
{
  Py_DECREF(fixatoms_name);
}

void VelocityVerlet::Run2(int nsteps, PyObject *observers, PyObject *self)
{
  //std::cerr << "Thread " << thread << ": steps = " << nsteps << std::endl;
  DEBUGPRINT;

  ParseObservers(observers);
  const vector<double> &invmasses = atoms->GetInverseMasses();
  DEBUGPRINT;
  const Vec *F = GetForces();
  Vec *r = atoms->GetPositions();
  Vec *p = atoms->GetMomenta();
  const asap_z_int *Z = atoms->GetAtomicNumbers();
  double *mult = atoms->GetDoubleDataMaybe(fixatoms_name);
  double halftimestep = 0.5 * timestep;
  DEBUGPRINT;
  for (int n = 0; n < nsteps; n++)  // nsteps-1 steps
    {
      CHECKNOASAPERROR; // Unnecessary paranoia.
      DEBUGPRINT;
      if (mult == NULL)
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
          for (int i = 0; i < nAtoms; i++)
            {
              p[i] += halftimestep * F[i];
              r[i] += timestep * invmasses[Z[i]] * p[i];  // Position update step n.
            }
        }
      else
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
          for (int i = 0; i < nAtoms; i++)
            {
              p[i] = mult[i] * (p[i] + halftimestep * F[i]);
              r[i] += timestep * invmasses[Z[i]] * p[i];  // Position update step n.
            }
        }
      F = GetForces();
      r = atoms->GetPositions();  // May have changed in parallel simulation
      p = atoms->GetMomenta();
      Z = atoms->GetAtomicNumbers();
      mult = atoms->GetDoubleDataMaybe(fixatoms_name);
      if (mult == NULL)
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
          for (int i = 0; i < nAtoms; i++)
            {
              p[i] += halftimestep * F[i];
            }
        }
      else
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
          for (int i = 0; i < nAtoms; i++)
            {
              p[i] += halftimestep * F[i] * mult[i];
            }
        }
      steps++;
      CallObservers(self);
    }
  CleanupObservers();
  UpdateStepsInPython(self);
}
