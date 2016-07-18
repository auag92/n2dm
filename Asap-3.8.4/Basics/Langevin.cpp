// Langevin.cpp  --  The Velocity Verlet molecular dynamics algorithm.
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

#include "Langevin.h"
#include "Potential.h"
#include "DynamicAtoms.h"
#include "RandomNumbers.h"
//#define ASAPDEBUG
#include "Debug.h"

#include <iostream>
#include <cmath>

Langevin::Langevin(PyObject *py_atoms, Potential *calc, double timestep,
                   PyObject *sdpos_name, PyObject *sdmom_name, PyObject *c1_name,
                   PyObject *c2_name, bool fixcm, unsigned int seed) :
    MolecularDynamics(py_atoms, calc, timestep)
{
  vectorconstants = false;
  npy_uint32 seed32 = seed;
  random = new AsapRandomThread(seed32);
  this->sdpos_name = sdpos_name;
  this->sdmom_name = sdmom_name;
  this->c1_name = c1_name;
  this->c2_name = c2_name;
  this->fixcm = fixcm;
  Py_INCREF(this->sdpos_name);
  Py_INCREF(this->sdmom_name);
  Py_INCREF(this->c1_name);
  Py_INCREF(this->c2_name);
}

Langevin::~Langevin()
{
  Py_DECREF(sdpos_name);
  Py_DECREF(sdmom_name);
  Py_DECREF(c1_name);
  Py_DECREF(c2_name);
  if (vectorconstants)
    ClearPyNames();
  delete random;
}

void Langevin::SetScalarConstants(double act0, double c3, double c4, double pmcor,
                                  double cnst)
{
  if (vectorconstants)
    ClearPyNames();
  vectorconstants = false;
  this->act0 = act0;
  this->c3 = c3;
  this->c4 = c4;
  this->pmcor = pmcor;
  this->cnst = cnst;
}

void Langevin::SetVectorConstants(PyObject *act0_name, PyObject *c3_name,
                                  PyObject *c4_name, PyObject *pmcor_name,
                                  PyObject *cnst_name)
{
  if (vectorconstants)
    ClearPyNames();
  vectorconstants = true;
  this->act0_name = act0_name;
  this->c3_name = c3_name;
  this->c4_name = c4_name;
  this->pmcor_name = pmcor_name;
  this->cnst_name = cnst_name;
  Py_INCREF(this->act0_name);
  Py_INCREF(this->c3_name);
  Py_INCREF(this->c4_name);
  Py_INCREF(this->pmcor_name);
  Py_INCREF(this->cnst_name);
}

void Langevin::ClearPyNames()
{
  Py_DECREF(act0_name);
  Py_DECREF(c3_name);
  Py_DECREF(c4_name);
  Py_DECREF(pmcor_name);
  Py_DECREF(cnst_name);
}

void Langevin::Run2(int nsteps, PyObject *observers, PyObject *self)
{
  DEBUGPRINT;

  ParseObservers(observers);
  DEBUGPRINT;
  const Vec *F = GetForces();
  Vec *r = atoms->GetPositions();
  Vec *p = atoms->GetMomenta();
  vector<Vec> rnd1;
  vector<Vec> rnd2;
  vector<Vec> rrnd;
  vector<Vec> prnd;
  rnd1.reserve(nAtoms + nAtoms/20);
  rnd2.reserve(nAtoms + nAtoms/20);
  rrnd.reserve(nAtoms + nAtoms/20);
  prnd.reserve(nAtoms + nAtoms/20);


  for (int n = 0; n < nsteps; n++)  // nsteps-1 steps
    {
      CHECKNOASAPERROR; // Unnecessary paranoia.
      // A buffer for the random numbers
      double *sdpos = atoms->GetDoubleData(sdpos_name);
      double *sdmom = atoms->GetDoubleData(sdmom_name);
      double *c1 = atoms->GetDoubleData(c1_name);
      double *c2 = atoms->GetDoubleData(c2_name);
      double *pmcor_v = NULL;
      double *cnst_v = NULL;
      double *act0_v = NULL;
      double *c3_v = NULL;
      if (vectorconstants)
        {
          pmcor_v = atoms->GetDoubleData(pmcor_name);
          cnst_v = atoms->GetDoubleData(cnst_name);
          act0_v = atoms->GetDoubleData(act0_name);
          c3_v = atoms->GetDoubleData(c3_name);
        }
      GetRandom(rnd1, rnd2);
      rrnd.resize(nAtoms);
      prnd.resize(nAtoms);
      if (vectorconstants)
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
          for (int i = 0; i < nAtoms; i++)
            {
              rrnd[i] = sdpos[i] * rnd1[i];
              prnd[i] = sdmom[i] * pmcor_v[i] * rnd1[i] +
                  sdmom[i] * cnst_v[i] * rnd2[i];
            }
        }
      else
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
          for (int i = 0; i < nAtoms; i++)
            {
              rrnd[i] = sdpos[i] * rnd1[i];
              prnd[i] = sdmom[i] * pmcor * rnd1[i] +
                  sdmom[i] * cnst * rnd2[i];
            }
        }
      if (fixcm)
        {
          double a,b,c,d,e,f;
          a = b = c = d = e = f = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : a,b,c,d,e,f)
#endif // _OPENMP
          for (int i = 0; i < nAtoms; i++)
            {
              // This sum is done element by element using tmp variables to allow OpenMP parallelization
              a += rrnd[i][0];
              b += rrnd[i][1];
              c += rrnd[i][2];
              d += prnd[i][0];
              e += prnd[i][1];
              f += prnd[i][2];
            }
          Vec sumrrnd(a,b,c);
          Vec sumprnd(d,e,f);
          sumrrnd *= 1.0/nAtoms;
          sumprnd *= 1.0/nAtoms;
          double factor = sqrt(nAtoms / (nAtoms - 1.0));
#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
          for (int i = 0; i < nAtoms; i++)
            {
              rrnd[i] = (rrnd[i] - sumrrnd) * factor;
              prnd[i] = (prnd[i] - sumprnd) * factor;
            }
        }

      // Now the positions and momenta are updated
      if (vectorconstants)
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
          for (int i = 0; i < nAtoms; i++)
            {
              r[i] += c1[i] * p[i] + c2[i] * F[i] + rrnd[i];
              p[i] = p[i] * act0_v[i] + c3_v[i] * F[i] + prnd[i];
            }
        }
      else
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
          for (int i = 0; i < nAtoms; i++)
            {
              r[i] += c1[i] * p[i] + c2[i] * F[i] + rrnd[i];
              p[i] = p[i] * act0 + c3 * F[i] + prnd[i];
            }
        }
      // Get forces, this may migrate atoms.
      F = GetForces();
      r = atoms->GetPositions();
      p = atoms->GetMomenta();
      if (vectorconstants)
        {
          double *c4_v = atoms->GetDoubleData(c4_name);
#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
          for (int i = 0; i < nAtoms; i++)
            p[i] += c4_v[i] * F[i];
        }
      else
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
          for (int i = 0; i < nAtoms; i++)
            p[i] += c4 * F[i];
        }

      steps++;
      CallObservers(self);
    }
  CleanupObservers();
  UpdateStepsInPython(self);
}

void Langevin::GetRandom(std::vector<Vec> &x1, std::vector<Vec> &x2, bool gaussian /*=true*/)
{
  x1.resize(nAtoms);
  x2.resize(nAtoms);
  random->RandomDoubles((double *) &x1[0], 3*nAtoms);
  random->RandomDoubles((double *) &x2[0], 3*nAtoms);
  if (!gaussian)
    return;
  // Make them Gaussian
  double u1, u2, f;
#ifdef _OPENMP
#pragma omp parallel for private(u1, u2, f)
#endif // _OPENMP
  for (int i = 0; i < nAtoms; i++)
    for (int j = 0; j < 3; j++)
      {
        u1 = x1[i][j];
	u2 = 2 * M_PI * x2[i][j];
	f = sqrt(-2.0 * log(1.0 - u1));
        x1[i][j] = f * sin(u2);
        x2[i][j] = f * cos(u2);
      }
}
