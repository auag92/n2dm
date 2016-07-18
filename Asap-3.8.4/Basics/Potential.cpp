// -*- C++ -*-
// Potential.cpp: Common functionality for all potentials.
//
// Copyright (C) 2001-2012 Jakob Schiotz and Center for Individual
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

#include "AsapPython.h"
#include "Asap.h"
#include "Atoms.h"
#include "Potential.h"
#include "Debug.h"

// Standard mapping of the six independent parts of the stress tensor to
// vector notation
const static int stresscomp[3][3] = {{0, 5, 4}, {5, 1, 3}, {4, 3, 2}};

void Potential::RecoverAfterException()
{
  DEBUGPRINT;
  if (atoms != NULL && atoms->IsActive())
    atoms->End();
}

bool Potential::CalcReq_Stress(PyObject *pyatoms)
{
  atoms->Begin(pyatoms, true);  // Expect reopen.
  bool result = (CalcReq_Virials(pyatoms) ||
                 (stresses_cnt != atoms->GetPositionsCounter()) ||
                 (stresses_mom_cnt != atoms->GetMomentaCounter()));
  atoms->End();
  return result;
}

const vector<SymTensor> &Potential::GetStresses(PyObject *a)
{
  DEBUGPRINT;
  assert(atoms != NULL);
  atoms->Begin(a, true);  // Expect reopen.
  int cnt = atoms->GetPositionsCounter();
  int momcnt = atoms->GetMomentaCounter();
  if ((cnt != stresses_cnt) || (momcnt != stresses_mom_cnt))
    {
      DEBUGPRINT;
      stress.clear();
      const vector<SymTensor> &virials = GetVirials(a);
      vector<double> vols;
      GetAtomicVolumes(vols);
      double invvol = atoms->GetTotalNumberOfAtoms() / atoms->GetVolume();
      const Vec *momenta = atoms->GetMomenta();
      const double *masses = atoms->GetMasses();
      stresses.resize(virials.size());
      for (int i = 0; i < stresses.size(); i++)
        {
          if (vols.size())
            invvol = 1.0 / vols[i];
          double invmass = 1.0 / masses[i];
          for (int alpha = 0; alpha < 3; alpha++)
            for (int beta = alpha; beta < 3; beta++)
              {
                int j = stresscomp[alpha][beta];
                double s = virials[i][j];
                if (momenta != NULL)
                  s -= momenta[i][alpha] * momenta[i][beta] * invmass;
                stress[j] += s;
                stresses[i][j] = s *invvol;
              }
        }
      stress *= 1.0 / atoms->GetVolume();
      stresses_cnt = cnt;
      stresses_mom_cnt = momcnt;
    }
  atoms->End();
  DEBUGPRINT;
  return stresses;
}

SymTensor Potential::GetVirial(PyObject *a)
{
  DEBUGPRINT;
  SymTensor result;
  for (int i = 0; i < 6; i++)
    result[i] = 0;
  const vector<SymTensor> &virials = GetVirials(a);
  for (int i = 0; i < virials.size(); i++)
    result += virials[i];
  DEBUGPRINT;
  return result;
}

SymTensor Potential::GetStress(PyObject *a)
{
  DEBUGPRINT;
  (void) GetStresses(a);
  DEBUGPRINT;
  return stress;
}

#ifdef ASAP_FOR_KIM
// When building as an OpenKIM model, there is no Python layer to pass through.
void Potential::SetAtoms_ThroughPython(PyObject *pyatoms, Atoms* accessobj /* = NULL */)
{
  SetAtoms(pyatoms, accessobj);
}
#else
// The normal version
void Potential::SetAtoms_ThroughPython(PyObject *pyatoms, Atoms* accessobj /* = NULL */)
{
  DEBUGPRINT;
  // Call self.set_atoms(...) as implemented in Python
  PyObject *py_accessobj;
  if (accessobj == NULL)
    {
      py_accessobj = Py_None;
      Py_INCREF(Py_None);
    }
  else
    py_accessobj = PyCObject_FromVoidPtr(accessobj, NULL);
  if (py_accessobj == NULL)
    throw AsapPythonError();
  PyObject *method = PyString_FromString("set_atoms");
  if (method == NULL)
    throw AsapPythonError();
  PyObject *result = PyObject_CallMethodObjArgs(self, method, pyatoms, py_accessobj, NULL);
  bool error = (result == NULL);
  Py_XDECREF(result);
  Py_DECREF(method);
  Py_DECREF(py_accessobj);
  if (error)
    throw AsapPythonError();
  DEBUGPRINT;
}
#endif
