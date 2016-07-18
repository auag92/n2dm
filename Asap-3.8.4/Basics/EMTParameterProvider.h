// EMTParameterProvider.h  --  Abstract class for objects providing EMT parameters.
//
// Copyright (C) 2001-2011 Jakob Schiotz and Center for Individual
// Nanoparticle Functionality, Department of Physics, Technical
// University of Denmark.  Email: schiotz@fysik.dtu.dk
//
// This file is part of Asap version 3.
// Asap is released under the GNU Lesser Public License (LGPL) version 3.
// However, the parts of Asap distributed within the OpenKIM project
// (including this file) are also released under the Common Development
// and Distribution License (CDDL) version 1.0.
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

#ifndef _EMTPARAMETERPROVIDER_H
#define _EMTPARAMETERPROVIDER_H

#include "AsapPython.h"
#include "Asap.h"
#include "AsapObject.h"
#include "TinyMatrix.h"
#include <string>
using std::string;

namespace ASAPSPACE {

// The stucture used to for store EMT parameters.
struct emt_parameters
{
  double e0, seq, neq, V0, eta2, kappa, lambda, mass, invmass;
  double gamma1, gamma2;
  // double pairA, pairD;
  double lengthscale;  // Not an EMT parameter.  1/sqrt(2) of the lattice constant.
  int Z;
  std::string name;
  int index;   // An index into various arrays.  Counts the elements in the simulation.
};

// A geometric constant.
// static const double Beta = 1.80939979; // ((16*pi/3)^(1/3))/sqrt(2)
static const double Beta = 1.809; // Preserve same rounding as in ARTwork.


/// An EMTParameterProvider provides the EMT parameters to the potential.
///
/// The GetParameters() method returns the parameters for a
/// given element.  If it is called multiple times with the same
/// element it must return a pointer to the same emt_parameters struct.
/// The Provider owns the parameters, they are deleted when the
/// provider is deleted.  After getting _all_ the elements,
/// CalcGammaEtc() should be called to calculate the cutoff radius, the
/// gammas and other quantities that cannot be calculated before it is
/// known which elements are in the simulation.  The gammas and the
/// cutoff are determined by the Provider since changing the way they
/// are determined corresponds to changing the EMT potential.
class EMTParameterProvider : public AsapObject
{
public:
  virtual ~EMTParameterProvider() {};
  virtual const emt_parameters *GetParameters(int element) = 0;
  virtual void CalcGammaEtc() = 0;
  virtual double GetCutoffDistance() = 0; // Can it be made element-dependent?
  virtual double GetCutoffSlope() = 0;
  virtual double GetListCutoffDistance() = 0;   // Cutoff for neighbor list.
  virtual double GetMaxListCutoffDistance() = 0;  // Max value, useful before initialization.
  virtual double GetLengthScale() = 0;    // The potential delegates this to the provider.
  virtual const TinyDoubleMatrix *GetChi() = 0;
  virtual int GetNumberOfElements() = 0;
};


/// The Python object corresponding to an EMTParameterProvider object.
typedef struct {
  PyObject_HEAD
  EMTParameterProvider *cobj;
  PyObject *weakrefs;
} PyAsap_EMTParamProvObject;

} // end namespace


#endif // ! _EMTPARAMETERPROVIDER_H
