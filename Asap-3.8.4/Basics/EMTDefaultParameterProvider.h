// EMTDefaultParameterProvider.h  --  provides the default EMT parameters.
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

#ifndef _EMTDEFAULTPARAMETERPROVIDER_H
#define _EMTDEFAULTPARAMETERPROVIDER_H

#include "EMTParameterProvider.h"
#include <vector>
using std::vector;

namespace ASAPSPACE {

class EMTDefaultParameterProvider : public EMTParameterProvider
{
public:
  EMTDefaultParameterProvider();
  virtual ~EMTDefaultParameterProvider();
  virtual string GetName() const {return "EMTDefaultParameterProvider";}
  virtual const emt_parameters *GetParameters(int element);
  virtual void CalcGammaEtc();
  // Can the cutoff be made element-dependent?
  virtual double GetCutoffDistance() {return cutoff;}
  virtual double GetListCutoffDistance() {return cutoff * listcutofffactor;}
  virtual double GetMaxListCutoffDistance();  // Max value, useful before initialization.
  virtual double GetCutoffSlope() {return cutslope;}
  // The potential delegates GetLengthScale to the provider.
  // What should be done for multiple elements?
  virtual double GetLengthScale() {return params[0]->lengthscale;}
  virtual int GetNumberOfElements() {return params.size();}
  virtual const TinyDoubleMatrix *GetChi() {return chi;}
  virtual void Debug();
	
protected:
  virtual emt_parameters *GetNewParameters(int element);
  virtual void calc_cutoff();
  virtual void calc_gammas();
  virtual void calc_chi();
	
  std::vector<emt_parameters *> params;
  TinyDoubleMatrix *chi;
  double maxseq;
  double cutoff;
  double cutslope;
  double listcutofffactor;

  // Two constants: the last shell included in interaction range, and in nb list
  static const int shell0;
  static const int shell1;

  // The length unit original used by EMT
  static const double bohr;
};

} // end namespace

#endif // ! _EMTDEFAULTPARAMETERPROVIDER_H
