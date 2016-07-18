// EMTDefaultParameterProvider.cpp  --  provides the default EMT parameters.

// Copyright (C) 2008-2011 Jakob Schiotz and Center for Individual
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


#include "EMTDefaultParameterProvider.h"
#include "Exception.h"
#include <string.h>
#include <cmath>
#include <iostream>

using namespace std;

const int EMTDefaultParameterProvider::shell0 = 3;
const int EMTDefaultParameterProvider::shell1 = 4;
const double EMTDefaultParameterProvider::bohr = 0.5291772; // Angstrom

typedef vector<emt_parameters *>::iterator param_p_iter;

EMTDefaultParameterProvider::EMTDefaultParameterProvider()
{
  chi = 0;
  listcutofffactor = 1.04500185048;
}

EMTDefaultParameterProvider::~EMTDefaultParameterProvider()
{
  if (chi != 0) delete chi;
  for (param_p_iter ep = params.begin(); ep != params.end(); ++ep) 
  {
      emt_parameters *e = *ep;
      delete e;
  }
}

void EMTDefaultParameterProvider::Debug()
{
  cerr << "EMTDefaultParameterProvider debug information:" << endl;
  cerr << "Length of params vector: " << params.size() << endl;
  for (param_p_iter e = params.begin(); e != params.end(); ++e)
    cerr << "   " << (*e)->name << " Z=" << (*e)->Z << endl;
  if (chi == 0) {
    cerr << "Chi matrix: NOT ALLOCATED." << endl;
  } else {
    int n = params.size();
    cerr << "Chi matrix: " << n << " x " << n << endl;
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) 
        cerr << "    chi[" << i << "][" << j << "] = " << (*chi)[i][j] << endl;
  }
}

const emt_parameters *EMTDefaultParameterProvider::GetParameters(int element)
{
  emt_parameters *newpar;
  // Check if we already have parameters for this element.  If so:
  // return them.
  for (param_p_iter i=params.begin(); i != params.end(); ++i) {
    if ((*i)->Z == element)
      return *i;
  }

  // Get the parameters and add them to the vector.
  newpar = GetNewParameters(element);
  newpar->index = params.size();  // Count the elements, starting at zero.
  params.push_back(newpar);
  return newpar;
}

emt_parameters *EMTDefaultParameterProvider::GetNewParameters(int element)
{
  // This version returns the parameters using the units Angstrom
  // and eV.  To prevent typos the parameters are not converted in
  // the if statement below, instead the conversion from bohr to
  // Angstrom is done when the values are copied to the
  // emt_parameters struct.
	
  double E0, S0, n0, V0, eta2, kappa, lambda, mass;
  double ls;
  int Z;
  string name;
  emt_parameters *p;
	
  if(element == 13){
    name = "Al";
    E0=-3.280; S0=3.000; V0=1.493; Z=13; eta2=1.240; kappa=2.000;
    lambda=1.169; mass=26.98 ; n0=0.007 ; ls = 7.54871784;
  } 
  else if(element == 29){
    name = "Cu";
    E0=-3.510; S0=2.67 ; V0=2.476; Z=29; eta2=1.652; kappa=2.74;
    lambda=1.906; mass=63.54 ; n0=0.0091 ;
    ls = 6.789382809;
  }
  else if(element == 47){
    name = "Ag";
    E0=-2.96 ; S0=3.01 ; V0=2.132; Z=47; eta2=1.652; kappa=2.790;
    lambda=1.892; mass=107.87; n0=0.00547;
    ls = 7.6790043;
  }
  else if(element == 79){
    name = "Au";
    E0=-3.80 ; S0=3.00 ; V0=2.321; Z=79; eta2=1.674; kappa=2.873;
    lambda=2.182; mass=196.97; n0=0.00703; 
    ls = 7.66504117182;
  }
  else if(element == 28){
    name = "Ni";
    E0=-4.44 ; S0=2.60 ; V0=3.673; Z=28; eta2=1.669; kappa=2.757;
    lambda=1.948; mass=58.71 ; n0=0.0103 ;
    ls = 6.598896;
  }
  else if(element == 46){
    name = "Pd";
    E0=-3.90 ; S0=2.87 ; V0=2.773; Z=46; eta2=1.818; kappa=3.107;
    lambda=2.155; mass=106.4 ; n0=0.00688;
    ls = 7.330378;
  }
  else if(element == 78){
    name = "Pt";
    E0=-5.85 ; S0=2.90 ; V0=4.067; Z=78; eta2=1.812; kappa=3.145;
    lambda=2.192; mass=195.09; n0=0.00802;
    ls = 7.41119853;
  }
  else if(element == 12){
    name = "Mg";
	// these Mg parameters have been generated by N. Bailey using properties
	// of Mg derived from DFT, as well as properties of Mg-Cu alloys (Cu
	// parameters were also refit, but the defaults will be left as is)
	// The fitting process used eV and Angstrom, so the reverse of the
	// scaling below needs to be done here for parameters involving length
	// units. The fitting was done in December 2002
    E0=-1.487 ; S0=1.7664 / bohr ; V0=2.2298; Z=12; eta2=2.5411*bohr; kappa=4.435 * bohr;
    lambda=3.2927 * bohr; mass=24.3050; n0=0.03554 * (bohr*bohr*bohr);
    ls = 4.52004 / bohr;
  }
  else 
    {
		throw AsapError("This element isn't defined in EMT.");
    }

  p = new emt_parameters;
  p->e0 = E0;
  p->seq = S0 * bohr;
  p->neq = n0 / (bohr*bohr*bohr);
  p->V0 = V0;
  p->eta2 = eta2 / bohr;
  p->kappa = kappa / bohr;
  p->lambda = lambda / bohr;
  p->mass = mass;
  p->invmass = 1.0/mass;
  p->gamma1 = 0.0;        // These are calculated later!
  p->gamma2 = 0.0;
  p->Z = Z;
  assert(element == Z);
  p->name = name;
  p->lengthscale = ls / sqrt(2.0) * bohr;

  return p;
}

void EMTDefaultParameterProvider::CalcGammaEtc()
{
  // These can be overridden on an individual basis
  calc_cutoff();
  calc_gammas();
  calc_chi();
}

double EMTDefaultParameterProvider::GetMaxListCutoffDistance()  // Max value, useful before initialization.
{
  maxseq = 3.01 * bohr; // Value for Ag
  double maxcutoff = 0.5 * maxseq * Beta * (sqrt((double) shell0) +
                                            sqrt((double) shell0 + 1));
  return maxcutoff * listcutofffactor;
}

void EMTDefaultParameterProvider::calc_cutoff()
{
  double r;
	
  maxseq = 0.0;
  for (param_p_iter i = params.begin(); i != params.end(); ++i) 
    if ((*i)->seq > maxseq)
      maxseq = (*i)->seq;

  cutoff = 0.5 * maxseq * Beta * (sqrt((double) shell0) +
                                  sqrt((double) shell0 + 1));
  // The cutoff slope is set (rather arbitrarily) to the same value
  // used in the original ARTwork code.  The idea of this
  // expression is that the slope is such that the Fermi function
  // has the value 1e-5 where the neighbor list is cut off.
  r = cutoff * 2.0 * sqrt((double) shell0 + 1) / (sqrt((double) shell0) + 
						  sqrt((double) shell0 + 1));
  cutslope = log(9999.0) / (r - cutoff);
#if 0
  cerr << "cutoff " << cutoff << endl;
  cerr << "maxseq " << maxseq << endl;
  cerr << "Beta " << Beta << endl;
  cerr << "cutslope " << cutslope << endl;
  cerr << "shell0 " << shell0 << endl;
  cerr << "shell1 " << shell1 << endl;
#endif
}

void EMTDefaultParameterProvider::calc_gammas()
{
  static int shellpop[5] = {12, 6, 24, 12, 24}; // Population of the first 5 shells (fcc).
  double w, d;
	
  for (param_p_iter ep = params.begin(); ep != params.end(); ++ep) 
    {
      emt_parameters *e = *ep;
      
      e->gamma1 = e->gamma2 = 0.0;
      // The neighborlist object does not return neighbors in th 4. shell!
//    for (int is = 0; is < shell1 && is < 5; is++)
      for (int is = 0; is < shell1 - 1 && is < 5; is++) // !!!!! xxxxx
        {
          d = sqrt((double) (is + 1)) * Beta * e->seq;
          w = 1. / (1. + exp(cutslope * (d - cutoff)));
          e->gamma1 += w * shellpop[is] * exp(-d * e->eta2);
          e->gamma2 += w * shellpop[is] * exp(-d * e->kappa / Beta);
        }
      e->gamma1 /= shellpop[0] * exp(-Beta * e->seq * e->eta2);
      e->gamma2 /= shellpop[0] * exp(-e->seq * e->kappa);
    }
}

void EMTDefaultParameterProvider::calc_chi()
{
  int n = params.size();

  chi = new TinyDoubleMatrix(n,n);

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      (*chi)[i][j] = params[j]->neq / params[i]->neq;
    }
}

