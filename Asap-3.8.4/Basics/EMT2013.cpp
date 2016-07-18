/*
 * EMT2013.cpp
 *
 *  Created on: Apr 29, 2011
 *      Author: s072162
 */

#include "EMT2013.h"
#include "Asap.h"
#include "EMTParameterProvider.h"
#include "Atoms.h"
#include "Vec.h"
#include "NeighborCellLocator.h"
#include "NeighborList2013.h"
#include "Exception.h"
#include "Timing.h"
#include "mass.h"
#include "Debug.h"
#include <math.h>
#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
using std::vector;
using std::set;
using std::cerr;
using std::endl;
using std::flush;
using std::sort;
using std::cout;

#undef INTERLEAVETHREADS

// Choose the neighbor list implementation.
// 2: NeighborCellLocator:  Simple, but cost almost a factor 2 in performance.
// 3: NeighborList:  New 'standard' neighborlist implementation in Asap-3.
#define NBLIST 3


#if 1
#define VERB(x) if (verbose == 1) cerr << x
#else
#define VERB(x)
#endif

using namespace std;

// We should ALWAYS use beta, not Beta
#define Beta ERROR_BETA_IS_NOT_ALLOWED

const float EMT2013::beta = 1.809399790563555;  // ((16*pi/3)^(1./3.))/sqrt(2)

EMT2013::EMT2013(PyObject *self, PyObject *params, bool no_new_elements) : EMT(self, NULL)
{
  if (!PyDict_Check(params))
    throw AsapError("EMT2013 parameters must be a dictionary.");
  param_obj = params;
  this->no_new_elements = no_new_elements;
  Py_INCREF(param_obj);
  chi = NULL;
}

EMT2013::~EMT2013()
{
  if (chi != NULL)
    delete chi;
  Py_DECREF(param_obj);
}

void EMT2013::InitParameters()
{
  DEBUGPRINT;
  // Extract the elements from the list of atoms, and set up an array
  // of EMT parameters.
  set<int> elements_set;
  vector<int> elements;

  // Unlike the original EMT we are able to initialize all elements
  // right away.  This may give a performance hit, and we should
  // perhaps add the posibility to only initialize existing elements
  // (but then a sanity check needs to be added).  Or perhaps the user
  // should then only initialize with the relevant parameters.
  if (no_new_elements)
    {
      // Extract the elements from the list of atoms.
      atoms->GetListOfElements(elements_set);
    }
  else
    {
      // Extract the elements from the parameter dictionary.
      GetListOfElements(elements_set);
    }
  for (set<int>::iterator i = elements_set.begin(); i != elements_set.end();
       ++i)
    elements.push_back(*i);
  nelements = elements.size();
  assert(nelements == elements_set.size());
  sort(elements.begin(), elements.end());

  // Get the EMT parameters
  parameters.clear();
  for (vector<int>::iterator i = elements.begin(); i != elements.end(); ++i)
    parameters.push_back(ExtractParameters(*i));

  /* DELETE */
  //assert(nelements == provider->GetNumberOfElements());
  /* DELETE */

  /* Calculate the quantities that can only be calculated, once the
     elements in the simulations are known. */

  // Cutoff distances in the system and maximum cutoff distance are calculated
  CalculateCutoffDistances();

  /* DELETE */
  // UNUSED PARAMETER IS INITIALIZED
  rFermi = cutoffslope = 0;   // Unused
  /* DELETE */

  // The relevant Chi-values are calculated
  CalculateChi();

  /* DELETE */
  //if (verbose)
  //    cerr << "EMT::InitParameters:  rFermi = " << rFermi << "  rNbCut = "
  //         << rNbCut << "  cutoffslope = " << cutoffslope << endl;
  /* DELETE */
  DEBUGPRINT;
}

void EMT2013::GetListOfElements(set<int> &elements)
{
  PyObject *key;
  PyObject *value;
  Py_ssize_t pos = 0;
  while (PyDict_Next(param_obj, &pos, &key, &value)) {
      long z = PyInt_AsLong(key);
      assert(z != -1);  // Either an error occurred or an element has Z=-1
      elements.insert((int) z);
  }
}

emt_parameters *EMT2013::ExtractParameters(int z)
{
  // Initialization
  PyObject *key = PyInt_FromLong(z);
  PyObject *param_this_element = PyDict_GetItem(param_obj, key);
  Py_DECREF(key);

  // Error handling
  if (param_this_element == NULL)
    throw AsapError("No EMT2013 parameters for element ") << z;
  if (!PyDict_Check(param_this_element))
    throw AsapError("EMT2013 parameters must be a dictionary. Z = ") << z;

  // Data extraction
  emt_parameters *result = new emt_parameters;
  result->eta2 = ExtractParam_double(param_this_element, "eta2");
  result->lambda = ExtractParam_double(param_this_element, "lambda");
  result->kappa = ExtractParam_double(param_this_element, "kappa");
  result->e0 = ExtractParam_double(param_this_element, "E0");
  result->V0 = ExtractParam_double(param_this_element, "V0");
  result->seq = ExtractParam_double(param_this_element, "S0");
  result->neq = ExtractParam_double(param_this_element, "n0");
  result->mass = ExtractParam_double(param_this_element, "mass");
  result->invmass = 1. / result->mass;
  result->Z = z;

  // Data calculations
  // Parameters
  double eta2 = result->eta2;
  double seq = result->seq;
  double kappa = result->kappa;
  double rcut = 0.5 * (sqrt(1.5) + sqrt(2.)) * sqrt(2)
                                   * beta * seq;
  double s1r_cut = exp(- eta2 * ( rcut - beta * seq) );
  double ds1dr_cut = -eta2 * s1r_cut;
  double s2r_cut = exp(- kappa / beta * ( rcut - beta * seq) );
  double ds2dr_cut = -kappa / beta * s2r_cut;

  // Gamma's are calculated and saved in the parameter arrays
  result->gamma1 = 12 * (1 - (ds1dr_cut * (beta * seq - rcut) + s1r_cut ) ) +
                    6 * (exp( eta2 * beta * seq * (1 - sqrt(2)) ) - (ds1dr_cut * (sqrt(2) * beta * seq - rcut) + s1r_cut) )  +
                   24 * (exp( eta2 * beta * seq * (1 - sqrt(3)) ) - (ds1dr_cut * (sqrt(3) * beta * seq - rcut) + s1r_cut) );

  result->gamma2 = 12 * (1 - (ds2dr_cut * (beta * seq - rcut) + s2r_cut) ) +
                    6 * (exp( kappa * seq * (1 - sqrt(2)) ) - (ds2dr_cut * (sqrt(2) * beta * seq - rcut) + s2r_cut) )  +
                   24 * (exp( kappa * seq * (1 - sqrt(3)) ) - (ds2dr_cut * (sqrt(3) * beta * seq - rcut) + s2r_cut) );

  return result;
}

double EMT2013::ExtractParam_double(PyObject *params_for_z, const char *param_name)
{
  if (!PyDict_Check(params_for_z))
    throw AsapError("EMT2013 parameters must be a dictionary");

  PyObject *Py_value = PyDict_GetItemString(params_for_z, param_name);
  if (Py_value == NULL)
    throw AsapError("EMT2013 parameter missing: ") << param_name;

  if (!PyFloat_Check(Py_value))
    throw AsapError("The EMT2013 parameter must be a double.");

  return PyFloat_AsDouble(Py_value);
};

void EMT2013::CalculateCutoffDistances()
{
  // Maximum atom number
  int zmax = 0;

  // cutoff matrix for use by EMT object initialized
  rcut2.Allocate(nelements,nelements);

  // cutoff matrix for use by Neighborlist initialized
  for (int i = 0; i<nelements; i++)
    {
      if (zmax < parameters[i]->Z)
        {
          zmax = parameters[i]->Z;
        }
    }
  rcut2_NB.Allocate(zmax+1,zmax+1);
  for (int i = 0; i<zmax+1; i++)
    {
    for (int j = 0; j<zmax+1; j++)
      rcut2_NB[i][j] = 0.0;
    }

  // rcut2 and rcut2_NB filled
  rNbCut = 0.0;
  for (int i = 0; i<nelements ;i++)
    {
    for (int j = 0; j<nelements ;j++)
      {
        double TemprCut = 0.5 * (sqrt(1.5) + sqrt(2.)) * sqrt(2.)
                               * beta * max(parameters[i]->seq,parameters[j]->seq);
        if (TemprCut > rNbCut)
          rNbCut = TemprCut;

        rcut2[i][j] = TemprCut * TemprCut;
        rcut2_NB[parameters[i]->Z][parameters[j]->Z] = TemprCut * TemprCut;
      }
    }
  return;
}

void EMT2013::CalculateChi()
{
  TinyDoubleMatrix *chi = new TinyDoubleMatrix(nelements,nelements);

  for (int i = 0; i<nelements; i++)
    for (int j = 0; j<nelements; j++)
        (*chi)[i][j] = parameters[j]->neq / parameters[i]->neq;
  this->chi = chi;
}

void EMT2013::sigma_batch(int *self, int *other, Vec rnb[],
                      double *sq_dist, int zs, int zo, int n,
                      bool calculatesigma2, bool partialupdate /* = false */)
{
  USETIMER("EMT::sigma_batch");
  double *temporary = new double[4*BUFLEN];
  double *dsigma1s = &temporary[0];
  double *dsigma2s = dsigma1s + BUFLEN;
  double *ds1o = dsigma2s + BUFLEN;
  double *ds2o = ds1o + BUFLEN;
  assert(ds2o - &temporary[0] == (4 - 1) * BUFLEN);
  double *dsigma1o = NULL;
  double *dsigma2o = NULL;

  // initialization of lsf parameters;
  /* TEMP */ double rcut, ds1sdrcut, ds2sdrcut, s1srcut, s2srcut,
                    ds1odrcut, ds2odrcut, s1orcut, s2orcut; /* TEMP */

  double other_eta2betaseq, self_eta2betaseq;
  double other_kappaoverbeta, self_kappaoverbeta;
  double other_kappaseq, self_kappaseq;
  double *s1s, *s1o, *s2s, *s2o;
  const emt_parameters *emtself, *emtother;

  assert(n <= BUFLEN);

  /* Get EMT parameters !!! REWRITE !!! XXXX */
  emtself = parameters[zs];
  emtother = parameters[zo];
  /* DELETE  cutslopecutdist = cutoffslope * rFermi;   DELETE */
  other_eta2betaseq = emtother->eta2 * beta * emtother->seq;
  self_eta2betaseq = emtself->eta2 * beta * emtself->seq;
  other_kappaoverbeta = emtother->kappa / beta;
  self_kappaoverbeta = emtself->kappa / beta;
  other_kappaseq = emtother->kappa * emtother->seq;
  self_kappaseq = emtself->kappa * emtself->seq;
  double other_eta2 = emtother->eta2;
  double self_eta2 = emtself->eta2;

  rcut = sqrt(rcut2[zs][zo]);
  s1srcut = exp(-other_eta2 * rcut + other_eta2betaseq);
  ds1sdrcut = -other_eta2 * s1srcut;
  if (zs != zo)
  {
	  s1orcut = exp(-self_eta2 * rcut + self_eta2betaseq);
	  ds1odrcut = -self_eta2 * s1orcut;
  }
  /* if sigma2 is to be calculated the parameters needed for the cutoff are initialized */
  if (calculatesigma2)
  	  {
	  	  s2srcut = exp(-other_kappaoverbeta * rcut + other_kappaseq);
	  	  ds2sdrcut = -other_kappaoverbeta * s2srcut;
	  	  if (zs != zo)
	  	  {
	  		  s2orcut = exp(-self_kappaoverbeta * rcut + self_kappaseq);
	  		  ds2odrcut = -self_kappaoverbeta * s2orcut;
	  	  }
 	  }

  assert(n <= BUFLEN);
  if ((zs == zo) && !calculatesigma2)
    {
      for (int i = 0; i < n; i++)
	{
	  /* dist = sqrt(sq_dist),  distances between atoms */
	  double dist = sqrt(sq_dist[i]);
	  /* Calculate Linear Subtraction Function (New cutoff function) */
	  double lsfs1 = ds1sdrcut * (dist - rcut) + s1srcut;

	  /* Contribution to sigma1 */
	  dsigma1s[i] = exp(-other_eta2 * dist + other_eta2betaseq) - lsfs1;
	}
      dsigma1o = dsigma1s;
    }
  else if ((zs != zo) && !calculatesigma2)
    {
      for (int i = 0; i < n; i++)
	{
      /* dist = sqrt(sq_dist),  distances between atoms */
      double dist = sqrt(sq_dist[i]);
      double distmrcut = dist - rcut;
      /* Calculate Linear Subtraction Functions (NEW cutoff function NEW) */
      double lsfs1 = ds1sdrcut * distmrcut + s1srcut;
      double lsfo1 = ds1odrcut * distmrcut + s1orcut;
      /* Contribution to sigma1 */
      dsigma1s[i] = exp(-other_eta2 * dist + other_eta2betaseq) - lsfs1;
	  ds1o[i] = exp(-self_eta2 * dist + self_eta2betaseq) - lsfo1;
	}
      dsigma1o = ds1o;
    }
  else if ((zs == zo) && calculatesigma2)
    {

      for (int i = 0; i < n; i++)
	{
   	  /* dist = sqrt(sq_dist),  distances between atoms */
   	  double dist = sqrt(sq_dist[i]);

  	  double distmrcut = dist - rcut;
   	  /* Calculate Linear Subtraction Functions (NEW cutoff function NEW) */
   	  double lsfs1 = ds1sdrcut * distmrcut + s1srcut;
   	  double lsfs2 = ds2sdrcut * distmrcut + s2srcut;

   	  /* Contribution to sigma1 */
   	  dsigma1s[i] = exp(-other_eta2 * dist + other_eta2betaseq) - lsfs1;
	  dsigma2s[i] = exp(-other_kappaoverbeta * dist + other_kappaseq) - lsfs2;

	}
      dsigma1o = dsigma1s;
      dsigma2o = dsigma2s;
    }
  else
    {
      // zs != zo && calculatesigma2
      for (int i = 0; i < n; i++)
	{
  	  /* dist = sqrt(sq_dist),  distances between atoms */
   	  double dist = sqrt(sq_dist[i]);
   	  double distmrcut = dist - rcut;
   	  /* Calculate Linear Subtraction Functions (NEW cutoff function NEW) */
   	  double lsfs1 = ds1sdrcut * distmrcut + s1srcut;
   	  double lsfs2 = ds2sdrcut * distmrcut + s2srcut;
   	  double lsfo1 = ds1odrcut * distmrcut + s1orcut;
   	  double lsfo2 = ds2odrcut * distmrcut + s2orcut;
   	  /* Contribution to sigma1 */
   	  dsigma1s[i] = exp(-other_eta2 * dist + other_eta2betaseq) - lsfs1;
	  ds1o[i] = exp(-self_eta2 * dist + self_eta2betaseq) - lsfo1;
	  dsigma2s[i] = exp(-other_kappaoverbeta * dist + other_kappaseq) - lsfs2;
	  ds2o[i] = exp(-self_kappaoverbeta * dist + self_kappaseq) - lsfo2;
	}
      dsigma1o = ds1o;
      dsigma2o = ds2o;
    }

  // Distribute the results to the atoms.
  if (!partialupdate)
    {
      // This branch is always taken by EMT.  MonteCarloEMT may take the other.
      if (calculatesigma2)
	{
	  /* Distribute contributions to sigma1 and sigma2 */
	  s1o = &sigma1[zo][0];
	  s1s = &sigma1[zs][0];
	  s2o = &sigma2[zo][0];
	  s2s = &sigma2[zs][0];
#ifdef _OPENMP
#pragma omp critical
#endif // _OPENMP
	  {
	    for (int i = 0; i < n; i++)
	      {
	        int s = self[i];
	        int o = other[i];
	        s1o[s] += dsigma1s[i];
	        s2o[s] += dsigma2s[i];
	        if (o < nAtoms)        // Dont add to ghost atoms
	          {
	            s1s[o] += dsigma1o[i];
	            s2s[o] += dsigma2o[i];
	          }
	      }
	  }
	}
      else
	{
	  /* Distribute contributions to sigma1. */
	  s1o = &sigma1[zo][0];
	  s1s = &sigma1[zs][0];
#ifdef _OPENMP
#pragma omp critical
#endif // _OPENMP
	  {
	    for (int i = 0; i < n; i++)
	      {
	        int s = self[i];
	        int o = other[i];
	        s1o[s] += dsigma1s[i];
	        if (o < nAtoms)
	          s1s[o] += dsigma1o[i];
	      }
	  }
	}
    }
  else  // partialupdate
    {
      // This branch may be taken by MonteCarloEMT.  Since it uses
      // full neighbor lists for partial updates, only update the
      // 'self' atom, not the 'other' atom.
      s1o = &sigma1[zo][0];
      s2o = &sigma2[zo][0];
#ifdef _OPENMP
#pragma omp critical
#endif // _OPENMP
      {
        for (int i = 0; i < n; i++)
          {
            int s = self[i];
            s1o[s] += dsigma1s[i];
            s2o[s] += dsigma2s[i];
          }
      }
    }
  delete[] temporary;
}

// Calculate the energies if Epot != 0, otherwise just calculate the
// derivatives needed for the forces.
void EMT2013::CalculateEnergiesAfterSigmas(bool calc_Epot)
{
  USETIMER("EMT::CalculateEnergiesAfterSigmas");
  DEBUGPRINT;

  //vector<double> sigma(nSize);
  double *sigma = &tmp_double[0];
  bool dosigmapart = (recalc.beforeforces || (calc_Epot && recalc.energies));
  bool doepotpart = (calc_Epot && recalc.energies);

  int zs, zo;
  double s;
  // Better performance if static ???
  assert(nelements < NMAXELEMENTS);
  double invgamma1[NMAXELEMENTS];
  double neginvbetaeta2[NMAXELEMENTS];
  double neglambda[NMAXELEMENTS];
  double lambdaseq[NMAXELEMENTS];
  double negkappa[NMAXELEMENTS];
  double kappaseq[NMAXELEMENTS];
  double nege0lambdalambda[NMAXELEMENTS];
  double e0lambdalambdaseq[NMAXELEMENTS];
  double neg6v0kappa[NMAXELEMENTS];
  double e0lambda[NMAXELEMENTS];
  double eccnst[NMAXELEMENTS];
  double sixv0[NMAXELEMENTS];
  double negsixv0overgamma2[NMAXELEMENTS];
  double seq[NMAXELEMENTS];
  int *id = &(this->id)[0];

  if (dosigmapart)
    {
#ifdef _OPENMP
#pragma omp master
#endif // _OPENMP
      {
        VERB("b");
        DEBUGPRINT;
      }
      /* Calculate total sigma1 */
#ifdef _OPENMP
#pragma omp for
#endif // _OPENMP
      for (int i = 0; i < nAtoms; i++)
        {
          double s = 0.0;
          int zs = id[i];
          for (int zo = 0; zo < nelements; zo++)
            s += (*chi)[zs][zo] * sigma1[zo][i];
          if (s < 1.0e-9)
            s = 1.0e-9;
          sigma[i] = s;
        }
      assert(nAtoms == radius.size() && nAtoms == Ec.size() &&
             nSize == dEds.size());
    }

  /* Calculate combinations of EMT parameters */
  for (int i = 0; i < nelements; i++)
    {
      invgamma1[i] = 1.0 / (parameters[i]->gamma1);
      neginvbetaeta2[i] = -1.0 / (beta * parameters[i]->eta2);
      neglambda[i] = - parameters[i]->lambda;
      lambdaseq[i] = parameters[i]->lambda * parameters[i]->seq;
      negkappa[i] = - parameters[i]->kappa;
      kappaseq[i] = parameters[i]->kappa * parameters[i]->seq;
      nege0lambdalambda[i] = - parameters[i]->e0 * parameters[i]->lambda *
        parameters[i]->lambda;
      e0lambdalambdaseq[i] = parameters[i]->e0 * parameters[i]->lambda *
        parameters[i]->lambda * parameters[i]->seq;
      neg6v0kappa[i] = - 6.0 * parameters[i]->V0 * parameters[i]->kappa;
      e0lambda[i] = parameters[i]->e0 * parameters[i]->lambda;
      eccnst[i] = parameters[i]->e0 * (1.0 - parameters[i]->lambda *
                                       parameters[i]->seq);
      sixv0[i] = 6.0 * parameters[i]->V0;
      negsixv0overgamma2[i] = -6 * parameters[i]->V0 /
        parameters[i]->gamma2;
      seq[i] = parameters[i]->seq;
    }

  if (dosigmapart)
  {
#ifdef _OPENMP
#pragma omp for
#endif // _OPENMP
    for (int i = 0; i < nAtoms; i++)
    {
      asap_z_int z = id[i];  // Actually not Z (atomic no) but ID no.
      double r = seq[z] + neginvbetaeta2[z] * log(sigma[i] * invgamma1[z]);
      radius[i] = r;
      double ex1 = exp(neglambda[z] * r + lambdaseq[z]);
      double e2 = exp(negkappa[z] * r + kappaseq[z]);
      ex2[i] = e2;
      double tmp = (nege0lambdalambda[z] * r + e0lambdalambdaseq[z]) * ex1
           + neg6v0kappa[z] * e2;
      dEds[i] = tmp * neginvbetaeta2[z] / sigma[i];
      Ec[i] = (e0lambda[z] * r + eccnst[z]) * ex1;
    }
#ifdef _OPENMP
#pragma omp for
#endif // _OPENMP
    for (int i = nAtoms; i < nSize; i++)
      dEds[i] = 0.0;
  }
  if (calc_Epot)
    {
      DEBUGPRINT;
      if (doepotpart)
        {
#ifdef _OPENMP
#pragma omp master
#endif // _OPENMP
          {
            VERB("e");
            DEBUGPRINT;
          }
          /* We also need Eas, but only for the real atoms */
          assert(sigma2isvalid);
          assert(counters.sigma2 == atoms->GetPositionsCounter());
          /* Calculate total sigma2 */
#ifdef _OPENMP
#pragma omp for
#endif // _OPENMP
          for (int i = 0; i < nAtoms; i++)
            {
              s = 0.0;
              zs = id[i];
              for (zo = 0; zo < nelements; zo++)
                s += (*chi)[zs][zo] * sigma2[zo][i];
              Eas[i] = sixv0[zs] * ex2[i] + negsixv0overgamma2[zs] * s;
            }
        }
      /* Add the energies */
      DEBUGPRINT;
      if(subtractE0)
        {
#ifdef _OPENMP
#pragma omp for
#endif // _OPENMP
          for (int i = 0; i < nAtoms; i++)
            Epot[i] = Ec[i] + Eas[i] - parameters[id[i]]->e0;
        }
      else
        {
#ifdef _OPENMP
#pragma omp for
#endif // _OPENMP
          for (int i = 0; i < nAtoms; i++)
            Epot[i] = Ec[i] + Eas[i];
        }
    } // if (calc_Epot)

  DEBUGPRINT;
}

void EMT2013::force_batch(const int *self, const int *other,
                      const Vec *rnb, const double *sq_dist,
                      const double *dEdss, const double *dEdso,
                      int zs, int zo, int n)
{
  USETIMER("EMT::force_batch");
  double *temporary = new double[BUFLEN];
  double *df = temporary;
  /* DELETE double cutslopecutdist;  DELETE */
  double other_eta2betaseq, other_kappaoverbeta;
  double other_kappaseq, self_eta2betaseq, self_kappaoverbeta;
  double self_kappaseq;
  const emt_parameters *emtself, *emtother;
  //double pairA, exprcut, pairD;
  //double *tmp;


  // initialization of lsf parameters;
  /* TEMP */ double rcut, ds1sdrcut, ds2sdrcut, s1srcut, s2srcut,
    	  	  	  	ds1odrcut, ds2odrcut, s1orcut, s2orcut; /* TEMP */


  assert(n <= BUFLEN);

  /* Get EMT parameters */
  emtself = parameters[zs];
  emtother = parameters[zo];
  // cutslopecutdist = cutoffslope * rFermi;
  other_eta2betaseq = emtother->eta2 * beta * emtother->seq;
  other_kappaoverbeta = emtother->kappa / beta;
  other_kappaseq = emtother->kappa * emtother->seq;
  self_eta2betaseq = emtself->eta2 * beta * emtself->seq;
  self_kappaoverbeta = emtself->kappa / beta;
  self_kappaseq = emtself->kappa * emtself->seq;
  double other_eta2 = emtother->eta2;
  double self_eta2 = emtself->eta2;


  /* weight function parameters, could/should be provided as Arrays to avoid calculations! */
  /* TEMP rcut = 0.5 * (sqrt(3. / 2.) + sqrt(2.)) * sqrt(2)
     	            * beta * max(emtother->seq,emtself->seq); TEMP */
  rcut = sqrt(rcut2[zs][zo]);
  s1srcut = exp(-other_eta2 * rcut + other_eta2betaseq);
  /* ds1sdrcut are equal to the derivative of the cutoff function */
  ds1sdrcut = -other_eta2 * s1srcut;
  if (zs != zo)
    {
     s1orcut = exp(-self_eta2 * rcut + self_eta2betaseq);
     /* ds1odrcut are equal to the derivative of the cutoff function */
     ds1odrcut = -self_eta2 * s1orcut;
    }
  s2srcut = exp(-other_kappaoverbeta * rcut + other_kappaseq);
  /* ds2sdrcut are equal to the derivative of the cutoff function */
  ds2sdrcut = -other_kappaoverbeta * s2srcut;
  if (zs != zo)
    {
      s2orcut = exp(-self_kappaoverbeta * rcut + self_kappaseq);
      /* ds2odrcut are equal to the derivative of the cutoff function */
   	  ds2odrcut = -self_kappaoverbeta * s2orcut;
   	}

  // chi for the two elements involved
    double chi_zs_zo = (*chi)[zs][zo];
    double chi_zo_zs = (*chi)[zo][zs];

  /* Derivative of AS energy with respect to sigma_2
   * (chi is multiplied here for performance increase) */
  double dEasds_s = -6 * emtself->V0 / emtself->gamma2 * chi_zs_zo;
  double dEasds_o = -6 * emtother->V0 / emtother->gamma2 * chi_zo_zs;

  if (zs == zo)
    {
      for (int i = 0; i < n; i++)
	  {
	    /* Get the distances from their squares */
	    double dist = sqrt(sq_dist[i]);
	    double inv_dist = 1.0 / dist;
	    /* Calculations of the derivatives of sigma_1 and sigma_2 */
	    double dsigma1dr = -other_eta2 *
	    		           exp(-other_eta2 * dist + other_eta2betaseq)
	                       - ds1sdrcut;
	    double dsigma2dr = -other_kappaoverbeta *
	                       exp(-other_kappaoverbeta * dist + other_kappaseq)
	                       - ds2sdrcut;
	    /* Force contribution */
	    df[i] = inv_dist * (dsigma1dr * dEdss[i] * chi_zs_zo
		  	        + dEasds_s * dsigma2dr
			        + dsigma1dr * dEdso[i] * chi_zo_zs
			        + dEasds_o * dsigma2dr * (other[i] < nAtoms));
	  }
    }
  else
    {
      for (int i = 0; i < n; i++)
	{
	  /* Get the distances from their squares */
	  double dist = sqrt(sq_dist[i]);
	  double inv_dist = 1.0 / dist;
      /* Calculations of the derivatives of sigma_1 and sigma_2 */
      double dsigma1dr_o = -other_eta2 *
	                       exp(-other_eta2 * dist + other_eta2betaseq)
	                       - ds1sdrcut;
	  double dsigma2dr_o = -other_kappaoverbeta *
	                       exp(-other_kappaoverbeta * dist + other_kappaseq)
	                       - ds2sdrcut;
      double dsigma1dr_s = -self_eta2 *
	                       exp(-self_eta2 * dist + self_eta2betaseq)
	                       - ds1odrcut;
	  double dsigma2dr_s = -self_kappaoverbeta *
	                       exp(-self_kappaoverbeta * dist + self_kappaseq)
	                       - ds2odrcut;
	/* Force contribution */
    df[i] = inv_dist * (dsigma1dr_o * dEdss[i] * chi_zs_zo
			      + dEasds_s * dsigma2dr_o
			      + dsigma1dr_s * dEdso[i] * chi_zo_zs
			      + dEasds_o * dsigma2dr_s * (other[i] < nAtoms));
	}
    }
  distribute_force(self, other, df, rnb, n);
  delete[] temporary;
}

// Will be replaced in MonteCarloEMT
void EMT2013::CreateNeighborList()
{
  MEMORY;
  USETIMER("EMT::CreateNeighborList");
  if (!initialized)
    THROW_RETURN( AsapError("EMT object has not been initialized!") );
#ifdef _OPENMP
#pragma omp single
#endif // _OPENMP
  {
      PyAsap_NeighborLocatorObject *nbl;
      nbl = PyAsap_NewNeighborList2013(atoms, rNbCut, driftfactor, rcut2_NB);
      nblist = nbl->cobj;
      nblist_obj = (PyObject *) nbl;
      MEMORY;
  }
  nblist->UpdateNeighborList();
  MEMORY;
}

/// IDs are a recoding of the atomic numbers used as indices into
/// various arrays.  Atomic numbers are not practical for this as they
/// are not contiguous.
void EMT2013::CalculateIDs()
{
  RETURNIFASAPERROR;
  USETIMER("EMT::CalculateIDs");
  DEBUGPRINT;
  if (!recalc.ids)
    return;
  VERB("i");
  DEBUGPRINT;
  // If there is only one element, all IDs are 0 and do not have to be
  // calculated.  Otherwise the ID is the offset into the list of elements.
  const asap_z_int *z = atoms->GetAtomicNumbers();
#ifdef _OPENMP
#pragma omp master
#endif // _OPENMP
  id_assigned = 0; // Member variable to permit OpenMP reduction.
#ifdef _OPENMP
#pragma omp barrier
#endif // _OPENMP
  
  int assigned = 0;
  for (int i = 0; i < nelements; i++)
    {
      int zcand = parameters[i]->Z;
#ifdef _OPENMP
#pragma omp for
#endif // _OPENMP
      for (int j = 0; j < nSize; j++)
        if (z[j] == zcand)
          {
            id[j] = i;
            assigned += 1;
          }
    }
#ifdef _OPENMP
#pragma omp critical
#endif // _OPENMP
  id_assigned += assigned; // Manual reduction needed inside function.

#ifdef _OPENMP
#pragma omp barrier
#endif // _OPENMP  
  if (id_assigned != nSize)
    THROW(AsapError("An unknown element was encountered."));
#ifdef _OPENMP
#pragma omp master
#endif // _OPENMP
  counters.ids = atoms->GetPositionsCounter();
  DEBUGPRINT;
}

/// Return a copy (new reference) of the parameter dictionary
PyObject *EMT2013::GetParameterDict() const
{
  // Perform a deep copy of the parameter object
  Py_ssize_t n = 0;
  PyObject *key;
  PyObject *value;
  PyObject *result = PyDict_New();
  assert(result != NULL);
  while (PyDict_Next(param_obj, &n, &key, &value))
    {
      if (!PyDict_Check(value))
        throw AsapError("EMT2013::GetParameterDict found non-dictionary in parameter dictionary.");
      PyObject *cp_value = PyDict_Copy(value);
      if (PyDict_SetItem(result, key, cp_value) < 0)
        throw AsapError("EMT2013::GetParameterDict failed to copy parameter dictionary");
      Py_DECREF(cp_value);
    }
  return result;
}
