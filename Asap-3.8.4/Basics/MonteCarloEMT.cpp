// Copyright (C) 2008-2011 Jakob Schiotz and Center for Individual
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

#include "MonteCarloEMT.h"
#include "NeighborList.h"
#include "MonteCarloAtoms.h"
#include "Timing.h"
#include "Debug.h"
#include <iostream>
#include <math.h>

using std::cerr;
using std::endl;
using std::flush;

extern int verbose;
#define VERB(x) if (verbose == 1) cerr << x

MonteCarloEMT::MonteCarloEMT(PyObject *self, PyObject *prov)
  : EMT(self, prov)
{
  mc_atoms = NULL;
}

MonteCarloEMT::~MonteCarloEMT()
{
  if (mc_atoms != NULL)
    AsapAtoms_DECREF(mc_atoms);
}

void MonteCarloEMT::SetAtoms(PyObject *pyatoms, Atoms* accessobj)
{
  assert(accessobj == NULL);
  if (atoms != NULL)
    throw AsapError("Cannot use the same MonteCarloEMT object for multiple Atoms objects.");
  mc_atoms = new MonteCarloAtoms();
  assert(mc_atoms != NULL);
  EMT::SetAtoms(pyatoms, mc_atoms);
}

void MonteCarloEMT::CreateNeighborList()
{
#ifdef _OPENMP
#pragma omp single
#endif // _OPENMP
  {
    PyAsap_NeighborLocatorObject *nbl;
    nbl = PyAsap_NewNeighborList(atoms, rNbCut, driftfactor);
    nblist = nbl->cobj;
    nblist_obj = (PyObject *) nbl;
    NeighborList *realnblist = dynamic_cast<NeighborList*>(nblist);
    assert(realnblist != NULL);
    realnblist->EnableFullNeighborLists();
  }
  nblist->UpdateNeighborList();
}

const vector<double> &MonteCarloEMT::GetPotentialEnergies(PyObject *pyatoms)
{
  DEBUGPRINT;
  atoms->Begin(pyatoms);
  recalc.energies = (counters.energies != atoms->GetPositionsCounter());
  bool mc_optim = (mc_atoms->GetMonteCarloRelevant());
  if (mc_optim)
    {
      // Monte Carlo optimizations ...
      USETIMER("MonteCarloEMT::GetPotentialEnergies");
      VERB(" MCEnergies[");
      DEBUGPRINT;
      if (counters.energies != atoms->GetPositionsCounter())
	{
          DEBUGPRINT;
	  const set<int> &modified_atoms = mc_atoms->GetModifiedAtoms();
	  set<int> affected_atoms;
	  assert(modified_atoms.size() > 0);
	  PartialCalculateIDs(modified_atoms);
	  NeighborList *nbl = dynamic_cast<NeighborList*>(nblist);
	  assert(nbl != NULL);
	  nbl->RemakeLists(modified_atoms, affected_atoms);
	  counters.nblist = atoms->GetPositionsCounter();
	  VERB("(" << modified_atoms.size() << "/" << affected_atoms.size()
	       << " atoms)");
	  PartialCalculateSigmas(affected_atoms);
	  PartialCalculateEnergiesAfterSigmas(affected_atoms);
	  // CheckNeighborLists();
	  // CalculateIDs();
	  // CalculateSigmas();
	  // CalculateEnergiesAfterSigmas(&(potentialenergy[0]));
	}
      else
	{
	  VERB("-");
	}
      VERB("]" << flush);
      DEBUGPRINT;
      mc_atoms->MonteCarloEnd();
      atoms->End();
      return Epot;
    }
  else
    {
      // Either the last update was not a single atom update, or the energy
      // was not up to date before the sequence of single atom updates.
      skip_begin = true;  // Atoms are already open
      // Make an ordinary energy calculation
      DEBUGPRINT;
      const vector<double> &result = EMT::GetPotentialEnergies(pyatoms);  
      // Atoms were closed in EMT::GetPotentialEnergies
      DEBUGPRINT;
      mc_atoms->MonteCarloEnd();
      return result;
    }
}

/// Update the identities of modified atoms, in case the modification
/// included changing the atomic number.  No new elements may be
/// introduced.
void MonteCarloEMT::PartialCalculateIDs(const set<int> &changedatoms)
{
  USETIMER("MonteCarloEMT::PartialCalculateIDs");
  if (counters.ids == atoms->GetPositionsCounter() || nelements == 1)
    return;
  counters.ids = atoms->GetPositionsCounter();
  if (nelements > 1)
    {
      // Make a map from atomic numbers to ids.  The map is only made
      // the first time.
      if (zmap.size() == 0)
	{
	  for (int i = 0; i < nelements; i++)
	      zmap[parameters[i]->Z] = i;
	}
      // Update the IDs.
      const asap_z_int *z = atoms->GetAtomicNumbers();
      for (set<int>::const_iterator i = changedatoms.begin();
	   i != changedatoms.end(); ++i)
	id[*i] = zmap[z[*i]];
    }
}

void MonteCarloEMT::PartialCalculateSigmas(const set<int> &changedatoms)
{
  USETIMER("MonteCarloEMT::PartialCalculateSigmas");
  int ctr = atoms->GetPositionsCounter();
  if ((counters.sigma1 == ctr) && (counters.sigma2 == ctr))
    return;
  counters.sigma1 = counters.sigma2 = ctr;
  assert(sigma2isvalid);
  
  int maxnblen = nblist->MaxNeighborListLength();
  if (maxnblen > BUFLEN)
    throw AsapError("Neighborlist overrun (did you squeeze your atoms?).  Longest neighbor list is ") << maxnblen;

  // Buffer data:
  TinyMatrix<int> nbatch(nelements,nelements);
  TinyMatrix<int[BUFLEN]> self(nelements,nelements);
  TinyMatrix<int[BUFLEN]> other(nelements,nelements);
  TinyMatrix<Vec[BUFLEN]> rnb(nelements,nelements);
  TinyMatrix<double[BUFLEN]> sqdist(nelements,nelements);
  int other_buf[BUFLEN];
  Vec rnb_buf[BUFLEN];
  double sqdist_buf[BUFLEN];

  /* Set sigmas of affected atoms to zero */
  for (set<int>::const_iterator i = changedatoms.begin();
       i != changedatoms.end(); ++i)
    for (int j = 0; j < nelements; ++j)
      sigma1[j][*i] = sigma2[j][*i] = 0.0;

  /* No atoms in batch pools */
  for (int i = 0; i < nelements; i++)
    for (int j = 0; j < nelements; j++)
      nbatch[i][j] = 0;

  // Loop over atoms
  for (set<int>::const_iterator atomptr = changedatoms.begin();
       atomptr != changedatoms.end(); ++atomptr) 
    {
      int atom = *atomptr;
      int zself = id[atom];
      assert(zself < nelements);
      // Get neighbors and loop over them.  Simplest if only one element
      if (nelements == 1) 
        {
          int nbat = nbatch[0][0];  // only element
          int size = BUFLEN-nbat;
          int n = nblist->GetFullNeighbors(atom, other[0][0]+nbat,
					   rnb[0][0]+nbat, sqdist[0][0]+nbat,
					   size);
          assert(size >= 0);    // REMOVE LATER !!!
          for (int i = nbat; i < nbat+n; i++)
            self[0][0][i] = atom;
          nbatch[0][0] += n;
        } 
      else 
        {
          int size = BUFLEN;
          int n = nblist->GetFullNeighbors(atom, other_buf, rnb_buf, sqdist_buf,
					   size);
          assert(size >= 0);     // REMOVE LATER !!!
          for (int i = 0; i < n; i++) 
            {
              int zother = id[other_buf[i]];
	      assert(zother < nelements);
              int nbat = nbatch[zself][zother]++;  // Count this atom
              self[zself][zother][nbat] = atom;
              other[zself][zother][nbat] = other_buf[i];
              rnb[zself][zother][nbat][0] = rnb_buf[i][0];
              rnb[zself][zother][nbat][1] = rnb_buf[i][1];
              rnb[zself][zother][nbat][2] = rnb_buf[i][2];
              sqdist[zself][zother][nbat] = sqdist_buf[i];
            }
        }
      // Now process any full batch
      for (int zo = 0; zo < nelements; zo++)
        if (nbatch[zself][zo] >= BUFLEN - maxnblen) 
          {
            sigma_batch(self[zself][zo], other[zself][zo], rnb[zself][zo],
                        sqdist[zself][zo], zself, zo, nbatch[zself][zo],
                        true, true);
            nbatch[zself][zo] = 0;
          }
    }  // Loop over atoms
  /* Process the remaining incomplete batches */
  for (int zs = 0; zs < nelements; zs++)
    for (int zo = 0; zo < nelements; zo++)
      if (nbatch[zs][zo])
        sigma_batch(self[zs][zo], other[zs][zo], rnb[zs][zo],
                    sqdist[zs][zo], zs, zo, nbatch[zs][zo], true, true);
}

void MonteCarloEMT::PartialCalculateEnergiesAfterSigmas(const set<int>
							&changedatoms)
{
  int i;
  int zo;
  // Better performance if static ???
  assert(nelements < NMAXELEMENTS);
  double inv12gamma1[NMAXELEMENTS];
  double neglambda[NMAXELEMENTS];
  double lambdaseq[NMAXELEMENTS];
  double negkappa[NMAXELEMENTS];
  double kappaseq[NMAXELEMENTS];
  double nege0lambdalambda[NMAXELEMENTS];
  double e0lambdalambdaseq[NMAXELEMENTS];
  double neg6v0kappa[NMAXELEMENTS];
  double invbetaeta2[NMAXELEMENTS];
  double e0lambda[NMAXELEMENTS];
  double eccnst[NMAXELEMENTS];
  double sixv0[NMAXELEMENTS];
  double neghalfv0overgamma2[NMAXELEMENTS];
  double seq[NMAXELEMENTS];
  int *id = &(this->id)[0];
    
  /* Calculate conbinations of EMT parameters */
  for (i = 0; i < nelements; i++)
    { 
      inv12gamma1[i] = 1.0 / (12.0 * parameters[i]->gamma1);
      neglambda[i] = - parameters[i]->lambda;
      lambdaseq[i] = parameters[i]->lambda * parameters[i]->seq;
      negkappa[i] = - parameters[i]->kappa;
      kappaseq[i] = parameters[i]->kappa * parameters[i]->seq;
      nege0lambdalambda[i] = - parameters[i]->e0 * parameters[i]->lambda *
        parameters[i]->lambda;
      e0lambdalambdaseq[i] = parameters[i]->e0 * parameters[i]->lambda *
        parameters[i]->lambda * parameters[i]->seq;
      neg6v0kappa[i] = - 6.0 * parameters[i]->V0 * parameters[i]->kappa;
      invbetaeta2[i] = 1.0 / (Beta * parameters[i]->eta2);
      e0lambda[i] = parameters[i]->e0 * parameters[i]->lambda;
      eccnst[i] = parameters[i]->e0 * (1.0 - parameters[i]->lambda *
                                       parameters[i]->seq);
      sixv0[i] = 6.0 * parameters[i]->V0;
      neghalfv0overgamma2[i] = -0.5 * parameters[i]->V0 /
        parameters[i]->gamma2;
      seq[i] = parameters[i]->seq;
    }

  assert(counters.beforeforces != atoms->GetPositionsCounter() ||
	 counters.energies != atoms->GetPositionsCounter());
    
  counters.beforeforces = counters.energies = atoms->GetPositionsCounter();
  VERB("E");
  /* Calculate total sigma1 */
  for (set<int>::const_iterator aptr = changedatoms.begin();
       aptr != changedatoms.end(); ++aptr)
    {
      int i = *aptr;  // The offset of the atom
      double sigma = 0.0;
      int z = id[i];
      for (zo = 0; zo < nelements; zo++)
	sigma += (*chi)[z][zo] * sigma1[zo][i];
      if (sigma < 1.0e-9)
	sigma = 1.0e-9;
      
      radius[i] = seq[z] - invbetaeta2[z] * log(sigma * inv12gamma1[z]);
      /* dEds */
      double ex1 = exp(neglambda[z] * radius[i] + lambdaseq[z]);
      double ex2 = exp(negkappa[z] * radius[i] + kappaseq[z]);
      dEds[i] = (nege0lambdalambda[z] * radius[i] + e0lambdalambdaseq[z])
	* ex1 + neg6v0kappa[z] * ex2;
      dEds[i] *= - invbetaeta2[z] / sigma;
      /* Cohesive energy */
      Ec[i] = (e0lambda[z] * radius[i] + eccnst[z]) * ex1;
      /* We also need Eas, but only for the real atoms */
      assert(sigma2isvalid);
      assert(counters.sigma2 == atoms->GetPositionsCounter());
      /* Calculate total sigma2 */
      sigma = 0.0;
      z = id[i];
      for (zo = 0; zo < nelements; zo++)
	sigma += (*chi)[z][zo] * sigma2[zo][i];
      if (sigma < 1.0e-9)
	sigma = 1.0e-9;
      /* Atomic-sphere energy */
      Eas[i] = sixv0[z] * ex2 + neghalfv0overgamma2[z] * sigma;
      if(subtractE0)
	Epot[i] = Ec[i] + Eas[i] - parameters[id[i]]->e0;
      else
	Epot[i] = Ec[i] + Eas[i];
    }
}
