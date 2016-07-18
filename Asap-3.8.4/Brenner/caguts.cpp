/*
 This code was released into the public domain by Peter McCluskey on 6/2/2000.
 It is adapted from code released into the public domain by Donald Brenner on
 2/29/2000.
*/

/* The code in this file should be "pure", in the sense that it should only
   access the physical state of the system through the API defined in api.h and
   neighborstate.h.  Then, when I want to make the code access Fungimol data
   structures, all I have to do is change the next two includes.  Tim Freeman
   5 Sep 2000. */

///////////////////////////////////////////////////////////////
//
// Modifications for ASAP (5. Feb 2010, Jakob Schiotz):
//
// * Converted to C++.
//
// * New header files.
//
// * bren_vec is now an Asap Vec.
//
// * caguts() made a method of BrennerPotential, and info about the
//   atoms are now data in BrennerPotential.
//
// * AtomPairInfoState now uses STL vectors instead of own memory management.
//
// * Use ASAP neighbor list instead of own list, as the latter
//   apparently contains memory corruption bugs at least on 64 bit
//   machines.
//
// * getKtype(i, info) replaced by getKtype(i) which is now a method.
//   It still returns a number in the range 1-4 instead of the more
//   natural 0-3.
//
// * Distance calculation and periodic boundaries now handled by the
//   NeighborList object.
//
///////////////////////////////////////////////////////////// 



//#include "Asap.h"
#include "BrennerPotential.h"
#include "AtomPairInfoState.h"
#include "asapbrenner.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#undef min
#define min(x, y) ((x) < (y) ? (x) : (y))


static Float drtable[4][4][NTAB];
/* ddtab is the amount we multiply a radius by to get an index into
   tabfc and tabdfc.  Note that it used to be an amount we *divide*
   by, but I changed this  5 Aug 2000 because multiplication is
   faster.  Tim Freeman */
static Float ddtab[4][4];
/* rtable holds values for V super R sub i j, defined in equation 4 on page
   7. */
static Float rtable[4][4][NTAB];
/* atable holds values for V super A sub i j, defined in equation 5 on page
   7. */
static Float atable[4][4][NTAB];
/* datable holds values for the derivative of V super A sub i j with respect to
   r. */
static Float datable[4][4][NTAB];
/* tabfc encodes f super c sub ij (r sub ij), mentioned in equations 4 and 5 on
   page 7, defined in equation 18 on page 16. The first coordinates are the
   atomic numbers of the two atoms, and the last coordinate is the radius,
   divided by a constant chosen to make the range right, then truncated to an
   integer. */ 
static Float tabfc[4][4][NTAB];
/* tabdfc has the derivative of f super c sub i j with respect to r. */
static Float tabdfc[4][4][NTAB];

/* Beware that this initializes rb2, which is only used in ljguts.  Thus you'd
   better call this before calling ljguts. */
void BrennerPotential::init_c()
{
  DEBUGPRINT;
  mtable(drtable, ddtab, rtable, atable, datable, tabfc, tabdfc);
}


/* We've already figured out neighbors.  Do the physics for Brenner's bond
   order potential.  This routine and all code that it calls should be "pure",
   in the sense that it should only access the state through the api defined in
   api.h, because that gives hope of making this into a Fungimol plugin.

   I can't remember what "ca" stands for; replacing this sentence
   with positive information would be of value.  Tim Freeman  5 Sep 2000*/
//Float
//caguts(BrennerMainInfo *info)

Float BrennerPotential::caguts()
{
  // The number of atoms
  const int num_atms = nAtoms;
  /* k is the position of the current neighbor in apis. */
  int k = 0;
  int l;
  // Allocate buffers for neighbor list data
  int nmaxnb = nblist->MaxNeighborListLength();
  vector<int> iNeighbors(nmaxnb);
  vector<Vec> distances(nmaxnb);
  vector<double> sq_distances(nmaxnb);
  
  /* tote is the total energy, for the return value of caguts. */
  Float tote = 0;
  
  allocateAtomPairAtom(apis, num_atms);  
  /* Loop over each bond.  Actually we'll see each bond twice, once from each
     end. */
  {
    int i;
    for (i = 0; i < num_atms; i++) {
      int sizemax = nmaxnb;
      const int iNeighborCount =
	nblist->GetFullNeighbors(i, &iNeighbors[0], &distances[0],
				 &sq_distances[0], sizemax);
      const int ki = getKtype(i) - 1;
      struct AtomPairInfo *pair_i;
      int jn;
      apis->ai_start[i] = k; 
      k += iNeighborCount;
      /* Early escape if there is no neighbours to this atoms. This also avoids a memory access error */
      if (iNeighborCount == 0) continue;	
      allocate_AtomPairInfo (apis, k);
      pair_i = pairs (i, apis);
      for (jn = 0; jn < iNeighborCount; jn++) {
        const int j = iNeighbors [jn];
        const int kj = getKtype(j) - 1; 
        int it, floor_rt;
        /* rsq will be the squared length of the bond. */
        Float rsq;
        /* rc will be the length of the bond. */
        Float rc, rt, vv, rp;
        /* We're going to fill in *pair_k. */
        struct AtomPairInfo *thisPair = &(pair_i[jn]);
        thisPair->lcheck = 0;
        //rsq = CALC_DIST(info, i, j, &thisPair->cor);
	rsq = sq_distances[jn];
	thisPair->cor = -distances[jn];
        /* Now rsq is the square of the distance between the two neighbors, and
           thisPair->cor is the vector from j to i. */
        /* Hmm, many pairs of atoms will be far apart by any measure, so we
           might want to have one number that's the maximum rmax to avoid some
           array indexing in this common case.  Small deal. */
        if(rsq > rmax[ki][kj]) continue;
        /* If both are hydrogen or carbon, set lcheck to 1. */
        if(kj <= 1 && ki <= 1) thisPair->lcheck = 1;
        /* If neither are hydrogen or carbon, set lcheck to 2. */
        if(kj >= 2 && ki >= 2) thisPair->lcheck = 2;
        /* otherwise lcheck remains 0. */
        /* set rc to the distance between them. */
        rc = sqrt(rsq);
        /* We used to divide by the value in ddtab here, but I changed the
           definition of ddtab so we multiply now.  Tim Freeman  5 Aug 2000. */
        rt = rc * ddtab[ki][kj];
        floor_rt = (int)floor(rt);
        it = min(floor_rt, NTAB-2);
        thisPair->rcor = rc;
        /* Linear interpolation to compute f super c sub i j (rc), defined in
           equation 18 on page 16. */
        thisPair->ww = tabfc[ki][kj][it]
          + (tabfc[ki][kj][it+1] - tabfc[ki][kj][it])*(rt-it);
        thisPair->dww = tabdfc[ki][kj][it]
          + (tabdfc[ki][kj][it+1] - tabdfc[ki][kj][it])*(rt-it);
        thisPair->exx1 = atable[ki][kj][it]
          + (atable[ki][kj][it+1] - atable[ki][kj][it])*(rt-it);
        thisPair->dexx1 = datable[ki][kj][it] +
          (datable[ki][kj][it+1] - datable[ki][kj][it])*(rt-it);
        if(i >= j) continue;
        /* For the rest of the loop we're only seeing each bond once. */
        vv = rtable[ki][kj][it]
          + (rtable[ki][kj][it+1] - rtable[ki][kj][it])*(rt-it);
        /* Now vv is the energy for the pair due to their pairwise repulsion. */
        rp = drtable[ki][kj][it]
          + (drtable[ki][kj][it+1] - drtable[ki][kj][it])*(rt-it);
        /* Now rp is the magnitude of the force from pairwise repulsion. */
        tote = tote + vv;
        
        /* rpp is the repulsive force between the pair. */
        thisPair->rpp = rp*thisPair->cor;
      } /* End loop over neighors. */
    } /* End loop over atoms. */
    assert (num_atms == i);
    apis->ai_start [i] = k;
  } /* Forget i. */
  
  { /* FIXME Is there any reason not to merge this loop with the previous one?
     */
    int i;
    for (i = 0; i < num_atms; i++) {
      // This extra call to GetFullNeighbors is costly!
      int sizemax = nmaxnb;
      const int iNeighborCount =
	nblist->GetFullNeighbors(i, &iNeighbors[0], &distances[0],
				 &sq_distances[0], sizemax);
      /* Early escape if there is no neighbours to this atoms. This also avoids a memory access error */
      if (iNeighborCount == 0) continue;				 
      const struct AtomPairInfo *const apiNeighbors = pairs (i, apis);
      int jn;
      assert (iNeighborCount == numPairs(i, apis));
      for (jn = 0; jn < iNeighborCount; jn++) {
        const int j = iNeighbors [jn];
        if(apiNeighbors[jn].lcheck == 0) continue;
        if(i >= j) continue;
        //transForcex (i, j, info, apiNeighbors[jn].rpp.x);
        //transForcey (i, j, info, apiNeighbors[jn].rpp.y);
        //transForcez (i, j, info, apiNeighbors[jn].rpp.z);
	transForce(i, j, apiNeighbors[jn].rpp);
      }
    }
  }
  /* If there are hydrogens or carbons, then do pibond. */
  if(getNoa(1) + getNoa(2) != 0)   // Bug?  This was 0 and 1 originally.
    tote += pibond();
  /* If there are silicons or germaniums, then do sili_germ. */
  if(getNoa(3) + getNoa(4) != 0)
    tote += sili_germ();
  return tote;
}
