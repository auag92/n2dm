// -*- C++ -*-

/*
 This code was released into the public domain by Tim Freeman on 29 Aug 2000.
*/

#ifndef __ATOM_PAIR_INFO_STATE_H__
#define __ATOM_PAIR_INFO_STATE_H__

#include "asapbrenner.h"

namespace ASAPSPACE {

/* This data structure holds miscellaneous information about the bonds.  It's
   similar to the CCNeighborState structureq.  The difference is that I intend
   to replace CCNeighborState with something else when this becomes a Fungimol
   plugin, but AtomPairInfoState will be permanently a part of the bond order
   force field.  Thus the force field can manipulate AtomPairInfoState
   directly, but it has to go through an API to deal with CCNeighborState.  Tim
   Freeman 28 Aug 2000 */
struct AtomPairInfoState {
public:
  vector<AtomPairInfo> ai_list;
  /* For atom number i, ai_list [ai_start [i]] is the first
     neighbor of atom i.  ai_list [ai_start [i+1] - 1] is the last
     neighbor of atom i.  Looking at the last atom, we'll need room for one
     more value in ai_start than we have atoms. */
  vector<int> ai_start;
};

/* Call this every time after a new neighbor is added.  AtomPairInfo_end should
   be the position of where we want to add the *next* neighbor, and we'll make
   room for it.  By using conservative high estimates for AtomPairInfo_end, we
   can call this only once per atom instead of once per neighbor.  Thus this
   routine doesn't need to be inline. */
inline void allocate_AtomPairInfo(struct AtomPairInfoState *ais,
				  const int AtomPairInfo_end)
{
  ais->ai_list.resize(AtomPairInfo_end);
}

inline int numPairs(int ai, const struct AtomPairInfoState *s)
{
    return s->ai_start[ai+1] - s->ai_start[ai];
}

inline struct AtomPairInfo *pairs(int ai, struct AtomPairInfoState *s)
{
  return &s->ai_list[s->ai_start[ai]];
}

inline const struct AtomPairInfo *pairsconst(int ai,
					     const struct AtomPairInfoState *s)
{
  return &s->ai_list[s->ai_start[ai]];
} 

inline void allocateAtomPairAtom(struct AtomPairInfoState *ljns, int max_atoms)
{
  ljns->ai_start.resize(max_atoms+1);
}

struct AtomPairInfoState *newAtomPairInfoState ();

void deallocateAtomPairInfoState (struct AtomPairInfoState *ljns);

} // end namespace

#endif

