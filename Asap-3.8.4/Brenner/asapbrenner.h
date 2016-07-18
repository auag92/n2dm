// Common header file for the ASAP Brenner potential.

#ifndef ASAPBRENNER_H
#define ASAPBRENNER_H

#include <Vec.h>
#include "BrennerPotential.h"
//#define ASAPDEBUG
#include "Debug.h"

typedef Vec bren_vector;
typedef double Float;
typedef double Double;

// NTAB is defined as BRENNERNTAB in BrennerPotential.h to avoid name clashes.
#define NTAB BRENNERNTAB

namespace ASAPSPACE {

struct AtomPairInfo
{
  /* lcheck is 0, 1, or 2.  1 means both ends are either hydrogen or carbon.
     2 means both ends are either silicon or germanium.  0 means some other
     combination.  We set this up in the main caguts routine. */
  short lcheck;
  /* cor is the bren_vector from the first atom to the second. */
  bren_vector cor;
  /* rpp is the pairwise repulsion between them. */
  bren_vector rpp;
  /* rcor is the length of cor. */
  Float rcor;
  /* ww is the value of f super c sub i j defined in equation 18 on page 16 for
     rcor. */
  Float ww;
  /* dww is the value of partial (f super c sub i j (r)) per partial r
     evaluated at rcor. */
  Float dww;
  /* exx1 is the attraction between the pair, V super A sub i j in equation 5
     on page 7. */
  Float exx1;
  /* dexx1 is the partial of exx1 with respect to rcor. */
  Float dexx1;
};

} // end namespace

#endif // ASAPBRENNER_H
