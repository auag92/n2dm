#include "expand.h"
#include "Exception.h"

/* Set values in ti such that ti[j] is the coefficient of x**j in the expansion
   of (x+(L-1))**i.  Only works for i < 4.
*/

namespace ASAPSPACE {

void expand (int L, int ti [4], int i) {
  int j;
  for (j = 0; j < 4; j++) ti[j] = 0;
  switch (i) {
  case 0:
    ti[0] = 1;
    break;
  case 1:
    ti[0] = L;
    ti[1] = 1;
    break;
  case 2:
    ti [0] = L * L;
    ti [1] = 2 * L;
    ti [2] = 1;
    break;
  case 3:
    ti [0] = L * L * L;
    ti [1] = 3 * L * L;
    ti [2] = 3 * L;
    ti [3] = 1;
    break;
  default:
    assert (0 && "i should be between 0 and 3 inclusive");
  }
}

} // end namespace
