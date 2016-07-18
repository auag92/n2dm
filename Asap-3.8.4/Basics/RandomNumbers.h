
#ifndef _RANDOMNUMBERS_H
#define _RANDOMNUMBERS_H

#include "AsapPython.h"
#include "Asap.h"
#include <vector>

namespace ASAPSPACE {

class AsapRandom
{
public:
  AsapRandom(npy_uint32 seed);

  void RandomDoubles(double *p, int n);
  npy_uint64 RandomInt64() {return RandRersResrResdra();} // More useful name.

private:
  void SeedRersResrResdra(npy_uint32 seed);
  npy_uint64 RandRersResrResdra();

private:
  npy_uint64 xx, yy, zz;
  npy_uint64 rand_mask;
  double rand_factor;
};

class AsapRandomThread
{
public:
  AsapRandomThread(npy_uint32 seed);

  void RandomDoubles(double *p, int n);

private:
  AsapRandom **generators;
  int nthreads;
};

} // end namespace

#endif // _RANDOMNUMBERS_H
