#include "RandomNumbers.h"
#include "Exception.h"
#include <cfloat>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

// http://www.drdobbs.com/tools/fast-high-quality-parallel-random-number/229625477?pgno=1

#define rotl(r,n) (((r)<<(n)) | ((r)>>((8*sizeof(r))-(n))))

AsapRandom::AsapRandom(npy_uint32 seed)
{
  SeedRersResrResdra(seed);
  npy_uint64 rand_max = (1ULL << DBL_MANT_DIG);
  rand_mask = rand_max - 1;
  rand_factor = 1.0 / rand_mask;
}

void AsapRandom::SeedRersResrResdra(npy_uint32 seed) {
  unsigned n;
  xx =    914489ULL;
  yy =   8675416ULL;
  zz = 439754684ULL;
  for (n=((seed>>22)&0x3ff)+20; n>0; n--) { xx = rotl(xx,8) - rotl(xx,29); }
  for (n=((seed>>11)&0x7ff)+20; n>0; n--) { yy = rotl(yy,21) - yy;  yy = rotl(yy,20); }
  for (n=((seed    )&0x7ff)+20; n>0; n--) { zz = rotl(zz,42) - zz;  zz = rotl(zz,14) + zz; }
}

npy_uint64 AsapRandom::RandRersResrResdra() {  // Combined period = 2^116.23
  xx = rotl(xx,8) - rotl(xx,29);                 //RERS,   period = 4758085248529 (prime)
  yy = rotl(yy,21) - yy;  yy = rotl(yy,20);      //RESR,   period = 3841428396121 (prime)
  zz = rotl(zz,42) - zz;  zz = zz + rotl(zz,14); //RESDRA, period = 5345004409 (prime)
  return xx ^ yy ^ zz;
}

void AsapRandom::RandomDoubles(double *p, int n)
{
  npy_uint64 xx = this->xx;
  npy_uint64 yy = this->yy;
  npy_uint64 zz = this->zz;
  for (double *pp = p; pp < p + n; pp++)
    {
      xx = rotl(xx,8) - rotl(xx,29);                 //RERS,   period = 4758085248529 (prime)
      yy = rotl(yy,21) - yy;  yy = rotl(yy,20);      //RESR,   period = 3841428396121 (prime)
      zz = rotl(zz,42) - zz;  zz = zz + rotl(zz,14); //RESDRA, period = 5345004409 (prime)
      npy_uint64 tmp = xx ^ yy ^ zz;
      *pp = rand_factor * (tmp & rand_mask);
    }
  this->xx = xx;
  this->yy = yy;
  this->zz = zz;
}

AsapRandomThread::AsapRandomThread(npy_uint32 seed)
{
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#else
  nthreads = 1;
#endif
  // Generate seeds for the RNGs.  Even if only one generator is needed,
  // it is not good to seed it with the user-provided seed, as it is often
  // a small number with many zero-bits, leading to parts of the RNG being
  // seeded with 0.
  AsapRandom seeder(seed);
  generators = new AsapRandom*[nthreads];
  int threadnum = 0;
  npy_uint32 s;
#ifdef _OPENMP
#pragma omp parallel private(threadnum, s)
  {
    threadnum = omp_get_thread_num();
#pragma omp critical
#endif
    {
      s = seeder.RandomInt64() & ((1ULL << 32) - 1);
      generators[threadnum] = new AsapRandom(s);
    }
#ifdef _OPENMP
  }
#endif
}

void AsapRandomThread::RandomDoubles(double *p, int n)
{
  if ((nthreads == 1) || (n < 250))
    {
      // It does not pay to multi-thread
      generators[0]->RandomDoubles(p, n);
      return;
    }
  else
    {
#ifdef _OPENMP
      int threadnum;
      int numthreads;
#pragma omp parallel private(threadnum, numthreads)
      {
        threadnum = omp_get_thread_num();
        numthreads = omp_get_num_threads();
        int blocksize = n / numthreads;
        int lastblock = n - (numthreads - 1) * blocksize;
        //std::cerr << "nthr=" << numthreads << " thr#=" << threadnum << " blocksize=" << blocksize << " lastblock=" << lastblock << std::endl;
        if (threadnum < numthreads - 1)
          generators[threadnum]->RandomDoubles(p + threadnum * blocksize, blocksize);
        else
          generators[threadnum]->RandomDoubles(p + threadnum * blocksize, lastblock);
      }
#if 0
      // Find number of blocks
      int nblocks = (int) floor(sqrt(n));
      int blocksize = n / nblocks;
      int lastblocksize = n % nblocks;
      int threadnum;
      if (lastblocksize < 5)
        {
          // Don't want a tiny block, include in previous.
          nblocks--;
          lastblocksize += blocksize;
        }

#pragma omp parallel private(threadnum)
      {
        threadnum = omp_get_thread_num();
#pragma omp for
        for (int i = 0; i <= nblocks; i++)
          {
            int size = blocksize;
            if (i == nblocks)
              size = lastblocksize;
            generators[threadnum].RandomDoubles(p + i * blocksize, size);
          }
      }
#endif // 0
#else // _OPENMP
      throw AsapError("Multithreading not supported in the RNG.");
#endif // _OPENMP
    }
}

