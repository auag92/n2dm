// Timing.h  --     -*- C++ -*-
// If compiled with -DTIMING the USETIMER("functionname") macro
// can be used to time functions.  If TIMING is not defined, USETIMER does
// nothing.
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


#ifndef _TIMING_H
#define _TIMING_H

#include "Asap.h"

#ifdef TIMING
#define USETIMER(x) Timing_administrator timer__(x)
#else
#define USETIMER(x)
#endif

#ifdef TIMING
#include <sys/times.h>
#include <unistd.h>
#include <map>
#include <stack>
#include <sstream>
#include <vector>
#include "Exception.h"
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace ASAPSPACE {

class Timing_timer {
public:
  inline Timing_timer(string name);
  inline void Start(Timing_timer *parent = 0);
  inline void Start(clock_t w, clock_t c);
  inline void Stop();
  inline void GoToChild(clock_t w, clock_t c);
  inline void ReturnFromChild(clock_t w, clock_t c);
  inline void Sync();
  inline void Report(string &name, double &wall, double &cpu,
                     double &wall_c, double &cpu_c, int &count);
public:
  string name;
  clock_t cputime;
  clock_t walltime;
  clock_t cputime_c;
  clock_t walltime_c;
  int count;
  int active;
  std::map<string, Timing_timer *> masters;
private:
  clock_t runningcpu;
  clock_t runningwall;
  Timing_timer *parenttimer;
  Timing_timer *activemaster;
};

class Timing_administrator {
public:
  inline Timing_administrator(string name, bool inuse = true);
  inline ~Timing_administrator();
private:
  Timing_timer *timer;
  int threadnum;
  bool inuse;
};
  
// Global variables
#ifdef ALLOCATE_TIMING
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN std::map<string, Timing_timer *> Timing_timerpool;
EXTERN std::vector<std::stack<Timing_timer *> > Timing_timerstack;
EXTERN Timing_timer *Timing_metatimer;
EXTERN clock_t Timing_Resolution;
#undef EXTERN

void Timing_assert_failed(const char *assertion, const char *file,
                          const unsigned int line);

// Inline functions

inline void Timing_init() {
  Timing_timerpool.clear();
#ifdef _OPENMP
  int ncpu = omp_get_num_procs();
  Timing_timerstack.resize(ncpu);
#else
  Timing_timerstack.resize(1);
#endif
  Timing_timer *global = new Timing_timer("Global");
  Timing_timerpool["Global"] = global;
  Timing_timerstack[0].push(global);
  global->Start();
  Timing_Resolution = sysconf(_SC_CLK_TCK);
  Timing_metatimer = new Timing_timer("Timing overhead");
  Timing_timerpool["Timing overhead"] = Timing_metatimer;
}

inline Timing_administrator::Timing_administrator(string name, bool inuse) {
  this->inuse = inuse;
  if (inuse)
    {
#ifdef _OPENMP
      threadnum = omp_get_thread_num();
      char buffer[100];
      if (threadnum == 0)
        Timing_metatimer->Start();
      sprintf(buffer, "-%d", threadnum);
      name += buffer;
#else
      threadnum = 0;
      Timing_metatimer->Start();
#endif
#ifdef _OPENMP
#pragma omp critical
#endif // _OPENMP
      {
        timer = Timing_timerpool[name];
        if (!timer)
          {
            timer = new Timing_timer(name);
            Timing_timerpool[name] = timer;
          }
      }
      Timing_timer *lastt = NULL;
      if (Timing_timerstack[threadnum].size())
        lastt = Timing_timerstack[threadnum].top();
      timer->Start(lastt);
      Timing_timerstack[threadnum].push(timer);
      if (threadnum == 0)
        Timing_metatimer->Stop();
    }
}

inline Timing_administrator::~Timing_administrator() {
  if (inuse)
    {
      if (threadnum == 0)
        Timing_metatimer->Start();
      Timing_timer *lastt = Timing_timerstack[threadnum].top();
      Timing_timerstack[threadnum].pop();
      if (lastt != timer)
	throw AsapError("INCONSISTENT TIMER: Stopping ") << timer->name
		       << " but last timer is " << lastt->name;
      timer->Stop();
      if (threadnum == 0)
        Timing_metatimer->Stop();
    }
}
       
inline void Timing_gettimes(clock_t *wall, clock_t *cpu) {
  struct tms t;
  *wall = times(&t);
  *cpu = t.tms_utime + t.tms_stime;
}

Timing_timer::Timing_timer(string name) : masters()
{
  this->name = name;
  cputime = 0;
  walltime = 0;
  cputime_c = 0;
  walltime_c = 0;
  active = 0;
  parenttimer = 0;
  activemaster = 0;
  count = 0;
}

void Timing_timer::Start(Timing_timer *parent)
{
  assert(active == 0);
  active = 1;
  count++;
  Timing_gettimes(&runningwall, &runningcpu);
  parenttimer = parent;
  if (parent) {
    parent->GoToChild(runningwall, runningcpu);
    activemaster = masters[parent->name];
    if (!activemaster) {
      activemaster = new Timing_timer(parent->name);
      masters[parent->name] = activemaster;
    }
    activemaster->Start(runningwall, runningcpu);
  }
  else
    activemaster = 0;
}

void Timing_timer::Start(clock_t w, clock_t c)
{
  assert(active == 0);
  active = 1;
  count++;
  runningwall = w;
  runningcpu = c;
  parenttimer = 0;
  activemaster = 0;
}

void Timing_timer::Stop()
{
  assert(active == 1);
  active = 0;
  clock_t w, c;
  Timing_gettimes(&w, &c);
  walltime += w - runningwall;
  cputime += c - runningcpu;
  if (parenttimer)
    parenttimer->ReturnFromChild(w, c);
  if (activemaster)
    activemaster->Stop();
}

void Timing_timer::GoToChild(clock_t w, clock_t c)
{
  assert(active == 1);
  active = 2;
  walltime += w - runningwall;
  cputime += c - runningcpu;
  runningwall = w;
  runningcpu = c;
  if (activemaster)
    activemaster->GoToChild(w, c);
}

void Timing_timer::ReturnFromChild(clock_t w, clock_t c)
{
  assert(active == 2);
  active = 1;
  walltime_c += w - runningwall;
  cputime_c += c - runningcpu;
  runningwall = w;
  runningcpu = c;
  if (activemaster)
    activemaster->ReturnFromChild(w, c);
}

void Timing_timer::Sync()
{
  if (active) {
    clock_t w, c;
    Timing_gettimes(&w, &c);
    if (active == 1)
      {
        walltime += w - runningwall;
        cputime += c - runningcpu;
      }
    if (active == 2)
      {
        walltime_c += w - runningwall;
        cputime_c += c - runningcpu;
      }
    runningwall = w;
    runningcpu = c;
  }
  if (activemaster)
    activemaster->Sync();
}

inline void Timing_timer::Report(string &name, double &wall, double &cpu,
                                 double &wall_c, double &cpu_c, int &count)
{
  Sync();
  name = this->name;
  wall = 1.0 * walltime / Timing_Resolution;
  cpu = 1.0 * cputime / Timing_Resolution;
  wall_c = 1.0 * walltime_c / Timing_Resolution;
  cpu_c = 1.0 * cputime_c / Timing_Resolution;
  count = this->count;
}

} // end namespace

#endif // TIMING
#endif // _TIMING_H
