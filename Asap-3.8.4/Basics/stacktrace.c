// -*- C++ -*-
//
// stacktrace.cpp: Provide a stack trace when Asap crashes.
//
// Copyright (C) 2002-2011 Jakob Schiotz and Center for Individual
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

#if defined(STACKTRACELINUX) || defined(STACKTRACEBFD) || defined(STACKTRACEALPHA) || defined(STACKTRACEGDB)
#define STACKTRACE
#endif

#include "stacktrace.h"

#ifdef STACKTRACE
#include <Python.h>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#ifdef STACKTRACEALPHA
#define USE_LADEBUG
#define USE_DEBUGGER
#endif

#ifdef STACKTRACEGDB
#define USE_DEBUGGER
#define USE_GDB
#endif

#ifdef STACKTRACEBFD
#define STACKTRACELINUX
#define BFD_BACKTRACE
#endif

#ifdef STACKTRACELINUX
#include <execinfo.h>
#ifndef __USE_GNU
#define __USE_GNU
#endif
// Note: Use sys/ucontext.h instead of ucontext.h for Max OS X compatibility.
#include <sys/ucontext.h>
#undef USE_DEBUGGER
/* Define HACK_STACK if glibc ruins the stack frame when a signal is called. */
#undef HACK_STACK
#define STACK_OFFSET 2
/* Define GLIBC_BACKTRACE to use glibc's backtrace_symbols to produce
   the stack trace, or BFD_BACKTRACE to use libbfd's backtrace
   (requires linking with -lbfd, produces more readable backtraces,
   but may fail if the heap is corrupted).  If both are enabled, both
   stack traces are produced. */
#define GLIBC_BACKTRACE
#endif /* STACKTRACE == linux */

#ifdef BFD_BACKTRACE
#include <dlfcn.h>
#include <bfd.h>

static void resolve_address(FILE *outfile, char *address);
#endif /* BFD_BACKTRACE */

static char mpipython[] = "asap-python";
static char *execname;
static char crashfilename[100];
static abortfunc *abortfunction;

static void sig_handler(int sig, siginfo_t *sip, void *extra);
static void sethandler(int signal);

void Asap_setSignalHandlers(int node, abortfunc *ab)
{
  abortfunction = ab;
  if (node >= 0)
    {
      sprintf(crashfilename, "CRASH.%d", node);
      execname = mpipython;  /* Py_GetProgramFullPath() is not realiable */
    }
  else
    {
      strcpy(crashfilename, "CRASH");
      execname = Py_GetProgramFullPath();
    }

  /* Set the signal handlers */
  struct sigaction sa;
  sa.sa_sigaction = sig_handler;
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = SA_RESTART | SA_SIGINFO;
  sigaction(SIGBUS, &sa, NULL);
  sigaction(SIGILL, &sa, NULL);
  sigaction(SIGFPE, &sa, NULL);
  sigaction(SIGSEGV, &sa, NULL);
  sigaction(SIGUSR1, &sa, NULL);
}

static void sethandler(int signal)
{
  struct sigaction sa;
  int ret;
  /* sa.sa_handler = sig_handler; */
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = SA_SIGINFO;
  ret = sigaction(signal, &sa, 0);
  if (ret) {
    perror("sigaction");
    exit(1);
  }
}

static void sig_handler(int sig, siginfo_t *sip, void *extra)
{
  char pid[100];
  int fpid;
  FILE *debug;
  FILE *crash;
  char debugname[100];
  char *name;
  char *r;
  
  /* Find which signal occurred */
  switch (sig) {
  case SIGBUS: name = "SIGBUS: Bus error"; break;
  case SIGCHLD: name = "SIGCHLD: Child process event"; break;
  case SIGILL: name = "SIGILL: Illegal operation"; break;
  case SIGFPE: name = "SIGFPE: Arithmetic error"; break;
  case SIGSEGV: name = "SIGSEGV: Segmentation fault"; break;
  default: name = "Unknown signal"; 
  }

  /* Print which signal occurred. */
  fprintf(stderr, "\n\nASAP signal handler: An unexpected error occurred.\n");
  fprintf(stderr, "Process %s got signal %d (%s).\n  si_code = %d.  si_addr = 0x%lx\n", execname, sig, name, sip->si_code, (long) sip->si_addr);

  crash = fopen(crashfilename, "w");
  fprintf(crash, "\nProcess %s got signal %d (%s).\n  si_code = %d.  si_addr = 0x%lx\n", execname, sig, name, sip->si_code, (long) sip->si_addr);
  

  /* Try to detect the reason for the signal, starting with generic reasons...*/
  switch(sip->si_code) {
  case SI_USER: r = "sent by kill or raise"; break;
#ifdef SI_KERNEL
  case SI_KERNEL: r = "sent by kernel"; break;
#endif
  default: r = "Unknown reason";
  }

  /* ... then look for signal-dependent reasons. */
  switch (sig) {
  case SIGBUS:
    switch(sip->si_code) {
    case BUS_ADRALN: r = "invalid address alignment"; break;
    case BUS_ADRERR: r = "non-existent physical address"; break;
    case BUS_OBJERR: r = "object specific hardware error"; break;
    }
    break;
  case SIGCHLD:
    switch(sip->si_code) {
    case CLD_EXITED: r ="child has exited"; break;
    case CLD_KILLED: r = "child was killed"; break;
    case CLD_DUMPED: r = "child terminated abnormally"; break;
    case CLD_TRAPPED: r = "traced child has trapped"; break;
    case CLD_STOPPED: r = "child has stopped"; break;
    case CLD_CONTINUED: r = "stopped child has continued"; break;
#ifdef STACKTRACEALPHA
    case CLD_SIGEXITING: r = "child is about to exit because it received a fatal signal"; break;
#endif
    }
    break;
  case SIGILL:
    switch(sip->si_code) {
    case ILL_ILLOPC: r = "illegal opcode"; break;
    case ILL_ILLOPN: r = "illegal operand"; break;
    case ILL_ILLADR: r = "illegal addressing mode"; break;
    case ILL_ILLTRP: r = "illegal trap"; break;
    case ILL_PRVOPC: r = "privileged opcode"; break;
    case ILL_PRVREG: r = "privileged register"; break;
    case ILL_COPROC: r = "coprocessor error"; break;
    case ILL_BADSTK: r = "internal stack error"; break;
    }
    break;
  case SIGFPE:
    switch(sip->si_code) {
    case FPE_INTDIV: r = "integer divide by zero"; break;
    case FPE_INTOVF: r = "integer overflow"; break;
    case FPE_FLTDIV: r = "floating point divide by zero"; break;
    case FPE_FLTOVF: r = "floating point overflow"; break;
    case FPE_FLTUND: r = "floating point underflow"; break;
    case FPE_FLTRES: r = "floating point inexact result"; break;
    case FPE_FLTINV: r = "invalid floating point operation"; break;
    case FPE_FLTSUB: r = "subscript out of range"; break;
    }
    break;
  case SIGSEGV:
    switch(sip->si_code) {
    case SEGV_MAPERR: r = "address not mapped to object"; break;
    case SEGV_ACCERR: r = "invalid permissions for mapped object"; break;
    }
    break;
  }
  fprintf(stderr, "Cause of signal: %s\n\n", r);
  fprintf(crash, "Cause of signal: %s\n\n", r);
  fprintf(stderr, "\nA traceback will be written to the file %s\n",
	  crashfilename);

  fprintf(stderr, "Please email this file to the programmers for debugging usage.\n\n");
  fflush(NULL);  /* Flush all file buffers */

#ifdef STACKTRACELINUX
  static void *trace[50];
  int trace_size = 0;
  
#ifdef GLIBC_BACKTRACE
  /*  This section generates a stack trace using the GLIBC functions
   *  backtrace() and backtrace_symbols_fd(), it is thus not generally
   *  portable but will only work on glibc systems.
   *
   *  HACK_STACK: On some systems, calling the signal handler
   *  overwrites the most important stack frame: the one creating the
   *  error.  To get a useful stack trace, that stack frame is
   *  restored using an UNDOCUMENTED and PROCESSOR-SPECIFIC hack.
   *  This will only work on Intel/Linux systems (the i686 and x86_64
   *  architectures).  Similar hacks may be possible on other
   *  architectures.  THIS HACK DOES NOT SEEM TO BE NEEDED!
   */

#ifdef HACK_STACK
  ucontext_t *uc = (ucontext_t *) extra;
#endif
  /* Get the backtrace as adresses */
  trace_size = backtrace(trace, 50);
#ifdef HACK_STACK
  /* Overwrite sigaction stack frame with caller's address */
#if __WORDSIZE == 64
  trace[1] = (void *) uc->uc_mcontext.gregs[REG_RIP];
#else /* __WORDSIZE == 32 */
  trace[1] = (void *) uc->uc_mcontext.gregs[REG_EIP];
#endif
#endif
  /* Print stack trace to crash file. */
  fprintf(crash, "\n\nStack trace (glibc):\n");
  fprintf(crash, "   If error is not in a .so file, addresses in [] can be translated\n");
  fprintf(crash, "   to line numbers with the command\n");
  fprintf(crash, "   addr2line -e EXECUTABLE ADDRESS\n\n");
  
  fflush(NULL);
  backtrace_symbols_fd(trace+STACK_OFFSET, trace_size-STACK_OFFSET, fileno(crash));
  fflush(NULL);
  
#endif /* GLIBC_BACKTRACE */
#ifdef BFD_BACKTRACE
  /*  This section generates a stack trace using the GLIBC function
   *  backtrace() and the functions in libbfd.
   */

  fprintf(crash, "\n\nStack trace (BFD):\n");
  fflush(NULL);
   /* Get the backtrace as adresses */
  trace_size = backtrace(trace, 50);
  void *approx_text_end = (void*) ((128+100) * 2<<20);
  int i;
  for (i = STACK_OFFSET; i < trace_size; i++)
    {
      backtrace_symbols_fd(trace+i, 1, fileno(crash));
      fflush(NULL);
      if (trace[i] < approx_text_end)
	resolve_address(crash, trace[i]);
      fflush(NULL);
    }
  
#endif /* BFD_BACKTRACE */
  
  /* Exit */
  if (abortfunction) {
    fprintf(stderr, "Exiting using abort function (MPI)....\n");
    abortfunction();
  } else {
    fprintf(stderr, "Exiting using exit....\n");
    exit(3);
  }
#endif /* STACKTRACE == linux */



  
#ifdef USE_DEBUGGER
  /*  This section is LEGACY CODE, left to be revived when backtrace.h
   *  methods cannot be used to generate a stack trace.  Instead, an
   *  external debugger is invoked, it appears to be less reliable.
   */
     
  sprintf(pid, "%d", getpid());
#ifdef USE_LADEBUG
  sprintf(debugname, "/tmp/dump.ladebug.%d", getpid());
#else
  sprintf(debugname, "/tmp/dump.gdb.%d", getpid());
#endif
  debug = fopen(debugname, "w");
  if (debug == 0)
    {
      perror("fopen failed");
      exit(3);
    }
  fprintf(debug, "where\ndetach\nquit\n");
  fclose(debug);
  fflush(crash);
  fpid = fork();
  if (fpid == -1)
    {
      perror("fork failed");
      exit(1);
    }
  if (fpid == 0)
    {
      /* The child */
      dup2(fileno(crash), fileno(stderr));
      dup2(fileno(crash), fileno(stdout));
#if defined(USE_LADEBUG)
      execlp("ladebug", "ladebug", "-pid", pid, "-c", debugname, execname,
	     (char *) 0);
#elif defined(USE_IDB_GDB)
      execlp("idb", "idb", "-gdb", "-batch", "-pid", pid, "-command", debugname, execname,
	     (char *) 0);
#else
      fprintf(stderr, "Exec'ing gdb...\n");
      fflush(stderr);
      execlp("gdb", "gdb", "-x", debugname, execname, pid, (char *) 0);
#endif
    }
  else
    {
      sleep(5);
#ifdef USE_LADEBUG
      fprintf(stderr, "Transferring control to debugger using SIGINT\n");
      kill(getpid(), SIGINT);
#endif
      sleep(3);
      unlink(debugname);
      if (abortfunction) {
	fprintf(stderr, "Exiting using abort function (MPI)....\n");
	abortfunction();
      } else {
	fprintf(stderr, "Exiting using exit....\n");
	exit(3);
      }
    }
#endif /* USE_DEBUGGER */


  /* Common fallback section. */
  fprintf(stderr, "\n\nReached end of signal handler: This should not happen!!!!\n");
  exit(42);
}


#ifdef BFD_BACKTRACE
/* globals retained across calls to resolve_address. */
static bfd* abfd = 0;
static asymbol **syms = 0;
static asection *text = 0;

static void resolve_address(FILE *outfile, char *address) {
  if (!abfd)
    {
      char ename[1024];
      int l = readlink("/proc/self/exe",ename,sizeof(ename));
      if (l == -1)
	{
	  perror("failed to find executable\n");
	  return;
 	}
      ename[l] = 0;
      
      bfd_init();
      
      abfd = bfd_openr(ename, 0);
      if (!abfd)
	{
	  perror("bfd_openr failed: ");
	  return;
 	}
      /* oddly, this is required for it to work... */
      bfd_check_format(abfd,bfd_object);
      
      unsigned storage_needed = bfd_get_symtab_upper_bound(abfd);
      syms = (asymbol **) malloc(storage_needed);
      unsigned cSymbols = bfd_canonicalize_symtab(abfd, syms);
      
      text = bfd_get_section_by_name(abfd, ".text");
    }
  long offset = ((long)address) - text->vma;
  if (offset > 0)
    {
      const char *file;
      const char *func;
      unsigned line;
      if (bfd_find_nearest_line(abfd, text, syms, offset, &file, &func, &line)
	  && file)
	{
	  fprintf(outfile, "   %s    +%u %s\n", func, line, file);
	}
    }
}
#endif /* BFD_BACKTRACE */

#endif /* STACKTRACE */
