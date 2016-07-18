// -*- C++ -*-
// AsapMPI.h: Interface to the MPI library.
//
// Copyright (C) 2001-2011 Jakob Schiotz and Center for Individual
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


#ifndef _ASAP_MPI_H
#define _ASAP_MPI_H

#include <mpi.h>
#include <vector>
using std::vector;
#include <string>
using std::string;

#include "Timing.h"

namespace ASAPSPACE {

/// The Communicator provides a simplified interface to the MPI protocol.

/// Technical details: The MPI protocol needs to be accessed through
/// the interface exported by Scientific Python's MPI module.
/// Unfortunately, importing said MPI module's header file from a C++
/// file fails, as it is a C header, causing mpi.h to be imported
/// inside an 'extern "C" {}' statement, which will fail since mpi.h
/// defines a C++ interface to MPI.  The solution is to encapsulate
/// the actual MPI calls in a C file.
class Communicator
{
public:
  inline Communicator() {
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    int ok = 0;
    waiting = false;
    recvwaiting = false;
    MPI_Initialized(&ok);
    assert(ok);
    nProcessors = 0;
    MPI_Comm_size(comm, &nProcessors);
    MPI_Comm_rank(comm, &nProcessor);
  }

  inline ~Communicator() {
    MPI_Comm_free(&comm);
  }
  
  /// Get the number of processors in the simulation.
  inline int GetNumberOfProcessors() const {return nProcessors;}
  /// Get the number of this processor (the "rank" in the MPI communicator).
  inline int GetProcessorNumber() const {return nProcessor;}

  /// Send the buffer to another processor,
  inline void Send(const vector<char> &buffer, int dest) {
    USETIMER("Communicator::Send");
#ifdef USESYNCSEND
    MPI_Ssend((void *) &buffer[0], buffer.size(), MPI_BYTE, dest,
	      7, comm);  
#else
    MPI_Send((void *) &buffer[0], buffer.size(), MPI_BYTE, dest,
	     7, comm);  
#endif /* USESYNCSEND */    
    //asap_mpi_send(&(buffer[0]), buffer.size(), nProcessor);
  }

  /// Send the buffer to another processor, but do not block sending process.
  inline void NonBlockingSend(const vector<char> &buffer, int dest) {
    USETIMER("Communicator::NonBlockingSend");
    assert(!waiting);
#ifdef USESYNCSEND
    MPI_Issend((void *) &buffer[0], buffer.size(), MPI_BYTE, dest, 7,
	       comm, &request);  
#else /* USESYNCSEND */
    MPI_Isend((void *) &buffer[0], buffer.size(), MPI_BYTE, dest, 7,
	      comm, &request);  
#endif /* USESYNCSEND  */
    //asap_mpi_nonblocking_send(&buffer[0], buffer.size(), nProcessor, request);
    waiting = true;
  }

  /// \brief Receive a message, appending it to the buffer (which is
  /// resized accordingly).
  inline void Receive(vector<char> &buffer, int src) {
    USETIMER("Communicator::Receive");
    MPI_Status status;
    int cnt;
    MPI_Probe(src, 7, comm, &status);
    MPI_Get_count(&status, MPI_BYTE, &cnt);
    //int cnt = asap_mpi_probe_and_count(nProcessor);
    int n = buffer.size();
    if (cnt)
      {
	// Got real data, read it.
	buffer.resize(n + cnt);
        MPI_Recv(&buffer[n], cnt, MPI_BYTE, src, 7, comm,
		 MPI_STATUS_IGNORE);
      }
    else
      {
	// Zero length message: Read it but discard it.  Apparently,
	// MPI_Recv must be passed a valid address, the end of
	// buffer[] is not OK.
	char dummy[128];  
	MPI_Recv(dummy, cnt, MPI_BYTE, src, 7, comm,
		 MPI_STATUS_IGNORE);
      }
  }

  /// Receive a message without blocking this process.

  /// Overwrites the buffer, which must be big enough to hold the
  /// message.  The corresponding WaitReceive call shrinks the buffer
  /// to the size of the actual message.
  inline void NonBlockingReceive(vector<char> &buffer, int src) {
    USETIMER("Communicator::NonBlockingReceive");
    assert(!recvwaiting);
    MPI_Irecv(&buffer[0], buffer.size(), MPI_BYTE, src, 7, comm,
	      &recvrequest);
    //asap_mpi_nonblocking_receive(&buffer[0], buffer.size(), nProcessor,
    //				 recvrequest);
    recvbuffer = &buffer;
    recvwaiting = true;
  }

#if 0  // Requires a patch to Scientific.MPI, and is not fastest.
  
  // Simultaneous send and receive.  The receive buffer is overwritten and
  // resized.
  void SendReceive(const vector<char> &sendbuf, int destination,
                   vector<char> &recvbuf, int source);
#endif

  /// \name Reduction operations.
  /// They work as expected. :-)

  //@{
  inline bool LogicalOr(bool boolean) {
    USETIMER("Communicator::LogicalOr");
    int send = (int)boolean;
    int receive;
    MPI_Allreduce(&send, &receive, 1, MPI_INT, MPI_SUM, comm);
    //    asap_mpi_allreduce_sum_int(&send, &receive, 1);
    return (receive > 0);
  }
    
  inline int Min(int x) {
    USETIMER("Communicator::Min");
    int min;
    MPI_Allreduce(&x, &min, 1, MPI_INT, MPI_MIN, comm);
    //    asap_mpi_allreduce_min_double(&x, &min, 1);
    return min;
  }

  inline double Min(double x) {
    USETIMER("Communicator::Min");
    double min;
    MPI_Allreduce(&x, &min, 1, MPI_DOUBLE, MPI_MIN, comm);
    //    asap_mpi_allreduce_min_double(&x, &min, 1);
    return min;
  }

  inline int Max(int x) {
    USETIMER("Communicator::Max");
    int max;
    MPI_Allreduce(&x, &max, 1, MPI_INT, MPI_MAX, comm);
    //    asap_mpi_allreduce_max_double(&x, &max, 1);
    return max;
  }
  
  inline double Max(double x) {
    USETIMER("Communicator::Max");
    double max;
    MPI_Allreduce(&x, &max, 1, MPI_DOUBLE, MPI_MAX, comm);
    //    asap_mpi_allreduce_max_double(&x, &max, 1);
    return max;
  }
  
  inline void Max(vector<int> &x, vector<int> &sum) {
    USETIMER("Communicator::Max(vector<int>)");
    sum.resize(x.size());
    MPI_Allreduce(&x[0], &sum[0], x.size(), MPI_INT, MPI_MAX, comm);
  }
    
  inline double Add(double x) {
    USETIMER("Communicator::Add(double)");
    double sum;
    MPI_Allreduce(&x, &sum, 1, MPI_DOUBLE, MPI_SUM, comm);
    //    asap_mpi_allreduce_sum_double(&x, &sum, 1);
    return sum;
  }

  inline int Add(int x) {
    USETIMER("Communicator::Add(int)");
    int sum;
    MPI_Allreduce(&x, &sum, 1, MPI_INT, MPI_SUM, comm);
    //    asap_mpi_allreduce_sum_int(&x, &sum, 1);
    return sum;
  }
  
  inline long Add(long x) {
    USETIMER("Communicator::Add(long)");
    long sum;
    MPI_Allreduce(&x, &sum, 1, MPI_LONG, MPI_SUM, comm);
    //    asap_mpi_allreduce_sum_long(&x, &sum, 1);
    return sum;
  }
  
  inline void Add(vector<int> &x, vector<int> &sum) {
    USETIMER("Communicator::Add(vector<int>)");
    sum.resize(x.size());
    MPI_Allreduce(&x[0], &sum[0], x.size(), MPI_INT, MPI_SUM, comm);
    //    asap_mpi_allreduce_sum_int(&x[0], &sum[0], x.size());
  }
  
  inline void Add(vector<long> &x, vector<long> &sum) {
    USETIMER("Communicator::Add(vector<long>)");
    sum.resize(x.size());
    MPI_Allreduce(&x[0], &sum[0], x.size(), MPI_LONG, MPI_SUM, comm);
    //    asap_mpi_allreduce_sum_int(&x[0], &sum[0], x.size());
  }
  
  inline void Add(vector<double> &x, vector<double> &sum) {
    USETIMER("Communicator::Add(vector<double>)");
    sum.resize(x.size());
    MPI_Allreduce(&x[0], &sum[0], x.size(), MPI_DOUBLE, MPI_SUM, comm);
    //    asap_mpi_allreduce_sum_double(&x[0], &sum[0], x.size());
  }
  //@}
  
  /// Wait for a nonblocking send to complete.

  ///
  /// When Wait() returns, the send buffer may be overwritten.
  inline void Wait() {
    USETIMER("Communicator::Wait");
    assert(waiting);
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    //    asap_mpi_wait(request);
    waiting = false;
  }

  /// Wait for a nonblocking receive to complete.

  ///
  /// At return, the receive buffer contains the incoming message.
  inline void WaitReceive() {
    USETIMER("Communicator::WaitReceive");
    assert(recvwaiting);
    int cnt;
    MPI_Status status;
    MPI_Wait(&recvrequest, &status);
    MPI_Get_count(&status, MPI_BYTE, &cnt);
    //int cnt = asap_mpi_wait_and_count(request);
    recvbuffer->resize(cnt);
    recvwaiting = false;
  }
  
  /// All to all communication.

  /// Send size integers to each process, receiving as many.  The size
  /// of the sendbuffer must be size * nProcessors, the receive buffer
  /// is resized to this same size.
  inline void AllToAll(vector<int> &sendbuf, vector<int> &recvbuf, int size) {
    USETIMER("Communicator::AllToAll");
    assert(sendbuf.size() == size * nProcessors);
    recvbuf.resize(size * nProcessors);
    MPI_Alltoall(&sendbuf[0], size, MPI_INT, &recvbuf[0], size, MPI_INT,
		 comm);
    //    asap_mpi_all_to_all_int(&sendbuf[0], &recvbuf[0], size);
  }

  /// AllGather.

  /// Send size integers to all the processors, receive as
  /// many from each processor, size*nProcessor in total.  The size of
  /// the send buffer must be size, the receive buffer is resized to
  /// size*nProcessors.
  inline void AllGather(vector<int> &sendbuf, vector<int> &recvbuf, int size) {
    USETIMER("Communicator::AllGather");
    assert(sendbuf.size() == size);
    recvbuf.resize(size * nProcessors);
    MPI_Allgather(&sendbuf[0], size, MPI_INT, &recvbuf[0], size, MPI_INT,
		  comm);
    //    asap_mpi_allgather_int(&sendbuf[0], &recvbuf[0], size);
  }
  
private:
  MPI_Comm comm;     ///< Communicator object
  bool waiting;      ///< Waiting for nonblocking send
  bool recvwaiting;  ///< Waiting for nonblocking receive
  int nProcessor;    ///< The number of this processor
  int nProcessors;   ///< The total number of processors in the simulation.
  MPI_Request request;     ///< For nonblocking send.
  MPI_Request recvrequest; ///< For nonblocking recv.
  vector<char> *recvbuffer; ///< The buffer used in nonblocking recv.
  // string filename;   // Not used?
};

} // end namespace

#endif//  _ASAP_MPI_H
