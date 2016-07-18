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

#define PARALLEL
#include "mpimodule.h"
#include "ParallelAtomsInterface.h"
#include "ParallelPotentialInterface.h"

// Now include the main code
#include "AsapModule.cpp"

//#include <mpi.h>

int main(int argc, char **argv)
{
  int status;
  
#ifdef _OPENMP
  int provided = 0;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
#if 0
  if (provided < MPI_THREAD_FUNNELED)
    {
      std::cerr << "Combined MPI / OpenMP parallelization not supported by MPI." << std::endl;
      std::cerr << "provided: " << provided << " requested: " << MPI_THREAD_FUNNELED << std::endl;
      return -1;
    }
#endif
#else
  MPI_Init(&argc, &argv);
#endif
  Py_Initialize();
  initasapparallel3();
  status = Py_Main(argc, argv);
  MPI_Finalize();
  return status;
}

  
