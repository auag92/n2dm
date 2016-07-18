#!/bin/bash
#PBS -N comparepar
#PBS -q medium
#PBS -l nodes=1:ppn=8:xeon8

mpirun --mca mpi_paffinity_alone 1 asap-python ParallelTiming.py

export PYTHONPATH=`echo $PYTHONPATH | sed 's/development/production/g'`
export PATH=`echo $PATH | sed 's/development/production/g'`

mpirun --mca mpi_paffinity_alone 1 asap-python ParallelTiming.py

