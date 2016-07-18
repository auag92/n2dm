#!/bin/sh

TMPFILE=timing.$$.pickle
NPROCS=2
echo "Running on `hostname` on $NPROC cores."

export OMP_NUM_THREADS=1
python ForceTiming.py A $TMPFILE
mpiexec -np $NPROCS asap-python ForceTiming.py S $TMPFILE
mpiexec -np $NPROCS asap-python ForceTiming.py M $TMPFILE
export OMP_NUM_THREADS=$NPROCS
python ForceTiming.py O $TMPFILE

python analyse.py $TMPFILE >> result.txt
rm $TMPFILE


