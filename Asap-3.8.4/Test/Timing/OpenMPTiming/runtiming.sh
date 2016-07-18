#!/bin/sh

#PBS -N ForceTiming
#PBS -q medium
#PBS -l nodes=1:xeon:ppn=8

TMPFILE=timing.$$.pickle
NPROC=`cat $PBS_NODEFILE | wc -l`
echo "Running on `hostname` on $NPROC cores."
cat /proc/cpuinfo | grep 'model name' | head -1

#export O64_OMP_SET_AFFINITY=FALSE
export OMP_NUM_THREADS=1
python ForceTiming.py A $TMPFILE > A.txt
mpiexec asap-python ForceTiming.py S $TMPFILE > S.txt
mpiexec asap-python ForceTiming.py M $TMPFILE > M.txt
export OMP_NUM_THREADS=$NPROC
python ForceTiming.py O $TMPFILE > O.txt

echo "" >> result.txt
cat /proc/cpuinfo | grep 'model name' | head -1 >> result.txt
echo "Default library" >> result.txt
python analyse.py $TMPFILE >> result.txt
rm $TMPFILE
numactl --show
printenv


