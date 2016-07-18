#!/bin/bash
#PBS -N compare
#PBS -q medium
#PBS -l nodes=1:ppn=4:opteron4

python OfficialTiming.py
asap-python OfficialTiming.py -t

export PYTHONPATH=`echo $PYTHONPATH | sed 's/development/production/g'`
export PATH=`echo $PATH | sed 's/development/production/g'`

python OfficialTiming.py
asap-python OfficialTiming.py -t

