#!/bin/sh

mpirun -np 2 `which asap-python` TestAll.py --parallel
