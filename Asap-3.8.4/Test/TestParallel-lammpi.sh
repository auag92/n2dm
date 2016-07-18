#!/bin/sh

lamboot -v localhost
mpirun -np 2 ../`arch`/asap-python TestAll.py --parallel
lamhalt
