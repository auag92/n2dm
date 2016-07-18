#! /bin/sh

PYVER=`python -c 'import sys,string; print string.join(string.split(string.split(sys.version, " ")[0], ".")[:2], ".")'`
export PYTHONPATH=${PYTHONPATH}:/usr/lib/python$PYVER/site-packages/Scientific/linux2

epydoc --output epydoc --inheritance included --name Asap --url http://www.fysik.dtu.dk/campos/Asap ../Python/Asap
