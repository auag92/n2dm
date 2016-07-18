#! /bin/sh

PYVER=`python -c 'import sys,string; print string.join(string.split(string.split(sys.version, " ")[0], ".")[:2], ".")'`

export PYTHONPATH=${PYTHONPATH}:/usr/lib/python$PYVER/site-packages/Scientific/linux2

set -e

EPYDOC=`which epydoc`
echo "Using this epydoc: $EPYDOC"
MPIPYT=`which mpipython`
echo "Using this mpipython: $MPIPYT"
rm -rf epydoc
lamboot -v epydoc.host
mpirun -O -np 1 $MPIPYT $EPYDOC --output epydoc --inheritance included --name Asap --url http://www.fysik.dtu.dk/campos/Asap ../Python/Asap
lamhalt
PUBDIR=/home/camp2/schiotz/WWW/comp/Asap/epydoc
echo "Copying files to $PUBDIR"
rm -rf ${PUBDIR}.old
[[ -d ${PUBDIR} ]] && mv ${PUBDIR} ${PUBDIR}.old
cp -pr epydoc ${PUBDIR}


