#!/bin/bash

if [[ $# -lt 1 ]]; then
cat <<EOF
This is a script for compiling Python packages with distutils on Gentoo Linux.

Usage: gentoo-setup.sh commandline

This replaces the usual: python setup.py commandline

Example: To install the package, use
  gentoo-setup.sh install


Rationale:

Due to a buglet in the Gentoo build system, the optimization flags for
Python were not registed in the Python Makefile, and will therefore
not be used when compiling packages.  This script reads the relevant
environment variables from /etc/make.conf before calling setup.py

EOF
exit 1
fi

source /etc/make.conf
export CFLAGS
export CXXFLAGS
python setup.py "$@"
