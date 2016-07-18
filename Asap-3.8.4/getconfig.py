#!/usr/bin/env python

"""Get Python configuration variables.

Usage: python getconfig.py VARIABLE

Used to extract variables from Python's Makefile, intended to be
called from Asap's Makefile.
"""

import sys
import distutils.sysconfig

if len(sys.argv) != 2:
    print >>sys.stderr, "\nERROR: Got %s arguments, expected 1.\n\n" % (len(sys.argv)-1,)
    print >>sys.stderr,  __doc__
    sys.exit(-1)

key = sys.argv[1]
if key == "SITEPACKAGES":
    print distutils.sysconfig.get_python_lib(plat_specific=True)
else:
    cfgDict = distutils.sysconfig.get_config_vars()
    print cfgDict[key]
