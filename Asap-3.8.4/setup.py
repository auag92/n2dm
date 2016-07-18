#!/usr/bin/env python

from distutils.core import setup, Extension
from glob import glob
import re
import os
from os.path import join
import sys
import numpy

folders = ['Basics', 'Interface', 'Brenner', 'Tools']
kim_folders = ['OpenKIMimport']
exclude_files = ['Basics/MoPotential.cpp', 'Interface/AsapModule.cpp']

# Get the version number from Python/Asap/version.py
locs = {}
execfile('Python/asap3/version.py', {}, locs)
version = locs['__version__']
print "Asap version number:", version

def runcmd(cmd):
    x = os.system(cmd)
    if x:
        raise RuntimeError, "Command failed: "+cmd

if not os.path.exists("Distutils-autogen"):
    os.mkdir("Distutils-autogen")
runcmd("python recordversion.py '"+version+"' 'serial' 'distutils' '' > Distutils-autogen/version.cpp")

if 'ASAP_KIM_DIR' in os.environ:
    akd = os.environ['ASAP_KIM_DIR']
    kim_inc = os.path.join(akd, 'src')
    kim_lib = kim_inc
    folders.extend(kim_folders)
elif 'KIM_HOME' in os.environ:
    kh = os.environ['KIM_HOME']
    kim_inc = os.path.join(kh, 'include', 'kim-api-v1')
    kim_lib = os.path.join(kh, 'lib')
    folders.extend(kim_folders)
else:
    kim_inc = kim_lib = None
        
incl_dirs = folders[:]
incl_dirs.append(numpy.get_include())
if kim_inc:
    incl_dirs.append(kim_inc)
    
if kim_lib:
    kim_args = {'library_dirs': [kim_lib], 
                'libraries': ['kim-api-v1'], 
                'define_macros': [('WITH_OPENKIM', '1')],
                }
    msg = "OpenKIM support enabled."
else:
    kim_args = {}
    msg = "OpenKIM support not enabled."
    
print "kim_inc", kim_inc
print "kim_lib", kim_lib
print "kim_args", kim_args

print "Reading source files"
src_files = []
for d in folders:
    for f in os.listdir(d):
        if f.endswith('.cpp'):
            fn = os.path.join(d,f)
            if fn not in exclude_files:
                src_files.append(fn)
                print "  ", fn
src_files.append("Distutils-autogen/version.cpp")

packages = ['asap3',
            'asap3.analysis',
            'asap3.md',
            'asap3.optimize',
            'asap3.io',
            'asap3.Setup',
            'asap3.Utilities',
            'asap3.Internal',
            'asap3.MonteCarlo',
            'asap3.Tools',
            'asap3.nanoparticle_mc',
            ]

setup(name="asap3",
      version=version,
      description="Atomic SimulAtion Program - As Soon As Possible",
      author="Jakob Schiotz et. al.",
      author_email="schiotz@fysik.dtu.dk",
      url="http://www.fysik.dtu.dk/Software",
      packages=packages,
      package_dir={'asap3': 'Python/asap3'},
      ext_modules=[Extension('asapserial3',
                             src_files,
                             include_dirs=incl_dirs,
                             **kim_args
                             )]
      )

print ""
print msg
