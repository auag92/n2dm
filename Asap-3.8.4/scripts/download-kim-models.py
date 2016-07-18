#!/usr/bin/env python

"""Download all missing KIM Models from openkim.org.

Usage: 
  download-kim-models.py kim-api-dir
or
  download-kim-models.py .
to use the current directory (note the dot!).
"""

import subprocess
import sys
import os
import json

if len(sys.argv) != 2:
    print __doc__
    raise RuntimeError("Invalid number of arguments, got %i expected 1" % (len(sys.argv)-1,))
                       
kimdir = sys.argv[1]

if os.path.isdir(kimdir):
    print "Working directory: ", kimdir
    os.chdir(kimdir)

if not os.path.isdir('src/models') or not os.path.isdir('src/model_drivers'):
    raise RuntimeError("Folder '%s' does not look like the KIM project" % (kimdir,))

def getdirs(d):
    return dict.fromkeys([x for x in os.listdir(d) if os.path.isdir(os.path.join(d,x))], True)

ok_models = getdirs('src/models')
ok_drivers = getdirs('src/model_drivers')

print "Found %d models and %d model drivers" % (len(ok_models), len(ok_drivers))
print ""
print "Downloading list of OpenKIM models."
cmd = """curl -d 'fields={"kimcode":1,"kim-api-version":1}' -d 'database=obj' -d 'query={"type":"mo"}' https://query.openkim.org/api"""
print "Running", cmd
pipe = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout
candidates_json = pipe.read()
pipe.close()
candidates_raw = json.loads(candidates_json)
candidates = [x['kimcode'] for x in candidates_raw if x['kim-api-version'] >= "1.6"]
print "Found %d relevant KIM models on the website" % (len(candidates))
to_dl = [x for x in candidates if x not in ok_models]
print "Found %d KIM models to download" % (len(to_dl))
for m in to_dl:
    cmd = ["make",  "add-%s" % (m,)]  # Don't pass downloaded data to the shell !
    print " ".join(cmd)
    subprocess.check_call(cmd)

