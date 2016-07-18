"""Export the files needed by the KIM api."""

import os
import glob
import shutil
import string
from copy import copy
import subprocess
if not hasattr(subprocess, 'check_output'):
    # Pre-2.7 hack
    print "WARNING: Fixing pre Python 2.7 compatibility"
    from subprocess_py27 import check_output
    subprocess.check_output = check_output
    
# List all models to be exported, and the files that are specific to that model.
# The keys are the models, the values the list of files.  If the list of files is
# an empty list, it is assumed that only the common files are to be exported.  If
# the list of files is False, no files except the directory is copied, this is for
# models depending on a model driver.

emtfilelist = ['EMT.cpp', 'EMTDefaultParameterProvider.cpp', 'EMTParameterProvider.h']
model_files = {
               'EMT_Asap__MD_128315414717_002': emtfilelist,     # Model driver
               'EMT_Asap_Standard_Jacobsen_Stoltze_Norskov_AlAgAuCuNiPdPt__MO_118428466217_002': False, # Based on EMT_Asap driver
               'EMT_Asap_MetalGlass_CuMgZr__MO_655725647552_002': False,       # Based on EMT_Asap driver
               }

# List all model drivers, so the script can distinguish between model drivers and models.
all_drivers = ['EMT_Asap__MD_128315414717_002']

# List all files from the Basics directory.  If a .cpp file is given, the corresponding
# .h file does not need to be specified, it is added automatically.
basic_files = ['NeighborLocator.h', 'NeighborCellLocator.cpp', 'Asap.h',
               'AsapObject.cpp', 'Vec.cpp', 'TinyMatrix.h', 'Potential.cpp', 'SymTensor.h', 
               'Exception.cpp', 'IVec.h', 'Timing.cpp', 'mass.h', 'Debug.h',
               'Matrix3x3.cpp', 'TimingResults.h']

fake_headers = ['AsapPython.h', 'Atoms.h', 'Templates.h']

fake_header_template = """
// Fake header replacing an Asap header file.

#include "Kim%s"
"""

svnrevision = None
asapversion = None

def main():
    global svnrevision, asapversion
    svnrevision = subprocess.check_output('svnversion .', shell=True).strip()
    asapversion = subprocess.check_output('python ../Python/asap3/version.py', shell=True).strip()
    models = model_files.keys()
    try:
        kimdir = os.environ['ASAP_KIM_DIR']
    except KeyError:
        if os.path.isdir('../../kim-api'):
            kimdir = '../../kim-api'
        else:
            raise RuntimeError("Cannot locate the openkim-api directory - please set $ASAP_KIM_DIR.")
    kimmodeldir = os.path.join(kimdir, 'src', 'models')
    kimdriverdir = os.path.join(kimdir, 'src', 'model_drivers')
    # Add header files
    local_files = glob.glob('*.h') + glob.glob('*.cpp')
    for model in models:
        isdriver = model in all_drivers
        dependent_model = model_files[model] is False
        if dependent_model:
            my_basic_files = copy(basic_files)
        else:
            my_basic_files = basic_files + model_files[model]
        hdr = []
        for f in my_basic_files:
            if f.endswith('.cpp'):
                hdr.append(f[:-4]+'.h')
        my_basic_files.extend(hdr)
        if isdriver:
            model_dir = os.path.join(kimdriverdir, model)
        else:
            model_dir = os.path.join(kimmodeldir, model)
        print "Populating", model_dir
        if not os.path.exists(model_dir):
            os.mkdir(model_dir)
        populate(model_dir, model, os.listdir(model))
        if not dependent_model:
            populate(model_dir, '../Basics', my_basic_files)
            populate(model_dir, '.', local_files)
            for fake in fake_headers:
                fn = os.path.join(model_dir, fake)
                if not os.path.exists(fn):
                    print "Creating", fn
                    f = open(fn, "w")
                    f.write(fake_header_template % (fake,))
                    f.close()
            shutil.copy('Makefiletemplate.mk', os.path.join(model_dir, 'Makefile'))
            # Create model definition include file
            incfilename = os.path.join(model_dir, 'Modeldefinition.mk')
            incfile = open(incfilename, 'w')
            if isdriver:
                incfile.write('MODEL_DRIVER_NAME := %s\n' % (model,))
                # There should be exactly one *.kim.tpl file
                template = glob.glob('%s/*.kim.tpl' % (model,))
                if len(template) != 1:
                    raise RuntimeError('Expected one template file, found %i: %s' % (len(template), str(template)))
                incfile.write('MODEL_DRIVER_KIM_FILE_TEMPLATE := %s\n' % (os.path.basename(template[0]),))
                incfile.write('MODEL_DRIVER_INIT_FUNCTION_NAME := model_driver_init\n')
            else:
                incfile.write('MODEL_NAME := %s\n' % (model,))
            allfiles = os.listdir(model_dir)
            cppfiles = []
            objfiles = []
            for f in allfiles:
                if f.endswith('.cpp'):
                    cppfiles.append(f)
            incfile.write('MODEL_SOURCES := %s\n' % (' '.join(cppfiles)),)
            incfile.close()
            # Make the file with the dependencies
            cmd = '''cd %s && bash -c "g++ -MM -I'%s/src' *.cpp > Depend.tmp"''' % (model_dir, kimdir)
            print cmd
            subprocess.check_call(cmd, shell=True)
            infile = os.path.join(model_dir, 'Depend.tmp')
            outfile = os.path.join(model_dir, 'Depend.mk')
            out = open(outfile, "w")
            for line in open(infile):
                words = line.split(' ')
                words = [x for x in words if '/' not in x]
                keep = False
                for w in words:
                    if w != '' and w != '\\' and w != '\\\n' and w != '\n':
                        keep = True
                if keep:
                    out.write(" ".join(words))
            out.close()
            os.remove(infile)
    print "\nAsap svn version:", svnrevision
    if 'M' in svnrevision or ':' in svnrevision:
        print "\nWARNING: SVN is not up to date (locally modified files or mixed revision)."
        print "This is OK for daily use, but not when submitting to OpenKIM.org"
       
def populate(targetdir, sourcedir, files):
    for f in files:
        if f == '.svn':
            continue            
        sourcefile = os.path.join(sourcedir, f)
        targetfile = os.path.join(targetdir, f)
        if f == 'README':
            # README file is a template that should always be updated
            readme = string.Template(open(sourcefile).read())
            ff = open(targetfile, 'w')
            ff.write(readme.substitute(svnrevision=svnrevision, version=asapversion))
            ff.close()
            print "Updated (and substituted)", targetfile        
        else:
            if not os.path.exists(targetfile):
                newer = True
            else:
                newer = (os.stat(sourcefile).st_mtime > os.stat(targetfile).st_mtime)
            if newer:
                shutil.copy(os.path.join(sourcedir, f), targetfile) 
                print "Updated", targetfile        

main()

            
            