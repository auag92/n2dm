"""detect.py - support for the ASAP makefile.

Usage:

python detect.py compiler

Prints which compiler it seems reasonable to use.  The output will be
"pathscale" if the PathScale compiler is installed (pathcc and pathCC),
"intel" if the Intel C++ compiler is installed (icc and icpc), or
"gnu" if none of the above are installed.

python detect.py processor

Prints which processor.  Can return "unknown", "Pentium4", "Core2", or
"Opteron".  More possibilities may come.

"""


import sys, os

def main():
    if len(sys.argv) != 2:
        print >> sys.stderr, __doc__
        sys.exit(2)
    action = sys.argv[1]
    if action == "compiler":
        checkCompiler()
    elif action == "processor":
        checkProcessor()
    else:
        print >> sys.stderr, __doc__
        sys.exit(2)


def checkCompiler():
    print findCompiler()
    
def findCompiler():
    if checkpath("icc") and checkpath("icpc"):
        return "intel"
    elif checkpath("pathcc") and checkpath("pathCC"):
        return "pathscale"
    elif checkpath("gcc") and checkpath("g++"):
        return "gnu"
    else:
        print >> sys.stderr, "No compiler found (not even gcc and g++)!"

def checkpath(program):
    "Check if a program is on the path"
    for d in os.environ["PATH"].split(":"):
        candidate = os.path.join(d, program)
        if os.access(candidate, os.X_OK):
            return True
    return False

def checkProcessor():
    print findProcessor()

def findProcessor():
    try:
        lines = open("/proc/cpuinfo").readlines()
    except:
        print >> sys.stderr, "WARNING: Could not open /proc/cpuinfo - unknown processor"
        return "unknown"
    for l in lines:
        if l.find("Pentium(R) 4") != -1:
            return "Pentium4"
        elif l.find("Opteron") != -1:
            return "Opteron"
        elif l.find("Core(TM)2") != -1:
            return "Core2"
        elif l.find("Xeon") != -1:
            return "Xeon"
    
    print >> sys.stderr, "WARNING: Did not recognize processor type"
    return "unknown"

if __name__ == "__main__":    
    main()
