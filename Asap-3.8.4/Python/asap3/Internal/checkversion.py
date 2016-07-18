"""Tests if the version numbers are consistent in Python and compiled code."""

from asap3 import _asap, __version__
from asap3 import __file__ as mainfile


versionerror = """


   OOPS  -  BAD ASAP INSTALLATION: INCONSISTENT VERSION NUMBERS

Version number of Python code: %s
Version number of compiled code: %s

Perhaps some modules are loaded from the wrong place.
  Python main module: %s
  Compiled module: %s

Or maybe you just forgot to compile the code after an upgrade.


"""


def check_version(verbose = False):
    "Check if the version numbers are consistent in Python and compiled code."
    try:
        compiled = _asap.get_short_version().strip("'")
    except AttributeError:
        compiled = "unknown (probably 3.0.0)"
    if compiled != __version__:
        compiledfile = _asap.__file__
        print versionerror % (__version__, compiled, mainfile, compiledfile)
        raise RuntimeError, "Inconsistent Asap version numbers (see message above)"

