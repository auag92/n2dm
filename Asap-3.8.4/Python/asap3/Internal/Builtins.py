"""Asap builtin objects.

This module loads the modules from the compiled parts of Asap.  It
chooses between the parallel and serial version based on whether MPI
is active at runtime.
"""

__docformat__ = "restructuredtext en"

import sys

parallelpossible = sys.modules.has_key('asapparallel3')

# Now import the relevant C module
if parallelpossible:
    import asapparallel3 as _asap
else:
    import asapserial3 as _asap

AsapError = _asap.AsapError
get_version = _asap.get_version
get_short_version = _asap.get_short_version
timing_results = _asap.timing_results
