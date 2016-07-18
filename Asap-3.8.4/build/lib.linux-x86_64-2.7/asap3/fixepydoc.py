"""Fix epydoc documentation.

The function fix is used to fix epydoc documentation by fooling
epydoc into documenting objects although they are really documented in
a submodule.
"""

import inspect, sys

if sys.modules.has_key('epydoc'):
    def fix(*objects):
        """Fix epydoc documentation for objects.
        
        The function fixdoc is used to fix epydoc documentation by
        fooling epydoc into documenting objects although they are
        really documented in a submodule.
        
        Usage:
    
        import Asap.fixepydoc
        from submodule import a, b
        Asap.fixepydoc.fix(a, b)
        
        This will cause epydoc to document objects (classes,
        functions, ...) a and b as if they were actually defined here.
        If epydoc is not running, nothing happens.
        
        Do not call fix on the same object in more than one module.
        """
        caller = inspect.getmodule(inspect.currentframe().f_back).__name__
        assert caller is not None
        for o in objects:
            assert hasattr(o, '__module__')
            o.__module__ = caller
else:
    def fix(*objects):
        "Fix epydoc documentation: Do nothing as epydoc is not running."
        pass
    
            
