"""WrappedPositions - a filter wrapping the positions of an Asap list of atoms.

The positions of an Asap list of atoms are wrapped within the
computational box.  Implemented by accessing the wrapped positions in
the underlying C++ implementation.
"""

__docformat__ = "restructuredtext en"

class WrappedPositions:
    """WrappedPositions - a filter wrapping the positions of an Asap list of atoms.

    The positions of a `ListOfAtoms` are wrapped within the
    computational box.  Implemented by accessing the wrapped positions
    in the underlying C++ implementation.
    """

    def __init__(self, atoms):
        self.atoms = atoms

    def GetCartesianPositions(self):
        """Get the positions, but wrapped into the computational box."""
        return self.atoms.GetWrappedPositions()

    def SetCartesianPositions(self, pos):
        """Set the positions from positions wrapped into the computational box.

        The actual positions of the underlying list of atoms are not
        necessarily within the computational box, since the offset
        subtracted from the positions by `GetCartesianPositions` is
        added back again.

        """
        self.atoms.SetWrappedPositions(pos)

    def __getattr__(self, attr):
        return getattr(self.atoms, attr)

    
