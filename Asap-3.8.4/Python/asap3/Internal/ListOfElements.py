"ListOfElements function returns which elements are present in an Atoms obj."

import numpy

def ListOfElements(atoms):
    """Return a list of the elements present in an Atoms object."""

    if getattr(atoms, "parallel", False):
        return atoms.get_list_of_elements()
    
    z = atoms.get_atomic_numbers()
    minz = min(z)
    maxz = max(z)
    loe = [minz]
    if minz != maxz:
        # More than one element
        for i in range (minz+1, maxz):
            if numpy.equal(z, i).any():
                loe.append(i)
        loe.append(maxz)
    return loe
