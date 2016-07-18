import sys
import numpy as np

from asap3.Tools.ParameterOptimization.Quantities import QuantityCalc

def ParameterPerformance(calculator, quantities, latticeconstants, txt=None,
                         debug=False, noprint=False):
    """
    Calculates the performance of a calculator based on some quantities.

    Parameters
    -----------
    calculator: Calculator object.

    quantities: List of quantities which should be calculated. The list should
    be on the form [(quantity_name, structure, elements, target), ...].

    txt: Output file. If none standard output is used.
    """
    quantitycalc = QuantityCalc(latticeconstants, calculator, debug=debug)

    if txt is None:
        f = sys.stdout
    else:
        f = open(txt, 'a')

    results = []
    for quan in quantities:
        name = quan[0]
        structure = quan[1]
        elements = quan[2]
        target = quan[3]
        valuefunc = getattr(quantitycalc, 'get_' + name)
        if isinstance(elements, (str, int)):
            elements = (elements,)
        else:
            elements = tuple(elements)
        value = valuefunc(structure, elements)

        results.append((name, structure, elements, target, value))

    if noprint:
        return results
    
    f.write("Potential Performance\n" +
            "%-45s  %-7s  " % ("Quantity, Structure:", "Target") +
            "%8s    %-9s\n" % ("Value", "Deviation") +
            "-" * 77 + "\n")
    devsum = 0
    for (name, structure, elements, target, value) in results:
        if len(elements) > 1:
            struct = "%s%s-%s" % (elements + (structure,))
        else:
            struct = "%s-%s" % (elements + (structure,))
        name = "%s, %s:" % (name, struct)
        dev = value / target * 100 - 100
        devsum += np.abs(dev)
        f.write("%-45s  %7.3f  " % (name, target) +
                "%10.6f  %8.2f%%\n" % (value, dev))
    f.write("%76.2f%%\n\n" % (devsum / len(results),))
    #f.close()


