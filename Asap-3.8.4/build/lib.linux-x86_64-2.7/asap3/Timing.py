from asap3.Internal.Builtins import timing_results
import sys, string
from numpy import array

def _reporttiming(file, name, data, all=1):
    if (data[0]+data[2]) != 0:
        ptot = "%.0f%%" % (100*(data[1]+data[3])/(data[0]+data[2]))
    else:
        ptot = "-"
    if data[0] != 0:
        pself = "%.0f%%" % (100*data[1]/data[0])
    else:
        pself = "-"
    if data[2] != 0:
        pch = "%.0f%%" % (100*data[3]/data[2])
    else:
        pch = "-"
    if all:
        xtotal = "Total: %.1f / %.1f sec (%s)." % (data[1]+data[3],
                                                   data[0]+data[2], ptot,)
    else:
        xtotal = ""
    xself = "Self: %.1f / %.1f sec (%s)." % (data[1], data[0], pself)
    if all:
        xch = "Children: %.1f / %.1f sec (%s)." % (data[3], data[2], pch)
    else:
        xch = ""
    if len(name) > 25:
        file.write("%s\n%-26s n=%-4d %-33s %-33s %s\n" % (name+":", "",
                                                          data[4], xtotal,
                                                          xself, xch))
    else:
        file.write("%-26s n=%-4d %-33s %-33s %s\n" % (name+":", data[4], 
                                                      xtotal,xself, xch))
    if all and data[5] is not None:
        keys = data[5].keys()
        keys.sort()
        for master in keys:
            mdat = data[5][master]
            file.write("     for %-45s n=%-4d Tot: %.1f / %.1f  Self: %.1f / %.1f  Ch: %.1f / %.1f\n"
                       % (master+":", mdat[4], mdat[1]+mdat[3],
                          mdat[0]+mdat[2], mdat[1], mdat[0], mdat[3], mdat[2]))
            
def report_timing(fil = sys.stdout, header = "\nASAP Timing results:\n"):
    if type(fil) == type("string"):
        fil = open(fil, "w")
    fil.write(header+"\n")
    t = timing_results()
    if t is None:
        fil.write("No timing data!  Compile with -DTIMING to enable it.\n")
    else:
        _reporttiming(fil, "Global", t["Global"])
        _reporttiming(fil, "Timing overhead", t["Timing overhead"])
        del t["Global"], t["Timing overhead"]
        fil.write("\nFunctions and methods:\n")
        keys = t.keys()
        keys.sort()
        for l in keys:
            _reporttiming(fil, l, t[l])
        fil.write("\nClasses:\n")
        cls = {}
        for l in keys:
            words = string.split(l, "::")
            if len(words) == 2:
                try:
                    cls[words[0]].append(t[l][:5])
                except KeyError:
                    cls[words[0]] = [t[l][:5]]
        keys = cls.keys()
        keys.sort()
        for l in keys:
            v = sum(array(cls[l]))
            _reporttiming(fil, l, v, all=0)

