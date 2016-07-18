import pickle
import sys

ncpu = {}
wall = {}

f = open(sys.argv[1])
try:
    while True:
        x, y, z = pickle.load(f)
        if x in ncpu:
            print "Warning: repeated data of type", x
        ncpu[x] = y
        wall[x] = z
except EOFError:
    pass

assert len(wall) == 4

print "Single CPU, running alone: %.3f" % (wall['A'],)
print "Serial, all CPUs active:   %.3f (%.0f%%)" % (wall['S'], 100*wall['A']/wall['S'])
print "Using the latter time as baseline for parallel performance."
base = wall["S"]
mfrac = base/wall['M']
print "Parallel using MPI:        %.3f (%.2f/%i, %.0f%%)" % (wall['M'], ncpu['M']*mfrac, ncpu['M'], 100*mfrac)
ofrac = base/wall['O']
print "Parallel using OpenMP:     %.3f (%.2f/%i, %.0f%%)" % (wall['O'], ncpu['O']*ofrac, ncpu['O'], 100*ofrac)
print
print "OpenMP versus MPI:         %.0f%%" % (100 * wall['M'] / wall['O'],)
