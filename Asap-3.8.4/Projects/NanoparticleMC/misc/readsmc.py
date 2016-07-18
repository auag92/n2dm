import pickle
import sys
f = open(sys.argv[1])
f.read(len("parametermontecarlo"))
p = pickle.load(f)
a = pickle.load(f)
print "Length of counts:" + str(len(a['counts']))
print "Sum of counts:" + str( sum(a['counts']) )
print a['counts']
#find the minimum number of counts and the number of configurations with this multiplicity
mn = min(a['counts'])
ma= [i for i in a['counts'] if i==mn]
mp = len(ma)
print "Minimum number of counts: "+str(mn)
print "Number of times this min. occours: "+str(mp)


f.close()
