from numpy import *
from asap3 import *
from ase.lattice.cubic import SimpleCubic
from asap3.testtools import *
import time
import commands
import re

optimize = True
element = "Cu"
systemsizes = [(5,5,5), (10,10,10), (15, 15, 15), (20, 20, 25), (46,46,47), (100,100,100)]
#systemsizes = [(5,5,5), (10,10,10)]
nsteps = (1000, 1000, 1000, 1000, 1000, 100)
kT = 0.1
logfilename = "montecarlotiming.log"

def MakeSystem(size, element, perturbation):
    lc = data.reference_states[data.atomic_numbers[element]]['a'] / 2.0**(1.0/3.0)
    
    atoms = MonteCarloAtoms(SimpleCubic(directions=((1,0,0),(0,1,0),(0,0,1)),
                                        size=size, symbol=element,
                                        latticeconstant=lc,
                                        pbc=True))
    rs = atoms.get_scaled_positions()
    atoms.set_scaled_positions(0.5 * (rs + rs*rs))
    r = atoms.get_positions()
    dr = perturbation * random.standard_normal(r.shape)
    atoms.set_positions(r + dr)
    return atoms


percentages = []
sizes = []
timings = []

print_version()

for i,size in enumerate(systemsizes):
    atoms = MakeSystem(size, element, 0.1)
    n = len(atoms)
    #p1 = RasMol(atoms.Repeat((2,2,2)))
    # Pregenerate the trial steps as not to time random
    trial_id = random.randint(0, len(atoms), (nsteps[i],))
    lc = data.reference_states[data.atomic_numbers[element]]['a'] / 2.0**(1.0/3.0)
    trial_move = random.normal(0, lc, (nsteps[i],3))
    if optimize:
        atoms.set_calculator(MonteCarloEMT())
    else:
        atoms.set_calculator(EMT())
    eorig = e = atoms.get_potential_energy()

    nrej = 0
    startcpu, startwall = time.clock(), time.time()
    for j in range(nsteps[i]):
        atom = atoms[trial_id[j]]
        oldpos = atom.get_position()
        atom.set_position(oldpos + trial_move[j])
        newe = atoms.get_potential_energy()
        de = newe - e
        #print de
        if de > 500.0*kT:
            de = 500.0*kT
        if (de > 0.0 and random.random() > exp(-de/kT)):
            atom.set_position(oldpos)
            nrej = nrej + 1
        else:
            e = newe
    tcpu, twall = time.clock() - startcpu, time.time() - startwall
    atoms.set_scaled_positions(atoms.get_scaled_positions())
    #p2 = RasMol(atoms.Repeat((2,2,2)))
    pc = 100.0 * tcpu / twall
    timing = 1e6 * tcpu / nsteps[i]

    print
    print "Acceptance ratio:", 100.0 * (nsteps[i]-nrej) / nsteps[i]
    print "Energy change:", atoms.get_potential_energy() - eorig, "eV"
    print "CPU usage:", pc, "%"
    print "Number of atoms:", n
    print "Timing:", timing, "usec/step"
    percentages.append(pc)
    timings.append(timing)
    sizes.append(n)
    
now=time.strftime("%Y/%m/%d %H:%M")
asapversion = get_version()
version = asapversion.split()[2]
compiler = asapversion.split("'")[1]
host = commands.getoutput("hostname")
if re.match("^n\d\d\d.dcsc.fysik.dtu.dk$", host):
    print "    This is a d512 node on Niflheim."
    fullhost = "niflheim-d512/%s" % (host.split(".")[0])
    host = "niflheim-d512"
elif re.match("^[stu]\d\d\d.dcsc.fysik.dtu.dk$", host):
    print "    This is an s50 node on Niflheim."
    fullhost = "niflheim-s50/%s" % (host.split(".")[0])
    host = "niflheim-s50"
else:
    fullhost = host


loginfo = """Running on %s at %s
  Asap version: %s
  Monte Carlo Optimization: %s
""" % (fullhost, now, version, str(optimize))
for i in range(len(percentages)):
    loginfo = loginfo + ("    %d atoms: %.1f usec/step  (CPU usage: %.1f)\n"
                         % (sizes[i], timings[i], percentages[i]))

print "\nLog info:\n"
print loginfo
logfile = open(logfilename, "a")
logfile.write(loginfo+"\n")
logfile.close()
