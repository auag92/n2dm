import os
import sys
import time

from string import Template

from mystuff.niflheim import submit_job

from asap3.Tools.MaterialProperties import MaterialPropertiesData

def use_template(target, template, args):
    f = open(template, 'r')
    temp = Template(f.read())
    f.close()

    f = open(target, 'w')
    f.write(temp.substitute(args))
    f.close()

#metals = ['Al', 'Ni', 'Cu', 'Pd', 'Ag', 'Pt', 'Au']
metals = ['Ni', 'Cu', 'Pd', 'Ag', 'Pt', 'Au']
#metals = ['Pt']
properties = [('a',    0.001),
              ('B',    0.01),
              ('A',    0.05),
              ('C11',  0.05),
              ('C12',  0.05),
              ('C44',  0.01),
              ('Ecoh', 0.001),
              ('E111', 0.02),
              ('E100', 0.02),
              ('E111_100', 0.01),
              ('fmatch', 0.03),
              ('Esf', 0.02),
              ]

mp = MaterialPropertiesData(['properties_metals.dat', 'properties_alloys.dat'])

oops = open('oops.sh', 'a')
oops.write('# ' + time.asctime() + '\n')

for m in metals:
    dirname = m + '_' + time.strftime('%d%m%y')
    if not os.path.exists(dirname):
        os.mkdir(dirname)

    args = {'symbol': m, 'lc': '%.3f' % (mp.get(m, 'a'),)}
    for p, w in properties:
        if p == 'E111_100':
            args[p] = '%.4f, %.5f' % (mp.get(m, 'E111') / mp.get(m, 'E100'), w)
        elif p == 'fmatch':
            args[p] = '%.5f' % w
        else:
            args[p] = '%.4f, %.5f' % (mp.get(m, p), w)

    use_template(dirname + '/fit.py', 'fit.py', args)

    if 'start' in sys.argv:
        os.chdir(dirname)
        jobid =  submit_job(['gpaw-qsub', 'fit.py'])
        print 'Job id:', jobid
        oops.write('qdel %s' % (jobid,))
        os.chdir('..')
