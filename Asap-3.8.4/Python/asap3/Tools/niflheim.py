import sys
import subprocess

def get_proces_output(command):
    qsubproc = subprocess.Popen(command, stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE)
    (out, err) = qsubproc.communicate()
    errcode = qsubproc.wait()
    if errcode:
        sys.stderr.write(err)
        #sys.stderr.write('Error %i in: %s\n' % (errcode, ' '.join(command)))
        sys.exit()
    return out

def submit_job(command):
    jobid = get_proces_output(command)
    #print 'Job id:', jobid
    if 'asap-qsub' in command:
        jobid = jobid[jobid.find('JOB ID') + 8:-1]
    return jobid

def submit_script(script, filename, command):
    scriptfile = open(filename, 'w')
    scriptfile.write(script)
    scriptfile.close()
    return submit_job(command + [filename])

