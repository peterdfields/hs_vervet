#!/usr/bin/env python
import sys, time, subprocess, socket


host = socket.gethostname()

if "login" in host or 'dmn' in host:
    call = ""
else: 
    call = "ssh login.mendel.gmi.oeaw.ac.at nohup "

def check_finished(pbs_id):
    p = subprocess.Popen("{}qstat {}".format(call,pbs_id),shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out, err = p.communicate()
    #if job is finished, the output goes to std_err
    if err:
        if "Job has finished" in err:
            return True
        else:
            raise Exception("Unexpected qstat error: "+err)
    else:
        return False  

def check_exit_status(pbs_id):
    p = subprocess.Popen("{}qstat -xf {}".format(call,pbs_id),shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out, err = p.communicate()
    info = [l.strip() for l in out.split("\n")]
    status_l = [l.split(" = ")[1].strip() for l in info if "job_state" in l.split(" = ")[0]]
    if len(status_l)==1:
        status=status_l[0]
    else:
        raise Exception('Job status could not be identified unambigously. Improve parsing.')

    if status != 'F':
        raise Exception("Job " + str(pbs_id) + " is not finished. PBS job status is "+ status)
    else:
        exit_status_l = [int(l.split(" = ")[1].strip()) for l in info if "Exit_status" in l.split(" = ")[0]]
        if len(exit_status_l)==1:
            exit_status = exit_status_l[0]
        else:
            raise Exception('Exit status could not be identified unambigously. Improve parsing.')
        return exit_status

def wait_for_job(pbs_id,poll_interval=30):
    finished = False
    time.sleep(5)
    while not finished:
        finished = check_finished(pbs_id)
        time.sleep(poll_interval)
    return check_exit_status(pbs_id)
        
    


if __name__ == '__main__':
    pass
#    import argparse
#    parser = argparse.ArgumentParser(description="Wait while job is running on cluster.")
#    parser.add_argument("job-id",description="PBS_ID of the job to wait for.")
    
    #test it
#    p = subprocess.Popen("qsub /home/GMI/hannes.svardal/test_pbs_polling/jobscript.sh",shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
#    out, err = p.communicate()
#    print out
#    id = out.strip()
#    print wait_for_job(id)
#    print "finished"

