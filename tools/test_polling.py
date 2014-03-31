#!/usr/bin/env python
import subprocess

def poll_stats():
    """
    poll pbs stats once job has finished
    """
    p = subprocess.Popen("qstat -xf 229666.login0",shell=True,stdout=subprocess.PIPE)
    o,_ = p.communicate()
    print o
    
    #parse qstat
    

poll_stats()
