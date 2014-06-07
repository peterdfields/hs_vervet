#!/usr/bin/env python
from hs_vervet.tools import analysis as ana

def submit_nb_server_job(dir,port,walltime,ncpus):
    step = ana.Step(name="")
    #usw



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Start cluster job with python notebook server and bind to port.")
    parser.add_argument("-d","--dir",default="/home/GMI/hannes.svardal/script/hs_vervet/notebooks",
                        help="Directory to start notebook server in.")
    parser.add_argument("-p","--port",default=7000,help="Port to bind notebook server to.")
    parser.add_argument("-w","--walltime",default="12:00:00")
    parser.add_argument("-n","--ncpus",default=1)

    args = parser.parse_args()
    

