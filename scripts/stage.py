#!/usr/bin/env python
import os, sys, socket, subprocess

"""
input: non-flag arguments should be single or multiple files
include option to read filenames from file

the script should support relative path

this could be run: on lws12 or on mendel
first, log into dmn
then, there should be 3 staging modes:
stage:
-non-exist
-newer
-all

for non-exist, first check locally what is there, only keep files in the staging list that are non-existant

usecases from lws12 and from mendel are quite different
on mendel it would be best to:

keep it simple:
on lws12, I only stage to and from labfolder
on mendel, I only launch staging to and from scratch


there are two parts of the staging process:
1) locally
here: if mode==non-exist, check which files exist already and remove them from the staging list, if the staging list is empty, the job is done.
2) on dmn
do the rest:
i.e. if mode==newer, check which files are more recent at the source, remove those which are not from the staging list,
for the remaining files, check size
- if many small files -> compress
- if very large files: use mcp
for little data: use rsync


"""

def local_prepare_staging(file_ls,partner,direction,mode,project='vervet'):
    """
    this is run locally and establishes ssh connection to dmn 
    where staging proceeds
    """
        
    # sanity check for input 
    partners = ['lab','scratch']
    if partner not in partners:
        raise ValueError('staging partner (source/destination other than /project) must be in {}'.format(stage_partners))
    host = socket.gethostname()

    directions = ['in','out']
    if direction not in directions:
        raise ValueError('direction must be in {}'.format(directions))

    modes = ['non-exist','newer','force']
    if mode not in modes:
        raise ValueError('stage_mode must be in {}'.format(modes))
    
    project_base = os.path.join("~/",project + '_project')
    scratch_base = os.path.join("~/",project + '_'  +partner)    

    if direction == 'in':
        source_base = project_base
        destination_base = scratch_base
    elif direction == 'out':
        source_base = scratch_base
        destination_base = project_base
        
    if partner == 'scratch' and 'login' not in host and 'dmn' not in host:
        raise Exception('staging to scratch only implemented when running on mendel. Your host name is {}'.format(host)) 

    if mode == "non-exist":
        #check wether files exist locally
        nonexist_file_ls = [file for file in file_ls if not os.path.exists(os.path.join(os.path.expanduser(destination_base),file))]
        file_ls = nonexist_file_ls
    
    
    command = "dmn_stage.py -m {mode} {source_base} {target_base} {files}".format(mode=mode,source_base=source_base,target_base=destination_base,files=' '.join(file_ls))
    if 'dmn' in host:
        p = subprocess.Popen(command, shell=True)
    else:
        p = subprocess.Popen("ssh dmn.mendel.gmi.oeaw.ac.at nohup {0}".format(command), shell=True)#, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    #out, err = p.communicate()
    #self.returncode = p.returncode        
               
       




if __name__ == '__main__':
    sys.path.insert(0, os.path.expanduser('~/lib/python'))
    sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
    import argparse

    parser=argparse.ArgumentParser(description="Stage from /project on mendel\
                                                 to lsw12 or scratch")
    parser.add_argument("path_or_filename",help="Path or filename to stage.", nargs="+")
    parser.add_argument("-p","--partner",choices=['scratch','lab','auto'],default='auto',help="staging partner (source/destination other than /project)")
    parser.add_argument('-d',"--direction",default="in",choices=['in','out'],help="direction of staging (from or to project folder)")
    parser.add_argument("-m","--mode",choices=['non-exist','newer','force'],default='newer',help="Staging mode. Which files should be staged (only non existing, only newer or all.)")
    parser.add_argument("-v","--verbose",action="store_true")
    args =  parser.parse_args()
    if args.partner == 'auto':
        host = socket.gethostname()
        if host == 'gmi-lws12':
            partner = 'lab'
        elif 'login' in host or dmn in host:
            partner = 'scratch'
        else:
            raise Exception('"auto" staging partner determination only implemented for lws12 or mendel')
    else:
        partner = args.partner 
    local_prepare_staging(args.path_or_filename,partner,args.direction,args.mode)


