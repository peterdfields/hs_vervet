#!/usr/bin/env python
"""
input: non-flag arguments should be single or multiple files
include option to read filenames from file
"""
from __future__ import print_function
import os, sys, socket, subprocess, datetime
from hs_vervet_basics import v_print


def local_prepare_staging(file_ls_or_fnfname,partner,direction,mode,run_type='auto',afterok=None,afterany=None,startonhold=False,job_fn=None,out_fn=None,job_name=None,verbose=0,project='vervet',file_to_print_to=None):
    """
    this is run anywhere (lws12, mendel login or dmn) and establishes ssh connection to dmn 
    where staging proceeds
    file_ls_or_fnfname is either a list of files to stage or the path to a file that contains the filenames to stage
    """
    #empty the file to print to (usually things are appended
    if file_to_print_to is not None:
        v_print("",file=file_to_print_to,append=False)
    #prints to stdout if file_to_print_to is None
    vprint = lambda text,min_verb: v_print(text,min_verb,verbose,file_to_print_to)   

    vprint(datetime.datetime.now().strftime("%Y%m%d-%H%M%S")+" - running local_prepare_staging in stage.py",1)
    
    if afterany is None:
        afterany = []
    if afterok is None:
        afterok = []

    if type(file_ls_or_fnfname) is str:
        with open(file_ls_or_fnfname,'r') as f:
            file_ls = [s.strip() for s in f.readlines()]
    elif type(file_ls_or_fnfname) is list:
        file_ls = file_ls_or_fnfname
    else:
        raise TypeError('First argument should be list of filenames or path to a file that contains filenames. But it is {0}:{1}'.format(type(file_ls_or_fnfname),file_ls_or_fnfname))

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
    
    if not file_ls:
        print("Nothing to stage in mode {0}".format(mode))
        return (None, None,0)


    command = "dmn_stage.py -m {mode} -t {run_type}  {source_base} {target_base} {files}".format(mode=mode,run_type=run_type,source_base=source_base,target_base=destination_base,files=' '.join(file_ls))

    command += " -v " + str(verbose) 

    if job_fn is not None:
        command += " -j " + job_fn

    if out_fn is not None:
        command += " -o " + out_fn

    if job_name is not None:
        command += " -n " + job_name

    if afterok:
        command += " --afterok {}".format(' '.join(afterok))
    if afterany:
        command += " --afterany {}".format(' '.join(afterany))
    if startonhold:
        command += " -H"
    if file_to_print_to is not None:
        command += " -l {}.dmn".format(file_to_print_to)
    
    #print(command)

    if 'dmn' in host:
        vprint('command:',1)
        p = subprocess.Popen(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    else:
        vprint('command submitted to dmn via ssh dmn.mendel.gmi.oeaw.ac.at nohup <command>:',1)
        p = subprocess.Popen("ssh dmn.mendel.gmi.oeaw.ac.at nohup {0}".format(command), shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out, err = p.communicate()
    rc = p.returncode        
    
    
    vprint(command,1)

    if run_type != 'submit':
        if out is not None:
            print('dmn_stage.py','out:',out, file=sys.stdout)
        if err is not None:    
            print('dmn_stage.py','err:',err, file=sys.stderr)

    return out, err, rc            
       




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
    parser.add_argument("-t","--run-type",choices=['direct','submit','auto','dry_run'],default='auto',help='Decides wheter staging should be run on data mover node directly or through qsub. Default (auto) make choice depending on file size and number.')
    parser.add_argument("-j","--job-fname",default=None,help="Full filename of jobfile. By default it is put in ~/vervet_project/staging/...")
    parser.add_argument("-n","--job-name",default=None,help="Name of the job used in the stage-script.")
    parser.add_argument("-o","--stdout-fname",default=None,help="Full filename of out/err filenames. By default it is put in ~/vervet_project/staging/...")
    parser.add_argument("-b","--project-base",default="vervet",type=str)
    parser.add_argument("--dry-run",action="store_true")
    parser.add_argument("-l","--local-print-file",default=None,type=str)
    parser.add_argument("-v","--verbose",type=int,default=0)
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
    
    
    possible_base_dirs = ['/projects/vervetmonkey/','/lustre/scratch/projects/vervetmonkey/','/net/gmi.oeaw.ac.at/nordborg/lab/Projects/vervetpopgen/']
    
    #print('input fn:',args.path_or_filename)
    #get filenames relative to project dir
    rel_fnames=[]
    for file in args.path_or_filename:
        real_path = os.path.realpath(file)
        print(real_path)
        rel_fn = real_path
        for bd in possible_base_dirs:
            if real_path.startswith(bd):
                    rel_fn = real_path[len(bd):]
                    break
        rel_fnames.append(rel_fn)
    if len(rel_fnames) != len(args.path_or_filename):
        raise Exception('File number not constistent after preparation. before: {0}, after: {1}'.format(args.path_or_filename,rel_fnames))
    #print('===========STAGE================')
    #print('rel_fnames:',rel_fnames)
    #sys.exit(0)
    


    local_prepare_staging(rel_fnames,partner,args.direction,args.mode,run_type=args.run_type,job_fn=args.job_fname,out_fn=args.stdout_fname,verbose=args.verbose,file_to_print_to=args.local_print_file,job_name=args.job_name,project=args.project_base)


