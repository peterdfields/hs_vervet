#!/usr/bin/env python
from __future__ import print_function
import socket, warnings, os, datetime, textwrap, subprocess, sys
import numpy as np
from hs_vervet_basics import v_print

#min total files size in bytes to use qsub for staging
min_mcp_size = 500000000 
min_qsub_size = 500000000

def stage(file_ls,source_base,target_base,mode,run_type='auto',ignore_ls=None,job_fn=None,out_fn=None,name=None,afterok=None,afterany=None,startonhold=False,verbose=0,file_to_print_to=None):
    #print("file_ls:",file_ls)
    # this function should only run on dmn,
    # where all file systems are seen
    #empty the file to print to (usually things are appended
    if file_to_print_to is not None:
        v_print("",file=file_to_print_to,append=False)
    #prints to stdout if file_to_print_to is None
    def vprint(*text,**kwa):
        """
        use similar to python3 print function
        additional argument mv or min_verbosity
        and append (to know whether append to file or not)
        """
        kwa.update({"file":file_to_print_to,"verbosity":verbose})
        v_print(*text,**kwa)

    def add_if_newer(file,file_ls):
        #print(file)
        source_time = os.path.getmtime(os.path.join(source_base,file))
        #vprint("source:",os.path.join(source_base,file), source_time,mv=1)
        try:
            target_time = os.path.getmtime(os.path.join(target_base,file))
            #vprint("target:",os.path.join(target_base,file),target_time,mv=1)
            #The int conversion is necessary, cause the 
            if int(source_time) > int(target_time):
                file_ls.append(file)
                vprint("Staging, source newer:"+file,mv=1)
            else:
                vprint("Not staging, not newer on source: "+file,mv=1)
        except OSError:
            file_ls.append(file)
            vprint("Staging, not exist on target:"+file,mv=1)
    
    def add_if_not_exist(file,file_ls):
        if not os.path.exists(os.path.join(os.path.expanduser(target_base),file)):
            file_ls.append(file)

    def add_if(file,file_ls):
        for ign in ignore_ls:
            #try:
            if ign in file:
                vprint("Not staging, matching ignore pattern {}: ".format(ign)+file,mv=1)
                return
            #except:
            #    print("ign, file:",ign, file)
            #    raise
        if mode == "newer":
            add_if_newer(file,file_ls)
        elif mode == "non-exist":
            add_if_not_exist(file,file_ls)
        elif mode == force:
            file_ls.append(file)

    vprint("output of",__file__,mv=10)


    #verboseprint = print if verbose else lambda *a, **k: None

    if afterany is None:
        afterany = []
    if afterok is None:
        afterok = []
    if ignore_ls is None:
        ignore_ls = []

    host = socket.gethostname()
    if 'dmn' not in host:
        warnings.warn("This script should be run on mendel data mover nodes where all file-systems are seen. However, filename is {}".format(host),UserWarning)

    source_base = os.path.expanduser(source_base)
    target_base = os.path.expanduser(target_base)    
    
    # attention, this is now duplicated (here and in the write_jobscript function), find a solution for this
    time = datetime.datetime.now().strftime("%Y%m%d-%H.%M%S")
    if job_fn is None:
        job_fn = os.path.expanduser("~/vervet_project/staging/jobscript/stage_{}.sh".format(time))
    if out_fn is None:
        out_fn = os.path.expanduser("~/vervet_project/staging/log/stage_{}".format(time))
     
    # all staging modes should preserve timestamp
    modes = ['non-exist','newer','force']
    if mode not in modes:
        raise ValueError('stage_mode must be in {}'.format(modes))
    
#    if mode == "non-exist":
#        #remove files from file-list that exist on target
#        nonexist_file_ls = [file for file in file_ls if not os.path.exists(os.path.join(os.path.expanduser(target_base),file))]
#        file_ls = nonexist_file_ls
#    
#    #print('before:',file_ls)    
#    vprint("staging form '" +source_base+"' to '" + target_base + "'",mv=1)
#    vprint("staging mode: "+mode,mv=1)
#    vprint("run type: "+run_type,mv=1)
#    n_files = 0
#    if mode == "newer":# and (run_type == 'direct' or run_type == 'auto'):
#        #remove files from file-list that are newer on target than on source
#        newer_on_source = []
#        for file in file_ls:
#            #print('from filelist:',file)
#            if os.path.isdir(os.path.join(source_base,file)):
#                #print('is a dir:',file)
#                for root, _, fs in os.walk(os.path.join(source_base,file)):
#                    #print('fs:',fs)
#                    for f in fs:
#                        #print('f:',f)
#                        n_files+=1
#                        add_if_newer(os.path.join(root[len(source_base)+1:],f),newer_on_source)
#            else:
#                n_files += 1
#                add_if_newer(file,newer_on_source)
#        vprint("Staging " + str(len(newer_on_source)) + " out of " + str(n_files) + " files." ,mv=1)
#        file_ls = newer_on_source
    

    #walk through the directories and decide whether to add the files
    #TODO: incorporate the size check into the function here! Then only one loop is necessary...
    n_files = 0
    retained_files = []
    for file in file_ls:
        #print("file:",file)
        if os.path.isdir(os.path.join(source_base,file)):   
            for root, _, fs in os.walk(os.path.join(source_base,file)):
                #print('fs:',fs)
                for f in fs:
                    #print('f:',f)
                    n_files+=1
                    rel_path = root[len(source_base)+1:]
                    add_if(os.path.join(rel_path,f),retained_files)
        else:
            n_files += 1
            add_if(file,retained_files)         
    vprint("Staging " + str(len(retained_files)) + " out of " + str(n_files) + " files." ,mv=1)
    file_ls = retained_files

    if not file_ls:
        vprint("Nothing to stage in mode {0}".format(mode))
        return (None, None, 0)
        

    sizes = []
    for file in file_ls:
        try:
            sizes.append(os.path.getsize(os.path.join(source_base,file)))
        except OSError, e:
            vprint(warnings.warn("Can't check size of file. Does not exist. Copy operation might not be optimised. "+str(e),UserWarning),mv=1)
            
            
    sizes=np.array(sizes)

    total = sizes.sum()
    number = sizes.shape[0]
    mean = sizes.mean()
    median = np.median(sizes)
    
    #implement a check for many small files (using number and median) and zip if necessary
    
    vprint("Total file size is " + str(total) + " bytes.",mv=1)  

    if total <= min_mcp_size:
        method = 'cp'
        vprint("using cp",mv=1)
    else:
        method = 'mcp'
        vprint("using mcp",mv=1)
    
    job_fn = write_jobscript(file_ls,source_base,target_base,mode,method,job_fn=job_fn,out_fn=out_fn,name=name,verbose=verbose)

    #print(job_fn)

    if run_type == 'dry_run':
        with open(os.path.expanduser(job_fn),'r') as f:
            for line in f:
                print(line.strip())
        sys.exit(0)
    

    if run_type == 'submit' or (run_type == 'auto' and total > min_qsub_size):
        depend_str=''
        if afterany or afterok:
            depend_str='-W depend='
            if afterany:
                depend_str += 'afterany'
                for pbs_id in afterany:
                    depend_str = depend_str + ':' + pbs_id
            if afterok:
                if afterany:
                    depend_str += ','
                depend_str += 'afterok'
                for pbs_id in afterok:
                    depend_str = depend_str + ':' + pbs_id 
        if startonhold:
            hold='-h '    
        else:
            hold=''
        command = "qsub {hold}{depend_str} {job_fn}".format(hold=hold,depend_str=depend_str,job_fn=job_fn)
        p = subprocess.Popen(command,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        vprint("submitted stage job with:",command,mv=1)
        out, err = p.communicate()    #ran_as = 'submit'
    else:
        vprint("running stage job at",socket.gethostname(),":",job_fn,mv=1)
        p = subprocess.Popen(job_fn,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out, err = p.communicate()
        with open(out_fn+'.o','w') as f:
            f.write(out)
        with open(out_fn+'.e','w') as f:
            f.write(err)
        #ran_as = 'direct'
    


    rc = p.returncode

    #if run_type != 'submit':
    #    if out is not None:
    #        print(job_fn,'out:',out,file=sys.stdout)
    #    if err is not None:    
    #        print(job_fn,'err:',err,file=sys.stderr)
    
    vprint('stage job out: ' + out.strip(),mv=1)
    vprint('stage job err: ' + err.strip(),mv=1)

    return (out, err, rc)

def write_jobscript(file_ls,source_base,destination_base,mode,method='mcp',job_fn=None,out_fn=None,name=None,verbose=0):
    """
    options (mcp options):
    p ... preserve timestamp
    r ... recursive
    """
    
    methods =['cp','mcp']
    if method not in methods:
        raise ValueError('stage_method must be in {}'.format(methods))    
    

    time = datetime.datetime.now().strftime("%Y%m%d-%H.%M%S")
    if job_fn is None:
        job_fn = os.path.expanduser("~/vervet_project/staging/jobscript/stage_{}.sh".format(time))
    if out_fn is None:
        out_fn = os.path.expanduser("~/vervet_project/staging/log/stage_{}".format(time))
    if name is None:
        name = 'stage_' + source_base.split('_')[-1][0] + destination_base.split('_')[-1][0]



    modes = ['non-exist','newer','force']
    if mode not in modes:
        raise ValueError('stage_mode must be in {}'.format(modes))

    if method == 'cp':
        options = 'a'    
        if mode == 'non-exist':
            options += 'n'
        elif mode == 'newer':
            options += 'u'
        elif mode == 'force':
            options += 'f'
        if verbose > 0:
            options += 'v'

        stage_script=textwrap.dedent("""\
        #!/bin/bash
        #PBS -N {name}
        #PBS -P vervetmonkey
        #PBS -q staging
        #PBS -o {out_fn}.o
        #PBS -e {out_fn}.e
        #PBS -l select=1:ncpus=1
        #PBS -l walltime=00:25:00
        cd {source_base}
        cp -{options} --parents {source} {destination_base}/
        """.format(name=name,out_fn=out_fn,source_base=source_base,options=options,source=' '.join([f for f in file_ls]),destination_base=destination_base))
    elif method == 'mcp':
        options = 'a'    
        if mode == 'non-exist':
            options += 'n'
        elif mode == 'newer':
            options += 'u'
        elif mode == 'force':
            options += 'f'
        if verbose > 0:
            options += 'v'

        stage_script=textwrap.dedent("""\
        #!/bin/bash
        #PBS -N {name}
        #PBS -P vervetmonkey
        #PBS -q staging
        #PBS -o {out_fn}.o
        #PBS -e {out_fn}.e
        #PBS -l select=2:ncpus=2:mpiprocs=4 -l place=scatter
        #PBS -l walltime=02:00:00
        module load mutil
        cd {source_base}
        mpirun mcp -{options} --parents --print-stats --direct-read --direct-write --threads=1 --double-buffer --mpi {source} {destination_base}/
        """.format(name=name,out_fn=out_fn,source_base=source_base,options=options,source=' '.join([f for f in file_ls]),destination_base=destination_base))
    
    with open(job_fn,'w') as f:
        f.write(stage_script)
    p = subprocess.call(['chmod','ug+x',os.path.expanduser(job_fn)])
    return job_fn

    

    
    

if __name__ == '__main__':
    sys.path.insert(0, os.path.expanduser('~/lib/python'))
    sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
    import argparse
    parser=argparse.ArgumentParser(description="Stage from /project on mendel\
                                                 to lsw12 or scratch")
    parser.add_argument("source_base",help="base directiory of source on dmn")
    parser.add_argument("target_base",help="base directiory of target on dmn")
    parser.add_argument("-m","--mode",choices=['non-exist','newer','force'],default='newer',help="Staging mode. Which files should be staged (only non existing, only newer or all.)")
    parser.add_argument("-t","--run-type",choices=['direct','submit','auto','dry_run'],default='auto',help='Decides wheter staging should be run on data mover node directly or through qsub. Default (auto) make choice depending on file size and number.')
    parser.add_argument("-j","--job-fname",default=None,help="Full filename of jobfile. By default it is put in ~/vervet_project/staging/...")
    parser.add_argument("-o","--stdout-fname",default=None,help="Full filename of out/err filenames. By default it is put in ~/vervet_project/staging/...")
    parser.add_argument("-n","--job-name",default=None,help="Name used if stage job is subnitted.")
    parser.add_argument("--afterok",nargs='+',help='IDs of jobs that this job depends on. Only runs on exit status 0. See PBS documentation.')
    parser.add_argument("--afterany",nargs='+',help='IDs of jobs that this job depends on. Runs on any exit status of depend job. See PBS documentation.')
    parser.add_argument("-H","--start-on-hold",action="store_true",help='Start the job on user hold. Only effective if run type is submit.')
    parser.add_argument("-v","--verbose",type=int,default=0)
    parser.add_argument("-l","--local-print-file",default=None,type=str)
    parser.add_argument("-i","--ignore",nargs="*",action="store",default=None)
    parser.add_argument("path_or_filename",default=None,help="Path or filename to stage.", nargs="*")
    parser.add_argument("-L","--file_list",default=None,help="Path to a list of files to stage.")
    
    args = parser.parse_args()
    
    if args.path_or_filename is None: 
        args.path_or_filename = []
    if args.file_list is None:
        args.file_list = ""

    if not args.path_or_filename and not args.file_list:
        raise parser.error("Either <path_or_filename> or --file_list <file_list> have to be specified.")
    
    if args.file_list:
        args.path_or_filename += [line.strip() for line in open(args.file_list, 'r')]

    #print("path_or_fn:",args.path_or_filename)

    #print("ignore",args.ignore)
    #sys.exit(1)
    (out, err, rc) = stage(args.path_or_filename,args.source_base,args.target_base,args.mode,run_type=args.run_type,ignore_ls=args.ignore,afterok=args.afterok,afterany=args.afterany,startonhold=args.start_on_hold,job_fn=args.job_fname,out_fn=args.stdout_fname,name=args.job_name,file_to_print_to=args.local_print_file,verbose=args.verbose)

    
    #print('in main of dmn_stage,rc:',rc)
    print(out,file=sys.stdout)
    print(err,file=sys.stderr)
    sys.exit(rc)

