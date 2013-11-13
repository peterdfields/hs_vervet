#!/usr/bin/env python
import socket, warnings, os, datetime, textwrap, subprocess
import numpy as np

#min total files size in bytes to use qsub for staging 
min_qsub_stage = 0#500000000

def stage(file_ls,source_base,target_base,mode,verbose=False):
    # this function should only run on dmn,
    # where all file systems are seen

    host = socket.gethostname()
    if 'dmn' not in host:
        warnings.warn("This script should be run on mendel data mover nodes where all file-systems are seen. However, filename is {}".format(host),UserWarning)

    source_base = os.path.expanduser(source_base)
    target_base = os.path.expanduser(target_base)    
    
    # all staging modes should preserve timestamp
    modes = ['non-exist','newer','force']
    if mode not in modes:
        raise ValueError('stage_mode must be in {}'.format(modes))
    
    if mode == "non-exist":
        #remove files from file-list that exist on target
        nonexist_file_ls = [file for file in file_ls if not os.path.exists(os.path.join(os.path.expanduser(target_base),file))]
        file_ls = nonexist_file_ls
    
    if mode == "newer":
        #remove files from file-list that are newer on target than on source
        def add_if_newer(file,file_ls):
            source_time = os.path.getmtime(os.path.join(source_base,file))
            try:
                target_time = os.path.getmtime(os.path.join(target_base,file))
                if source_time > target_time:
                    file_ls.append(file)
            except OSError:
                file_ls.append(file)
        newer_on_source = []
        for file in file_ls:
            if os.path.isdir(file):
                for root, _, fs in os.walk(os.path.join(source_base,file)):
                    for f in fs:
                        add_if_newer(os.path.join(root[len(source_base)+1:],f),newer_on_source)
            else:
                add_if_newer(file,newer_on_source)
        file_ls = newer_on_source
    
    if not file_ls:
        print "Nothing to stage in mode {0}".format(mode)
        return

    sizes = []
    for file in file_ls:
        sizes.append(os.path.getsize(os.path.join(source_base,file)))
    sizes=np.array(sizes)
    
    total = sizes.sum()
    number = sizes.shape[0]
    mean = sizes.mean()
    median = np.median(sizes)
    
    #implement a check for many small files (using number and median) and zip if necessary
    job_fn = write_jobscript(file_ls,source_base,target_base,mode)

    

    if total > min_qsub_stage:
        p = subprocess.Popen("qsub {0}".format(job_fn),shell=True)
    else:
        p = subprocess.Popen("{0}".format(job_fn),shell=True)


def write_jobscript(file_ls,source_base,destination_base,mode,job_fn=None,out_fn=None):
    """
    options (mcp options):
    p ... preserve timestamp
    r ... recursive
    """
    options = 'rp'    
    time = datetime.datetime.now().strftime("%Y%m%d-%H.%M%S")
    if job_fn is None:
        job_fn = os.path.expanduser("~/vervet_project/staging/jobscript/stage_{}.sh".format(time))
    if out_fn is None:
        out_fn = os.path.expanduser("~/vervet_project/staging/log/stage_{}".format(time))

    modes = ['non-exist','newer','force']
    if mode not in modes:
        raise ValueError('stage_mode must be in {}'.format(modes))
    if mode == 'non-exist':
        options += 'n'
    elif mode == 'newer':
        options += 'u'
    elif mode == 'force':
        options += 'f'


    stage_script=textwrap.dedent("""\
    #!/bin/bash
    #PBS -N stage
    #PBS -P vervetmonkey
    #PBS -q staging
    #PBS -o {out_fn}.o
    #PBS -e {out_fn}.e
    #PBS -l select=2:ncpus=2:mpiprocs=4 -l place=scatter
    #PBS -l walltime=02:00:00
    module load mutil
    cd {source_base}
    mpirun mcp -{options} --parents --print-stats --direct-read --direct-write --threads=1 --double-buffer --mpi {source} {destination_base}/
    """.format(out_fn=out_fn,source_base=source_base,options=options,source=' '.join([f for f in file_ls]),destination_base=destination_base))
    
    with open(job_fn,'w') as f:
        f.write(stage_script)
    p = subprocess.call(['chmod','ug+x',os.path.expanduser(job_fn)])
    return job_fn

    

    
    

if __name__ == '__main__':
    import sys
    sys.path.insert(0, os.path.expanduser('~/lib/python'))
    sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
    import argparse
    parser=argparse.ArgumentParser(description="Stage from /project on mendel\
                                                 to lsw12 or scratch")
    parser.add_argument("source_base",help="base directiory of source on dmn")
    parser.add_argument("target_base",help="base directiory of target on dmn")
    parser.add_argument("-m","--mode",choices=['non-exist','newer','force'],default='newer',help="Staging mode. Which files should be staged (only non existing, only newer or all.)")
    parser.add_argument("path_or_filename",help="Path or filename to stage.", nargs="+")
    
    args = parser.parse_args()

    stage(args.path_or_filename,args.source_base,args.target_base,args.mode)

