#!/usr/bin/env python
import os, filecmp, subprocess

def sync_analysis():
    
    unison -batch ~/vervet_data/analysis ~/vervet_data_project/analysis
    unison -batch ~/vervet_data_lab/analysis ~/vervet_data_project/analysis
    unison -batch ~/vervet_data_lab/analysis ~/vervet_data_project/analysis

def stage_in(in_fname,target_basedir):
    """
    in_fname        filename of  input_file (containing one name of a file to stage in per line)
                    file names in this file should be relative to the data directory, i.e. analysis/...
    target_basedir  e.g. '~/vervet_data_scratch'
    """
    #make sure that file names are relative to the data directory, i.e. analysis/...
    #for file in files:
    #    if file[0] == '/' or file[0] == '~':
    #        file = file[file.index('vervet_data/')+12:]
    
    
    target_basedir = os.path.expanduser(target_basedir)
    #check file size
    size = 0
    with open(in_fname,'r') as ff:
        for file in ff:
            fn_target = os.path.join(target_basedir, file)
            fn_project = os.path.join(os.path.expanduser(~/vervet_data), file)
            if os.path.exits(fn_target):
                if filecmp.cmp(fn_target,fn_project):
                    pass
                else:
                    size += os.path.getsize(fn_project)
            else:
                size += os.path.getsize(fn_project)

    if size>10^9:
        staging_qsub(fname,source_basedir,target_basedir)
    else:
        staging_direct()

def staging_direct(in_fname,source_basedir,target_basedir):
    """
    base filename is the path to the vervet_data directory
    on mendel, sym-links should assure that:
    ~/vervet_data ... scratch
    ~/vervet_data_project ... project
    ~/vervet_data_lab ... lab folder
    """
    
    #mem_limit=30*1024**2
    #command = "ulimit -v {mem_limit}; module load parallel; parallel -j -2 rsync -Rau :::  " 
    command = "rsync -au --files_from={in_fname} {source_basedir} {target_basedir}".format(**locals()) 
    #rsync -au {}
    p = subprocess.Popen(command.split())
    out, err = p.comunicate()
    rc = p.returncode
    return (rc, out, err)
    
def staging_job(in_fname,source_basedir,target_basedir):
    pass

            
    
