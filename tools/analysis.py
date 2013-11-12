#!/usr/bin/env python

import sys, os, datetime, subprocess, socket, filecmp
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from hs_vervet.tools import hs_vervet_basics as hvb


class Analysis(object):
    def __init__(self, name=None,project="~/vervet_project",scratch="~/vervet_scratch"):
        self.name = name
        self.host = socket.gethostname()
        self.scratch = scratch
        self.project = project
        self.ana_dir = os.path.join('analyses/',self.name)
        self.submit_log_fn = os.path.join(self.project,self.ana_dir,'log/', self.name+"_submit_log.txt") 
        for direc in ["_data","log","script","jobscript","io"]:
            hvb.try_make_dirs(os.path.join(self.project,self.ana_dir,direc))
        self.steps=[]

    #def ana_dir(self,loc):
    #    if loc == 'scratch':
    #        return os.path.join(self.scratch_dir,'analyses/',self.name)
    #    elif loc == 'project':
    #        return os.path.join(self.project_dir,'analyses/',self.name)
    #    else:
    #        raise Exception('location must be "scratch" or "project" but is {}'.format(loc))
        

    def append_step(self, step):
        self.steps.append(step)
        step.bind_to_analysis(self)
        
    def append_submit_log(self, jobfile, out, err):
        with open(os.path.expanduser(self.submit_log_fn),'a') as lf:
            date=datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
            lf.write(date+'\t'+jobfile+'\t'+out+'\t'+err+'\n')
    def write_and_qsub_all(self):
        for step in self.steps:
            for job in step.jobs:
                job.write_jobscript()
                job.qsub_jobscript()
            

class AnalysisStep(object):

    def __init__(self, analysis=None, name=None, jobs=None, append_to_ana=True, depend_steps = None, stagein=True, stageout=True, default_run='qsub',  check_date_input = True, verbose = False):
        #add options stagein=True and stageout=True that triger automatic staging
        # the next to lines are just for the following methods to work
        # the attributes will be set again by their constructors
        self.name = name
        self.jobs = ([] if jobs is None else jobs)
        self.stagein = stagein
        self.stageout = stageout
        self.default_run = default_run
        self.stagein_job = None
        self.stageout_job = None
        if analysis is None:
            self.analysis = analysis
        else: 
            if append_to_ana:
                analysis.append_step(self)
            else:
                self.bind_to_analysis(analysis)
        self.add_name(name)
        self.depend_steps = ([] if depend_steps is None else depend_steps)
        self.jobs = []
        if jobs != None:
            for job in jobs:
                self.append_job(job)
        self.check_date_input = check_date_input 
        
    def add_name(self,name):
        if name is None:
            self.name = 'no_name'
        else:
            self.name = str(name)
        self.short_name = ''.join([s[0].upper() for s in self.name.split('_')])
        #if name is None:
        #    self.in_fname = None
        #    self.out_fname = None
        #else:
        #    self.in_fname = os.path.join(self.analysis.ana_dir,'io/',name+'.in')    
        #    self.out_fname = os.path.join(self.analysis.ana_dir,'io/',name+'.out')    
        #    #self.stage_job os.path.join(self.scratch_dir,'analyses/',self.name)
    
    def bind_to_analysis(self, analysis):
        self.analysis = analysis
        #self.in_fname = os.path.join(self.analysis.ana_dir,'io/',self.name+'.in')    
        #self.out_fname = os.path.join(self.analysis.ana_dir,'io/',self.name+'.out')    
        for job in self.jobs:
            job.analysis = analysis
        
    def append_job(self, job):
        #only add the job if it is not added yet
        if job not in self.jobs:
            self.jobs.append(job)
            job.bind_to_analysis_step(self)
            for depend_step in self.depend_steps:
                for depend_job in depend_step.jobs:
                    job.depends.append(depend_job)
        else:
            UserWarning("Wont append job {0} to step {1}, because it is already in {1}'s job-list".format(job.name,self.name))
                  
    def append_jobs(self,jobs):
        for job in jobs:
            self.append_job(job)


    def run(self,mode=None):
        modes=['qsub','local','write_jobscripts']
        if mode is None:
            mode = self.default_run
        if mode not in modes:
            raise Exception("mode must be in {0} but is {1}".format(modes,mode))
        self.prepare_staging()
        if mode == 'qsub':
            self.qsub()
        elif mode == 'local':
            self.local_run()
        elif mode == 'write_jobscripts':
            self.write_jobscripts()
    
    def prepare_staging(self):
        if self.stagein: #here we could also check if stagein_job already exists
            self.add_stagein()
        if self.stageout:
            self.add_stageout()
                
    def write_jobscripts(self):    
        if self.stagein_job is not None:
            self.stagein_job.write_jobscript()
        for job in self.jobs:
            job.write_jobscript()
        if self.stageout_job is not None:
            self.stageout_job.write_jobscript()

    def local_run(self):
        if self.stagein_job is not None:
            self.stagein_job.write_jobscript()
            self.stagein_job.run_jobscript()
        for job in self.jobs:
            job.write_jobscript()
            job.run_jobscript()
        if self.stageout_job is not None:
            self.stageout_job.write_jobscript()
            self.stageout_job.run_jobscript()

    def qsub(self):
        if self.stagein_job is not None:
            self.stagein_job.write_jobscript()
            self.stagein_job.qsub_jobscript()
        for job in self.jobs:
            job.write_jobscript()
            job.qsub_jobscript()
        if self.stageout_job is not None:
            self.stageout_job.write_jobscript()
            self.stageout_job.qsub_jobscript()    
    
    def print_summary(self,job_summary=True):
        print "="*60
        print "analysis:", (None if self.analysis is None else self.analysis.name)
        print "analysis_step:", ('unnamed' if self.name is None else self.name)
        print "stagein_job:", (None if self.stagein_job is None else self.stagein_job.name)
        if job_summary and self.stagein_job is not None:
            self.stagein_job.print_summary()
        for i,job in enumerate(self.jobs):
            print 'job {0}:'.format(i), ('unnamed' if job.name is None else job.name)
            if job_summary:
                job.print_summary()
        print "stageout_job:", (None if self.stageout_job is None else self.stageout_job.name)
        if job_summary and self.stageout_job is not None:
            self.stageout_job.print_summary()
        print "="*60

    def add_stagein(self, stage_analysis_dir=True):
        in_files = []
        for job in self.jobs:
            for file in job.input_files:
                in_files.append(file)
        in_files=list(set(in_files)) #only unique files
        #remove a possible pre-existing stagein job from the jobs
        if self.stagein_job is not None:
            for job in self.jobs:
                while self.stagein_job in job.depends:
                    job.depends.remove(self.stagein_job)
        stagein_job = StageJob(direction='in',files=in_files,analysis_step=self,stage_analysis_dir=stage_analysis_dir)
        for job in self.jobs:
            job.depends.append(stagein_job)
        self.stagein_job = stagein_job
  

    def add_stageout(self, stage_analysis_dir=True, add_to_jobs=False):
        out_files = []
        for job in self.jobs:
            for file in job.output_files:
                out_files.append(file)
        out_files=list(set(out_files)) #only unique files
        stageout_job = StageJob(direction='out',files=out_files,analysis_step=self,stage_analysis_dir=stage_analysis_dir,depends=[job for job in self.jobs])
        self.stageout_job = stageout_job       
        

"""
    def write_stage_in_file(self):
        with open(self.in_fname) as f:
            for fn in self.input_files:
                f.write(fn+'\n')

    def write_stage_out_file(self):
        with open(self.out_fname) as f:
            for fn in self.output_files:
                f.write(fn+'\n')
    
    def write_stage_script(self,mode='direct'):
        pass

    def stage_in(self):
        # check total filesize for things to stage in
        # limit for considering a file as small in (bytes):
        #small_limit = 2*2**20
        # starting form how many small files should we start zipping (not implemented yet)
        #max_small = 500
        host = socket.gethostname()
        if host == 'gmi-lws12':
            target_basedir = self.analysis.base_dir+'_lab'
        elif 'login' in host or 'dmn' in host:
            target_basedir = self.analysis.base_dir + '_scratch'  
        else:
            raise Exception('host must be a login node, dnm or gmi-lws12, but it is {}'.format(host))
        source_basedir = self.analysis.base_dir        
        
        self.write_stage_in_file()
        
        analysis_sync
        
        check_size(in_fname)
        
        if 'dmn' not in host:
            ssh dmn run stagin.py self.in_fname
        else:
            run stagin.py self.in_fname
                

        size = 0
        for file in files:
            fn_target = os.path.join(target_basedir, file)
            fn_project = os.path.join(os.path.expanduser(~/vervet_data_project), file)
            if os.path.exits(fn_target):
                if filecmp.cmp(fn_target,fn_project):
                    pass
                else:
                    size += os.path.getsize(fn_project)
            else:
                size += os.path.getsize(fn_project)
"""

"""
    if size>10^9:
        qsub_staging
    else:
        local stagingtage_in_direct(self):
        if socket.gethostname()=='gmi-lws12':
            pass
        else:
            pass
"""
 
class Job(object):
    def __init__(self,commands=None, modules=None, analysis_step=None,append_to_ana=True, id='', depends=None, input_files=None, output_files=None, walltime='08:00:00',ncpus=1, mem=None, debug=False):
        self.commands = ([] if commands is None else commands)
        self.modules = ([] if modules is None else modules)
        self.depends = ([] if depends is None else depends)
        self.input_files = ([] if input_files is None else input_files)
        self.output_files = ([] if output_files is None else output_files)
        self.walltime = walltime
        self.ncpus = ncpus
        self.debug = debug
        self.pbs_id = None
        self.returncode = None
        if mem is None:
            self.mem = str(ncpus * 3825) + 'mb'
        self.id = str(id)
        if analysis_step is not None:
            self.bind_to_analysis_step(analysis_step)
        else:
            self.name = None
            self.analysis_step = None
            self.file_name = None
        if self.analysis_step is not None and append_to_ana:
            self.analysis_step.append_job(self)

    def bind_to_analysis_step(self,analysis_step):
        self.analysis_step = analysis_step
        self.analysis = analysis_step.analysis
        id = self.id
        self.name = str(analysis_step.name) + ('_' if len(id)>0 else '') + id
        self.file_name = os.path.join(self.analysis.project,self.analysis.ana_dir,"jobscript/",self.name+".sh")   
        self.oe_fn=os.path.join(self.analysis.ana_dir,"log/",self.name)
        if self.debug:
            self.oe_fn+='_debug'

    def write_jobscript(self):
        commands = self.commands
        modules = self.modules
        id = self.id
        walltime = self.walltime
        ncpus = self.ncpus
        debug = self.debug
        mem = self.mem

        with open(os.path.expanduser(self.file_name), "w") as jf:
            jf.write("#!/bin/bash\n")
            sn = self.analysis_step.short_name
            name = (sn if len(sn)<=3 else sn[:3]) + ('_' if len(id)>0 else '') + (id if len(id) <=7 else id[:7])
            #print 'name:',name
            jf.write("#PBS -N {}\n".format(name))
            jf.write("#PBS -P vervetmonkey\n")
            if debug:
                jf.write("#PBS -q debug\n")
            else:
                jf.write("#PBS -l mem={}\n".format(mem))
                jf.write("#PBS -l ncpus={}\n".format(ncpus))
                jf.write("#PBS -l walltime={}\n".format(walltime))
            jf.write("#PBS -o {}.o\n".format(os.path.join(os.path.expanduser(self.analysis.project),self.oe_fn)))
            jf.write("#PBS -e {}.e\n".format(os.path.join(os.path.expanduser(self.analysis.project),self.oe_fn)))

            #load modules:
            for module in modules:
                jf.write("module load {}\n".format(module))
            #command:
            for command in commands:
                jf.write(command)
                jf.write('\n')
        self.chmod_jobscript()
        #return jfn
    
    def chmod_jobscript(self):
        p = subprocess.call(['chmod','ug+x',os.path.expanduser(self.file_name)])
        #out, err = p.communicate()
        

    def qsub_jobscript(self):
        depend_str=''
        if self.depends:
            for depend in self.depends:
                if type(depend) == str:
                    pbs_id = depend.strip()
                else:
                    if depend.pbs_id == None:
                        raise Exception("{0} depending on {1}. Qsub {0} before submitting {1}".format(depend.analysis_step.name+'_'+depend.id,self.analysis_step.name+'_'+self.id))
                    pbs_id = depend.pbs_id.strip()     
                depend_str = depend_str + (':' if len(depend_str)>0 else '') + pbs_id 
            command = 'qsub -W depend=afterok:{0} {1}'.format(depend_str,os.path.expanduser(self.file_name))
            P = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        else:
            command = 'qsub {0}'.format(os.path.expanduser(self.file_name))
            p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out, err = p.communicate()
        self.pbs_id = out.strip() 
        self.analysis.append_submit_log(command, out, err)
        return self.pbs_id        

    def run_jobscript(self):
        command = os.path.expanduser(self.file_name)
        self.execute(command)
        
    def execute(self,command):
        for depend_job in self.depends:
            if depend_job.returncode != 0:
                raise Exception("{0} finished with exit code {1}, won't start {2}".format(depend_job.name,depend_job.returncode,self.name))
        p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out, err = p.communicate()
        self.returncode = p.returncode
        for name, content in zip(['o','e'],[out,err]): 
            with open(os.path.expanduser(os.path.join(self.analysis.scratch,self.oe_fn+'.'+name)),'w') as oef:
                oef.write(content)
        

    def print_summary(self):
        print "-"*50
        print self.name
        print "depends on:", [job.name for job in self.depends]
        print "input files:", self.input_files
        print "output_files:", self.output_files 
        print "-"*50
        

class StageJob(Job):
    def __init__(self,  direction, files=None, analysis_step=None, stage_analysis_dir=True, depends=None,  verbose=False, debug=False):
        if direction not in ['in','out']:
            raise ValueError('direction should be "in" or "out" but is {}'.format(direction))
        self.direction = direction
        self.files = ([] if files is None else files)
        self.depends = ([] if depends is None else depends)
        self.debug = debug
        self.stage_analysis_dir = stage_analysis_dir
        if analysis_step is not None:
            self.bind_to_analysis_step(analysis_step)
        else:
            self.name = None
            self.analysis_step = None    
            self.file_name = None
            self.id = None
            self.stage_fn = None
        self.verbose = verbose
        self.input_files = None
        self.output_files = None    
        
    def bind_to_analysis_step(self,analysis_step):
        if analysis_step.analysis.host == 'gmi-lws12':
            hid = 'lws12'
        elif 'login' in  analysis_step.analysis.host or 'dmn' in analysis_step.analysis.host:
            hid = 'mendel'
        self.id = 'stage' + self.direction + '_'  + hid
        self.stage_fn = os.path.join(analysis_step.analysis.project,analysis_step.analysis.ana_dir,'io/',analysis_step.name+'.'+self.direction)
        Job.bind_to_analysis_step(self,analysis_step)
        self.scratch = ('~/vervet_lab' if hid == 'lws12' else '~/vervet_scratch')
        if self.stage_analysis_dir:
            self.files.insert(0,'analyses/'+analysis_step.analysis.name+'/')


    def write_stage_fn_file(self):
        with open(os.path.expanduser(self.stage_fn),'w') as f:
            for fn in self.files:
                f.write(fn+'\n')    

    def write_jobscript(self, mode='rsync',write_fn_file=True):
        if write_fn_file:
            self.write_stage_fn_file()
        # how should this function know where the dirs are -> check self.analysis.host !
        modes = ['rsync']
        source_dir = (self.analysis.project if self.direction == 'in' else self.scratch)
        target_dir = (self.scratch if self.direction == 'in' else self.analysis.project)

        mem = '3825mb'
        ncpus = 1
        walltime = '04:00:00'
        if mode not in modes:
            raise ValueError('mode must be in {0} but is {1}'.format(modes, mode))
        jfn = os.path.expanduser(self.file_name)
        with open(os.path.expanduser(self.file_name),'w') as jf: 

            jf.write("#!/bin/bash\n")
            id = self.id
            sn = self.analysis_step.short_name
            name = (sn if len(sn)<=3 else sn[:3]) + ('_' if len(id)>0 else '') + (id if len(id) <=7 else id[:7])
            jf.write("#PBS -N {0}\n".format(name))
            jf.write("#PBS -q staging\n")
            jf.write("#PBS -P vervetmonkey\n")
            jf.write("#PBS -o {0}.o\n".format(os.path.join(os.path.expanduser(self.analysis.project),self.oe_fn)))
            jf.write("#PBS -e {0}.e\n".format(os.path.join(os.path.expanduser(self.analysis.project),self.oe_fn)))
            jf.write("#PBS -l mem={0}\n".format(mem))
            jf.write("#PBS -l ncpus={0}\n".format(ncpus))
            jf.write("#PBS -l walltime={0}\n".format(walltime))
            commands = ['echo start','date','date1=$(date +"%s")',
                        'stage_fn={0}'.format(self.stage_fn),"rsync --files-from=$stage_fn -aur {source} {target}".format(source=source_dir, target=target_dir),
                        'date2=$(date +"%s")',
                        'diff=$(($date2-$date1))',
                         """if [ "$diff" -lt 10 ]
then
sleep 10   
fi""",'echo end','date']
            for command in commands:
                jf.write(command)
                jf.write('\n')
        self.chmod_jobscript()
        if self.verbose:
            print jfn, 'written'
            
    def local_check_if_staging_needed(self):
        """
        check only whether files exist on target
        don't compare modification dates
        """    
        pass
    

    def run_jobscript(self):
        if 'dmn' not in self.analysis.host:
            if 'login' not in self.analysis.host:
                raise Exception('Staging from non-mendel host not yet implemented')
            if self.verbose:
                print 'executing', self.file_name, 'at dmn trough ssh'
            command = 'ssh dmn.mendel.gmi.oeaw.ac.at nohup {}'.format(self.file_name)
        else:
            command = '{}'.format(self.file_name)
        self.execute(command)
                   
                    
                      
