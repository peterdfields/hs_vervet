#!/usr/bin/env python
"""
This python module provides classes for a very basic work flow system 
that allows to run analysis on a cluster and on local linux workstation.
It is developped for GMI's mendel cluster with PBS scheduling, but the code can easily
be adapted for for other clusters. 
Contact: hannes@svardal.at

Objects: 
--------------------
analysis
    ---> step0
            --> job0
            --> job1
            --> job2
                ...
    ---> step2
            --> job0
            --> job1
                ...
        ...
-------------------
The module provides 3 classes: Analysis, Step, Job
Each analysis is expected to contain several steps and each step one or more jobs.
One can think of the analysis as the workflow of an experiment that logically belongs together.
For instance, an analyis could be variant discovery for next generation sequencing. 
Steps are then discrete steps of this experiment. In example of variant discovery,
steps could be: read_mapping, initial_variant_calling, local_realignment, variant_calling, filtering, ...
Each step consists of several jobs (for instance for each individual or chromosme) and each job can have several commands.
A job can be seen as the unit that will be submitted to the cluster (qsub) or run locally.

Basic usage:
ana = Analysis(name="test_analysis")
first_step = Step(name="hello_world",analysis=ana)
first_step_first_job = Job(commands=['echo Hello World'],step=first_step)
first_step_second_job = Job(commands=['cd {}'.format(ana.ana_dir),'echo Hello world to file > testfile.txt'],step=first_step)
first_step.run(mode='qsub')

See each class docstring for more information on usage.
Note: By default the jobs are automatically appended to the step if the step keyword argument is given.
It is actually not necessary to assign the jobs to variables.



Features:
-- Creates a file hierarchy for the analysis on project and scratch.
-- Automatic staging between project and scratch.
-- Can make use of parallel execution locally.
-- Same code can run on cluster and linux workstations.
-- Creates jobscripts in <ana_dir>/jobscripts for easy inspection what is going on


Note: Get newest version of this script by pulling from the repository (if you have cloned the hs_vervet_repository before): 
cd <hs_vervet-dir> && git pull origin master
Or from python: subprocess.Popen("cd <hs_vervet_dir> && git pull origin master".format(self.scratch),shell=True)


Todo:
Make it possible to merge analysis steps.
Fix stageout to first submit a diagnostic job that then submits the real stage jobs. 
More sophisticated staging methods for stage out.

"""
import sys, os, datetime, subprocess, socket, filecmp, shutil
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from hs_vervet.tools import hs_vervet_basics as hvb
from hs_vervet.tools.stage import local_prepare_staging




class Analysis(object):
    """
    This is the top class in the hirarchy. (See also module docstring.)

    Usage:
    analysis = Analysis(name,project="~/vervet_project",scratch="~/vervet_scratch")

    Where name is the name of the current analysis/experiment and project and scratch
    are the root directiories for the project on the project- and on the scratch-
    file-systems.
    Instantiation creates a folder hierarchy
    <project>/analysis/<name>/_data
                              script
                              jobscript
                              log
                              output
                              io
    and the same where <project> is replaced by <scratch>. The script which 
    instantiates this class is automatically copied to the "script"-folder.
    We recommend for analysis-naming to follow the convention name=
    "YYYYMMSS_easily_identifiable_name".
    
    Important methods:
    (see methods' docstring for more info)
    analysis.append_step(step)
        ... Appends a step to the analysis. Note that this is not necessary if
            the step is instantiated with the keyword "analysis=analysis".
            Then the step is automatically appended to the analysis.

    Todo:
    analysis.run_all()
    analysis.join_steps()
    
    """
    def __init__(self, name=None,description=None,project_dir_prefix="vervet",project_name="vervetmonkey",verbose=0):
        callingframe = sys._getframe(1)
        c = callingframe.f_locals["__file__"]
        self.calling_fn = os.path.join(os.getcwdu(),
                (c if c[:2]!='./' else c[2:]))
        self.name = name
        self.host = socket.gethostname()
        self.scratch = os.path.join("~", project_dir_prefix + "_scratch")
        self.project =  os.path.join("~", project_dir_prefix + "_project")
        self.dir_prefix = project_dir_prefix
        self.project_name = project_name 
        self.description = description
        self.ana_dir = os.path.join('analyses/',self.name)
        self.submit_log_fn = os.path.join(self.project,self.ana_dir,'log/', self.name+"_submit_log.txt") 
        for direc in ["_data","log","script","jobscript","io","output"]:
            hvb.try_make_dirs(os.path.join(self.project,self.ana_dir,direc))
            hvb.try_make_dirs(os.path.join(self.scratch,self.ana_dir,direc))
        self.steps=[]
        try:
            shutil.copy(self.calling_fn,os.path.join(os.path.expanduser(self.project),
                        self.ana_dir,"script",self.calling_fn.split("/")[-1]))
        except shutil.Error, e:
            if "are the same file" in str(e.message): 
                pass
            else:
                raise
        self.verbose = verbose
        self.vprint = lambda text, min_verb=0: v_print(text,min_verb,self.verbose)
    
            
            
        #subprocess.Popen("cd {}/script/hs_vervet && git pull origin master".format(self.scratch),shell=True)

    def append_step(self, step):
        self.steps.append(step)
        step.bind_to_analysis(self)
        
    #def join_steps(self,join_jobs=False,join_jobs_on='id'):
    #    pass

    def append_submit_log(self, jobfile, out, err):
        with open(os.path.expanduser(self.submit_log_fn),'a') as lf:
            date = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
            lf.write(date+'\t'+jobfile+'\t'+out+'\t'+err+'\n')
    def write_and_qsub_all(self):
        for step in self.steps:
            for job in step.jobs:
                job.write_jobscript()
                job.qsub_jobscript()
    def join_steps(steps=None,all=False):
        pass
            
    


class Step(object):

    def __init__(self, analysis=None, name=None, jobs=None, append_to_ana=True, description=None, depend_steps=None, stagein=True, stageout=True, default_run='qsub',  check_date_input = True, verbose = None):
        # the next to lines are just for the following methods to work
        # the attributes will be set again by their constructors
        self.name = name
        self.jobs = ([] if jobs is None else jobs)
        self.stagein = stagein
        self.stageout = stageout
        self.default_run = default_run
        self.stagein_job = None
        self.stageout_job = None
        self.description = None
        self.verbose = verbose
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
        if self.verbose is None:
            self.verbose=analysis.verbose
        self.vprint = lambda text, min_verb=0: v_print(text,min_verb,self.verbose)
        
    def append_job(self, job):
        #only add the job if it is not added yet
        if job not in self.jobs:
            self.jobs.append(job)
            job.bind_to_step(self)
            for depend_step in self.depend_steps:
                for depend_job in depend_step.jobs:
                    job.depends.append(depend_job)
        else:
            UserWarning("Wont append job {0} to step {1}, because it is already in {1}'s job-list".format(job.name,self.name))
                  
    def append_jobs(self,jobs):
        for job in jobs:
            self.append_job(job)


    def join(self,step,mode='append',join_jobs=False,join_jobs_on='id'):
        """
        step will be appended to self
        """
        modes = ['append','prepend']
        if mode not in modes:
            raise ValueError('mode shoud be in {0} but is {1}'.format(modes,modes))
        if mode == 'append':
            self.append_jobs(step.jobs)
        elif mode == 'prepend':
            raise UserException('Not implemented.')

        if join_jobs:
            self.join_jobs(join_on=join_jobs_on)
        if self.analysis is not None:
            self.analysis.steps.remove(step)

    def join_jobs(self,join_on='id'):
        """
        jobs will be joined and commands joined from left to right
        """
        join_ons=['id','all']
        if join_on not in join_ons:
            raise ValueError('join_on shoud be in {0} but is {1}'.format(join_ons,join_on))
        if join_on == 'all':
            for job in self.jobs[1:]:
                self.jobs[0].join(job)
        elif join_on == 'id':
            ids = [job.id for job in self.jobs]
            unique_ids = list(set(ids))
            job_sets = [[job for job in self.jobs if job.id == id] for id in unique_ids]
            for js in job_sets:
                for job in js[1:]:
                    js[0].join(job)
                
            
        

    def run(self,mode=None,print_summary=False,staging=True,parallel=False,nprocs='auto'):
        modes=['qsub','scratch_run','write_jobscripts','project_run']
        if mode is None:
            mode = self.default_run
        if mode not in modes:
            raise Exception("mode must be in {0} but is {1}".format(modes,mode))
        if mode == 'project_run':
            staging = False
        if staging == True:
            #print "stage_job_created"
            self.prepare_staging()
        if print_summary:
            self.vprint_summary()
        if mode == 'qsub':
            self.qsub()
        elif mode == 'scratch_run':
            self.scratch_run(parallel=parallel,nprocs=nprocs)
        elif mode == 'write_jobscripts':
            self.write_jobscripts()
        elif mode == 'project_run':
            self.project_run(parallel=parallel,nprocs=nprocs)
    
    def prepare_staging(self):
        if self.stagein: #here we could also check if stagein_job already exists
            self.add_stagein()
        if self.stageout:
            self.add_stageout()
    
    #combine all the following methods to avoid parallelism!!!
            
    def write_jobscripts(self):    
        if self.stagein_job is not None:
            self.stagein_job.stage(run_type='dry_run')
        for job in self.jobs:
            job.commands.insert(0,'PROJECT_HOME=' + self.analysis.scratch)
            job.write_jobscript()
        if self.stageout_job is not None:
            self.stageout_job.stage(run_type='dry_run')

    def scratch_run(self,parallel=False,nprocs='auto'):
        if self.stagein_job is not None:
            self.stagein_job.stage(run_type='direct')
        for job in self.jobs:
            job.commands.insert(0,'PROJECT_HOME=' + self.analysis.scratch)
            job.write_jobscript()
        if parallel:
            if nprocs == 'auto':
                hvb.parmap(lambda j: j.run_jobscript(),self.jobs)
            else:
                hvb.parmap(lambda j: j.run_jobscript(),self.jobs,nprocs=nprocs)
        else:
            for job in self.jobs:
                job.run_jobscript()
        if self.stageout_job is not None:
            self.stageout_job.stage(run_type='direct')

    def project_run(self,parallel=False,nprocs='auto'):   
        for job in self.jobs:
            job.commands.insert(0,'PROJECT_HOME='+self.analysis.project)
            job.write_jobscript()
        if parallel:
            if nprocs == 'auto':
                hvb.parmap(lambda j: j.run_jobscript(),self.jobs)
            else:
                hvb.parmap(lambda j: j.run_jobscript(),self.jobs,nprocs=nprocs)
        else:
            for job in self.jobs:
                job.run_jobscript()
        

    def qsub(self):
        #hold makes sure that the depend jobs are seen by the dependent jobs
        
        if self.stagein_job is not None:
            self.stagein_job.stage(run_type='submit')
        for job in self.jobs:
            job.commands.insert(0,'PROJECT_HOME=' + self.analysis.scratch)
            job.write_jobscript()
            job.qsub_jobscript()
        if self.stageout_job is not None:
            self.stageout_job.write_prepare_jobscript()
            self.stageout_job.qsub_prepare_jobscript()
            #self.stageout_job.release()
        for job in self.jobs:
            job.release()
        if self.stagein_job is not None:
            self.stagein_job.release()
    
    def print_summary(self,job_summary=True):
        print "="*60
        print "analysis:", (None if self.analysis is None else self.analysis.name)
        print "step:", ('unnamed' if self.name is None else self.name)
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
        stagein_job = StageJob(direction='in',files=in_files,step=self,stage_analysis_dir=stage_analysis_dir)
        for job in self.jobs:
            job.depends.append(stagein_job)
        self.stagein_job = stagein_job
  

    def add_stageout(self, stage_analysis_dir=True, add_to_jobs=False):
        out_files = []
        for job in self.jobs:
            for file in job.output_files:
                out_files.append(file)
        out_files=list(set(out_files)) #only unique files
        stageout_job = StageJob(direction='out',files=out_files,step=self,stage_analysis_dir=stage_analysis_dir,depends=[job for job in self.jobs])
        self.stageout_job = stageout_job       
        
    def remove_correctly_finished_jobs(self):
        self.jobs = [job for job in self.jobs if not job.ran_noerror()]


#for backwards compatability
AnalysisStep = Step

class JoinedStep(Step):
    def __init__(self):
        pass

 
class Job(object):
    def __init__(self,commands=None, modules=None,cluster_modules=None,local_modules=None,cluster_commands=None, local_commands=None, step = None, analysis_step = None,append_to_ana=True, id='', depends=None, input_files=None, output_files=None, walltime='08:00:00',ncpus=1, mem=None, exit_on_error=True, description=None, debug=False):
        self.depends = ([] if depends is None else depends)
        self.input_files = ([] if input_files is None else input_files)
        self.output_files = ([] if output_files is None else output_files)
        self.walltime = walltime
        self.ncpus = ncpus
        self.debug = debug
        self.pbs_id = None
        self.returncode = None
        self.exit_on_error = exit_on_error
        self.description = description
        if mem is None:
            self.mem = str(ncpus * 3825) + 'mb'
        else:
            self.mem = mem
        self.id = str(id)
        

        if step is not None:
            self.bind_to_step(step)
            if self.step.analysis is not None:
                host = self.analysis.host
            else:
                host = socket.gethostname()
        #for backwards comatability
        #analysis_step is depriciated, use step
        elif analysis_step is not None:
            self.bind_to_step(analysis_step)
            if self.step.analysis is not None:
                host = self.analysis.host
            else:
                host = socket.gethostname()
        else:
            self.name = None
            self.step = None
            self.file_name = None
            host = socket.gethostname()
        cluster_commands = ([] if cluster_commands is None else cluster_commands)
        local_commands = ([] if local_commands is None else local_commands)
        self.commands = ["cd $PROJECT_HOME"] + (local_commands if 'lws12' in host else cluster_commands) +  ([] if commands is None else commands)
        cluster_modules = ([] if cluster_modules is None else cluster_modules)
        local_modules = ([] if local_modules is None else local_modules)
        self.modules = ([] if modules is None else modules) + (local_modules if 'lws12' in host else cluster_modules)
        if self.step is not None and append_to_ana:
            self.step.append_job(self)

    def bind_to_step(self,step):
        self.step = step
        self.analysis = step.analysis
        id = self.id
        self.name = str(step.name) + ('_' if len(id)>0 else '') + id
        self.file_name = os.path.join(self.analysis.project,self.analysis.ana_dir,"jobscript/",self.name+".sh")   
        self.oe_fn=os.path.join(self.analysis.ana_dir,"log/",self.name)
        if self.debug:
            self.oe_fn+='_debug'

    def join(self,job):
        """
        appends job commands and modules to self
        """
        self.modules += job.modules
        self.modules = list(set(self.modules))
        self.commands += job.commands
        self.input_files += job.input_files
        self.output_files += job.output_files
        if self.step is not None:
            self.step.jobs.remove(job)

    
    def ran_noerror(self):
        f = os.path.join(os.path.expanduser(self.analysis.project),self.oe_fn+'.e')
        if os.path.isfile(f) and os.path.getsize(f) == 0:
            ran_noerror = True
        else:
            ran_noerror = False
        return ran_noerror


    def write_jobscript(self):
        
        commands = self.commands
        modules = self.modules 
        id = self.id
        walltime = self.walltime
        ncpus = self.ncpus
        debug = self.debug
        try:
            mem = self.mem
        except:
            print self.name
            raise

        with open(os.path.expanduser(self.file_name), "w") as jf:
            jf.write("#!/bin/bash\n")
            sn = self.step.short_name
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

            if self.exit_on_error:
                jf.write("set -e\n")
            #load modules:
            for module in modules:
                jf.write("module load {}\n".format(module))
            #command:
            for command in commands:
                #try:
                jf.write(command)
                #except:
                #    print 'command:',  command
                #    print 'commands:', commands
                #    raise
                jf.write('\n')
        self.chmod_jobscript()
        #return jfn
    
    def chmod_jobscript(self,file='auto'):
        if file=='auto':
            file = self.file_name
        p = subprocess.call(['chmod','ug+x',os.path.expanduser(file)])
        #out, err = p.communicate()
        

    def qsub_jobscript(self):
        depend_str=''
        if self.depends:
            for depend in self.depends:
                if type(depend) == str:
                    pbs_id = depend.strip()
                else:
                    if depend.pbs_id == None:
                        raise Exception("{0} depending on {1}. Qsub {0} before submitting {1}".format(depend.step.name+'_'+depend.id,self.step.name+'_'+self.id))
                    pbs_id = depend.pbs_id.strip()     
                depend_str = depend_str + (':' if len(depend_str)>0 else '') + pbs_id 
            command = 'qsub -h  -W depend=afterok:{0} {1}'.format(depend_str,os.path.expanduser(self.file_name))
            p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        else:
            command = 'qsub -h {0}'.format(os.path.expanduser(self.file_name))
            p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out, err = p.communicate()
        self.pbs_id = out.strip() 
        self.analysis.append_submit_log(command, out, err)
        return self.pbs_id        

    def release(self):
        command = "qrls {0}".format(self.pbs_id)
        #print "releasing", self.name, command
        p = subprocess.Popen(command, shell=True)

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
        print "job filename:", self.file_name
        print "oe filenames:", self.oe_fn
        print "input files:", self.input_files
        print "output_files:", self.output_files 
        print "-"*50
        

class StageJob(Job):
    def __init__(self,  direction, files=None, step=None, stage_analysis_dir=True, mode='newer',depends=None,  verbose=None, debug=False):
        if direction not in ['in','out']:
            raise ValueError('direction should be "in" or "out" but is {}'.format(direction))
        self.direction = direction
        self.files = ([] if files is None else files)
        self.depends = ([] if depends is None else depends)
        self.debug = debug
        self.stage_analysis_dir = stage_analysis_dir
        self.mode = mode
        self.verbose = verbose
        if step is not None:
            self.bind_to_step(step)
        else:
            self.name = None
            self.step = None    
            self.file_name = None
            self.id = None
            self.stage_fn = None
        self.input_files = None
        self.output_files = None    
        
    def bind_to_step(self,step):
        if step.analysis.host == 'gmi-lws12':
            hid = 'lws12'
        elif 'login' in  step.analysis.host or 'dmn' in step.analysis.host:
            hid = 'mendel'
        self.id = 'stage' + self.direction + '_'  + hid
        self.stage_fn = os.path.join(step.analysis.project,step.analysis.ana_dir,'io/',step.name+'.'+self.direction)
        Job.bind_to_step(self,step)
        self.scratch = ('~/vervet_lab' if hid == 'lws12' else '~/vervet_scratch')
        if self.stage_analysis_dir:
            self.files.insert(0,'analyses/'+step.analysis.name+'/')
        if self.verbose is None:
            self.verbose = step.verbose
        self.vprint = lambda text, min_verb=0: v_print(text,min_verb,self.verbose)
        #choose appropriate name and use this, attention with different file systems...
        #-> different files for local and remote ...
        #print self.oe_fn
        #print self.analysis.project
        self.local_output = os.path.join(self.analysis.project, self.oe_fn+'_local.o')

    def stage(self,run_type='auto'):
        # todo: incorporate depends!
        if self.analysis.host == 'gmi-lws12':
            partner = 'lab'
        else:
            partner = 'scratch'
        id = self.id
        sn = self.step.short_name
        name = (sn if len(sn)<=3 else sn[:3]) + ('_' if len(id)>0 else '') + (id if len(id) <=7 else id[:7])
        if run_type=='dry_run' or run_type=='direct':
            depends = None
        else:
            depends = [job.pbs_id.strip() for job in self.depends]
        #print self.local_output
        (out, err, rc) = local_prepare_staging(self.files,partner,self.direction,self.mode,run_type=run_type,afterok=depends,startonhold=True,job_fn=self.file_name,out_fn=os.path.join(self.analysis.project,self.oe_fn),job_name=name,project=self.analysis.dir_prefix,verbose=self.verbose,file_to_print_to=self.local_output)
        #if out is not None and out:
        #    print >>sys.stdout, 'stage.local_prepare_staging','out:',out 
        #if err is not None:
        #    print >>sys.stderr, 'stage.local_prepare_staging','err:',err
        self.returncode = rc
        #print "ana_stage, out",out
        #print "ana_stage,err", err
        if run_type == 'submit':
            #print self.name,'submitted with pbs_id', out
            #print self.name
            self.pbs_id = out.strip()
        
    #the following two functions are only used by the stage-out job
    #this is basically a dummy job that tests first how much files really need to be staged and
    #then submits or directly runs a job accordingly
    def write_prepare_jobscript(self):
        mem = '3825mb'
        ncpus = 1
        walltime = '00:20:00'
        id = self.id
        sn = self.step.short_name
        name = (sn if len(sn)<=3 else sn[:3]) + ('_' if len(id)>0 else '') + (id if len(id) <=7 else id[:7])
#        if mode not in modes:
#            raise ValueError('mode must be in {0} but is {1}'.format(modes, mode))
#        jfn = os.path.expanduser(self.file_name)
        fn = os.path.splitext(self.file_name)[0] + "_prep" + os.path.splitext(self.file_name)[1]
        with open(os.path.expanduser(fn),'w') as jf: 
            jf.write("#!/bin/bash\n")
            id = self.id
            sn = self.step.short_name
            name = (sn if len(sn)<=3 else sn[:3]) + ('_' if len(id)>0 else '') + (id if len(id) <=7 else id[:7])
            jf.write("#PBS -N {0}\n".format(name))
            jf.write("#PBS -q staging\n")
            jf.write("#PBS -P vervetmonkey\n")
            jf.write("#PBS -o {0}_prep.o\n".format(os.path.join(os.path.expanduser(self.analysis.project),self.oe_fn)))
            jf.write("#PBS -e {0}_prep.e\n".format(os.path.join(os.path.expanduser(self.analysis.project),self.oe_fn)))
            jf.write("#PBS -l mem={0}\n".format(mem))
            jf.write("#PBS -l ncpus={0}\n".format(ncpus))
            jf.write("#PBS -l walltime={0}\n".format(walltime))
            jf.write("module load Python/2.7.3-goolf-1.4.10\n")
            commands = ["stage -p scratch -d out -m {mode} -b {project} -v {verbose} -l {print_to} -j {job_fn} -o {out_fn} -n {job_name} {files}".format(
                mode=self.mode,
                project=self.analysis.dir_prefix, verbose=self.verbose, print_to=self.local_output,job_fn=self.file_name,
                out_fn=os.path.join(self.analysis.project,self.oe_fn),job_name=name,files=" ".join(self.files))]            
            for command in commands:
                jf.write(command)
                jf.write('\n')
        #implement the following...
        self.chmod_jobscript(fn)
        self.prep_file_name = fn
        return fn

        
    def qsub_prepare_jobscript(self):
        depend_str=''
        if self.depends:
            for depend in self.depends:
                if type(depend) == str:
                    pbs_id = depend.strip()
                else:
                    if depend.pbs_id == None:
                        raise Exception("{0} depending on {1}. Qsub {0} before submitting {1}".format(depend.step.name+'_'+depend.id,self.step.name+'_'+self.id))
                    pbs_id = depend.pbs_id.strip()
                depend_str = depend_str + (':' if len(depend_str)>0 else '') + pbs_id
            command = 'qsub  -W depend=afterok:{0} {1}'.format(depend_str,os.path.expanduser(self.prep_file_name))
            p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        else:
            command = 'qsub {0}'.format(os.path.expanduser(self.prep_file_name))
            p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out, err = p.communicate()
        self.pbs_id = out.strip()
        self.analysis.append_submit_log(command, out, err)
        return self.pbs_id

        


    def write_stage_fn_file(self):
        with open(os.path.expanduser(self.stage_fn),'w') as f:
            for fn in self.files:
                f.write(fn+'\n')    

#    def write_jobscript(self, mode='rsync',write_fn_file=True):
#        if write_fn_file:
#            self.write_stage_fn_file()
#        # how should this function know where the dirs are -> check self.analysis.host !
#        modes = ['rsync']
#        source_dir = (self.analysis.project if self.direction == 'in' else self.scratch)
#        target_dir = (self.scratch if self.direction == 'in' else self.analysis.project)
#
#        mem = '3825mb'
#        ncpus = 1
#        walltime = '04:00:00'
#        if mode not in modes:
#            raise ValueError('mode must be in {0} but is {1}'.format(modes, mode))
#        jfn = os.path.expanduser(self.file_name)
#        with open(os.path.expanduser(self.file_name),'w') as jf: 
#
#            jf.write("#!/bin/bash\n")
#            id = self.id
#            sn = self.step.short_name
#            name = (sn if len(sn)<=3 else sn[:3]) + ('_' if len(id)>0 else '') + (id if len(id) <=7 else id[:7])
#            jf.write("#PBS -N {0}\n".format(name))
#            jf.write("#PBS -q staging\n")
#            jf.write("#PBS -P vervetmonkey\n")
#            jf.write("#PBS -o {0}.o\n".format(os.path.join(os.path.expanduser(self.analysis.project),self.oe_fn)))
#            jf.write("#PBS -e {0}.e\n".format(os.path.join(os.path.expanduser(self.analysis.project),self.oe_fn)))
#            jf.write("#PBS -l mem={0}\n".format(mem))
#            jf.write("#PBS -l ncpus={0}\n".format(ncpus))
#            jf.write("#PBS -l walltime={0}\n".format(walltime))
#            commands = ['echo start','date','date1=$(date +"%s")',
#                        'stage_fn={0}'.format(self.stage_fn),"rsync --files-from=$stage_fn -aur {source} {target}".format(source=source_dir, target=target_dir),
#                        'date2=$(date +"%s")',
#                        'diff=$(($date2-$date1))',
#                         """if [ "$diff" -lt 10 ]
#then
#sleep 10   
#fi""",'echo end','date']
#            for command in commands:
#                jf.write(command)
#                jf.write('\n')
#        self.chmod_jobscript()
#        if self.verbose:
#            print jfn, 'written'
            
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

class Command(object):
    def __init__(self,command,description):
        self.command = command
        self.description = description
                   
                    
                      
