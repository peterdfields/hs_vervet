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
For instance, an analysis could be variant discovery for next generation sequencing. 
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
More sophisticated staging methods for stage out.

Generic use for different options:
scenarios:
----
1)
host = mendel
run = qsub
staging = yes
status = works
------------
2)
host = mendel
run = direct_run_mendel
staging = yes
status = test it!
------------
3)
host = mendel
run = direct_run_mendel
staging = no
status = test it!
4)
host = mendel
run = run on lws
staging = yes
status = implement
5)
host = lws
run = run on lws
staging = yes
status = this is the most tricky, because staging should happen from mendel; the project has to be created on mendel somehow
6)
host = lws
run = run on lws
staging = no
status = this should be easy to implement, everythin happend on lws, no staging at all

Note: 4 and 5 should do exactly the same thing, only that it is executed from mendel once and on lws the other time.
Maybe we should implement 4 and not offer 5, as a first approach? Or only 5?

only offer 1,3,5,6?

#implement check whether input exists for no-staging mode

The analysis.config should be read by the analysis directly at runtime!??

"""
import sys, os, datetime, subprocess, socket, filecmp, shutil
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from hs_vervet.tools import hs_vervet_basics as hvb




#some global variables
#max walltime in seconds
max_walltime = 48*3600



def pretty_comment(comment,title,sep='-'):
    fill_to = 60
    #offset = 5
    head = "#" + sep*(fill_to-1) + '\n' + "#" + title + '\n' #sep*offset + title + sep*(fill_to-1-offset-len(title)) + '\n'
    body = ''
    for line in comment.split("\n"):
        if line:
            body += "#" + line + '\n' # sep*(fill_to-1-len(line)) + '\n'
    body += "#" + sep*(fill_to-1) + '\n'
    return head+body

def make_ana_dirs(base_dir,ana_dir):
    for direc in ["_data","log","script","jobscript","io","output"]:
        hvb.try_make_dirs(os.path.join(base_dir,ana_dir,direc))

def is_cluster():
    #check whether on mendel
    try:
        os.environ["BC_HPC_SYSTEM"]
        return True
    except KeyError:
        return False

def value_check(el,lst):
    if el not in lst:
        raise ValueError("Argument must be in {0} but is {1}".format(lst,el))

def create_cluster_ana_folder(project_dir,ana_dir):
            command = '''ssh dmn.mendel.gmi.oeaw.ac.at nohup "python -c 'from hs_vervet.tools.analysis import make_ana_dirs; make_ana_dirs(\\"{}\\",\\"{}\\")'"'''.format(project_dir,ana_dir)
            #print 'bla:',command
            subprocess.Popen(command,shell=True)

class BaseClass(object):
    """
    This template is used as parent of Analysis, Step and Job.
    """
    def vprint(self,*text,**kwa):
        """
        use similar to python3 print function
        additional argument mv or min_verbosity
        and append (to know whether append to file or not)
        """
        kwa.update({"verbosity":self.verbose})
        hvb.v_print(*text,**kwa)



class Analysis(BaseClass):
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
    consistent local and mendel usage with optionall staging...
    """
    def __init__(self, name, default_run = None, project_dir="~/vervet_project",scratch_dir="~/vervet_scratch",lab_dir="~/vervet_lab",project_name="vervetmonkey",description=None,verbose=0):

        callingframe = sys._getframe(1)
        c = callingframe.f_locals["__file__"]
        self.calling_fn = os.path.join(os.getcwdu(),
                (c if c[:2]!='./' else c[2:]))
        self.name = name
        self.host = 'cluster' if is_cluster() else 'workstation'
        if default_run is None:
            if self.host == 'cluster':
                default_run = 'qsub'
            elif self.host == 'workstation':
                default_run = 'scratch_run'
        default_runs = ['qsub','scratch_run','project_run','write_jobscripts']
        value_check(default_run, default_runs)
        self.default_run = default_run
        #print self.default_run
        self.scratch = scratch_dir
        self.project =  project_dir
        self.lab_dir = lab_dir
        self.project_name = project_name 
        self.description = (description if description is not None else '')
        self.ana_dir = os.path.join('analyses/',self.name)
        self.submit_log_fn = os.path.join(self.project,self.ana_dir,'log/', self.name+"_submit_log.txt") 
        make_ana_dirs(self.project,self.ana_dir)
        make_ana_dirs(self.scratch,self.ana_dir)
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
            
            
            
        #subprocess.Popen("cd {}/script/hs_vervet && git pull origin master".format(self.scratch),shell=True)

    def append_step(self, step):
        self.vprint("appending step",step.name,"to analysis",self.name,mv=1)
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

    def join_steps(self,steps,name=None,join_jobs_on='id'):
        #check whether steps list is a sublist of analysis.steps
        def position_sublist(lst, sublst):
            n = len(sublst)
            return [(sublst == lst[i:i+n]) for i in xrange(len(lst)-n+1)].index(True)
        try:
            pos = position_sublist(self.steps,steps)
        except ValueError:
            raise Exception("The steps to join must be a sublist of analysis.steps")
        if name == None:
            name = 'joined' + '_' + '_'.join([step.name for step in steps])
        js = JoinedStep(steps,name,self,join_jobs_on=join_jobs_on)
        self.steps = self.steps[:pos] + [js] + self.steps[pos+len(steps):]       
        return js
            
            
    


class Step(BaseClass):

    def __init__(self, analysis=None, name=None, jobs=None, append_to_ana=True, description=None, depend_steps=None, stagein=True, stageout=True, default_run=None,  check_date_input = True, verbose = None):
        # the next to lines are just for the following methods to work
        # the attributes will be set again by their constructors
        self.name = name
        self.jobs = ([] if jobs is None else jobs)
        self.stagein = stagein
        self.stageout = stageout
        self.default_run = default_run
        self.stagein_job = None
        self.stageout_job = None
        self.description = (description if description is not None else '')
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
        if self.default_run is None:
            self.default_run = analysis.default_run
        
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

    def run(self,mode=None,print_summary=False,parallel=False,nprocs='auto'):
        modes=['qsub','scratch_run','write_jobscripts','project_run']
        if mode is None:
            mode = self.default_run
        if mode not in modes:
            raise Exception("mode must be in {0} but is {1}".format(modes,mode))
        if mode !=  'project_run':
            self.prepare_staging()
        if print_summary:
            self.print_summary()
        if mode == 'qsub':
            if self.analysis.host == 'workstation':
                raise Exception("mode qsub not available on workstation. Run on the cluster.")
            self.run_type = 'qsub'
            self.qsub()
        elif mode == 'scratch_run':
            if self.analysis.host == 'workstation':
                # if on local workstation and staging enabled, create analysis folder hierarchy on cluster
                self.vprint("Creating analysis folder hierarchy on cluster.",mv=1)       
                create_cluster_ana_folder(self.analysis.project,self.analysis.ana_dir)
            self.run_type = 'run'
            self.scratch_run(parallel=parallel,nprocs=nprocs)
        elif mode == 'write_jobscripts':
            self.run_type = None
            self.write_jobscripts()
        elif mode == 'project_run':
            self.run_type = 'run'
            self.project_run(parallel=parallel,nprocs=nprocs)
    
    def run2(self,mode=None,print_summary=False,parallel=False,nprocs='auto'):
        modes=['qsub','scratch_run','write_jobscripts','project_run']
        if mode is None:
            mode = self.default_run
        if mode not in modes:
            raise Exception("mode must be in {0} but is {1}".format(modes,mode))
        if mode !=  'project_run':
            self.prepare_staging()
        if print_summary:
            self.print_summary()
        if mode == 'qsub':
            if self.analysis.host == 'workstation':
                raise Exception("mode qsub not available on workstation. Run on the cluster.")
            self.run_type = 'qsub'
            self.qsub()
        elif mode == 'scratch_run':
            if self.analysis.host == 'workstation':
                # if on local workstation and staging enabled, create analysis folder hierarchy on cluster
                self.vprint("Creating analysis folder hierarchy on cluster.",mv=1)       
                create_cluster_ana_folder(self.analysis.project,self.analysis.ana_dir)
            self.run_type = 'run'
            self.scratch_run(parallel=parallel,nprocs=nprocs)
        elif mode == 'write_jobscripts':
            self.run_type = None
            self.write_jobscripts()
        elif mode == 'project_run':
            self.run_type = 'run'
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
            job.write_jobscript()
        if self.stageout_job is not None:
            self.stageout_job.stage(run_type='dry_run')

    def scratch_run(self,parallel=False,nprocs='auto'):
        if self.stagein_job is not None:
            self.stagein_job.stage(run_type='auto')
        for job in self.jobs:
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
            self.stageout_job.stage(run_type='auto')
            #self.vprint("Attention stagout might be submitted to the cluster. Check whether it finished.",mv=1)

    def project_run(self,parallel=False,nprocs='auto'):   
        for job in self.jobs:
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
        
        if self.stageout_job is not None:
            self.stageout_job.write_prepare_jobscript()
        if self.stagein_job is not None:
            self.stagein_job.stage(run_type='submit')
        for job in self.jobs:
            job.write_jobscript()
            job.qsub_jobscript()
        if self.stageout_job is not None:
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
        self.vprint("adding stagein job to",self.name,mv=1)
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
        self.vprint("adding stageout job to",self.name,mv=1)
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
    def __init__(self,steps, name, analysis, join_jobs=True,join_jobs_on='id', append_to_ana=True, description=None, stagein=True, stageout=True, default_run='qsub', verbose = None):
        self.jobs = []
        self.name = name
        self.stagein = stagein
        self.stageout = stageout
        self.default_run = default_run
        self.stagein_job = None
        self.stageout_job = None
        self.description = (description if description is not None else '')
        self.verbose = verbose
        self.depend_steps = []

        if analysis is None:
            self.analysis = analysis
        else: 
            if append_to_ana:
                analysis.append_step(self)
            else:
                self.bind_to_analysis(analysis)
        self.add_name(name)
        if join_jobs:
            join_ons=['id','all']
            if join_jobs_on not in join_ons:
                raise ValueError('join_on shoud be in {0} but is {1}'.format(join_ons,join_on))
            if join_jobs_on == 'id':
                #check whether all steps have jobs with identical ids:
                ids = lambda step: [job.id for job in step.jobs]
                if not all(ids(s) == ids(steps[0]) for s in steps):
                    raise Exception("To join jobs on ids, all steps must have jobs with the same ids.")
                #join the jobs on the id
                jobs = [JoinedJob([[job for job in step.jobs if job.id == id][0] for step in steps]) for id in ids(steps[0])]
            elif join_jobs_on == 'all':
                jobs = [JoinedJob([job for step in steps for job in step.jobs])]
        else:
            raise Exception('Not joining jobs in a Joined Step is not implemented. Use join_jobs=True.')
            
        for job in jobs:
            self.append_job(job)



 
class Job(BaseClass):
    def __init__(self,commands=None, modules=None,cluster_modules=None,local_modules=None,cluster_commands=None, local_commands=None, step = None, analysis_step = None,append_to_ana=True, id='', depends=None, input_files=None, output_files=None, walltime='04:00:00',ncpus=1, mem=None, exit_on_error=True, description=None, debug=False, verbose=None):
        self.depends = ([] if depends is None else depends)
        self.input_files = ([] if input_files is None else input_files)
        self.output_files = ([] if output_files is None else output_files)
        self.walltime = walltime
        self.ncpus = ncpus
        self.debug = debug
        self.pbs_id = None
        self.returncode = None
        self.exit_on_error = exit_on_error
        self.description = (description if description is not None else '')
        self.verbose = verbose
        if mem is None:
            self.mem = str(ncpus * 3825) + 'mb'
        else:
            self.mem = mem
        self.id = str(id)
        

        if step is not None:
            self.bind_to_step(step)
        #for backwards comatability
        #analysis_step is depriciated, use step
        elif analysis_step is not None:
            UserWarning('analysis_step is depriciated, use step')
            self.bind_to_step(analysis_step)
        else:
            self.name = None
            self.step = None
            self.file_name = None
        
        d = {"commands":commands,"local_commands":local_commands,"cluster_commands":cluster_commands,"modules":modules,"cluster_modules":cluster_modules,"local_modules":local_modules}

        #embed single string commands in lists:
        for k,v in d.items():
            if type(v) == str or type(v) == Command:
                d[k]=[v]
                #print "changing", d,k,v

        #print d
        #print commands, cluster_commands,local_commands

        self.cluster_commands = ([] if d['cluster_commands'] is None else d['cluster_commands'])
        self.local_commands = ([] if d['local_commands'] is None else d['local_commands'])
        self.commands =  ([] if d['commands'] is None else d['commands'])
        self.cluster_modules = ([] if d['cluster_modules'] is None else d['cluster_modules'])
        self.local_modules = ([] if d['local_modules'] is None else d['local_modules'])
        self.modules = ([] if d['modules'] is None else d['modules'])

        #for backwards compatability
        for commands in [self.cluster_commands,self.local_commands,self.commands]:
            for i,c in enumerate(commands):
                if type(c) != Command:
                    try:
                        commands[i] = Command(c)
                    except:
                        print commands[i]
                        raise
        

        if self.step is not None and append_to_ana:
            self.step.append_job(self)

    def bind_to_step(self,step):
        self.step = step
        self.analysis = step.analysis
        id = self.id
        self.name = str(step.name) + ('_' if len(id)>0 else '') + id
        sn = self.step.short_name
        self.job_name = (sn if len(sn)<=3 else sn[:3]) + ('_' if len(id)>0 else '') + (id if len(id) <=8 else id[:8])
        self.file_name = os.path.join(self.analysis.project,self.analysis.ana_dir,"jobscript/",self.name+".sh")   
        self.oe_fn=os.path.join(self.analysis.project,self.analysis.ana_dir,"log/",self.name)
        if self.verbose is None:
            self.verbose = step.verbose
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
        f = os.path.expanduser(self.oe_fn+'.e')
        if os.path.isfile(f) and os.path.getsize(f) == 0:
            ran_noerror = True
        else:
            ran_noerror = False
        return ran_noerror
    
    def pbs_header(self):
        header = ''  
        header += "#!/usr/bin/env bash\n"
        header += "#PBS -N {}\n".format(self.job_name)
        header += "#PBS -P {}\n".format(self.analysis.project_name)
        if self.debug:
            header += "#PBS -q debug\n"
        else:
            header += "#PBS -l mem={}\n".format(self.mem)
            header += "#PBS -l ncpus={}\n".format(self.ncpus)
            header += "#PBS -l walltime={}\n".format(self.walltime)
            header += "#PBS -o {}.o\n".format(os.path.expanduser(self.oe_fn))
            header += "#PBS -e {}.e\n".format(os.path.expanduser(self.oe_fn))
        if self.exit_on_error:
            header += "#exit on error\n"
            header += "set -e\n"
        header += pretty_comment(self.analysis.description,"analysis: "+self.analysis.name,'=')
        header += 'if [ -n "$PBS_ENVIRONMENT" ]; then\n'
        header += 'PROJECT_HOME=' + self.analysis.scratch + '\n'
        header += 'else\n'
        header += 'PROJECT_HOME=' + self.analysis.project + '\n'
        header += 'fi\n'
        return header

    def command_string(self):
        cmd_str = ""
        cmd_str += pretty_comment(self.step.description,"step: "+self.step.name,'+')
        cmd_str += pretty_comment(self.description,"job: "+self.name,'-')
            
        #local and cluster modules:
        if self.cluster_modules:
            cmd_str += 'if [ -n "$BC_HPC_SYSTEM" ]; then\n'
            for mod in self.cluster_modules:
                cmd_str += "module load {}\n".format(mod)
            cmd_str += "fi\n"
        if self.local_modules:
            cmd_str += 'if [ -z "$BC_HPC_SYSTEM" ]; then\n'
            for mod in self.local_modules:
                cmd_str += "module load {}\n".format(mod)
            cmd_str += "fi\n"       
        #load modules:
        for module in self.modules:
            cmd_str += "module load {}\n".format(module)

        cmd_str += "cd $PROJECT_HOME\n"       

        if self.cluster_commands:
            cmd_str += 'if [ -n "$BC_HPC_SYSTEM" ]; then\n'
            for cmd in self.cluster_commands:
                cmd_str += cmd.cmd_str()
            cmd_str += "fi\n"
        if self.local_commands:
            cmd_str += 'if [ -z "$BC_HPC_SYSTEM" ]; then\n'
            for cmd in self.local_commands:
                cmd_str += cmd.cmd_str()
            cmd_str += "fi\n"
        #command:
        #print "the commands to write:",commands
        for command in self.commands:
            cmd_str += command.cmd_str()
        return  cmd_str

    def write_to_jobscript(self,strings):
        with open(os.path.expanduser(self.file_name), "w") as jf:
            for string in strings:
                jf.write(string)

    #depriciated:
    def write_jobscript(self):
        self.write_to_jobscript([self.pbs_header(),self.command_string()])
        self.chmod_jobscript()
        """
        self.vprint("writing jobscript for",self.name,"to",os.path.expanduser(self.file_name),mv=1)
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
            #print 'name:',name
            jf.write("#PBS -N {}\n".format(self.job_name))
            jf.write("#PBS -P vervetmonkey\n")
            if debug:
                jf.write("#PBS -q debug\n")
            else:
                jf.write("#PBS -l mem={}\n".format(mem))
                jf.write("#PBS -l ncpus={}\n".format(ncpus))
                jf.write("#PBS -l walltime={}\n".format(walltime))
            jf.write("#PBS -o {}.o\n".format(os.path.expanduser(self.oe_fn)))
            jf.write("#PBS -e {}.e\n".format(os.path.expanduser(self.oe_fn)))
            #basic description of the job
            jf.write(pretty_comment(self.analysis.description,"analysis: "+self.analysis.name,'='))
            #jf.write("\n")
            #if self.analysis.description:
            #    descr = '#' + self.analysis.description.replace("\n","\n#")
            #    jf.write(descr)
            #    jf.write("\n")

            jf.write(pretty_comment(self.step.description,"step: "+self.step.name,'+'))
            jf.write(pretty_comment(self.description,"job: "+self.name,'-'))
            
            if self.exit_on_error:
                jf.write("#exit on error\n")
                jf.write("set -e\n")
            #local and cluster modules:
            if self.cluster_modules:
                jf.write('if [ -n "$BC_HPC_SYSTEM" ]; then\n')
                for mod in self.cluster_modules:
                    jf.write("module load {}\n".format(mod))
                jf.write("fi\n")
            if self.local_modules:
                jf.write('if [ -z "$BC_HPC_SYSTEM" ]; then\n')
                for mod in self.local_modules:
                    jf.write("module load {}\n".format(mod))
                jf.write("fi\n")           
            #load modules:
            for module in modules:
                jf.write("module load {}\n".format(module))


            if self.cluster_commands:
                jf.write('if [ -n "$BC_HPC_SYSTEM" ]; then\n')
                for cmd in self.cluster_commands:
                    cmd.write(jf)
                    jf.write("\n")
                jf.write("fi\n")
            if self.local_commands:
                jf.write('if [ -z "$BC_HPC_SYSTEM" ]; then\n')
                for cmd in self.local_commands:
                    cmd.write(jf)
                    jf.write("\n")
                jf.write("fi\n")
            #command:
            #print "the commands to write:",commands
            for command in commands:
                command.write(jf)
                jf.write('\n')
        self.chmod_jobscript()
        #return jfn
    """


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
        self.vprint("submitted",self.name,"with job ID",self.pbs_id,mv=1)
        self.vprint("submit command was",command,mv=2)
        self.vprint("out =",out,mv=2)
        self.vprint("err =",err,mv=2)
        return self.pbs_id        

    def release(self):
        command = "qrls {0}".format(self.pbs_id)
        p = subprocess.Popen(command, shell=True)
        self.vprint("released job",self.name,"with",command,mv=1)

    def run_jobscript(self):
        command = os.path.expanduser(self.file_name)
        self.vprint("running locally",self.name)
        self.execute(command)
        
    def execute(self,command):
        import select
        for depend_job in self.depends:
            if depend_job.returncode != 0:
                raise Exception("{0} finished with exit code {1}, won't start {2}".format(depend_job.name,depend_job.returncode,self.name))
        print "running", self.name
        p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        if self.verbose >= 1:
            stdout = []
            stderr = []
            while True:
                reads = [p.stdout.fileno(), p.stderr.fileno()]
                ret = select.select(reads, [], [])

                for fd in ret[0]:
                    if fd == p.stdout.fileno():
                        read = p.stdout.readline()
                        #only print if the line is non-empty
                        if read.strip():
                            sys.stdout.write(read)
                        stdout.append(read)
                    if fd == p.stderr.fileno():
                        read = p.stderr.readline()
                        if read.strip():
                            sys.stderr.write(read)
                        stderr.append(read)
                if p.poll() != None:
                    break
            out = ''.join(stdout)
            err = ''.join(stderr)
        else:
            out, err = p.communicate()

        for name, content in zip(['o','e'],[out,err]):
            with open(os.path.expanduser(self.oe_fn+'.'+name),'w') as oef:
                oef.write(content)
        self.returncode = p.returncode

    def print_summary(self):
        print "-"*50
        print self.name
        print "depends on:", [job.name for job in self.depends]
        print "job filename:", self.file_name
        print "oe filenames:", self.oe_fn
        print "input files:", self.input_files
        print "output_files:", self.output_files 
        print "-"*50
        
class JoinedJob(Job):
    """
    This class holds several job objects to create a single job.
    Per default it takes the sum of walltimes, the maximim of cores
    and of memory requested.
    """
    def __init__(self,jobs, step = None, append_to_ana=True, id='', depends=None, ncpus=None, mem=None, walltime=None, exit_on_error=True, description=None, debug=False, verbose=None):
        def mem_in_mb(mem_str):
            if mem_str[-2:] == 'gb':
                mem = int(mem_str[:-2])*1024
            elif mem_str[-2:] == 'mb':
                mem = int(mem_str[:-2])
            else:
                raise Exception('Cannot parse memory string: '+mem_str)
            return mem
    
        def time_delta(time_str):
            t = datetime.datetime.strptime(time_str,"%H:%M:%S")
            return datetime.timedelta(hours=t.hour,minutes=t.minute,seconds=t.second)
 
        #implement time comparison
        def seconds(time_str):
            times = map(int,time_str.split(":"))
            if len(times) == 3:
                return times[2] + 60*(times[1] +times[0]*60)   
            else:
                raise ValueError("walltime should have format 'HH:MM:SS'")
        
        def time_str(seconds):
            hours, rest = divmod(seconds, 3600)
            minutes, seconds = divmod(rest, 60)
            return "{}:{}:{}".format(hours,minutes,seconds)

        self.jobs = jobs

        if ncpus is None:
            self.ncpus = max([job.ncpus for job in jobs])
        else:
            self.ncpus = ncpus

        if mem is None:
            self.mem = str(max([mem_in_mb(job.mem) for job in jobs])) + 'mb'
        else:
            self.mem = mem

        if walltime is None:
            #add up the individual job walltimes
            sum_walltime = reduce(lambda x,y: x+y, map(lambda job: seconds(job.walltime), jobs))
            if sum_walltime > max_walltime:
                raise Exception("Sum of walltimes exceeds maximum walltime. Consider setting it with keyword walltime.") 
            self.walltime = time_str(sum_walltime)
        else:
            self.walltime = walltime

        if not id:
            if not all(job.id == jobs[0].id for job in jobs):
                raise Exception("Jobs don't have same id. Provide id for JoinedStep in this case! IDs: "+str([job.id for job in jobs]))
            else:
                self.id = jobs[0].id
        else:
            self.id = id
        
        #output and input files are the union of those of the jobs to be joined
        self.input_files = list(set(reduce(lambda x,y: x+y,[job.input_files for job in jobs])))    
        self.output_files = list(set(reduce(lambda x,y: x+y,[job.output_files for job in jobs])))    
        
        #only add 
        self.depends = list(set((depends if depends is not None else []) +
                    reduce(lambda x,y: x+y,[[j for j in job.depends if j not in jobs] for job in jobs])))

        self.debug = debug
        self.pbs_id = None
        self.returncode = None
        self.exit_on_error = exit_on_error
        self.description = (description if description is not None else '')
        self.verbose = verbose
        

        if step is not None:
            self.bind_to_step(step)
        else:
            self.name = None
            self.step = None
            self.file_name = None
       
        if self.step is not None and append_to_ana:
            self.step.append_job(self)
    
   
    def write_jobscript(self):
        self.write_to_jobscript([self.pbs_header()]+[job.command_string() for job in self.jobs])
        self.chmod_jobscript()
    
    


class StageJob(Job):
    def __init__(self,  direction, files=None, step=None, stage_analysis_dir=True, mode='newer',depends=None, description=None,  verbose=None, debug=False):
        if direction not in ['in','out']:
            raise ValueError('direction should be "in" or "out" but is {}'.format(direction))
        self.direction = direction
        self.files = ([] if files is None else files)
        self.depends = ([] if depends is None else depends)
        self.debug = debug
        self.stage_analysis_dir = stage_analysis_dir
        self.mode = mode
        self.verbose = verbose
        self.description = (description if description is not None else '')
        if step is not None:
            self.bind_to_step(step)
        else:
            self.name = None
            self.step = None    
            self.file_name = None
            self.id = None
        self.input_files = None
        self.output_files = None    
        
    def bind_to_step(self,step):
        if step.analysis.host == 'workstation':
            hid = 'lab'
            if self.direction == 'in':
                self.source = step.analysis.project
                self.target = step.analysis.lab_dir
            elif self.direction == 'out':
                self.source = step.analysis.lab_dir
                self.target = step.analysis.project
        elif step.analysis.host == 'cluster':
            hid = 'scratch'
            if self.direction == 'in':
                self.source = step.analysis.project
                self.target = step.analysis.scratch
            elif self.direction == 'out':
                self.source = step.analysis.scratch
                self.target = step.analysis.project
            
        self.id = 'stg' + self.direction + '_'  + hid
        Job.bind_to_step(self,step)
        #self.scratch = ('~/vervet_lab' if hid == 'lws12' else '~/vervet_scratch')
        if self.stage_analysis_dir:
            self.files.insert(0,'analyses/'+step.analysis.name+'/')
        #if self.verbose is None:
        #    self.verbose = step.verbose
        #choose appropriate name and use this, attention with different file systems...
        #-> different files for local and remote ...    
    
    def stage(self,run_type='auto'):
        if self.step.run_type=='run':
            #depends = []
            wait = True
            start_on_hold = False
        elif self.step.run_type=='qsub':
            #depends = [job.pbs_id.strip() for job in self.depends]
            wait = False
            start_on_hold =True
        elif self.step.run_type is None:
            wait = False
            start_on_hold = False        

        if self.depends:
            raise Exception('Dependencies not implemented for stage job in "stage".')        

        self.vprint("staging",self.name,"in mode",run_type,mv=1)
        
        stage_command = self.stage_command(run_type,start_on_hold=start_on_hold)
        if "dmn" not in self.analysis.host:
            stage_command = "ssh dmn.mendel.gmi.oeaw.ac.at nohup " + stage_command
        
        p = subprocess.Popen(stage_command, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        
        out, err = p.communicate()
        self.returncode = p.returncode

        out = out if out.strip() else None
        err = err if err.strip() else None

        self.vprint("running stage command:",stage_command,mv=2)
        self.vprint("out is:",out,mv=2)
        self.vprint("err is:",err,mv=2)
        
        #naive check whether a job was submitted
        if 'login' in out:
            self.pbs_id = out.strip()
            self.vprint("Submitted",self.name,"with job ID",self.pbs_id,mv=1)
            if wait:
                import poll_pbs
                self.vprint("Waiting for job {} to finish...".format(self.pbs_id))
                self.returncode = poll_pbs.wait_for_job(out.strip())
        else:
            self.vprint("Staging job {} directly.".format(self.name),mv=1)


        
    def stage_command(self,run_type='auto',start_on_hold=False):
        #todo: implement a separate staging out that happens in any case, even if the jobs failed...
        run_types =  ['direct','submit','auto','dry_run']
        if run_type not in run_types:
            raise ValueError("run_type should be in {} but is {}".format(run_types,run_type))        
        if start_on_hold:
            hold = "-H "
        else:
            hold = ""
        stage_command = "dmn_stage.py -m {mode} -t {run_type} -v {verbose} -j {job_fn} -o {oe_fn}" \
                        " -n {job_name} -l {print_to} {hold}{source} {target} {files}".format(
                        mode=self.mode,run_type=run_type,verbose=self.verbose,job_fn=self.file_name,
                        oe_fn=self.oe_fn,
                        job_name=self.job_name,print_to=self.oe_fn+"_prep_dmn.o", hold=hold,
                        source=self.source,target=self.target,
                        files=" ".join(self.files))

        return stage_command
        


    #the following two functions are only used by the stage-out job
    #this is basically a dummy job that tests first how much files really need to be staged and
    #then submits or directly runs a job accordingly
    def write_prepare_jobscript(self):
        mem = '3825mb'
        ncpus = 1
        walltime = '00:20:00'
#        if mode not in modes:
#            raise ValueError('mode must be in {0} but is {1}'.format(modes, mode))
#        jfn = os.path.expanduser(self.file_name)
        fn = os.path.splitext(self.file_name)[0] + "_prep" + os.path.splitext(self.file_name)[1]

        self.vprint("Writing prepare-jobscript for",self.name,fn,mv=1)
        with open(os.path.expanduser(fn),'w') as jf: 
            jf.write("#!/bin/bash\n")
            jf.write("#PBS -N {0}\n".format(self.job_name))
            jf.write("#PBS -q staging\n")
            jf.write("#PBS -P vervetmonkey\n")
            jf.write("#PBS -o {0}_prep_job.o\n".format(os.path.expanduser(self.oe_fn)))
            jf.write("#PBS -e {0}_prep_job.e\n".format(os.path.expanduser(self.oe_fn)))
            jf.write("#PBS -l mem={0}\n".format(mem))
            jf.write("#PBS -l ncpus={0}\n".format(ncpus))
            jf.write("#PBS -l walltime={0}\n".format(walltime))
            jf.write("module load Python/2.7.3-goolf-1.4.10\n")
            stage_command = self.stage_command()
                        
            commands = ["if [[ `hostname -s` = dmn* ]]; then","echo already at dmn",
                        stage_command,
                        "else","echo ssh to dmn",
                        "ssh dmn.mendel.gmi.oeaw.ac.at nohup {}".format(stage_command),
                        "fi"]
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
        #print self.name, command
        self.pbs_id = out.strip()
        self.analysis.append_submit_log(command, out, err)
        self.vprint("submitted",os.path.expanduser(self.prep_file_name),"with job ID",self.pbs_id,mv=1)
        return self.pbs_id

        
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
    def __init__(self,command,description=None,job=None,format=True):
        self.description = description if description is not None else ""
        self.command = command
        if job is not None:
            self.bind_to_job(job)
        #One could add an exception that format is not used when Command is
        #called from within this module (analysis.py)
        #e.g. by checking whether callingframe.f_locals["__file__"] exists
        if format:
            callingframe = sys._getframe(1).f_locals
            try:
                self.command = self.command.format(**callingframe)
            except KeyError:
                print "The dict used for formating:", callingframe
                raise
    
    def bind_to_job(self,job):
        self.job = job
        if self in job.commands:
            raise Exception("Trying to add command to job.commands that is already added.")
        job.commands.append(self)

    def cmd_str(self):
        cmd_str = ''
        if self.description:
            descr = '#' + self.description.replace("\n","\n#") + '\n'
            cmd_str += descr
        cmd_str += self.command
        cmd_str += '\n' 
        return cmd_str

    def write(self,file):
        if self.description:
            descr = '#' + self.description.replace("\n","\n#")
            file.write(descr)
            file.write("\n")
        file.write(self.command)

if __name__ == '__main__':
    import argparse, json

    config_example = {"default_project": "vervet", "projects": {"test": {"scratch_dir": "test_scratch", "lab_dir": "test_lab", "project_name": "vervetmonkey", "project_dir": "test_project"}, "vervet": {"scratch_dir": "~/vervet_scratch", "lab_dir": "~/vervet_lab", "project_name": "vervetmonkey", "project_dir": "~/vervet_project"}}}    

    try:
        with open("analysis.config",'rb') as conf:
            default_dict = json.load(conf)
        print "defaults loaded updated from analysis.config."
    except IOError,e:
        raise Exception("analysis.config could not be loaded, make sure that the file exists in the same " \
                        "directory as analysis.py and that is it of the form {}. " \
                        "The error is: {}".format(config_example,e.message)) 
      
    


    parser=argparse.ArgumentParser(description="Initialise an new analysis. " \
                                            "Default arguments can be supplied " \
                                            "in a file analysis.config, "\
                                            'which is a json dict. '\
                                            "Use the option --config-example to get an example " \
                                            "of such a dict.")        
    parser.add_argument("name",help="Name of the analysis (not including date).")
    parser.add_argument("--date",default=datetime.datetime.strftime(datetime.datetime.now(),"%Y%m%d"),help="This date is added to the analysis name and folder. Should be YYYYMMDD. Default is today.")        
    parser.add_argument('-p','--project',default=None,help="Project to which this analysis belongs. If not specified the default from analysis.config is used.")
    #parser.add_argument("--project_dir",default=None,help="The base directory in which the project will be created. To have the same path on cluster and workstation, you can use a symbolic link of the form '~/project_dir' that points to the real project location.")
    #parser.add_argument("--scratch_dir",default=None,help="Equivalent of <project_dir>. Location on the lustre file system to which analysis is staged.")
    #parser.add_argument("--project_name",default=None,help="Name of the project. This must be the name used in the PBS accounting system as supplied with qsub -N <project_name>.")       
    parser.add_argument("--config-example",help="Print an example of the content of the analysis.config file and exit.",action="store_true")       
    #parser.add_argument("--library_dir",default=None,help="Path to the analysis python module. Usually this is the directory in which this script resides. Avoid tilde expansion with single quotes to keep the generic ~, e.g. --library_dir '~/script'. (This applies if you want to run the analysis on different systems, such as cluster and workstation.)")
                    
    args = parser.parse_args() 
    
    if args.config_example:
        print config_example
        sys.exit(0)


    if args.project is not None:
        d =  default_dict['projects'][args.project]
    else:
        d = default_dict['projects'][default_dict["default_project"]] 

    ana_dir = os.path.join('analyses',args.date + '_' + args.name)

    #create analysis directories:     
    make_ana_dirs(d['project_dir'],ana_dir)
    print "creating project dir:"
    print os.path.join(d['project_dir'],ana_dir)
    make_ana_dirs(d['scratch_dir'],ana_dir)
    print "creating scratch dir:"
    print os.path.join(d['scratch_dir'],ana_dir)

    #create the anaylsis file
    fn = os.path.expanduser(os.path.join(d['project_dir'],ana_dir,'script',args.name+'.py'))
    ana_str = "#!/usr/bin/env python\n"
    ana_str += "#imports\n"
    #ana_str += "#add the location of the analysis module to the python path\n"
    #ana_str += "import sys, os\n"
    #ana_str += "sys.path.insert(0,os.path.expanduser('{}'))\n".format(args.library_dir)
    ana_str += "from hs_vervet.tools import analysis as ana\n"
    ana_str += "\n"
    ana_str += "#optional: set some global varialbes here\n"
    ana_str += "\n"
    ana_str += "#create the analysis object\n"
    ana_str += "analysis = ana.Analysis(name='{}',project_dir='{}',scratch_dir='{}',project_name='{}')\n".format(args.date + '_' + args.name,d['project_dir'],d['scratch_dir'],d['project_name'])
    #make an additional arguments that allows to write some job and step examples to the file ..
    ana_str += '\n'
    ana_str += '#define your steps and jobs here\n'
    ana_str += '\n'
    ana_str += "if __name__ == '__main__':\n"
    ana_str += "    #run your analysis here\n"
    ana_str += "    pass\n"
    if os.path.exists(fn):
        raise Exception('Analysis file {} exists. Delete it before running analysis.py.'.format(fn))
    with open(fn,"w") as f:
        f.write(ana_str)
    
    print "analysis script template written to:"
    print fn
        

                    


