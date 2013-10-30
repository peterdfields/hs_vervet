#!/usr/bin/env python
import os
import analysis as ana

#todo: decide wether binding a job to a step should append that step to the job list (i think no) and make the behaviour (also for init) consitent across Job and StageJob

analysis = ana.Analysis(name='20131024_test_analysis')

step1 = ana.AnalysisStep(name='test_step1',analysis=analysis)

test_files = ['test/test_dir/file_in_test_dir.txt','test/test2.test']

commands=['echo "added something to file_in_test_dir.txt" &> {0}'.format(os.path.join(analysis.scratch_dir,test_files[0]))]

job1 = ana.Job(commands,analysis_step=step1)

stagein = ana.StageJob('in',test_files,analysis_step=step1)
stageout = ana.StageJob('out',test_files,analysis_step=step1)

step1.jobs.insert(0,stagein)
step1.jobs.append(stageout)



#print analysis.ana_dir['project']


for job in step1.jobs:
    job.write_jobscript()
    job.run_jobscript()

