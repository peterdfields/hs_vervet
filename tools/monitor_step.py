#!/usr/bin/env python
import subprocess, os, time
import pandas as pd

def monitor_step(pbs_ids,ana_job_ids,summary_fn,poll_interval=300):
    """
    Runs a deamon (execute on login node),
    that polls the jobs with qstat, runs
    qstat -xf on finsished jobs and parses the output.
    The individual qstat -xf outputs are stored in the 
    same location/filename as the jobs output file but with
    extension .qstat"
    A summary file for all jobs is continuously updated.
    """
    # stats to add to the summary file (use keywords from qstat_dict)
    summary_stats = ["Job Id","Exit_status","qtime","Resource_List.walltime","resources_used.walltime",
                    "Resource_List.mem","resources_used.mem","Resource_List.ncpus",
                        "resources_used.cpupercent"]
        
    try:
        stat_df = pd.read_csv(summary_fn,sep="\t")
        stat_df.set_index("ana_job_id",drop=False,inplace=True)
    except IOError:
        stat_df = pd.DataFrame(columns=["ana_job_id"]+summary_stats,index=ana_job_ids)
    
    #with open(summary_fn,"w") as sf:
    #    sf.write("\t".join(["ana_job_id"]+summary_stats)+"\n")
    finished = False
    time.sleep(10)
    while pbs_ids:
        p = subprocess.Popen(["qstat"] + pbs_ids,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        o, e = p.communicate()
        #handle ended jobs
        if e.strip():
            for line in e.strip().split("\n"):
                id = line.split()[1]
                p2 = subprocess.Popen(["qstat","-xf",id],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                qstat, e2 = p2.communicate()
                assert not e2, "qstat -xf " + id + "threw error: " + e2
                qstat_dic = parse_qstat_f(qstat)
                #save qstat to file
                base_fn = os.path.splitext(qstat_dic["Output_Path"])[0].split(":")[1]
                fn =  base_fn + ".qstat"
                with open(fn,"w") as f:
                    f.write(qstat)
                #if ana_job_ids is not None:
                ana_job_id = ana_job_ids[pbs_ids.index(id)]
                #else:
                #    ana_job_id = base_fn.split("_")[-1] 
                #append summary
                
                #report_dic = {k:qstat_dic[k] for k in summary_stats}
                report_dic = {}
                for k in summary_stats:
                    try:
                        report_dic.update({k:qstat_dic[k]})
                    except KeyError:
                        print k, "not in qstat_dic:"
                        print qstat_dic
                report_dic.update({"ana_job_id":ana_job_id})
                stat_series = pd.Series(report_dic)
                try:
                    stat_df.ix[ana_job_id] = stat_series
                except KeyError:
                    stat_series.name = ana_job_id
                    stat_df = stat_df.append(stat_series)
                pbs_ids.remove(id)
                ana_job_ids.remove(ana_job_id)
            stat_df.to_csv(summary_fn,index=False,sep="\t",float_format="%i")
        if pbs_ids:
            time.sleep(poll_interval)

def parse_qstat_f(qstat_f):
    """
    parse the output of qstat -f
    into a python dictionary
    """
    #remove line continuations
    qstat_f = qstat_f.strip().replace("\n\t","")
    #split entries, giving special treatment to first line
    stat_ls = qstat_f.split("\n")
    id_entry = stat_ls[0].split(": ")
    stat_ls = [s.strip().split(" = ") for s in stat_ls[1:]]
    stat_dict = {s[0]:s[1] for s in stat_ls}
    stat_dict.update({id_entry[0]:id_entry[1]})
    #for k,v in stat_dict.iteritems():
    #    print k,":",v
    return stat_dict


if __name__ == '__main__':
    import argparse
    parser=argparse.ArgumentParser(description="Monitor PBS jobs. Once a job is finished the\n"
                                               "output of qstat -xf is written to a file\n"
                                                "and a summary file for all jobs is updated.",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("pbs_ids",nargs="+",help="PBS ids of the jobs to monitor.")
    parser.add_argument("-f","--stats_fname",required=True,help="Filename to store the stats in.")
    parser.add_argument("-a","--ana_job_ids",required=True,nargs="+",help="Job identifiers as they are used in the analysis script.")
    parser.add_argument("-i","--polling_interval",type=int,default=300,help="Time between two polls in seconds.")
    args=parser.parse_args()

    if args.ana_job_ids is not None:
        assert len(args.ana_job_ids) == len(args.pbs_ids)

    monitor_step(args.pbs_ids,args.ana_job_ids,args.stats_fname,poll_interval=args.polling_interval)
    #monitor_step(["235293.login0","235294.login0","235295.login0","235296.login0"],os.getcwd()+"/test_pollstat.txt",poll_intervall=10)

