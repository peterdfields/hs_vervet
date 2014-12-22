"""
Functions that return specific analysis steps that are reused.
"""
from hs_vervet.tools import analysis as ana


def test_enrichment(name, out_suffix,
                    runs_per_job, n_jobs, top_mode, top_val, max_dist,
                    gene_to_cat_fn, in_rod_fn, columns, 
                    cat_to_name_fn, gene_df_fn, analysis, 
                    sort_mode="descending", mem = "3870mb"):
    out_suffix = ("_"+out_suffix if out_suffix else '')
    top_modes = ["top_n","top_q","thresh"]
    assert top_mode in top_modes, "n=top_mode must be in {} but is."\
                                            .format(top_modes,top_mode)
    log_fn = "{}/_data/{}_{}{}_d{}{}.log".format(analysis.ana_dir,name,top_mode,top_val,max_dist,out_suffix)
    #with open
    log_str = "\\\n  --log_fn " + log_fn
    ana_dir = analysis.ana_dir
    columns = [str(i) for i in columns]
    cols = " ".join(columns)
    step = ana.Step(name=name,analysis=analysis)
    cluster_modules = ["numpy/1.6.2-goolf-1.4.10-Python-2.7.3",
                    "pandas/0.14.0-goolf-1.4.10-Python-2.7.3"]
    total_permut = runs_per_job * n_jobs
    jobs = []
    out_fns = []
    ra_job = ana.Job(id="real_assoc_"+top_mode+str(top_val)+out_suffix,step=step,cluster_modules=cluster_modules,
                            walltime="00:10:00")
    peaks_per_gene_fn = "{}/_data/{}_genes_{}{}_d{}{}.tsv".format(analysis.ana_dir,name,top_mode,top_val,max_dist,out_suffix)
    top_peaks_fn = "{}/_data/{}_peaks_{}{}_d{}{}.tsv".format(analysis.ana_dir,name,top_mode,top_val,max_dist,out_suffix)
    assoc_fn = "{}/_data/{}_assoc0_{}{}_d{}{}.tsv".format(analysis.ana_dir,name,top_mode,top_val,max_dist,out_suffix)
    ana.Command("test_enrichment.py --gene_to_cat_fn {gene_to_cat_fn} \\\n"
                "  real_assoc --in_rod {in_rod_fn} \\\n"
                "  --cols {cols} \\\n"
                "  --{top_mode} {top_val} --max_dist {max_dist} \\\n"
                "  --gene_df_fn {gene_df_fn} \\\n"
                "  --{sort_mode} \\\n"
                "  --top_peaks_fn {top_peaks_fn} \\\n"
                "  --peaks_per_gene_fn {peaks_per_gene_fn}"
                "  --assoc_fn {assoc_fn} {log_str}",job=ra_job)
    for n in range(n_jobs):


        job = ana.Job(id=top_mode+str(top_val)+"_"+str(n)+out_suffix,
                                step=step,cluster_modules=cluster_modules,
                                ncpus=1,walltime="08:00:00",afterok=[ra_job],mem=mem)
        out_fn = "{}/_data/temp.{}_go_enrichment_{}{}_d{}_{}{}.tsv"\
                                            .format(analysis.ana_dir,name,
                                                    top_mode,top_val,max_dist,
                                                                    n,out_suffix)
        out_fns.append(out_fn)
        ana.Command("test_enrichment.py --gene_to_cat_fn {gene_to_cat_fn} \\\n"
                    "  permutation --in_rod {in_rod_fn}  \\\n"
                    "  --cols {cols} \\\n"
                    "  --{top_mode} {top_val} --max_dist {max_dist} \\\n"
                    "  --gene_df_fn {gene_df_fn} \\\n"
                    "  --{sort_mode} \\\n"
                    "  --real_assoc_fn {assoc_fn} --n_runs {runs_per_job} \\\n"
                    " --out_fn {out_fn} {log_str}",job=job)
        jobs.append(job)
    permut_fns = " \\\n  ".join(out_fns)
    reduce_job = ana.Job(id="reduce_"+top_mode+str(top_val)+out_suffix,
                         step=step,afterany=jobs,cluster_modules=cluster_modules)
    ana.Command("test_enrichment.py --gene_to_cat_fn {gene_to_cat_fn} \\\n"
                "  reduce  \\\n"
                "  --pval_out {ana_dir}/_data/"
                "{name}_go_enrichment_{top_mode}{top_val}_d{max_dist}"
                "_permut{total_permut}{out_suffix}.tsv \\\n"
                "  --pval_sign_out {ana_dir}/_data/"
                "{name}_go_enrichment_{top_mode}{top_val}_d{max_dist}"
                "_permut{total_permut}{out_suffix}_sign0.05.tsv \\\n"
                "  --pval_sign_thresh 0.05 \\\n"
                "  --peaks_per_gene_fn {peaks_per_gene_fn} \\\n"
                "  --cat_to_name_fn {cat_to_name_fn} \\\n"
                "  --permut_fns {permut_fns} {log_str}",job=reduce_job)
    return step

