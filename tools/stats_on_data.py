#!/usr/bin/env python
"""
This module contains functions to calculate statistics over the genome
and supporting functions to do this.
"""

import sys, os
import numpy as np
import pandas as pd
import scipy
#matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from scipy.stats import nanmean

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import hs_vervet_basics as hvb

#-----------supporting functions------------

def rolling_window_ovl_old(a,window,overlap):
    """
    fast way to separate array "a" into list of rolling windows with overlap
    """
    #array with rolling windows along axis one
    shape = (int(round(a.shape[-1]/(window-overlap))), window)
    strides = (a.strides[0]*(window-overlap),a.strides[0])
    print shape,strides,a.shape
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

def rolling_window_ovl(a,window,overlap):
    """
    fast way to separate array "a" into list of rolling windows with overlap
    """
    #array with rolling windows along axis one
    shape = (int(round(a.shape[0]/(window-overlap))), window)+a.shape[1:]
    strides = (a.strides[0]*(window-overlap),a.strides[0]) +a.strides[1:]
    #print shape,strides,a.shape
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

def rolling_window(a, window,overlap):
    print a.shape[:-1], (a.shape[-1] - window + 1, window)
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    print a.strides,(a.strides[-1],)
    strides = a.strides + (a.strides[-1],)
    print shape,strides,a.shape
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)

def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """
    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
#    in our case of pairwise comparisons: len(arrays)==2
        out = np.zeros([n, len(arrays)], dtype=dtype)
    print out,arrays
    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

def combinations2(array1,array2):
    """
	array of all combinations of elements from the two arrays
    """
    array1 = np.asarray(array1)
    array2 = np.asarray(array2)
    dtype = array1.dtype

    n = array1.size*array2.size
    out = np.zeros([n, 2], dtype=dtype)
    out[:,0]=np.repeat(array1,array2.size)
    out[:,1]=np.tile(array2,array1.size)
    return out.tolist()

def zipC(array1,array2):
    """
    array of all tuples of elements from the two arrays
    """
    array1 = np.asarray(array1)
    array2 = np.asarray(array2)

    t1=np.repeat(array1,array2.size)
    t2=np.tile(array2,array1.size)
    return zip(t1,t2)


   
def combination_map(f,array1,array2):
    """
	applies function <f> to each element in 
	an array of all combinations of elements from the two arrays
    """
    array1 = np.asarray(array1)
    array2 = np.asarray(array2)
    dtype = array1.dtype

    n = array1.size*array2.shape[0]
    out = np.zeros([n, 2], dtype=dtype)
    out[:,0]=np.repeat(array1,array2.shape[0])
    out[:,1]=np.tile(array2,array1.shape[0])
    res=np.apply_along_axis(lambda x:f(*x),1,out)
    return res


#-----------------functions to chop a contig or chromosome in smaller pieces-----

def chopContBP(contig_data_dict,window_size,overlap=0):
    """
    Chops the data in a SNP data dictionary of a contig or chromosome into equal chunks
    and returns a list of data dictionaries, one for each window.

    input:
    contig_data_dict ...  a dictionary containing the SNP data (as produced by hsVCFEXtrationTools)
    (hint: only pass data in the dictionary that you actually need)
    window_size ... size of the window in bp


    output:
    list of dictionaries, one for each window

    """
    sp=contig_data_dict['snp_pos_a1']
    ws=window_size

    cut_left=range(0,sp[-1]+1,ws-overlap)
    cut_right=range(ws,sp[-1]+1,ws-overlap)
    split_left=sp.searchsorted(cut_left)
    split_right=sp.searchsorted(cut_right)
    split_dict={key:value for (key,value) in contig_data_dict.iteritems()}


    for key,value in split_dict.iteritems():
        if isinstance(value, (list,np.ndarray)) and len(value)==len(sp):
            split_dict[key]=[value[left:right] for left,right in zip(split_left,split_right)]

    split_ls=[]
    for i in xrange(len(split_dict['snp_pos_a1'])):
        dic={}
        for key,value in split_dict.iteritems():
            #print key, len(value),len(split_at+1),isinstance(value, (list,np.ndarray)), type(value)
            if isinstance(value, (list,np.ndarray)) and len(value)==len(split_dict['snp_pos_a1']):
                dic[key]=value[i]
            else:
                dic[key]=value
        split_ls.append(dic)

    return split_ls

def getSplitTuples(snp_pos,window,overlap):
    """
    Returns tuples of indices where a list of snp-data has to be split
    to get the snp-information for each window in base-pair units.
    """
    cut_left=range(0,snp_pos[-1]+1,window-overlap)
    cut_right=range(window,snp_pos[-1]+1+window-overlap,window-overlap)
    split_left=snp_pos.searchsorted(cut_left)
    split_right=snp_pos.searchsorted(cut_right)
    return zip(split_left,split_right)

def chopContSNP(contig_data_dict,window_size,overlap=0):
    """
    Takes the data dictionary and chops it into smaller dictionaries, 
    each with <window_size> SNPs.
    """
    sp=contig_data_dict['snp_pos_a1']
    #copy the contig data
    split_dict={key:value for (key,value) in contig_data_dict.iteritems()}
    for key,value in split_dict.iteritems():
        if isinstance(value,(list,np.ndarray)) and len(value)==len(sp):
            split_dict[key]=rolling_window_ovl(value,window_size,overlap)
    split_ls=[]
    for i in xrange(len(split_dict['snp_pos_a1'])):
        dic={}
        for key,value in split_dict.iteritems():

            if isinstance(value, (list,np.ndarray)) and len(value)==len(split_dict['snp_pos_a1']):
                dic[key]=value[i]
            else:
                dic[key]=value
        split_ls.append(dic)

    return split_ls


def momentOverWindow(arr,snp_pos,window,overlap,statistic=np.mean,mode='BP'):
    """
    Calculate a given statistic (identity, mean, variance,...) over windows over the data.
    """
    if mode=='BP':
        lst=[arr[left:right] for left,right in getSplitTuples(snp_pos,window,overlap)]
    elif mode=='SNP':
        lst=rolling_window_ovl(arr,window,overlap)
    else:
        print 'please specify mode "BP" (basepair-windows) or "SNP" (SNP-windows)'
    return map(statistic,lst)



def rowApply(df,stat_dict={'Mean':nanmean}):
    """
    Applies functions in stat_dict to each row. 
    Function is expected to return a scalar. 
    """
    return pd.DataFrame({name:df.apply(stat,axis=1) for name,stat in stat_dict.items()})


def windowSummaryNoOverlap(arr,snp_pos,window,stat_dict={'Mean':nanmean}):
    snp_pos=np.array(snp_pos)
    cut_val=np.arange(window,snp_pos[-1]+1,window)
    split_at=snp_pos.searchsorted(cut_val)
    split_data=np.split(arr,split_at)
    mid_pos=np.arange(window/2,snp_pos[-1]+1+window/2,window)
    moving_stat={statn:map(statv,split_data) for statn,statv in stat_dict.items()}
    moving_stat.update({'n_SNPs':map(len,split_data)})
    return pd.DataFrame(moving_stat,index=mid_pos)

def nSNPwindow(df,window):
    return windowSummaryNoOverlap(df.index,df.index,window,{})

def afsOverWindow(df,window):
    """
    takes SNP data frame
    returns df of histogram arrays at certain positions
    """
    sample_size=df.shape[1]
    bins=np.arange(0.5,2*sample_size-0.499,1)
    af_ls=df.mean(axis=1)
    histfun=lambda x : np.histogram(x*sample_size,bins=bins)[0]
    return windowSummaryNoOverlap(af_ls,df.index,window,{'afs':histfun})


def watersonFromAFS(afsDF):
    window_size=afsDF.index[1]-afsDF.index[0]
    n=afsDF['afs'].iloc[0].shape[0]+1
    wt=afsDF['afs'].apply(sum)/(window_size*np.sum(1./np.arange(1,n-1,1)))
    wt.name="thetaW"
    return wt

def piFromAFS(afsDF):
    """
    assuming that there is a single window size
    """
    window_size=afsDF.index[1]-afsDF.index[0]
    n=afsDF['afs'].iloc[0].shape[0]+1
    fact=scipy.misc.comb(n,2)
    def f(row):
        return (np.arange(1,n,1)*(n-np.arange(1,n,1))*row).sum()/fact
    pi=afsDF['afs'].apply(f)/window_size
    pi.name="thetaPi" 
    return pi

def colWindowSummaryDF(df,window,stat_dict={'Mean':nanmean}):
    """
    df.. DataFrame or Series
    each statistic is applied to each column of the DataFrame separately
    Non overlapping windows.
    Window should be specified in terms of the integer index.
    """
    df1=pd.DataFrame(index=[],columns=[])
    for col,series in df.iteritems():
        df1=df1.combine_first(windowSummaryNoOverlap(\
            series.values,series.index,window,\
            {(str(col)+'_'+key):val for key,val in stat_dict.items()}))
    return df1



#------ pandas 

def windowApplyDF(ser,window,overlap=0,f=nanmean):
    """
    non-overlapping window per base pair
    """
    lst=[ser[left:right] for left,right in getSplitTuples(ser.index,window,overlap)]
    


#-------------functions that calculate statistics --------------

def calcAFS():
    """
    How to handle missing data? (work with frequency!?)
    allele frequency spectrum
    """
    pass
    

def r2(s1,s2):
    """
    Calculate the squared correlation coefficient between two
    rows of the genotype matrix. (i.e. 2 snps)
    """
    D=(s1*s2).mean()-s1.mean()*s2.mean()
    return D**2/(s1.mean()*s2.mean()*(1-s1.mean())*(1-s2.mean()))  

def snpwiseR2(data,snp_pos,distance):
    """
    Calculate r2 between each snp and the next snp that is >=<distance> bp to the right.
    
    Output:
    dictionary with
    'r2'... list of r2 values
    'mid_pos'... list of average position of each pair of SNPs
    'real_dist'... list of real distance between each two SNPs 
    """
    snp_pos=np.array(snp_pos)
    snp_left=np.arange(snp_pos.shape[0])
    right_pos=snp_pos+distance
    right_pos=right_pos[:right_pos.searchsorted(snp_pos[-1])]
    snp_right=snp_pos.searchsorted(right_pos)
    snp_left=snp_left[:snp_right.shape[0]]
    real_dist=snp_pos[snp_right]-snp_pos[snp_left]
    mean_pos=np.around(snp_pos[snp_left]+real_dist/2)
    #print data[0,:]
    #print vecr2(data[0,:],data[5,:])
    r2_ls=[r2(data[l,:],data[r,:]) for l,r in zip(snp_left,snp_right)]
    return pd.DataFrame({'r2':r2_ls,\
                  'real_dist':real_dist},index=mean_pos)
    

def calcR2(data1,data2):
    """
    Calculate the squared correlation coefficient over the chromosome
    between the SNPs in two windows separated by a distance. Uses all pairwise
    combinations of SNPs from the two windows. Quickly runs out of memory.
    """
    def r2(s1,s2):
        D=(s1*s2).mean()-s1.mean()*s2.mean()
        return D**2/(s1.mean()*s2.mean()*(1-s1.mean())*(1-s2.mean()))  
    r2_ls=[r2(data1[k[0],:],data2[k[1],:]) for k in combinations2(range(data1.shape[0]),range(data2.shape[0]))]
    return sum(r2_ls)/len(r2_ls) if len(r2_ls)>0 else np.nan

def r2overChunk(data,window,n_distance):
    """
    ---don't use, memory problems, rather use snpwiseR2()---
    Calculate r2 between 2 windows that are a distance <n_distance> apart.
    Calculates r2 for all combinations of SNP from the two windows.
    Runs out of memory quickly.
	"""
    data_chunks=momentOverWindow(data,window,0,statistic=lambda x:x)
    return [calcR2(a,b) for a,b in zip(data_chunks[0:-n_distance],data_chunks[n_distance:])]

def calcAF(data):
    """
    Calculate array of allele frequencies.

    data ... raw 0,1 genotype matrix

    out ... array of allele frequencies for each snp

    """
    return nanmean(data,axis=1)#data.sum(axis=1)/(1.*data.shape[1]) 

def calcPi(data):
    """
    Calculate nucleotide diversity (pi) along the data.
    This is quite slow, but does not run out of memory.
    """
    diff=np.zeros(data.shape[0])
    cp=0
    for j in range(0,data.shape[1]):
        for k in range(j+1,data.shape[1]):
            cp+=1
            diff+=np.absolute(data[:,j]-data[:,k])
    pi_ls=diff/cp
    return pi_ls

def calcWaterson(data):
    pass

def calcAvgPi_np(datamat1):
    comparisons1=datamat1.shape[1]*(datamat1.shape[1]-1)/2
    pi_arr=np.empty((datamat1.shape[0],comparisons1))
    pi_arr[:]=np.NAN
    indx=0;
    for j in range(0,datamat1.shape[1]):
        for k in range(j+1,datamat1.shape[1]):
            pi_arr[:,indx]=np.absolute(datamat1[:,j]-datamat1[:,k])
            indx=indx+1;
    return np.average(np.average(pi_arr,axis=1))


def calcFST(datamat1,datamat2):
    """
    input are numpy 0 1 matricesq
    both matrices must be for the same contig (datamat1.shape[0]==datamat2.shape[0])
    
    why does this function produce nans, even when there are only segregating snps?
    """
    #calculate pairwise difference for each SNP
    #for each pair of individuals in datamat 1
    comparisons1=datamat1.shape[1]*(datamat1.shape[1]-1)/2
    pi_arr=np.empty((datamat1.shape[0],comparisons1))
    pi_arr[:]=np.NAN
    indx=0;
    for j in range(0,datamat1.shape[1]):
        for k in range(j+1,datamat1.shape[1]):
            pi_arr[:,indx]=np.absolute(datamat1[:,j]-datamat1[:,k])
            indx=indx+1;          
    pi_within1_ls=np.average(pi_arr,axis=1)
    #same for datamat2
    comparisons2=datamat2.shape[1]*(datamat2.shape[1]-1)/2
    pi_arr=np.empty((datamat2.shape[0],comparisons2))
    pi_arr[:]=np.NAN
    indx=0;
    for j in range(0,datamat2.shape[1]):
        for k in range(j+1,datamat2.shape[1]):
            pi_arr[:,indx]=np.absolute(datamat2[:,j]-datamat2[:,k])
            indx=indx+1;          
    pi_within2_ls=np.average(pi_arr,axis=1)
    
    #calculate the weighted average
    pi_within_ls=(pi_within1_ls*comparisons1+pi_within2_ls*comparisons2)/(comparisons1+comparisons2)
    
    #pairwise comparison between populations
    pi_arr=np.empty((datamat1.shape[0],datamat1.shape[1]*datamat2.shape[1]))
    pi_arr[:]=np.NAN
    indx=0;
    for j in range(0,datamat1.shape[1]):
        for k in range(0,datamat2.shape[1]):
            pi_arr[:,indx]=np.absolute(datamat1[:,j]-datamat2[:,k])
            indx=indx+1;          
    pi_between_ls=np.average(pi_arr,axis=1)
    
    fst_ls=(pi_between_ls-pi_within_ls)/pi_between_ls
    
    return fst_ls

def plotAlongChunk(stat,snp_pos,bp_per_line=2000000,fname='tmp.png',plot_dic={}):
    """
    Plot the vector stat at positions snp_pos along the chromosome. 
    One sub-plot for each <bp_per_line> base-pairs is produced in rows in one plot.
    """
    default_plot_dic={'ylabel':'stat ?','title':'contig ?'}
    
    for key in default_plot_dic.keys():
        if key not in plot_dic.keys():
            plot_dic[key]=default_plot_dic[key]
            
    max_dat=np.nanmax(stat)
    min_dat=np.nanmin(stat)
    print [min_dat-(max_dat-min_dat)*0.01,max_dat+(max_dat-min_dat)*0.01]      
    chunk_len=snp_pos[-1]
    n_subpl=int(chunk_len/bp_per_line+1)
    fig_l, fig_h=32., 18.
    aspect=fig_h/(fig_l*n_subpl)
    fig = plt.figure(num=None, figsize=(fig_l, fig_h))
    stat_ls=momentOverWindow(stat,snp_pos,window=bp_per_line,overlap=0,statistic=lambda x: x)
    snp_pos_ls=momentOverWindow(snp_pos,snp_pos,window=bp_per_line,overlap=0,statistic=lambda x: x)  
    for i in range(1,n_subpl+1):
        dat=stat_ls[i-1]
        snp_pos1=snp_pos_ls[i-1]
        bp_start=(i-1)*bp_per_line
        bp_end=i*bp_per_line
        ax = fig.add_subplot(n_subpl,1,i,aspect=aspect*chunk_len/n_subpl)
        ax.plot(snp_pos1,dat,'o', markersize=4)
        plt.xlim([bp_start,bp_end])
        plt.ylim([min_dat-(max_dat-min_dat)*0.01,max_dat+(max_dat-min_dat)*0.01])
        if i==1:
            plt.title(plot_dic['title'],fontdict={'fontsize':40})
        ax.set_xticks([round(bp_start+bp_per_line*j/5) for j in range(6)])
        ax.set_xticklabels(["{0:.1f}".format((bp_start+j/5.*bp_per_line)*0.000001) for j in range(6)],fontdict={'fontsize':20})
        ax.tick_params(axis='both', which='major', labelsize=20)
        plt.ylabel(plot_dic['ylabel'],fontdict={'fontsize':30})
    plt.xlabel('position (mb)',fontdict={'fontsize':30})   
    #plt.subplots_adjust(wspace=0,hspace=0)
    plt.tight_layout()
    plt.savefig(fname) 

