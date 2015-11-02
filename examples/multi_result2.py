"""
Feb 23 2014
Author: Dazhe Meng

Script to plot stuff
- New data structure: Incorporating numpy functions to make everything easier
- Being made to be compatible with other genomes

Note:
new definition for chr_offset[i]: the actual offset if chr is i, whereas it was i-1 before
"""
from tools import progress_display, sys_write
import math
import numpy as np

class Multi_result(object):
    " A result class that can handle multiple result files and plot them "
    def __init__(self, transform = "nlog"):
        self.num_dif_results = 0
        self.res = [] # creating null array
        self.transform_type = transform
        self.G = None # Token for storing gene information, if needed
        self.scoremax = None # Token for storing the maximum score, if needed
        self.allscore = None
    
    def load_result(self, result_file_name, delim = ",", guess_chr_size=False, filtermac=0):
        " new method for new data structure, assuming bjarni output "
        f = open(result_file_name)
        chr = []
        pos = []
        score = []
        plotpos = []
        f.readline()
        for l in f:
            o = l.strip().split(',')
            if filtermac:
                if int(o[4])<filtermac:
                    continue
            chr.append(int(o[0]))
            pos.append(int(o[1]))
            score.append(float(o[2]))
            plotpos.append(pos[-1]+self.chr_offsets[chr[-1]])
        new_res = np.array(zip(chr,pos,score,plotpos),dtype=[('chr',int),('pos',int),('score',float),('plotpos',int)])
        self.res.append(new_res)
        '''
        if guess_chr_size:
            schr = set(chr)
            maxchr = max(schr)
            minchr = min(schr)
            assert(len(schr)==maxchr-minchr+1,"There is chromosome without data")
            chr_sizes = []
            for i in xrange(minchr,maxchr+1):
                chr_sizes.append(pos[chr.index(i)-1])
            chr_offsets = [0]*(len(schr)+1)
            for i in xrange(1, 6):
                chr_offsets[i] = chr_offsets[i-1] + tair9_sizes[i-1]
            # not being finished due to availability of data
        '''
        
        
    def set_chr_sizes(self, species='lj'):
        " attempt to guess chromosome sizes from data "
        if species=='lj':
            self.chr_sizes = [192322135,62285374,43247325,45610869,42341900,34192293,26985540]
            self.chr_range = range(0,7)
            self.chr_offsets = {0:0}
            for i in self.chr_range[1:]:
                self.chr_offsets[i] = self.chr_offsets[i-1] + self.chr_sizes[i-1]
            self.chr_offsets[-1]=(sum(self.chr_sizes))
            self.chr_offsets[self.chr_range[-1]+1]=(sum(self.chr_sizes))
        elif species=='at':
            self.chr_sizes = [0,30427671,19698289,23459830,18585056,26975502]
            self.chr_range = range(1,6)
            self.chr_offsets = {1:0}
            for i in xrange(2, 6):
                self.chr_offsets[i] = self.chr_offsets[i-1] + self.chr_sizes[i-1]
            self.chr_offsets[-1]=sum(self.chr_sizes)
            self.chr_offsets[self.chr_range[-1]+1]=(sum(self.chr_sizes))
            
    
    def add_result(self, res_entries, res_pos):
        " function to be called from other modules "
        assert res_entries.shape[0] == res_pos.shape[0]
        self.num_dif_results += 1
        new_plotpos = res_pos[:,1].copy()
        for chr_i in self.chr_range: # only start adding in chromosome 2 onwards
            new_plotpos[res_pos[:,0]==(chr_i)]+=self.chr_offsets[chr_i]
        #new_res = np.hstack([res_pos,res_entries[:,np.newaxis],new_plotpos[:,np.newaxis]]) # stack pos and score together
        new_res = np.array(zip(res_pos[:,0],res_pos[:,1],res_entries,new_plotpos),dtype=[('chr',int),('pos',int),('score',float),('plotpos',int)])
        self.res.append(new_res)
    
    def transform(self):
        sys_write("Tranforming score values by ")
        if self.transform_type == "nlog":
            sys_write("Negative log... ")
            for r in self.res:
                r['score'] = -np.log10(r['score'])
        sys_write("Done.\n")
    
    def sort(self, whichcolumn = 2):
        " sort all the result entries once all results are loaded, not working "
        sys_write("Sorting through result entries... ")
        self.res = self.res[self.res[:,2].argsort()]
        sys_write("Done.\n")

    def find_threshold(self, percentile, universal = True):
        " find the threshold so only top #percentile are kept; Universal indicates whether the threshold is found for all result sets or individually for each "
        sys_write("Preparing data for plotting... ")
        self.allscore = np.hstack([r['score'] for r in self.res])
        self.allscore.sort()
        self.threshold = np.percentile(self.allscore, 100-percentile)
        self.max_score = self.allscore[-1]
        self.min_score = self.allscore[0]
        self.plot_padding = 0.05*(self.max_score - self.min_score)
        sys_write("Done.\n")

    def manhattan_plot(self, savefigname = "trial.png", percentile = 2, full_plot = True, plot_snp_list = None, tick_gap=10):
        """ gwa plot, manhattan style 
        percentile : the top #% to be included in the final plot
        full_plot : toggles whether the entire span of 5 chromosomes are plotted regardless of what is inside the data
        """
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        self.find_threshold(percentile)
        
        # initializing plot
        plt.figure(figsize=(12,2.8))
        plt.axes([0.045,0.15,0.95,0.71])
        
        # plotting data points
        sys_write("Plotting data points... ")
        for r in self.res:
            indices_to_plot = (r['score']>=self.threshold)
            plt.plot(r['plotpos'][indices_to_plot],r['score'][indices_to_plot],".",markersize=3,alpha=0.7)
        sys_write("Done.\n")
        
        # plotting axis
        sys_write("Plotting axes and extras... ")
        for chr_ind in self.chr_range[1:]:
            plt.plot([self.chr_offsets[chr_ind], self.chr_offsets[chr_ind]],[self.min_score - self.plot_padding, self.max_score + self.plot_padding],"k-", linewidth = 0.5)
        # bonferroni
        bonf = -np.log10(0.05/sum([len(r) for r in self.res]))
        if bonf<self.max_score+self.plot_padding: # only plot if relevant
            plt.plot([0,self.chr_offsets[-1]],[bonf,bonf],'r--')
        # ticks
        if full_plot == True:
            # This part is independent of data
            # Mostly copied from Bjarni's
            ticklist = []
            ticklabels = []
            for chr_ind in self.chr_range:
                for chr_tickpos in xrange(self.chr_offsets[chr_ind], self.chr_offsets[chr_ind+1], 5000000):
                    ticklist.append(chr_tickpos)
                for chr_ticklabel in xrange(0, self.chr_sizes[chr_ind], 5000000): # no tick for 0
                    if chr_ticklabel % (tick_gap*1000000) == 0 and chr_ticklabel < self.chr_offsets[chr_ind+1]-(tick_gap*300000) and chr_ticklabel>0:
                        ticklabels.append(chr_ticklabel/1000000)
                    else:
                        ticklabels.append("")
        plt.axis([0,self.chr_offsets[-1],self.min_score - self.plot_padding, self.max_score + self.plot_padding])
        plt.xticks(ticklist, ticklabels)
        plt.ylabel('$-log(p-$value$)$',size="large")
        # Plots custom snp list
        if plot_snp_list != None:
            for custom_snp in plot_snp_list:
                plt.plot([custom_snp.plotpos, custom_snp.plotpos], [self.min_score - self.plot_padding, self.max_score + self.plot_padding], "b-", linewidth = 0.5)
        sys_write("Done\n")
        
        # saving figure
        sys_write("Outputting figure to file... ")
        plt.savefig(savefigname,format="png",dpi=300,bbox_inches='tight')
        sys_write("Done.\n")

    def region_plot(self, chr, posbeg, posend, savefigname = "trial.png", percentile = 2, plot_snp_list = None, gene_file_name = "/Volumes/samaras/Data/TAIR9/TAIR9_GFF3_genes.gff"):
        """
        OBSELETE 
        gwa plot, regional
        percentile : the top #% to be included in the final plot
        full_plot : toggles whether the entire span of 5 chromosomes are plotted regardless of what is inside the data
        """
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        self.find_threshold(percentile)
        
        # initializing plot
        plt.figure(figsize=(12,7.5))
        plt.axes([0,0.25,1,0.75])
        
        # plotting data points
        sys_write("Plotting data points... ")
        for data_indicator in xrange(1, self.num_dif_results+1):
            scores = []
            positions = []
            for res in self.res:
                if res.chr != chr:
                    continue
                if res.pos > posend or res.pos < posbeg:
                    continue
                if res.indicator == data_indicator: # This is so not optimal! 
                    if res.score >= self.thresholds[data_indicator-1]:
                        scores.append(res.score)
                        positions.append(res.plotpos)
            plt.plot(positions,scores,".",markersize=3,alpha=0.7)
        sys_write("Done.\n")
        
        # plotting axis
        sys_write("Plotting axes and extras... ")
        # Creating a ticklist for the region; there should be separate scale lists
        ticklist = []
        ticklabels = []
        # dynamic labeling of intervals
        pos_range = posend - posbeg
        log10scale = int(math.log10(pos_range))
        log10residue = float(pos_range) / math.pow(10,log10scale)
        if log10residue < 1.2:
            log10multiplier = 0.2
        elif log10residue >=1.2 and log10residue < 2.6:
            log10multiplier = 0.5
        elif log10residue >= 2.6 and log10residue < 5.8:
            log10multiplier = 1
        else:
            log10multiplier = 2
        scale_interval = int(log10multiplier * math.pow(10,log10scale))
        if scale_interval>=1000000:
            label_scale_interval = scale_interval / 1000000
            label_scale_alphabet = "Mb"
        elif scale_interval>=1000:
            label_scale_interval = scale_interval / 1000
            label_scale_alphabet = "kb"
        else:
            label_scale_interval = scale_interval
            label_scale_alphabet = "bp"
        for i in xrange(int((posbeg-1)/scale_interval)+1, int((posend)/scale_interval)+1):
            ticklabels.append(i * label_scale_interval)
            ticklist.append(self.chr_offsets[chr-1]+i*scale_interval)
        plt.axis([self.chr_offsets[chr]+posbeg, self.chr_offsets[chr]+posend,self.min_score - self.plot_padding, self.max_score + self.plot_padding])
        plt.xticks(ticklist, ticklabels)
        plt.ylabel('$-log(p-$value$)$',size="large")
        plt.xlabel(label_scale_alphabet, size="large")
        # Plots custom snp list
        if plot_snp_list != None:
            for custom_snp in plot_snp_list:
                if custom_snp.chr == chr:
                    if custom_snp.pos >= posbeg and custom_snp.pos <= posend:
                        plt.plot([custom_snp.plotpos, custom_snp.plotpos], [self.min_score - self.plot_padding, self.max_score + self.plot_padding], "b-", linewidth = 0.5)
        sys_write("Done\n")

        # plotting genomes
        if not self.G:
            import gene_info
            sys_write("Loading gene model information... ")
            self.G = gene_info.genes(gene_file_name)        
            sys_write("Done.\n")

        l_genes = self.G.get_genes_in_range(chr, posbeg, posend)
        if len(l_genes)>50:
            print "Skipping gene plots: too many genes in range: %s. "%len(l_genes)
        else:
            sys_write("Plotting gene models... ")
            plt.axes([0,0,1,0.18])
            plt.axis([self.chr_offsets[chr-1]+posbeg, self.chr_offsets[chr-1]+posend, -3, 0])
            plt.axis('off')
            broken_barh_xranges = []
            broken_barh_yranges = (-0.5,0.5)
            annotate_y = -1
            for gene in l_genes:
                broken_barh_xranges.append((self.chr_offsets[chr-1]+gene.posbeg, gene.posend-gene.posbeg))
                plt.annotate(gene.id, (self.chr_offsets[chr-1]+gene.posbeg, annotate_y), rotation = 270, size="small") 
            plt.broken_barh(broken_barh_xranges, broken_barh_yranges, facecolors = 'yellow', alpha=0.5)
            sys_write("Done.\n")
        
        # saving figure
        sys_write("Outputting figure to file... ")
        plt.savefig(savefigname,format="png",dpi=300,bbox_inches='tight')
        sys_write("Done.\n")

    def get_prob_distribution(self, chr, posbeg, posend, savefigname):
        " Plot the distribution of SNP scores within a certain range (Sort of non-optimal implementation) "
        l_target_snp = []
        sys_write("Getting list of SNP within range... ")
        for snp in self.res:
            if snp.chr == chr and snp.pos > posbeg and snp.pos < posend:
                l_target_snp.append(snp)
        sys_write("Done.\n")
        
        sys_write("Plotting histogram... ")
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        l_scores = [snp.score for snp in l_target_snp]
        plt.hist(l_scores, bins=25)
        plt.savefig(savefigname,format="png",dpi=300,bbox_inches='tight')
        plt.cla()
        sys_write("Done.\n")     

    def debug(self, outfilename = "", distance_allow = 200000, percentile = 0.02, restricted = False):
        """
        debug function
        """
        sys_write("Ranking result entries... ")
        self.res.sort(key=lambda d:-d.score)
        sys_write("Done.\n")
        
        l_top = []
        for res in self.res[:int(len(self.res)*percentile/100)]:
            l_top.append(res)
        print "Selected top %s SNPs, with a minimum score of %s"%(len(l_top), l_top[-1].score)

        # Sort the top SNPs by position
        sys_write("Sorting top SNPs by position... ")
        l_top.sort()
        sys_write("Done.\n")
        
        return l_top
 
    def find_peaks(self, outfilename = "", distance_allow = 200000, percentile = 1, restricted = False, replace_with_top = False):
        """
        A different approach to the problem of peaking finding, here I just find everything that we can
        - Bad repetitive codes! Should be merged somewhere into a shared private function
        """
        sys_write("Ranking result entries... ")
        self.res.sort(key=lambda d:-d.score)
        sys_write("Done.\n")
        
        l_top = []
        for res in self.res[:int(len(self.res)*percentile/100)]:
            l_top.append(res)
        print "Selected top %s SNPs, with a minimum score of %s"%(len(l_top), l_top[-1].score)

        # Sort the top SNPs by position
        sys_write("Sorting top SNPs by position... ")
        l_top.sort()
        sys_write("Done.\n")
        
        sys_write("Calculating number of peaks... ")
        if replace_with_top == False:
            peaks = []
            current_chr = 0
            new_peak_left_bound = None
            for snp_i, snp in enumerate(l_top):
                if restricted == True:
                    if snp.indicator == 1:
                        continue
                if snp.chr > current_chr:
                    current_chr = snp.chr
                    if new_peak_left_bound != None:
                        peaks.append(l_top[snp_i-1])
                    new_peak_left_bound = None
                    last_peak_right_limit = 0
                if new_peak_left_bound != None:
                    if snp.pos > new_peak_left_bound + distance_allow:
                        new_peak_left_bound = None
                        peaks.append(l_top[snp_i-1])
                        last_peak_right_limit = l_top[snp_i-1].pos + distance_allow
                if new_peak_left_bound == None and snp.pos > last_peak_right_limit:
                    new_peak_left_bound = snp.pos
            # process far right: if reach end, simply add the last snp
            if new_peak_left_bound != None:
                peaks.append(snp)
        
        # Additional step: replace the peaks found with top snps in the peak range. A larger window is allowed in this case
        if replace_with_top == True:
            peaks = []
            current_chr = l_top[0].chr
            current_peak = l_top[0] 
            left_limit = l_top[0].pos
            for snp_i, snp in enumerate(l_top[1:]):
                if restricted == True:
                    if snp.indicator == 1:
                        continue
                if snp.chr > current_chr:
                    peaks.append(current_peak)
                    current_peak = snp
                    current_chr = snp.chr
                    left_limit = snp.pos
                else:
                    if snp.pos > current_peak.pos + distance_allow:
                        peaks.append(current_peak)
                        current_peak = snp
                        left_limit = snp.pos
                    elif snp.score > current_peak.score:
                        if snp.pos > left_limit + distance_allow:
                            # a higher peak, but does not cover the leftmost snp
                            peaks.append(current_peak)
                            current_peak = snp
                            left_limit = snp.pos
                        else:
                            current_peak = snp
            peaks.append(current_peak)
        sys_write("Done.\n")
        return peaks        
    
    '''
    # Old obselete stuff
    
    def compare_peaks(self, outfilename = "", type = "rank", distance_allow_old = 50000, distance_allow_peak = 50000, percentile=1):
        """
        A script to compare peaks between different dataset, calls respective slave functions
        - ONLY WORKS for 2 result sets
        - "New set" should have indicator 0 (Loaded first)
        """
        self.outfilename = outfilename
        if type == "rank":
            return self.__compare_peak_rank__(distance_allow_old, distance_allow_peak, percentile)
        
    def __compare_peak_rank__(self, distance_allow_old = 50000, distance_allow_peak = 50000, percentile = 1):
        """
        Script to compare peaks, defined by rank
        """
        # Sorting through results 
        sys_write("Ranking result entries... ")
        self.res.sort(key=lambda d:-d.score)
        sys_write("Done.\n")
        
        # Selects top percent of results, putting their pointers into a list
        l_top = []
        for res in self.res[:int(len(self.res)*percentile/100)]:
            l_top.append(res)
        print "Selected top %s SNPs, with a minimum score of %s"%(len(l_top), l_top[-1].score)
            
        # Sort the top SNPs by position
        sys_write("Sorting top SNPs by position... ")
        l_top.sort()
        sys_write("Done.\n")
        
        # Now goes through the list, and try to detect new top-snp that is not near old top-snp
        sys_write("Removing top SNPs close to old peaks... ")
        current_chr = 0
        new_top_snp = [] # pointer list of new-top-snps that are not near any old-top-snp
        last_new_topsnp_index = None
        
        new_snp_num = 0
        old_snp_num = 0
        for snp_i, snp in enumerate(l_top):
            if snp.chr > current_chr:
                # process end of chromosomes
                if current_chr != 0 and last_new_topsnp_index != None:
                    for snp_ii in xrange(last_new_topsnp_index, snp_i):
                        if l_top[snp_ii].indicator == 1:
                            new_top_snp.append(l_top[snp_ii])
                last_old_topsnp_end = 0 # This tracks the last position that is within an old-top-snp
                last_new_topsnp_index = None # This tracks the last new-top-snp that is not in range of an old-top-snp
                current_chr = snp.chr
            if snp.indicator != 1: # old set
                old_snp_num += 1
                last_old_topsnp_end = snp.pos + distance_allow_old
                if last_new_topsnp_index != None:
                    for snp_ii in xrange(last_new_topsnp_index, snp_i):
                        if l_top[snp_ii].pos < snp.pos - distance_allow_old:
                            if l_top[snp_ii].indicator == 1:
                                new_top_snp.append(l_top[snp_ii])
                        else:
                            break
                    last_new_topsnp_index = None
            elif snp.indicator == 1: # new set
                new_snp_num += 1
                if last_new_topsnp_index == None:
                    if snp.pos > last_old_topsnp_end:
                        last_new_topsnp_index = snp_i
            # process end of last chromosome
            if snp_i == len(l_top)-1 and last_new_topsnp_index != None:
                for snp_ii in xrange(last_new_topsnp_index, snp_i + 1):
                    if l_top[snp_ii].indicator == 1:
                        new_top_snp.append(l_top[snp_ii])
        sys_write("Done.\n")
        print "%s SNPs left in top SNP list (%s, %s)"%(len(new_top_snp), new_snp_num, old_snp_num)
        
        # Removing solitude?
                        
        # Lastly, calculate the number of new peaks from the new top snp list]
        # The number of peaks is defined here as the minimum number of snps so that every snp in the list are within a certain range of chosen snps
        # This is achieved by a greedy algorithm that proceed from left to right
        sys_write("Calculating number of new peaks... ")
        new_peaks = []
        current_chr = 0
        new_peak_left_bound = None
        for snp_i, snp in enumerate(new_top_snp):
            if snp.chr > current_chr:
                current_chr = snp.chr
                if new_peak_left_bound != None:
                    new_peaks.append(l_top[snp_i-1])
                new_peak_left_bound = None
                last_peak_right_limit = 0
            if new_peak_left_bound != None:
                if snp.pos > new_peak_left_bound + distance_allow_peak:
                    new_peak_left_bound = None
                    new_peaks.append(new_top_snp[snp_i-1])
                    last_peak_right_limit = new_top_snp[snp_i-1].pos + distance_allow_peak
            if new_peak_left_bound == None and snp.pos > last_peak_right_limit:
                new_peak_left_bound = snp.pos
        if new_peak_left_bound != None:
            new_peaks.append(snp)
        sys_write("Done.\n")
        
        return new_peaks
    
    def __compare_peak_hmm__(self, stickiness = 1):
        """
        Complex HMM model to detect peaks
        - Iterated to find the max prob peak vs non-peak distribution
        """
    
    def enrichment_SNS(self, percentile = 1, sns_file = "/Volumes/samaras/Results/Imputed_SNP_list_SNS.csv", snp_list = None):
        """
        Do an enrichment study for synonymous/nonsynonymous SNP
        """
        # get sns information
        sys_write("Importing sns information... ")
        import sns_class
        S = sns_class.SNS(sns_file)
        sys_write("Done.\n")
        
        # Creating the top snp list if none is supplied
            # Sorting through results 
        sys_write("Ranking result entries... ")
        self.res.sort(key=lambda d:-d.score)
        sys_write("Done.\n")
            
        if self.scoremax == None:
            self.scoremax = self.res[0].score
            
        if snp_list == None:
            # Selects top percent of results, putting their pointers into a list
            l_top = []
            for res in self.res[:int(len(self.res)*percentile/100)]:
                l_top.append(res)
            print "Selected top %s SNPs, with a minimum score of %s"%(len(l_top), l_top[-1].score)
            snp_list = l_top

        # Count the number of S and N
        S_count = 0
        N_count = 0
        O_count = 0
        for snp in snp_list:
            tmp_snp_sns = S.get(snp.chr, snp.pos)
            if tmp_snp_sns == "S":
                S_count +=1
            elif tmp_snp_sns == "N":
                N_count +=1
            else:
                O_count +=1
        
        return (S_count, N_count, O_count)
    
    def enrichment_NS_linked(self, percentile = 1, ns_linked_file = "/Volumes/samaras/Results/Link_to_NS/Chr%s.csv", snp_list = None, ns_linked_class = None):
        # get ns-linked information
        if ns_linked_class == None:
            sys_write("Importing ns-linked information... ")
            import sns_class
            S = sns_class.link_NS(ns_linked_file)
            sys_write("Done.\n")
        else:
            S = ns_linked_class
        
        # Creating the top snp list if none is supplied
            # Sorting through results 
        if percentile < 100:
            sys_write("Ranking result entries... ")
            self.res.sort(key=lambda d:-d.score)
            sys_write("Done.\n")
            
        if self.scoremax == None:
            self.scoremax = self.res[0].score
            
        if snp_list == None:
            # Selects top percent of results, putting their pointers into a list
            if percentile == 100:
                snp_list = self.res
            else:
                l_top = []
                for res in self.res[:int(len(self.res)*percentile/100)]:
                    l_top.append(res)
                print "Selected top %s SNPs, with a minimum score of %s"%(len(l_top), l_top[-1].score)
                snp_list = l_top

        # Count the number of 0 and 1
        zero_count = 0
        one_count = 0
        bad_count = 0
        for snp in snp_list:
            tmp_snp_sns = S.get(snp.chr, snp.pos)
            if tmp_snp_sns == 1:
                one_count +=1
            elif tmp_snp_sns == 0:
                zero_count +=1
            else:
                #print "%s,%s"%(snp.chr, snp.pos)
                one_count +=1
        return (zero_count, one_count)
    '''
    
    def plot_pval_dist(self, outfilename = "pval_dist.png"):
        "Plots the pvalue distribution for the different sets"
        pval_lists = []
        for i in xrange(self.num_dif_results):
            pval_lists.append([])
        for res in self.res:
            pval_lists[res.indicator-1]+=[res.score]
            
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        
        plt.hist(pval_lists, normed=True, bins=20)
        plt.savefig(outfilename, format="png")
        
    def plot_pval_qq(self, outfilename = "pval_qq.png"):
        " Plot the qq plot for pval distribution for the different sets"
        # first attempt at qq plot
        pval_lists = []
        for i in xrange(self.num_dif_results):
            pval_lists.append([])
        for res in self.res:
            pval_lists[res.indicator-1]+=[res.score]
        
        pval_lists[0].sort()    
        #import matplotlib
        #matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import scipy.stats as st
        
        #st.probplot(pval_lists[0], dist='uniform', plot=plt)       
        plt.plot(range(len(pval_lists[0])),pval_lists[0], ".")
        plt.savefig(outfilename, format='png')
        