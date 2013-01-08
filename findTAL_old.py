#!/usr/bin/python

# Usage: findTAL.py [options]

# Options:
#   -h, --help
#       show this help message and exit
#   -f FASTA, --fasta=FASTA
#       FASTA file containing promoter sequence(s).
#
#   -j JOB, --job=OUTPUT
#               the job name, output files will have the job name as a prefix
#   -m MIN, --min=MIN
#       minimum of the range of spacer sizes to try
#   -x MAX, --max=MAX
#       maximimum of the range of spacer sizes to try

#findTAL.py
#The same as find_engineered_TAL_sites_v3.py except that all potential binding 
#sites and halfsites are output, not just the best scoring ones at each position

#A script to find possible TAL binding half sites that meet the following
#criteria:

#T at position -1
#15 -21 bp in length
#Appropriate %A, C, G, T
#No T at position 1 or A at position 2
#No G at last or next to last position

#Each half site has ACGT % composition within 2 standard deviations of the mean
#for known binding sites

#Half-sites on opposite strands separated by an spacer of length 12-30 bp are
#considered candidate binding sites

#Modified 10/21/2010 so that the reported binding site of TAL2 is the position
#at which the first RVD binds.  Previously, it listed as the position just 
#after the spacer. Also added TAL length (#RVDs) to the output.

#Modified 11/9/2010 for longer half-sites (up to 30 bases long).  This is to 
#reflect the reality that longer TALs can be quickly assembled with the
#GoldenGate technique.

from Bio.Alphabet import generic_dna

from talconfig import BASE_DIR
from talutil import validate_options_handler, OptParser, FastaIterator, create_logger, check_fasta_pasta, OptionObject, TaskError

celery_found = True
try:
    from celery.task import task
    from celery.registry import tasks
    from talutil import BaseTask
except ImportError:
    celery_found = False

import re

import math
import pickle

import string

#import code
#import signal

#Define a binding site object
class Binding_site:
    def __init__(self, perfectTAL1, perfectTAL2, start1, seq1, start2, seq2_plus, spacer, spacerseq, upstream):
        self.perfectTAL1 = perfectTAL1
        self.perfectTAL2 = perfectTAL2
        self.start1 = start1
        self.seq1 = seq1
        self.start2 = start2
        self.seq2_plus = seq2_plus
        self.seq2_minus = self.seq2_plus.reverse_complement()
        self.spacer = spacer
        self.spacerseq = spacerseq
        self.upstream = upstream
        self.re_sites = []

#DNA list and dictionary
DNA = ['A', 'C', 'G', 'T']
DNA_dict = {'A':0, 'C':1, 'G':2, 'T':3}

#Average percent composition of known TAL binding sites
avg_percents = {'A':0.31, 'C':0.37, 'G':0.09, 'T':0.22}
stdev = {'A':0.16, 'C':0.13, 'G':0.08, 'T':0.10}

#Calculate ranges based on mean +/- 2*stdev
percent_comp_range_top = {}
percent_comp_range_bottom = {}

for nt in DNA:
    percent_comp_range_top[nt] = avg_percents[nt] + 2*stdev[nt]
    percent_comp_range_bottom[nt] = avg_percents[nt] - 2*stdev[nt]

with open(BASE_DIR + "/talent/re_dict_dump", "rb") as re_dict_file:
    NEB_RE_sites = pickle.load(re_dict_file)
    
if celery_found:
    @task(base=BaseTask)
    def FindTALOldTask(*args, **kwargs):
        RunFindTALOldTask(OptionObject(**kwargs))

def validateOptions(options):
    
    if options.fasta == 'NA':
        raise TaskError("FASTA file required.")
    
    if options.cupstream not in [0, 1, 2]:
        raise TaskError("Invalid cupstream value provided")
    
    if options.arraymin < 10 or options.arraymin > 35:
        raise TaskError("Minimum repeat array length must be between 10 and 35")
    
    if options.arraymax < 10 or options.arraymax > 35:
        raise TaskError("Maximum repeat array length must be between 10 and 35")
    
    if options.arraymax < options.arraymin:
        raise TaskError("Maximum repeat array length must be greater than the minimum repeat array length")
    
    if options.min < 10 or options.min > 35:
        raise TaskError("Minimum spacer length must be between 10 and 35")
    
    if options.max < 10 or options.max > 35:
        raise TaskError("Maximum spacer length must be between 10 and 35")
    
    if options.max < options.min:
        raise TaskError("Maximum spacer length cannot be less than the minimum spacer length")
    
    with open(options.fasta, 'r') as seq_file:
        check_fasta_pasta(seq_file)

def RunFindTALOldTask(options):
    
    logger = create_logger(options.logFilepath)
    
    seq_file = open(options.fasta, 'r')
    
    logger("Beginning")
    
    #Set other parameters
    if options.arraymin is None or options.arraymax is None:
        half_site_size = range(15, 31)
    else:
        half_site_size = range(options.arraymin, options.arraymax + 1)
    
    if options.min is None or options.max is None:
        spacer_size = [15, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
    else:
        spacer_size = range(options.min, options.max + 1)
    
    #Initialize half site data structures:
    gene_binding_sites = {}
    
    #Open and read FASTA sequence file
    genes = []
    
    for gene in FastaIterator(seq_file, alphabet=generic_dna):
        genes.append(gene)
    
    seq_file.close()
    
    for gene in genes:
        gene.seq = gene.seq.upper()
    
    #Set up binding site counter for each gene
    binding_site_count_genes = {}
    
    #Scan each gene sequence:
    for gene in genes: #Scan sequence based on above criteria:
        logger("Scanning %s for binding sites" % (gene.id)) 
        binding_site_count_genes[gene] = 0
        sequence = gene.seq
        
        #Check each position along the sequence for possible binding sites using all combinations of binding site lengths and spacer lengths
        for size1 in half_site_size:
            for spacer_len in spacer_size:
                for size2 in half_site_size:
                    for sindex in range(1, len(sequence)-(size1+spacer_len+size2)):
                        
                        #Check for T at -1 for each half_site (A on plus strand of second halfsite)
    
                        if ((options.cupstream != 1 and sequence[sindex-1] == 'T' and sequence[sindex+size1+spacer_len+size2] == 'A') or (options.cupstream != 0 and sequence[sindex-1] == 'C' and sequence[sindex+size1+spacer_len+size2] == 'G')) and len(set(DNA) | set(sequence[sindex:sindex+size1+spacer_len+size2])) ==4:
                            half_site1 = sequence[sindex:sindex+size1]
                            half_site2_plus = sequence[sindex+size1+spacer_len:sindex+size1+spacer_len+size2]
                            half_site2_minus = sequence[sindex+size1+spacer_len:sindex+size1+spacer_len+size2].reverse_complement()
                            spacer_seq = sequence[sindex+size1:sindex+size1+spacer_len]
                            
                            Binding_site_flag = True
    
                            #Check for not T at 1 (A at 1 for second halfsite) on both sites
                            if Binding_site_flag==True and options.t1==True:
                                if sequence[sindex] != 'T' and sequence[sindex+size1+spacer_len+size2-1]!= 'A':
                                    Binding_site_flag=True
                                else:
                                    Binding_site_flag=False
                        
                            #Check not A at 2 (T on plus for second halfsite) on both sites
                            if Binding_site_flag==True and options.a2==True:
                                if sequence[sindex+1] !='A' and sequence[sindex+size1+spacer_len+size2-2]!='T':
                                    Binding_site_flag=True
                                else:
                                    Binding_site_flag=False
    
                            #Require T at end of both sites so bound by NG (A on plus for second half site)
                            if Binding_site_flag==True and options.tn==True:
                                if sequence[sindex+size1-1] == 'T' and sequence[sindex+size1+spacer_len] =='A':
                                    Binding_site_flag=True
                                else:
                                    Binding_site_flag=False
                            
                            #Require last position to not be G's
                            if Binding_site_flag==True and options.gn==True:
                                if sequence[sindex+size1-1] != 'G' and sequence[sindex+size1+spacer_len] != 'C':
                                    Binding_site_flag=True
                                else:
                                    Binding_site_flag=False
    
                            #Check nucleotide composition of the binding site
                            if Binding_site_flag==True and options.comp==True:
                                A1 = half_site1.count('A')/float(len(half_site1))
                                C1 = half_site1.count('C')/float(len(half_site1))
                                G1 = half_site1.count('G')/float(len(half_site1))
                                T1 = half_site1.count('T')/float(len(half_site1))
                            
                                A2 = half_site2_minus.count('A')/float(len(half_site2_minus))
                                C2 = half_site2_minus.count('C')/float(len(half_site2_minus))
                                G2 = half_site2_minus.count('G')/float(len(half_site2_minus))
                                T2 = half_site2_minus.count('T')/float(len(half_site2_minus))
                
                                if A1<=percent_comp_range_top['A'] and A1>=percent_comp_range_bottom['A'] and C1<=percent_comp_range_top['C'] and C1>=percent_comp_range_bottom['C'] and G1<=percent_comp_range_top['G'] and G1>=percent_comp_range_bottom['G'] and T1<=percent_comp_range_top['T'] and T1>=percent_comp_range_bottom['T']:
                                    if A2<=percent_comp_range_top['A'] and A2>=percent_comp_range_bottom['A'] and C2<=percent_comp_range_top['C'] and C2>=percent_comp_range_bottom['C'] and G2<=percent_comp_range_top['G'] and G2>=percent_comp_range_bottom['G'] and T2<=percent_comp_range_top['T'] and T2>=percent_comp_range_bottom['T']:
                                        Binding_site_flag=True
                                    else:
                                        Binding_site_flag=False
                                else:
                                    Binding_site_flag=False
                                    
                            #Create a binding site if all enforced rules have been met
                            if Binding_site_flag==True:
                                binding_site = Binding_site(perfectTAL1 = 'none', perfectTAL2 = 'none', start1 = sindex, start2=sindex+size1+spacer_len+size2-1, seq1 = half_site1, seq2_plus=half_site2_plus, spacer=spacer_len, spacerseq=spacer_seq, upstream = sequence[sindex-1])
                                if gene not in gene_binding_sites.keys():
                                    gene_binding_sites[gene] = {}
    
                                if sindex not in gene_binding_sites[gene].keys():
                                    gene_binding_sites[gene][sindex] = []
                                        
                                gene_binding_sites[gene][sindex].append(binding_site)
                                binding_site_count_genes[gene] +=1
    
    #Compute TALs for each gene, using "strong-binding" RVDs for each nucleotide (binds the nucleotide more than half the time and we have more than 10 observations)
    logger('Designing best scoring perfect TALs for each potential site...')
    strong_binding_RVDs = {'A':'NI', 'C':'HD', 'G':'NN', 'T':'NG'}
    RVD_pairs = {} #dictionary of RVD pair counts indexed by gene, RVD1 and RVD2
    for gene in gene_binding_sites.keys():
        RVD_pairs[gene] = {}
        for start in gene_binding_sites[gene].keys():
            
            #Find the perfect RVD sequence from each potential plus strand start site
            for binding_site in gene_binding_sites[gene][start]:
                TAL_1 = []
                TAL_2 = []
    
                for bindex in range(0, len(binding_site.seq1)):
                    TAL_1.append(strong_binding_RVDs[binding_site.seq1[bindex]])
    
                TAL_1 = ' '.join(TAL_1)
                
                for bindex in range(0, len(binding_site.seq2_minus)):
                    TAL_2.append(strong_binding_RVDs[binding_site.seq2_minus[bindex]])
    
                TAL_2 = ' '.join(TAL_2)
    
                binding_site.perfectTAL1 = TAL_1
                binding_site.perfectTAL2 = TAL_2
                
                if TAL_1 not in RVD_pairs[gene].keys():
                    RVD_pairs[gene][TAL_1] = {}
                if TAL_2 not in RVD_pairs[gene][TAL_1].keys():
                    RVD_pairs[gene][TAL_1][TAL_2] = 0
                RVD_pairs[gene][TAL_1][TAL_2] += 1
    
    #Find lists of binding sites in each gene with unique (for that gene)  perfect RVD sequences
    binding_sites_unique_plus_minus_pairs = {} #binding sites whose RVD sequences don't make up other binding sites in the gene
    
    for gene in gene_binding_sites.keys():
        binding_sites_unique_plus_minus_pairs[gene] = {}
        for start in gene_binding_sites[gene].keys():
            for binding_site in gene_binding_sites[gene][start]:
                RVD1 = binding_site.perfectTAL1
                RVD2 = binding_site.perfectTAL2
                
                if (RVD1 not in RVD_pairs[gene].keys() or RVD1 not in RVD_pairs[gene][RVD1].keys() or RVD_pairs[gene][RVD1][RVD1] == 1) and (RVD1 not in RVD_pairs[gene].keys() or RVD2 not in RVD_pairs[gene][RVD1].keys() or RVD_pairs[gene][RVD1][RVD2] == 1) and (RVD2 not in RVD_pairs[gene].keys() or RVD1 not in RVD_pairs[gene][RVD2].keys() or RVD_pairs[gene][RVD2][RVD1] == 1) and (RVD2 not in RVD_pairs[gene].keys() or RVD2 not in RVD_pairs[gene][RVD2].keys() or RVD_pairs[gene][RVD2][RVD2] == 1):
                    if start not in binding_sites_unique_plus_minus_pairs[gene].keys():
                        binding_sites_unique_plus_minus_pairs[gene][start] = []
                    
                    binding_sites_unique_plus_minus_pairs[gene][start].append(binding_site)
    
    
    #Check binding sites for unique Restriction endonuclease sites within the spacer.
    #Unique sites are those that occur once in the spacer and do not occur in the 250 bases on either side of the spacer.
    logger('Searching for restriction enzymes sites within each spacer...')
    for gene in binding_sites_unique_plus_minus_pairs.keys():
        for start_site in binding_sites_unique_plus_minus_pairs[gene].keys():
            for binding_site in binding_sites_unique_plus_minus_pairs[gene][start_site]:
    
                enzymes_in_spacer = []
            
                spacer_start = binding_site.start1 + len(binding_site.seq1)
                spacer_end = binding_site.start1 + len(binding_site.seq1) + binding_site.spacer-1
                
                #identify sequence to check around the spacer for unique-ness
                if spacer_start >= 250:
                    seq_check_start = spacer_start - 250
                else:
                    seq_check_start = 0
                
                if len(gene.seq) - spacer_end >= 250:
                    seq_check_end = spacer_end + 250 + 1
                else:
                    seq_check_end = len(gene.seq)
                
                #For each enzyme check if it occurs once in the spacer:
                for enzyme in NEB_RE_sites:
                    if len(NEB_RE_sites[enzyme]["compiled"].findall(str(binding_site.spacerseq))) == 1:
                
                        #If unique in spacer, check that it doesen't occur in  the flanking sequence                
                        if len(NEB_RE_sites[enzyme]["compiled"].findall(str(gene.seq[seq_check_start : seq_check_end]))) == 1:                  
                            enzymes_in_spacer.append(enzyme)
                
                #Create a string listing the enzymes and their sequences that can printed in the output
                enzyme_string = ' '.join(["%s:%s" % (enzyme, NEB_RE_sites[enzyme]["short"]) for enzyme in enzymes_in_spacer])
                
                if len(enzyme_string) == 0:
                    enzyme_string = 'none'
                
                #append enzyme string to binding site object
                binding_site.re_sites = enzyme_string
                
    #Print output results to file: binding sites
    if options.outpath == 'NA':
      filename = options.outdir + options.job + options.outfile
    else:
      filename = options.outpath
    
    out = open(filename, 'w')
    table_ignores = ["TAL1 length", "TAL2 length", "Spacer length"]
    out.write("table_ignores:" + string.join(table_ignores, ",") + "\n")
    
    u_bases = []
    
    if options.cupstream != 1:
        u_bases.append("T")
    
    if options.cupstream != 0:
        u_bases.append("C")
        
    out.write("options_used:" + ', '.join([
        "array_min = " + str(options.arraymin),
        "array_max = " + str(options.arraymax),
        "spacer_min = " + str(options.min),
        "spacer_max = " + str(options.max),
        "upstream_base = " + (" or ".join(u_bases)),
        ("No T at position 1" if options.t1 else ""),
        ("No A at position 1" if options.a2 else ""),
        ("Sites must end in a T" if options.tn else ""),
        ("Sites may not end in G/NN" if options.gn else ""),
        ("Base composition rules enforced" if options.comp else ""),
    ]) + "\n")
    
    out.write('Sequence Name\tTAL1 start\tTAL2 start\tTAL1 length\tTAL2 length\tSpacer length\tSpacer range\tTAL1 RVDs\tTAL2 RVDs\tPlus strand sequence\tUnique_RE_sites_in_spacer\n')
    if len(binding_sites_unique_plus_minus_pairs.keys()) > 0:
        for gene in binding_sites_unique_plus_minus_pairs.keys():
            for start_site in binding_sites_unique_plus_minus_pairs[gene].keys():
                for binding_site in binding_sites_unique_plus_minus_pairs[gene][start_site]:
                    out.write(gene.id + '\t' + str(binding_site.start1) + '\t' + str(binding_site.start2) + '\t' + str(len(binding_site.seq1)) + '\t' + str(len(binding_site.seq2_plus)) + '\t' + str(binding_site.spacer) + '\t' + str(binding_site.start1+len(binding_site.seq1)) + '-' + str(binding_site.start1+len(binding_site.seq1) + binding_site.spacer-1) + '\t' + binding_site.perfectTAL1 + '\t' +  binding_site.perfectTAL2 + '\t' + binding_site.upstream.upper() + " " + str(binding_site.seq1) + ' ' + str(binding_site.spacerseq).lower() + ' ' + str(binding_site.seq2_plus) + " " + ("A" if binding_site.upstream == "T" else "G") + '\t' + binding_site.re_sites + '\n')
    
    out.close()
    logger('Finished')
    
if __name__ == '__main__':
    
    #signal.signal(signal.SIGUSR2, lambda sig, frame: code.interact())
    
    # import arguments and options
    usage = 'usage: %prog [options]'
    parser = OptParser(usage=usage)
    parser.add_option('-f', '--fasta', dest='fasta', type='string', default='NA', help='FASTA file containing promoter sequence(s).')
    parser.add_option('-j', '--job', dest='job', type='string', default='output', help='the job name, output files will have the job name as a prefix.')
    parser.add_option('-m', '--min', dest='min', type='int', default=None, help='the minimum spacer size to try')
    parser.add_option('-x', '--max', dest='max', type='int', default=None, help='the maximum spacer size to try')
    parser.add_option('-t', '--t1', dest='t1', action='store_false', default = True, help='Do not enforce rule "No T at position 1."')
    parser.add_option('-a', '--a2', dest='a2', action='store_false', default = True, help='Do not enforce rule "No A at position 1."')
    parser.add_option('-n', '--tn', dest='tn', action='store_false', default = True, help='Do not enforce rule "Sites must end in a T."')
    parser.add_option('-g', '--gn', dest='gn', action='store_false', default = True, help='Do not enforce rule "Sites may not end in G/NN."')
    parser.add_option('-c', '--comp', dest='comp', action='store_false', default = True, help='Do not enforce base composition rules.')
    parser.add_option('-d', '--outdir', dest='outdir', type='string', default = 'upload/', help='Directory in which to place output files.')
    parser.add_option('-o', '--outfile', dest='outfile', type='string', default = '_TALEN_pairs_all.txt', help='Optional filename for output file.')
    parser.add_option('-r', '--arraymin', dest='arraymin', type='int', default=None, help='the minimum repeat array length to try')
    parser.add_option('-y', '--arraymax', dest='arraymax', type='int', default=None, help='the maximum repeat array length to try')
    parser.add_option('-u', '--cupstream', dest='cupstream', type='int', default = 0, help='1 to look for C instead of T, 2 to look for either')
    #Drupal options
    parser.add_option('-p', '--outpath', dest='outpath', type='string', default = 'NA', help='Optional full path for output file; if set --job, --outdir and --outfile are ignored.')
    parser.add_option('-l', '--logpath', dest='logFilepath', type='string', default = 'NA', help='Process log file path')
    parser.add_option('-z', '--nodeid', dest='nodeID', type='int', default = '-1', help='Optional node id if this script was called from Drupal.')
    parser.add_option('-k', '--ipaddr', dest='ip_address', type='string', default = '', help='IP address of job submitter')
    (options, args) = parser.parse_args()
    
    validate_options_handler(validateOptions, options)
    
    if options.nodeID != -1:
        
        # if run from drupal then it should be queued as a task
        
        if not celery_found:
            raise TaskError("nodeID option was provided but Celery is not installed")
        
        logger = create_logger(options.logFilepath)
        logger("Your task has been queued and will be processed when a worker node becomes available")
        
        # This is ensures that the task is queued with the same module name
        # that the Celery workers are expecting 
        from findTAL_old import FindTALOldTask
        
        FindTALOldTask.apply_async(kwargs=vars(options), queue="findtal_old")
        
    else:
        
        RunFindTALOldTask(options)
