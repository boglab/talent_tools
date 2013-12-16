#!/usr/bin/python

# Usage: findSingleTALsite.py [options]

# Options:
#    -h, --help
#        show this help message and exit
#    -f FASTA, --fasta=FASTA
#        FASTA file containing promoter sequence(s).
#    -j JOB, --job=OUTPUT
#        the job name, output files will have the job name as a prefix
#    -r revcomp, --revcomp
#        A flag indicating that the reverse compliment of the supplied FASTA sequence should be searched

#findSingleTALsite.py
#Created May 25, 2011 by Erin Doyle

#A script to find possible TAL binding sites for single TALs
#Based on the script findTAL.py

#T at position -1
#15-21 bp in length
#Appropriate %A, C, G, T
#No T at position 1 or A at position 2
#No G at last or next to last position

#Each half site has ACGT % composition within 2 standard deviations of the mean
#for known binding sites

from Bio.Alphabet import generic_dna

from talconfig import GENOME_FILE, PROMOTEROME_FILE, VALID_GENOME_ORGANISMS, VALID_PROMOTEROME_ORGANISMS, OFFTARGET_COUNTING_SIZE_LIMIT
from talutil import validate_options_handler, OptParser, FastaIterator, create_logger, check_fasta_pasta, OptionObject, TaskError, Conditional
from entrez_cache import CachedEntrezFile

celery_found = True
try:
    from celery.task import task
    from celery.registry import tasks
    from talutil import BaseTask
except ImportError:
    celery_found = False

tfcount_found = True
try:
    from btfcount import TargetFinderCountTask
except ImportError:
    tfcount_found = False

import os
import math
import string

#import code
#import signal

#Define a binding site object
class Binding_site:
    def __init__(self, perfectTAL1, start1, seq1, is_plus, upstream, gene_id = ""):
        self.gene_id = gene_id
        self.perfectTAL1 = perfectTAL1
        self.start1 = start1
        self.seq1 = seq1
        self.is_plus = is_plus
        self.upstream = upstream
        self.offtarget_count = 0

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

if celery_found:
    @task(base=BaseTask)
    def FindSingleTALSiteTask(*args, **kwargs):
        RunFindSingleTALSiteTask(OptionObject(**kwargs))

def validateOptions(options):
    
    if options.fasta == 'NA':
        raise TaskError('FASTA file required.')
    
    if options.cupstream not in [0, 1, 2]:
        raise TaskError("Invalid cupstream value provided")
    
    if options.arraymin < 10 or options.arraymin > 35:
        raise TaskError("Minimum repeat array length must be between 10 and 35")
    
    if options.arraymax < 10 or options.arraymax > 35:
        raise TaskError("Maximum repeat array length must be between 10 and 35")
    
    if options.arraymax < options.arraymin:
        raise TaskError("Maximum repeat array length must be greater than the minimum repeat array length")
    
    if options.offtargets_ncbi != "NA":
        
        if options.offtargets_fasta != "NA" or options.genome or options.promoterome:
            raise TaskError("--offtargets-fasta, --genome and --promoterome options cannot be combined with --offtargets-ncbi")
        
        # NCBI sequence validation is performed after the task has started instead of here to avoid having to download large files more than once
        
        options.check_offtargets = True
    
    if options.offtargets_fasta != "NA":
        
        if options.offtargets_ncbi != "NA" or options.genome or options.promoterome:
            raise TaskError("--offtargets-ncbi, --genome and --promoterome options cannot be combined with --offtargets-fasta")
        
        if (not os.path.exists(options.offtargets_fasta) or os.path.getsize(options.offtargets_fasta) <= 2):
            raise TaskError("Off-target FASTA file must exist and be non-empty.")
        
        options.check_offtargets = True
    
    if options.check_offtargets:
        
        if ((options.genome and options.organism not in VALID_GENOME_ORGANISMS) or (options.promoterome and options.organism not in VALID_PROMOTEROME_ORGANISMS)):
            raise TaskError("Invalid organism specified.")
    
    with open(options.fasta, 'r') as seq_file:
        
        check_fasta_pasta(seq_file)
        
        if options.check_offtargets:
            
            gene_length_total = 0
            
            for gene in FastaIterator(seq_file, alphabet=generic_dna):
                gene_length_total += len(gene.seq)
            
            if gene_length_total > 1000:
                raise TaskError("Off-target counting is only available when designing TALEs for sequences that are 1000 bases or less")

def RunFindSingleTALSiteTask(options):

    logger = create_logger(options.logFilepath)
    
    logger("Beginning")
    
    if options.check_offtargets and options.offtargets_ncbi != "NA":
        logger("Retrieving NCBI off-target sequence. This could take a while if this sequence hasn't been used recently and needs to be downloaded from NCBI.")
    
    with Conditional(options.check_offtargets and options.offtargets_ncbi != "NA", CachedEntrezFile(logger, options.offtargets_ncbi)) as maybe_entrez_file:
        
        if options.check_offtargets:
            
            if not tfcount_found:
                raise TaskError("Non off-target counting worker attempted to process off-target counting task.")
            
            if options.offtargets_ncbi != "NA":
                
                logger("Finished retrieving NCBI off-target sequence.")
                
                # Validate downloaded sequence
                
                check_fasta_pasta(maybe_entrez_file.file)
                
                for record in FastaIterator(maybe_entrez_file.file, alphabet=generic_dna):
                    if len(record.seq) > OFFTARGET_COUNTING_SIZE_LIMIT:
                        raise TaskError("Off-Target counting is only supported for NCBI records where all individual sequences are under %d megabases in size" % (OFFTARGET_COUNTING_SIZE_LIMIT / 1000000))
            
            offtarget_seq_filename = ""
            
            if options.offtargets_fasta != "NA":
                offtarget_seq_filename = options.offtargets_fasta
            elif options.offtargets_ncbi != "NA":
                offtarget_seq_filename = maybe_entrez_file.filepath
            elif options.genome:
                offtarget_seq_filename = GENOME_FILE % options.organism
            elif options.promoterome:
                offtarget_seq_filename = PROMOTEROME_FILE % options.organism
            else:
                offtarget_seq_filename = options.fasta
        
        strong_binding_RVDs = {
            'A':'NI',
            'C':'HD',
            'G':'NN',
            'T':'NG'
        }
        
        if options.gspec:
            strong_binding_RVDs['G'] = 'NH'
        
        seq_file = open(options.fasta, 'r')
        
        #Set other parameters
        if options.arraymin is None or options.arraymax is None:
            half_site_size = range(15, 31)
        else:
            half_site_size = range(options.arraymin, options.arraymax + 1)
        
        #Initialize half site data structures:
        gene_binding_sites = {}
        
        #Open and read FASTA sequence file
        genes = []
        
        for gene in FastaIterator(seq_file, alphabet=generic_dna):
            genes.append(gene)
        
        seq_file.close()
        
        for gene in genes:
            gene.seq = gene.seq.upper()
        
        #Scan each gene sequence:
        for gene in genes: #Scan sequence based on above criteria:
            logger("Scanning %s for binding sites" % (gene.id))
            sequence = gene.seq
            
            #Check each position along the sequence for possible binding sites using all combinations of binding site lengths and spacer lengths
            for size1 in half_site_size:
                for sindex in range(1, len(sequence)-size1):
            
                    #Check for T at -1
                    if ((options.cupstream != 1 and sequence[sindex-1] == 'T') or (options.cupstream != 0 and sequence[sindex-1] == 'C')) and len(set(DNA) | set(sequence[sindex:sindex+size1])) ==4:
                        half_site1 = sequence[sindex:sindex+size1]
                        Binding_site_flag = True
        
                        #Check for not T at 1
                        if Binding_site_flag==True and options.t1==True:
                            if sequence[sindex] != 'T':
                                Binding_site_flag=True
                            else:
                                Binding_site_flag=False
                            
        
                        #Check not A at 2
                        if Binding_site_flag==True and options.a2==True:
                            if sequence[sindex+1] !='A':
                                Binding_site_flag=True
                            else:
                                Binding_site_flag=False
        
                        #Require T at end
                        if Binding_site_flag==True and options.tn==True:
                            if sequence[sindex+size1-1] == 'T':
                                Binding_site_flag=True
                            else:
                                Binding_site_flag=False
                                
                        #Require last position to not be G's
                        if Binding_site_flag==True and options.gn==True:
                            if sequence[sindex+size1-1] != 'G':
                                Binding_site_flag=True
                            else:
                                Binding_site_flag=False
        
                        #Check nucleotide composition of the binding site
                        if Binding_site_flag==True and options.comp==True:
                            A1 = half_site1.count('A')/float(len(half_site1))
                            C1 = half_site1.count('C')/float(len(half_site1))
                            G1 = half_site1.count('G')/float(len(half_site1))
                            T1 = half_site1.count('T')/float(len(half_site1))
                    
                            if A1<=percent_comp_range_top['A'] and A1>=percent_comp_range_bottom['A'] and C1<=percent_comp_range_top['C'] and C1>=percent_comp_range_bottom['C'] and G1<=percent_comp_range_top['G'] and G1>=percent_comp_range_bottom['G'] and T1<=percent_comp_range_top['T'] and T1>=percent_comp_range_bottom['T']:
                                Binding_site_flag=True
                            else:
                                Binding_site_flag=False
        
                                        
                        #Create a binding site if all enforced rules have been met
                        if Binding_site_flag==True:
                            binding_site = Binding_site(perfectTAL1 = 'none', start1 = sindex, seq1 = half_site1, is_plus=True, upstream=sequence[sindex-1])
                            if gene not in gene_binding_sites.keys():
                                gene_binding_sites[gene] = {}
        
                            if sindex not in gene_binding_sites[gene].keys():
                                gene_binding_sites[gene][sindex] = []
                                            
                            gene_binding_sites[gene][sindex].append(binding_site)
        
        
                if options.revcomp==True: #Search for binding sites on the reverse complement strand
                    for sindex in range(size1-1, len(sequence)-1):
                            
                        #Check for T at -1 for each half_site (A on plus strand)
                        
                        if ((options.cupstream != 1 and sequence[sindex+1] == 'A') or (options.cupstream != 0 and sequence[sindex+1] == 'G')) and len(set(DNA) | set(sequence[sindex-size1+1:sindex+1])) == 4:
                            half_site1 = sequence[sindex-size1+1:sindex+1]
                            Binding_site_flag = True
        
                            #Check for not T at 1 (A at 1 on plus strand)
                            if Binding_site_flag==True and options.t1==True:
                                if sequence[sindex] != 'A':
                                    Binding_site_flag=True
                                else:
                                    Binding_site_flag=False
                            
                            #Check not A at 2 (T on plus strand)
                            if Binding_site_flag==True and options.a2==True:
                                if sequence[sindex-1] !='T':
                                    Binding_site_flag=True
                                else:
                                    Binding_site_flag=False
        
                            #Require T at end so bound by NG (A on plus)
                            if Binding_site_flag==True and options.tn==True:
                                if sequence[sindex-size1+1] =='A':
                                    Binding_site_flag=True
                                else:
                                    Binding_site_flag=False
                                
                            #Require last position to not be G (C on plus)
                            if Binding_site_flag==True and options.gn==True:
                                if sequence[sindex-size1+1] != 'C':
                                    Binding_site_flag=True
                                else:
                                    Binding_site_flag=False
        
                            #Check nucleotide composition of the binding site
                            if Binding_site_flag==True and options.comp==True:
                                
                                A2 = half_site1.count('T')/float(len(half_site1))
                                C2 = half_site1.count('G')/float(len(half_site1))
                                G2 = half_site1.count('C')/float(len(half_site1))
                                T2 = half_site1.count('A')/float(len(half_site1))
                    
                                if A2<=percent_comp_range_top['A'] and A2>=percent_comp_range_bottom['A'] and C2<=percent_comp_range_top['C'] and C2>=percent_comp_range_bottom['C'] and G2<=percent_comp_range_top['G'] and G2>=percent_comp_range_bottom['G'] and T2<=percent_comp_range_top['T'] and T2>=percent_comp_range_bottom['T']:
                                    Binding_site_flag=True
                                else:
                                    Binding_site_flag=False
                                        
                            #Create a binding site if all enforced rules have been met
                            if Binding_site_flag==True:
                                binding_site = Binding_site(perfectTAL1 = 'none', start1 = sindex, seq1 = half_site1, is_plus=False, upstream=sequence[sindex+1])
                                if gene not in gene_binding_sites.keys():
                                    gene_binding_sites[gene] = {}
        
                                if sindex not in gene_binding_sites[gene].keys():
                                    gene_binding_sites[gene][sindex] = []
                                        
                                gene_binding_sites[gene][sindex].append(binding_site)
        
        
        #Compute TALs for each gene, using "strong-binding" RVDs for each nucleotide (binds the nucleotide more than half the time and we have more than 10 observations)
        logger('Designing best scoring perfect TALs for each potential site...')
        
        for gene in gene_binding_sites.keys():
            for start in gene_binding_sites[gene].keys():
                #Find the perfect RVD sequence from each potential plus strand start site
                for binding_site in gene_binding_sites[gene][start]:
                    TAL_1 = []
                    if binding_site.is_plus:
                        for bindex in range(0, len(binding_site.seq1)):
                            TAL_1.append(strong_binding_RVDs[binding_site.seq1[bindex]])
        
                        TAL_1 = ' '.join(TAL_1)
                    else:
                        rev_comp_seq = binding_site.seq1.reverse_complement()
                        for bindex in range(0, len(rev_comp_seq)):
                            TAL_1.append(strong_binding_RVDs[rev_comp_seq[bindex]])
                        TAL_1 = ' '.join(TAL_1)
        
                    binding_site.perfectTAL1 = TAL_1
        
        #Print output results to file: binding sites
        
        #filename = 'upload/'+ options.job + '_TALEN_pairs_all.txt'
        
        if options.outpath == 'NA':
          filename = options.outdir + options.job + options.outfile
        else:
          filename = options.outpath
        
        binding_sites = []
        if len(gene_binding_sites.keys()) > 0:
            for gene in sorted(gene_binding_sites.keys()):
                for start_site in gene_binding_sites[gene].keys():
                    for binding_site in gene_binding_sites[gene][start_site]:
                        binding_site.gene_id = gene.id
                        binding_sites.append(binding_site)
        
        if options.check_offtargets:
            
            if len(binding_sites) > 0:
                
                off_target_seqs = []
                
                for i, binding_site in enumerate(binding_sites):
                    off_target_seqs.append(binding_site.perfectTAL1)
                
                off_target_counts = TargetFinderCountTask(offtarget_seq_filename, options.logFilepath, options.cupstream, 3.0, off_target_seqs)
                
                for i, binding_site in enumerate(binding_sites):
                    binding_site.offtarget_count = off_target_counts[i]
        
        out = open(filename, 'w')
        table_ignores = []
        if not options.revcomp:
            table_ignores.append("Plus strand sequence")
        if len(table_ignores) > 0:
            out.write("table_ignores:" + string.join(table_ignores, ",") + "\n")
        
        u_bases = []
        
        if options.cupstream != 1:
            u_bases.append("T")
        
        if options.cupstream != 0:
            u_bases.append("C")
            
        out.write("options_used:" + ', '.join([
            "array_min = " + str(options.arraymin),
            "array_max = " + str(options.arraymax),
            "upstream_base = " + (" or ".join(u_bases)),
            ("No T at position 1" if options.t1 else ""),
            ("No A at position 1" if options.a2 else ""),
            ("Sites must end in a T" if options.tn else ""),
            ("Sites may not end in G/NN" if options.gn else ""),
            ("Base composition rules enforced" if options.comp else ""),
            ("Search reverse complement" if options.revcomp else ""),
        ]) + "\n")
        
        offtarget_header = "\tOff-Target Counts" if options.check_offtargets else ""
        
        out.write('Sequence Name\tTAL start\tTAL length\tRVD sequence\tStrand\tTarget sequence\tPlus strand sequence' + offtarget_header + '\n')
        
        for i, binding_site in enumerate(binding_sites):
            
            offtarget_string = ""
            
            if options.check_offtargets:
                offtarget_string = "\t%d" % binding_site.offtarget_count
            
            if binding_site.is_plus:
                out.write(binding_site.gene_id + '\t' + str(binding_site.start1) + '\t' + str(len(binding_site.seq1)) + '\t' + binding_site.perfectTAL1 + '\t' +  'Plus' + '\t' + binding_site.upstream + " " + str(binding_site.seq1) + '\t' + binding_site.upstream + " " + str(binding_site.seq1) + offtarget_string + '\n')
            else:
                out.write(binding_site.gene_id + '\t' + str(binding_site.start1) + '\t' + str(len(binding_site.seq1)) + '\t' + binding_site.perfectTAL1 + '\t' +  'Minus' + '\t' + ("T" if binding_site.upstream == "A" else "C") + " " + str(binding_site.seq1.reverse_complement()) + '\t' + str(binding_site.seq1) + " " + binding_site.upstream + offtarget_string + '\n')
        
        out.close()
        
        logger('Finished')
    
if __name__ == '__main__':
    
    #signal.signal(signal.SIGUSR2, lambda sig, frame: code.interact())
    
    # import arguments and options
    usage = 'usage: %prog [options]'
    parser = OptParser(usage=usage)
    parser.add_option('-f', '--fasta', dest='fasta', type='string', default='NA', help='FASTA file containing promoter sequence(s).')
    parser.add_option('-j', '--job', dest='job', type='string', default='output', help='the job name, output files will have the job name as a prefix.')
    parser.add_option('-t', '--t1', dest='t1', action='store_false', default = True, help='Do not enforce rule "No T at position 1."')
    parser.add_option('-a', '--a2', dest='a2', action='store_false', default = True, help='Do not enforce rule "No A at position 1."')
    parser.add_option('-n', '--tn', dest='tn', action='store_false', default = True, help='Do not enforce rule "Sites must end in a T."')
    parser.add_option('-g', '--gn', dest='gn', action='store_false', default = True, help='Do not enforce rule "Sites may not end in G/NN."')
    parser.add_option('-c', '--comp', dest='comp', action='store_false', default = True, help='Do not enforce base composition rules.')
    parser.add_option('-r', '--revcomp', dest='revcomp', action = 'store_true', default = False, help='Search the reverse complement of the input FASTA sequences')
    parser.add_option('-d', '--outdir', dest='outdir', type='string', default = 'upload/', help='Directory in which to place output files.')
    parser.add_option('-o', '--outfile', dest='outfile', type='string', default = '_TALEN_pairs_all.txt', help='Optional filename for output file.')
    parser.add_option('-b', '--arraymin', dest='arraymin', type='int', default=None, help='the minimum repeat array length to try')
    parser.add_option('-y', '--arraymax', dest='arraymax', type='int', default=None, help='the maximum repeat array length to try')
    parser.add_option('-u', '--cupstream', dest='cupstream', type='int', default = 0, help='1 to look for C instead of T, 2 to look for either')
    parser.add_option('-s', '--gspec', dest='gspec', action='store_true', default = False, help='If true, use NH instead of NN for G')
    # Offtarget Options
    parser.add_option('--offtargets', dest='check_offtargets', action = 'store_true', default = False, help='Check offtargets')
    parser.add_option('--offtargets-fasta', dest='offtargets_fasta', type='string', default='NA', help='FASTA file containing to search for off-targets')
    parser.add_option('--offtargets-ncbi', dest='offtargets_ncbi', type='string', default='NA', help='NCBI nucleotide sequence ID to search for off-targets')
    parser.add_option('--genome', dest='genome', action = 'store_true', default = False, help='Input is a genome file')
    parser.add_option('--promoterome', dest='promoterome', action = 'store_true', default = False, help='Input is a promoterome file')
    parser.add_option('--organism', dest='organism', type = 'string', default='NA', help='Name of organism for the genome to be searched.')
    #Drupal options
    parser.add_option('-p', '--outpath', dest='outpath', type='string', default = 'NA', help='Optional full path for output file; if set --job, --outdir and --outfile are ignored.')
    parser.add_option('-l', '--logpath', dest='logFilepath', type='string', default = 'NA', help='Process log file path')
    parser.add_option('-z', '--nodeid', dest='nodeID', type='int', default = '-1', help='Optional node id if this script was called from Drupal.')
    parser.add_option('-k', '--ipaddr', dest='ip_address', type='string', default = '', help='IP address of job submitter')
    (options, args) = parser.parse_args()
    
    validate_options_handler(validateOptions, options)
    
    if options.nodeID != -1:
        
        if not celery_found:
            raise TaskError("nodeID option was provided but Celery is not installed")
        
        logger = create_logger(options.logFilepath)
        logger("Your task has been queued and will be processed when a worker node becomes available")
        
        from findSingleTALsite import FindSingleTALSiteTask
        
        if options.check_offtargets:
            job_queue = "findtal_offtargets"
        else:
            job_queue = "findsingletal"
        
        FindSingleTALSiteTask.apply_async(kwargs=vars(options), queue=job_queue)
        
    else:
        
        RunFindSingleTALSiteTask(options)
