#!/usr/bin/python

from Bio.Alphabet import generic_dna

from talconfig import BASE_DIR, GENOME_FILE, PROMOTEROME_FILE, VALID_GENOME_ORGANISMS, VALID_PROMOTEROME_ORGANISMS, OFFTARGET_COUNTING_SIZE_LIMIT
from talutil import validate_options_handler, OptParser, FastaIterator, create_logger, check_fasta_pasta, OptionObject, TaskError, reverseComplement, Conditional
from entrez_cache import CachedEntrezFile

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
import os
from itertools import ifilterfalse

tfcount_found = True
try:
    from btfcount import PairedTargetFinderCountTask
except ImportError:
    tfcount_found = False

#Define a binding site object
class BindingSite:
    
    def __init__(self, **kwargs):
        
        self.seq_id = kwargs.pop("seq_id", 0)
        
        self.cutsite = kwargs.pop("cutsite", 0)
        self.seq1_start = kwargs.pop("seq1_start", 0)
        self.seq1_end = kwargs.pop("seq1_end", 0)
        self.seq1_seq = kwargs.pop("seq1_seq", "")
        self.seq1_rvd = kwargs.pop("seq1_rvd", "")
        
        self.spacer_start = kwargs.pop("spacer_start", 0)
        self.spacer_end = kwargs.pop("spacer_end", 0)
        self.spacer_seq = kwargs.pop("spacer_seq", "")
        
        self.seq2_start = kwargs.pop("seq2_start", 0)
        self.seq2_end = kwargs.pop("seq2_end", 0)
        self.seq2_seq = kwargs.pop("seq2_seq", "")
        self.seq2_rvd = kwargs.pop("seq2_rvd", "")
        
        self.upstream = kwargs.pop("upstream", "")
        
        self.cg_percent = kwargs.pop("cg_percent", 0)
        
        self.offtarget_counts = kwargs.pop("offtarget_counts", [0, 0, 0, 0, 0])
        
        self.re_sites = ""

with open(BASE_DIR + "/talent/re_dict_dump", "rb") as re_dict_file:
    NEB_RE_sites = pickle.load(re_dict_file)

streubel_at_streak_re = re.compile('[AT]{6,}')

if celery_found:
    @task(base=BaseTask)
    def FindTALTask(*args, **kwargs):
        RunFindTALTask(OptionObject(**kwargs))

def findRESitesInSpacer(sequence, binding_site):
    
    enzymes_in_spacer = []

    #identify sequence to check around the spacer for unique-ness
    
    if binding_site.spacer_start >= 250:
        seq_check_start = binding_site.spacer_start - 250
    else:
        seq_check_start = 0
    
    if len(sequence) - binding_site.spacer_end >= 250:
        seq_check_end = binding_site.spacer_end + 250 + 1
    else:
        seq_check_end = len(sequence)
    
    #For each enzyme check if it occurs once in the spacer:
    for enzyme in NEB_RE_sites:
        if len(NEB_RE_sites[enzyme]["compiled"].findall(binding_site.spacer_seq)) == 1:
    
            #If unique in spacer, check that it doesen't occur in  the flanking sequence                
            if len(NEB_RE_sites[enzyme]["compiled"].findall(sequence[seq_check_start : seq_check_end])) == 1:                   
                enzymes_in_spacer.append(enzyme)
    
    #Create a string listing the enzymes and their sequences that can printed in the output
    enzyme_string = ' '.join(["%s:%s" % (enzyme, NEB_RE_sites[enzyme]["short"]) for enzyme in enzymes_in_spacer])
    
    if len(enzyme_string) == 0:
        enzyme_string = 'none'
    
    #append enzyme string to binding site object
    binding_site.re_sites = enzyme_string

def filterByTALSize(x, y):
    x_average_tal_len = float(len(x.seq1_seq) + len(x.seq2_seq)) / 2
    y_average_tal_len = float(len(y.seq1_seq) + len(y.seq2_seq)) / 2
    
    if x_average_tal_len == y_average_tal_len:
        return y if len(y.spacer_seq) < len(x.spacer_seq) else x
    else:
        return y if y_average_tal_len > x_average_tal_len else x

def filterStreubel(binding_site):
    
    seq2 = reverseComplement(binding_site.seq2_seq)
    
    if float(binding_site.seq1_seq.count('C') + binding_site.seq1_seq.count('G')) / len(binding_site.seq1_seq) < 0.25:
        return True
    
    if len(streubel_at_streak_re.findall(binding_site.seq1_seq)) > 0:
        return True
    
    if float(seq2.count('C') + seq2.count('G')) / len(seq2) < 0.25:
        return True
    
    if len(streubel_at_streak_re.findall(seq2)) > 0:
        return True
    
    return False

def validateOptions(options):
    
    if options.fasta == 'NA':
        raise TaskError("FASTA file required.")
    
    if options.cupstream not in [0, 1, 2]:
        raise TaskError("Invalid cupstream value provided")
        
    if options.filter == 1 and options.filterbase == -1:
        raise TaskError("Filter by cut site selected but no cut site was provided")
    
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
    
    if options.offtargets_ncbi != "NA":
        
        options.offtargets_ncbi = options.offtargets_ncbi.strip()
        
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
        
        if options.filter == 2:
            raise TaskError("Off-target counting is not allowed for unfiltered queries.")
    
    with open(options.fasta, 'r') as seq_file:
        
        check_fasta_pasta(seq_file)
        
        if options.check_offtargets:
            
            gene_length_total = 0
            
            for gene in FastaIterator(seq_file, alphabet=generic_dna):
                gene_length_total += len(gene.seq)
            
            if gene_length_total > 1000:
                raise TaskError("Off-target counting is only available when designing TALENs for sequences that are 1000 bases or less")

def RunFindTALTask(options):
    
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
        
        if options.outpath == 'NA':
            output_filepath = options.outdir + options.job + options.outfile
        else:
            output_filepath = options.outpath
    
        out = open(output_filepath, 'w')
        
        table_ignores = ["TAL1 length", "TAL2 length", "Spacer length"]
        
        out.write("table_ignores:" + ','.join(table_ignores) + "\n")
    
        strand_min = 15 if options.arraymin is None else options.arraymin
        strand_max = 20 if options.arraymax is None else options.arraymax
        
        spacer_min = 15 if options.min is None else options.min
        spacer_max = 30 if options.max is None else options.max
        
        u_bases = []
        
        if options.cupstream != 1:
            u_bases.append("T")
        
        if options.cupstream != 0:
            u_bases.append("C")
            
        out.write("options_used:" + ', '.join([
            "array_min = " + str(strand_min),
            "array_max = " + str(strand_max),
            "spacer_min = " + str(spacer_min),
            "spacer_max = " + str(spacer_max),
            "upstream_base = " + (" or ".join(u_bases))
        ]) + "\n")
        
        offtarget_header = "\tOff-Target Counts" if options.check_offtargets else ""
        
        out.write('Sequence Name\tCut Site\tTAL1 start\tTAL2 start\tTAL1 length\tTAL2 length\tSpacer length\tSpacer range\tTAL1 RVDs\tTAL2 RVDs\tPlus strand sequence\tUnique RE sites in spacer\t% RVDs HD or NN/NH' + offtarget_header + '\n')
        
        binding_sites = []
        
        for gene in FastaIterator(seq_file, alphabet=generic_dna):
            
            sequence = str(gene.seq).upper()
            
            site_entry_counts = {}
            
            if options.filter == 1:
                if options.filterbase > len(sequence):
                    logger("Skipped %s as the provided cut site was greater than the sequence length" % (gene.id))
                    continue
                cut_site_positions = [options.filterbase]
            else:
                cut_site_positions = range(len(sequence))
            
            logger("Scanning %s for binding sites" % (gene.id))
            
            for i in cut_site_positions:
                
                cut_site_potential_sites = []
                
                for spacer_size in range(spacer_min, spacer_max + 1):
                    
                    spacer_potential_sites = []
                    
                    spacer_size_left = int(math.floor(float(spacer_size) / 2))
                    spacer_size_right = int(math.ceil(float(spacer_size) / 2))
                    
                    if i < (strand_min + spacer_size_left + 1) or i > (len(sequence) - (strand_min + spacer_size_right) - 1):
                        continue
                    
                    for u_base in u_bases:
                        
                        if u_base == "T":
                            d_base = "A"
                        elif u_base == "C":
                            d_base = "G"
                            
                        u_pos_search_start = i - (strand_max + spacer_size_left) - 1
                        
                        if u_pos_search_start < 0:
                            u_pos_search_start = 0
                            
                        u_pos_search_end = i - (strand_min + spacer_size_left)
                        
                        d_pos_search_start = i + (strand_min + spacer_size_right)
                        d_pos_search_end = i + (strand_max + spacer_size_right) + 1
                        
                        u_positions = []
        
                        u_pos = 0
                        
                        while True:
                            
                            u_pos = sequence.rfind(u_base, u_pos_search_start, u_pos_search_end)
                            
                            if u_pos == -1:
                                break
                            else:
                                u_pos_search_end = u_pos
                                u_positions.append(u_pos)
                        
                        d_positions = []
                        
                        d_pos = 0
                        
                        while True:
                            
                            d_pos = sequence.find(d_base, d_pos_search_start, d_pos_search_end)
                            
                            if d_pos == -1:
                                break
                            else:
                                d_pos_search_start = d_pos + 1
                                d_positions.append(d_pos)
                        
                        break_out = False
                        
                        for u_pos in reversed(u_positions):
                            
                            for d_pos in reversed(d_positions):
                                
                                #uses inclusive start, exclusive end
                                tal1_start = u_pos + 1
                                tal1_end = i - spacer_size_left
                                tal1_seq = sequence[tal1_start : tal1_end]
                                tal2_start = i + spacer_size_right
                                tal2_end = d_pos
                                tal2_seq = sequence[tal2_start : tal2_end]
                                
                                if not ((tal1_seq in site_entry_counts and tal2_seq in site_entry_counts[tal1_seq]) or \
                                (tal1_seq in site_entry_counts and tal1_seq in site_entry_counts[tal1_seq]) or \
                                (tal2_seq in site_entry_counts and tal1_seq in site_entry_counts[tal2_seq]) or \
                                (tal2_seq in site_entry_counts and tal2_seq in site_entry_counts[tal2_seq])):
                                    
                                    bad_site = False
                                    
                                    cg_count = 0
                                    
                                    tal1_rvd = []
                                    
                                    for c in tal1_seq:
                                        
                                        if c not in strong_binding_RVDs:
                                            bad_site = True
                                            break
                                        
                                        if c == 'C' or c == 'G':
                                            cg_count += 1
                                        
                                        tal1_rvd.append(strong_binding_RVDs[c])
                                    
                                    if bad_site:
                                        continue
                                    
                                    tal1_rvd = ' '.join(tal1_rvd)
                                    
                                    tal2_rvd = []
                                    
                                    for c in reverseComplement(tal2_seq):
                                        
                                        if c not in strong_binding_RVDs:
                                            bad_site = True
                                            break
                                        
                                        if c == 'C' or c == 'G':
                                            cg_count += 1
                                            
                                        tal2_rvd.append(strong_binding_RVDs[c])
                                    
                                    if bad_site:
                                        continue
                                    
                                    tal2_rvd = ' '.join(tal2_rvd)
                                    
                                    if options.filter == 0:
                                        break_out = True
                                    
                                    binding_site = BindingSite(seq_id = gene.id,
                                                   cutsite = i,
                                                   seq1_start = tal1_start,
                                                   seq1_end = tal1_end,
                                                   seq1_seq = tal1_seq,
                                                   seq1_rvd = tal1_rvd,
                                                   spacer_start = tal1_end,
                                                   spacer_end = tal2_start,
                                                   spacer_seq = sequence[tal1_end : tal2_start],
                                                   seq2_start = tal2_start,
                                                   seq2_end = tal2_end,
                                                   seq2_seq = tal2_seq,
                                                   seq2_rvd = tal2_rvd,
                                                   upstream = u_base,
                                                   cg_percent = int(round(float(cg_count) / (len(tal1_seq) + len(tal2_seq)), 2) * 100))
                                    
                                    findRESitesInSpacer(sequence, binding_site)
                                    
                                    if binding_site.seq1_seq not in site_entry_counts:
                                        site_entry_counts[binding_site.seq1_seq] = {}
                                        
                                    if binding_site.seq2_seq not in site_entry_counts[tal1_seq]:
                                        site_entry_counts[binding_site.seq1_seq][binding_site.seq2_seq] = []
                                        
                                    site_entry_counts[binding_site.seq1_seq][binding_site.seq2_seq].append(binding_site)
                                    spacer_potential_sites.append(binding_site)
                                    
                                if break_out:
                                    break
                                
                            if break_out:
                                break
                    
                    if len(spacer_potential_sites) > 0:
                        if options.filter == 0:
                            cut_site_potential_sites.append(reduce(filterByTALSize, spacer_potential_sites))
                        else:
                            cut_site_potential_sites.extend(spacer_potential_sites)
                    
                
                if len(cut_site_potential_sites) > 0:
                    if options.filter == 0:
                        binding_sites.append(reduce(filterByTALSize, cut_site_potential_sites))
                    else:
                        binding_sites.extend(cut_site_potential_sites)
        
        
        if options.streubel:
            binding_sites[:] = list(ifilterfalse(filterStreubel, binding_sites))
        
        if options.check_offtargets:
            
            if len(binding_sites) > 0:
                
                off_target_pairs = []
                
                for i, binding_site in enumerate(binding_sites):
                    off_target_pairs.append([binding_site.seq1_rvd, binding_site.seq2_rvd])
                
                
                off_target_counts = PairedTargetFinderCountTask(offtarget_seq_filename, options.logFilepath, options.cupstream, 3.0, spacer_min, spacer_max, off_target_pairs)
                
                for i, binding_site in enumerate(binding_sites):
                    binding_site.offtarget_counts = off_target_counts[i]
        
        for i, binding_site in enumerate(binding_sites):
            
            output_items = [
                str(binding_site.seq_id),
                str(binding_site.cutsite),
                str(binding_site.seq1_start),
                str(binding_site.seq2_end - 1),
                str(binding_site.seq1_end - binding_site.seq1_start),
                str(binding_site.seq2_end - binding_site.seq2_start),
                str(binding_site.spacer_end - binding_site.spacer_start),
                str(binding_site.spacer_start) + '-' + str(binding_site.spacer_end - 1),
                binding_site.seq1_rvd,
                binding_site.seq2_rvd,
                binding_site.upstream + ' ' + binding_site.seq1_seq + ' ' + binding_site.spacer_seq.lower() + ' ' + binding_site.seq2_seq + ' ' + ("A" if binding_site.upstream == "T" else "G"),
                binding_site.re_sites,
                str(binding_site.cg_percent)
            ]
            
            if options.check_offtargets:
                output_items.append(' '.join(str(binding_site.offtarget_counts[x]) for x in range(5)))
            
            out.write("\t".join(output_items) + "\n")
        
        out.close()
        seq_file.close()
        
        logger('Finished')

if __name__ == '__main__':
    
    # import arguments and options
    usage = 'usage: %prog [options]'
    parser = OptParser(usage=usage)
    parser.add_option('--fasta', dest='fasta', type='string', default='NA', help='FASTA file containing promoter sequence(s).')
    parser.add_option('--min', dest='min', type='int', default=None, help='the minimum spacer size to try')
    parser.add_option('--max', dest='max', type='int', default=None, help='the maximum spacer size to try')
    parser.add_option('--arraymin', dest='arraymin', type='int', default=None, help='the minimum repeat array length to try')
    parser.add_option('--arraymax', dest='arraymax', type='int', default=None, help='the maximum repeat array length to try')
    parser.add_option('--cupstream', dest='cupstream', type='int', default = 0, help='1 to look for C instead of T, 2 to look for either')
    parser.add_option('--filter', dest='filter', type='int', default = 0, help='0 for smallest at each cut site, 1 for each everything targetting a specific site, 2 for unfiltered')
    parser.add_option('--filterbase', dest='filterbase', type='int', default = -1, help='if filter is 1 this gives the cutpos')
    parser.add_option('--gspec', dest='gspec', action='store_true', default = False, help='If true, use NH instead of NN for G')
    parser.add_option('--streubel', dest='streubel', action='store_true', default = False, help='If true, filter out TALENs that don\'t mean Streubel et. al. design guidelines')
    # Offtarget Options
    parser.add_option('--offtargets', dest='check_offtargets', action = 'store_true', default = False, help='Check offtargets')
    parser.add_option('--offtargets-fasta', dest='offtargets_fasta', type='string', default='NA', help='FASTA file containing to search for off-targets')
    parser.add_option('--offtargets-ncbi', dest='offtargets_ncbi', type='string', default='NA', help='NCBI nucleotide sequence ID to search for off-targets')
    parser.add_option('--genome', dest='genome', action = 'store_true', default = False, help='Input is a genome file')
    parser.add_option('--promoterome', dest='promoterome', action = 'store_true', default = False, help='Input is a promoterome file')
    parser.add_option('--organism', dest='organism', type = 'string', default='NA', help='Name of organism for the genome to be searched.')
    #Drupal options
    parser.add_option('--outpath', dest='outpath', type='string', default = 'NA', help='Optional full path for output file; if set --job, --outdir and --outfile are ignored.')
    parser.add_option('--logpath', dest='logFilepath', type='string', default = 'NA', help='Process log file path')
    parser.add_option('--nodeid', dest='nodeID', type='int', default = '-1', help='Optional node id if this script was called from Drupal.')
    parser.add_option('--ipaddr', dest='ip_address', type='string', default = '', help='IP address of job submitter')
    (options, args) = parser.parse_args()
    
    validate_options_handler(validateOptions, options)
    
    if options.nodeID != -1:
        
        # if run from drupal then it should be queued as a task
        
        if not celery_found:
            raise TaskError("nodeID option was provided but Celery is not installed")
        
        logger = create_logger(options.logFilepath)
        logger("Your task has been queued and will be processed when a worker node becomes available")
        
        # This ensures that the task is queued with the same module name
        # that the Celery workers are expecting 
        from findTAL import FindTALTask
        
        if options.check_offtargets:
            job_queue = "findtal_offtargets"
        else:
            job_queue = "findtal"
        
        FindTALTask.apply_async(kwargs=vars(options), queue=job_queue)
        
    else:
        
        RunFindTALTask(options)
