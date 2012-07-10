#!/usr/bin/python

# Import needed BioPython modules
from Bio import SeqIO
from Bio import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqUtils

from talconfig import BASE_DIR, GENOME_FILE, PROMOTEROME_FILE, VALID_GENOME_ORGANISMS, VALID_PROMOTEROME_ORGANISMS
from talutil import OptParser, FastaIterator, create_logger, check_fasta_pasta, OptionObject, TaskError, reverseComplement

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

from btfcount import TargetFinderCountTask

#Define a binding site object
class BindingSite:
	
	def __init__(self, **kwargs):
		
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
		
		self.offtarget_counts = kwargs.pop("offtarget_counts", [])
		
		self.re_sites = ""

#DNA list and dictionary
DNA = ['A', 'C', 'G', 'T']
DNA_dict = {'A':0, 'C':1, 'G':2, 'T':3}


#Average percent composition of known TAL binding sites
avg_percents = {'A':0.31, 'C':0.37, 'G':0.09, 'T':0.22}
stdev = {'A':0.16, 'C':0.13, 'G':0.08, 'T':0.10}

strong_binding_RVDs = {'A':'NI', 'C':'HD', 'G':'NN', 'T':'NG'}

with open(BASE_DIR + "/re_dict_dump", "rb") as re_dict_file:
	NEB_RE_sites = pickle.load(re_dict_file)

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
		return y if y_average_tal_len < x_average_tal_len else x

def filterByOfftargetCount(x, y):
	x_average_tal_len = float(len(x.seq1_seq) + len(x.seq2_seq)) / 2
	y_average_tal_len = float(len(y.seq1_seq) + len(y.seq2_seq)) / 2
	
	if x_average_tal_len == y_average_tal_len:
		return y if len(y.spacer_seq) < len(x.spacer_seq) else x
	else:
		return y if y_average_tal_len > x_average_tal_len else x

def RunFindTALTask(options):
	
	logger = create_logger(options.logFilepath)
	
	logger("Beginning")
	
	if options.fasta == 'NA':
		raise TaskError("FASTA file required.")
	
	if options.cupstream not in [0, 1, 2]:
		raise TaskError("Invalid cupstream value provided")
		
	if options.filter == 1 and options.filterbase == -1:
		raise TaskError("Filter by cut site selected but no cut site was provided")

	if options.check_offtargets:

		if ((options.genome and options.organism not in VALID_GENOME_ORGANISMS) or (options.promoterome and options.organism not in VALID_PROMOTEROME_ORGANISMS)):
			raise TaskError("Invalid organism specified.")
		
		if options.filter == 2:
			raise TaskError("Off-target counting is not allowed for unfiltered queries.")
		
		offtarget_seq_filename = ""
		
		if options.genome:
			offtarget_seq_filename = GENOME_FILE % options.organism
		elif options.promoterome:
			offtarget_seq_filename = PROMOTEROME_FILE % options.organism
		else:
			offtarget_seq_filename = options.fasta

	seq_file = open(options.fasta, 'r')
	
	#Prescreen for FASTA pasta
	check_fasta_pasta(seq_file)
	
	if options.outpath == 'NA':
		output_filepath = options.outdir + options.job + options.outfile
	else:
		output_filepath = options.outpath

	out = open(output_filepath, 'w')
	
	table_ignores = ["TAL1 length", "TAL2 length", "Spacer length"]
	
	out.write("table_ignores:" + ','.join(table_ignores) + "\n")

	strand_min = 15 if (options.arraymin is None or options.arraymin < 0) else options.arraymin
	strand_max = 20 if (options.arraymax is None or options.arraymax < 0) else options.arraymax
	
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
	
	out.write('Sequence Name\tCut Site\tTAL1 start\tTAL2 start\tTAL1 length\tTAL2 length\tSpacer length\tSpacer range\tTAL1 RVDs\tTAL2 RVDs\tPlus strand sequence\tUnique_RE_sites_in_spacer' + offtarget_header + '\n')
	
	found_something = False
	
	for gene in FastaIterator(seq_file, alphabet=generic_dna):
		
		sequence = str(gene.seq).upper()
		
		site_entry_counts = {}
		binding_sites = []
		
		if options.filter == 1:
			if options.filterbase > len(sequence):
				logger("Skipped %s as the provided cut site was greater than the sequence length" % (gene.id))
				continue
			cut_site_positions = [options.filterbase]
		else:
			cut_site_positions = range(len(sequence))
		
		logger("Scanning %s for binding sites" % (gene.id))
		
		off_target_pairs = []
		
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
								
								tal1_rvd = []
								
								for c in tal1_seq:
									
									if c not in strong_binding_RVDs:
										bad_site = True
										break
									
									tal1_rvd.append(strong_binding_RVDs[c])
									
								if bad_site:
									continue
								
								tal1_rvd = ' '.join(tal1_rvd)
								
								tal2_rvd = []
								
								for c in reverseComplement(tal2_seq):
									
									if c not in strong_binding_RVDs:
										bad_site = True
										break
									
									tal2_rvd.append(strong_binding_RVDs[c])
									
								if bad_site:
									continue
								
								tal2_rvd = ' '.join(tal2_rvd)
								
								found_something = True
								
								if options.filter == 0:
									break_out = True
								
								binding_site = BindingSite(cutsite = i,
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
											   upstream = u_base)
								
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
						if options.check_offtargets:
							cut_site_potential_sites.append(reduce(filterByOfftargetCount, spacer_potential_sites))
						else:
							cut_site_potential_sites.append(reduce(filterByTALSize, spacer_potential_sites))
					else:
						cut_site_potential_sites.extend(spacer_potential_sites)
				
			
			if len(cut_site_potential_sites) > 0:
				if options.filter == 0:
					if options.check_offtargets:
						binding_sites.append(reduce(filterByOfftargetCount, cut_site_potential_sites))
					else:
						binding_sites.append(reduce(filterByTALSize, cut_site_potential_sites))
				else:
					binding_sites.extend(cut_site_potential_sites)
		
		if options.check_offtargets:
			
			for i, binding_site in enumerate(binding_sites):
				off_target_pairs.append([binding_site.seq1_rvd, binding_site.seq2_rvd])
			
			off_target_counts = TargetFinderCountTask(offtarget_seq_filename, options.cupstream, 3.0, spacer_min, spacer_max, off_target_pairs)
			
			for i, binding_site in enumerate(binding_sites):
				binding_site.offtarget_counts = off_target_counts[i]
		
		for i, binding_site in enumerate(binding_sites):
			
			output_items = [
				gene.id,
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
			]
			
			if options.check_offtargets:
				output_items.append(' '.join(str(binding_site.offtarget_counts[x]) for x in range(5)))
			
			out.write("\t".join(output_items) + "\n")

	if not found_something:
		out.write('No TALEN pairs matching your criteria were found.')
	
	out.close()
	seq_file.close()

	logger('Finished')
	
if __name__ == '__main__':
	
	# import arguments and options
	usage = 'usage: %prog [options]'
	parser = OptParser(usage=usage)
	parser.add_option('-f', '--fasta', dest='fasta', type='string', default='NA', help='FASTA file containing promoter sequence(s).')
	parser.add_option('-m', '--min', dest='min', type='int', default=None, help='the minimum spacer size to try')
	parser.add_option('-x', '--max', dest='max', type='int', default=None, help='the maximum spacer size to try')
	parser.add_option('-r', '--arraymin', dest='arraymin', type='int', default=None, help='the minimum repeat array length to try')
	parser.add_option('-y', '--arraymax', dest='arraymax', type='int', default=None, help='the maximum repeat array length to try')
	parser.add_option('-u', '--cupstream', dest='cupstream', type='int', default = 0, help='1 to look for C instead of T, 2 to look for either')
	parser.add_option('-i', '--filter', dest='filter', type='int', default = 0, help='0 for smallest at each cut site, 1 for each everything targetting a specific site, 2 for unfiltered')
	parser.add_option('-b', '--filterbase', dest='filterbase', type='int', default = -1, help='if filter is 1 this gives the cutpos')
	# Offtarget Options
	parser.add_option('-e', '--offtargets', dest='check_offtargets', action = 'store_true', default = False, help='Check offtargets')
	parser.add_option('-v', '--genome', dest='genome', action = 'store_true', default = False, help='Input is a genome file')
	parser.add_option('-w', '--promoterome', dest='promoterome', action = 'store_true', default = False, help='Input is a promoterome file')
	parser.add_option('-s', '--organism', dest='organism', type = 'string', default='NA', help='Name of organism for the genome to be searched.')
	#Legacy
	parser.add_option('-j', '--job', dest='job', type='string', default='output', help='the job name, output files will have the job name as a prefix.')
	parser.add_option('-d', '--outdir', dest='outdir', type='string', default = 'upload/', help='Directory in which to place output files.')
	parser.add_option('-o', '--outfile', dest='outfile', type='string', default = '_TALEN_pairs_all.txt', help='Optional filename for output file.')
	parser.add_option('-t', '--t1', dest='t1', action='store_false', default = True, help='Do not enforce rule "No T at position 1."')
	parser.add_option('-a', '--a2', dest='a2', action='store_false', default = True, help='Do not enforce rule "No A at position 1."')
	parser.add_option('-n', '--tn', dest='tn', action='store_false', default = True, help='Do not enforce rule "Sites must end in a T."')
	parser.add_option('-g', '--gn', dest='gn', action='store_false', default = True, help='Do not enforce rule "Sites may not end in G/NN."')
	parser.add_option('-c', '--comp', dest='comp', action='store_false', default = True, help='Do not enforce base composition rules.')
	#Drupal options
	parser.add_option('-p', '--outpath', dest='outpath', type='string', default = 'NA', help='Optional full path for output file; if set --job, --outdir and --outfile are ignored.')
	parser.add_option('-l', '--logpath', dest='logFilepath', type='string', default = 'NA', help='Process log file path')
	parser.add_option('-z', '--nodeid', dest='nodeID', type='int', default = '-1', help='Optional node id if this script was called from Drupal.')
	(options, args) = parser.parse_args()

	
	if options.nodeID != -1:
		
		# if run from drupal then it should be queued as a task
		
		if not celery_found:
			raise TaskError("nodeID option was provided but Celery is not installed")
		
		logger = create_logger(options.logFilepath)
		logger("Your task has been queued and will be processed when a worker node becomes available")
		
		# This is ensures that the task is queued with the same module name
		# that the Celery workers are expecting 
		from findTAL import FindTALTask
		
		FindTALTask.apply_async(kwargs=vars(options), queue="findtal")
		
	else:
		
		RunFindTALTask(options)
