#!/usr/bin/python

# Usage: findSingleTALsite.py [options]

# Options:
#	 -h, --help
#		 show this help message and exit
#	 -f FASTA, --fasta=FASTA
#		 FASTA file containing promoter sequence(s).
#	 -j JOB, --job=OUTPUT
#		 the job name, output files will have the job name as a prefix
#	 -r revcomp, --revcomp
#		 A flag indicating that the reverse compliment of the supplied FASTA sequence should be searched

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

from talutil import validate_options_handler, OptParser, FastaIterator, create_logger, check_fasta_pasta, OptionObject, TaskError

celery_found = True
try:
	from celery.task import task
	from celery.registry import tasks
	from talutil import BaseTask
except ImportError:
	celery_found = False

import math

import string

#import code
#import signal

#Define a binding site object
class Binding_site:
	def __init__(self, perfectTAL1, start1, seq1, is_plus, upstream):
		self.perfectTAL1 = perfectTAL1
		self.start1 = start1
		self.seq1 = seq1
		self.is_plus = is_plus
		self.upstream = upstream

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
	
	with open(options.fasta, 'r') as seq_file:
		check_fasta_pasta(seq_file)

def RunFindSingleTALSiteTask(options):

	logger = create_logger(options.logFilepath)
	
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
	
	out.write('Sequence Name\tTAL start\tTAL length\tRVD sequence\tStrand\tTarget sequence\tPlus strand sequence\n')
	if len(gene_binding_sites.keys()) > 0:
		for gene in sorted(gene_binding_sites.keys()):
			for start_site in gene_binding_sites[gene].keys():
				for binding_site in gene_binding_sites[gene][start_site]:
					if binding_site.is_plus:
						out.write(gene.id + '\t' + str(binding_site.start1) + '\t' + str(len(binding_site.seq1)) + '\t' + binding_site.perfectTAL1 + '\t' +  'Plus' + '\t' + binding_site.upstream + " " + str(binding_site.seq1) + '\t' + binding_site.upstream + " " + str(binding_site.seq1) + '\n')
					else:
						out.write(gene.id + '\t' + str(binding_site.start1) + '\t' + str(len(binding_site.seq1)) + '\t' + binding_site.perfectTAL1 + '\t' +  'Minus' + '\t' + ("T" if binding_site.upstream == "A" else "C") + " " + str(binding_site.seq1.reverse_complement()) + '\t' + str(binding_site.seq1) + " " + binding_site.upstream + '\n')
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
		# if run from drupal then it should be queued as a task
		FindSingleTALSiteTask.apply_async(kwargs=vars(options), queue="findsingletal")
		
	else:
		
		RunFindSingleTALSiteTask(options)
