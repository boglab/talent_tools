from paired_talesf import ScorePairedTalesfTask

from talconfig import RVD_SEQ_REGEX, GENOME_FILE, PROMOTEROME_FILE

from talutil import OptParser, create_logger, OptionObject, TaskError

celery_found = True
try:
	from celery.task import task
	from celery.registry import tasks
	from talutil import BaseTask
except ImportError:
	celery_found = False

import re

if celery_found:
	@task(base=BaseTask)
	def PairedTalesfTask(*args, **kwargs):
		RunPairedTalesfTask(OptionObject(**kwargs))

def RunPairedTalesfTask(options):
	
	logger = create_logger(options.logFilepath)
	
	logger("Beginning")
	
	if options.cupstream not in [0, 1, 2]:
		raise TaskError("Invalid cupstream value provided")
	
	rvdString = options.rvdString.strip().upper()
	rvdString2 = options.rvdString2.strip().upper()
	
	RVD_re = re.compile(RVD_SEQ_REGEX, re.IGNORECASE | re.MULTILINE)
	if not RVD_re.match(rvdString):
		raise TaskError("RVD sequence is not in the correct format.  Enter between 12 and 31 RVDs using the standard single letter amino acid abbreviations.")
	
	valid_genome_organisms = ['drosophila_melanogaster', 'arabidopsis_thaliana', 'mus_musculus', 'oryza_sativa', 'caenorhabditis_elegans', 'danio_rerio', 'homo_sapiens']
	valid_promoterome_organisms = ['drosophila_melanogaster', 'arabidopsis_thaliana', 'mus_musculus', 'oryza_sativa', 'caenorhabditis_elegans', 'danio_rerio', 'homo_sapiens']
	
	if ((options.genome and options.organism not in valid_genome_organisms) or (options.promoterome and options.organism not in valid_promoterome_organisms)):
		raise TaskError("Invalid organism specified.")
	
	if options.genome:
		seqFilename = GENOME_FILE % options.organism
	elif options.promoterome:
		seqFilename = PROMOTEROME_FILE % options.organism
	else:
		seqFilename = options.fasta
	
	result = ScorePairedTalesfTask(seqFilename, rvdString, rvdString2, options.outputFilepath, options.logFilepath, forwardOnly, options.cupstream, options.cutoff, 1, options.organism if options.genome else "")
	
	if(result == 1):
		raise TaskError()

if __name__ == '__main__':
	
	# import arguments and options
	usage = 'usage: %prog [options]'
	parser = OptParser(usage=usage)
	# input options
	parser.add_option('-f', '--fasta', dest='fasta', type='string', default='NA', help='Path to input file if input is not a genome or promoterome')
	parser.add_option('-y', '--genome', dest='genome', action = 'store_true', default = False, help='Input is a genome file')
	parser.add_option('-x', '--promoterome', dest='promoterome', action = 'store_true', default = False, help='Input is a promoterome file')
	parser.add_option('-o', '--organism', dest='organism', type = 'string', default='NA', help='Name of organism for the genome to be searched.')
	# output options
	parser.add_option('-p', '--outpath', dest='outputFilepath', type='string', default = 'NA', help='Template file path for output file')
	parser.add_option('-l', '--logpath', dest='logFilepath', type='string', default = 'NA', help='Process log file path')
	parser.add_option('-z', '--nodeid', dest='nodeID', type='int', default = '-1', help='Drupal node ID')
	# program options
	parser.add_option('-u', '--cupstream', dest='cupstream', type='int', default = 0, help='1 to look for C instead of T, 2 to look for either')
	parser.add_option('-t', '--cutoff', dest='cutoff', type='float', default = 3.0, help='The threshold score that results must meet')
	parser.add_option('-c', '--revcomp', dest='revcomp', action = 'store_true', default = False, help='Search the reverse complement of the input FASTA sequences')
	parser.add_option('-r', '--rvds', dest='rvdString', type = 'string', default='NA', help='RVD sequence seperated by spaces or underscores.')
	(options, args) = parser.parse_args()
	
	if options.genome:
		queue_name = 'talesf_genome'
	elif options.promoterome:
		queue_name = 'talesf_promoterome'
	else:
		queue_name = 'talesf_other'
	
	if options.nodeID != -1:
		
		if not celery_found:
			raise TaskError("nodeID option was provided but Celery is not installed")
			
		logger = create_logger(options.logFilepath)
		logger("Your task has been queued and will be processed when a worker node becomes available")
		
		from findRvdTAL import PairedTalesfTask
		#if run from drupal then it should be queued as a task
		PairedTalesfTask.apply_async(kwargs=vars(options), queue=queue_name)
		
	else:
		
		RunPairedTalesfTask(options)
