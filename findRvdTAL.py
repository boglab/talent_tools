from talesf import ScoreTalesfTask

from talconfig import RVD_SEQ_REGEX, GENOME_FILE, PROMOTEROME_FILE, VALID_GENOME_ORGANISMS, VALID_PROMOTEROME_ORGANISMS
from talutil import validate_options_handler, OptParser, create_logger, OptionObject, TaskError, check_fasta_pasta, Conditional
from entrez_cache import CachedEntrezFile

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
    def TalesfTask(*args, **kwargs):
        RunTalesfTask(OptionObject(**kwargs))

def validateOptions(options):
    
    if options.cupstream not in [0, 1, 2]:
        raise TaskError("Invalid cupstream value provided")
    
    if options.cutoff not in [3, 3.5, 4]:
        raise TaskError("Invalid cutoff value provided")
    
    RVD_re = re.compile(RVD_SEQ_REGEX, re.IGNORECASE | re.MULTILINE)
    
    if not RVD_re.match(options.rvdString):
        raise TaskError("RVD sequence is not in the correct format.  Enter between 12 and 31 RVDs using the standard single letter amino acid abbreviations.")
    
    if options.ncbi != "NA":
        
        if options.genome or options.promoterome:
            raise TaskError("--genome and --promoterome options cannot be combined with --ncbi")
        
        with CachedEntrezFile(options.ncbi) as ncbi_file:
            check_fasta_pasta(ncbi_file.file)
        
    else:
        
        if ((options.genome and options.organism not in VALID_GENOME_ORGANISMS) or (options.promoterome and options.organism not in VALID_PROMOTEROME_ORGANISMS)):
            raise TaskError("Invalid organism specified.")
        
        if not options.genome and not options.promoterome:
            with open(options.fasta, 'r') as seq_file:
                check_fasta_pasta(seq_file)

def RunTalesfTask(options):
    
    logger = create_logger(options.logFilepath)
    
    logger("Beginning")
    
    if options.revcomp:
        forwardOnly = False
    else:
        forwardOnly = True
    
    if options.ncbi != "NA":
        logger("Retrieving NCBI sequence. This could take a while if this sequence hasn't been used recently and needs to be downloaded from NCBI.")
    
    with Conditional(options.ncbi != "NA", CachedEntrezFile(options.ncbi)) as maybe_entrez_file:
        
        if options.ncbi != "NA":
            seqFilename = maybe_entrez_file.filepath
        elif options.genome:
            seqFilename = GENOME_FILE % options.organism
        elif options.promoterome:
            seqFilename = PROMOTEROME_FILE % options.organism
        else:
            seqFilename = options.fasta
        
        result = ScoreTalesfTask(seqFilename, options.rvdString, options.outputFilepath, options.logFilepath, forwardOnly, options.cupstream, options.cutoff, 4, options.organism if options.genome else "")
        
        if(result == 1):
            raise TaskError()

if __name__ == '__main__':
    
    # import arguments and options
    usage = 'usage: %prog [options]'
    parser = OptParser(usage=usage)
    # input options
    parser.add_option('-f', '--fasta', dest='fasta', type='string', default='NA', help='Path to input file if input is not a genome or promoterome')
    parser.add_option('--ncbi', dest='ncbi', type='string', default='NA', help='NCBI nucleotide sequence ID to search')
    parser.add_option('-y', '--genome', dest='genome', action = 'store_true', default = False, help='Input is a genome file')
    parser.add_option('-x', '--promoterome', dest='promoterome', action = 'store_true', default = False, help='Input is a promoterome file')
    parser.add_option('-o', '--organism', dest='organism', type = 'string', default='NA', help='Name of organism for the genome to be searched.')
    # program options
    parser.add_option('-u', '--cupstream', dest='cupstream', type='int', default = 0, help='1 to look for C instead of T, 2 to look for either')
    parser.add_option('-t', '--cutoff', dest='cutoff', type='float', default = 3.0, help='The threshold score that results must meet')
    parser.add_option('-c', '--revcomp', dest='revcomp', action = 'store_true', default = False, help='Search the reverse complement of the input FASTA sequences')
    parser.add_option('-r', '--rvds', dest='rvdString', type = 'string', default='NA', help='RVD sequence seperated by spaces or underscores.')
    #Drupal options
    parser.add_option('-p', '--outpath', dest='outputFilepath', type='string', default = 'NA', help='Template file path for output file')
    parser.add_option('-l', '--logpath', dest='logFilepath', type='string', default = 'NA', help='Process log file path')
    parser.add_option('-z', '--nodeid', dest='nodeID', type='int', default = '-1', help='Drupal node ID')
    parser.add_option('-k', '--ipaddr', dest='ip_address', type='string', default = '', help='IP address of job submitter')
    (options, args) = parser.parse_args()
    
    options.rvdString = options.rvdString.strip().upper()
    
    validate_options_handler(validateOptions, options)
    
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
        
        from findRvdTAL import TalesfTask
        #if run from drupal then it should be queued as a task
        TalesfTask.apply_async(kwargs=vars(options), queue=queue_name)
        
    else:
        
        RunTalesfTask(options)
