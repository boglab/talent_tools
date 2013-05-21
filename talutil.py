from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import single_letter_alphabet, generic_dna
from Bio.Seq import Seq

import string
import urllib
import sys
from datetime import datetime
from optparse import OptionParser

celery_found = True
try:
    from celery.task import task, Task
    from celery.registry import tasks
except ImportError:
    celery_found = False

from talconfig import DRUPAL_CALLBACK_URL
from entrez_cache import CachedEntrezFile

class TaskError(ValueError):
    pass

def task_on_success(kwargs):
    if "nodeID" in kwargs:
        urllib.urlopen(DRUPAL_CALLBACK_URL + str(kwargs['nodeID']) + "/0")

def task_on_error(exc, kwargs):
    logger = create_logger(kwargs['logFilepath'])
    
    errorType = type(exc).__name__
    
    if errorType == "TaskError":
        logger("%s: %s" % (errorType, str(exc)))
    else:
        logger("Internal error occurred. Please try again later.")
    
    if "nodeID" in kwargs:
        urllib.urlopen(DRUPAL_CALLBACK_URL + str(kwargs['nodeID']) + "/1")

if celery_found:
    class BaseTask(Task):
        
        abstract = True
        
        def on_success(self, retval, task_id, args, kwargs):
            task_on_success(kwargs)
        
        def on_failure(self, exc, task_id, args, kwargs, einfo):
            task_on_error(exc, kwargs)


def validate_options_handler(validator, options):
    try:
        validator(options)
    except TaskError as exc:
        task_on_error(exc, options.__dict__)
        raise

# Hack so that scripts can continue pretending they have an OptionParser object
# http://stackoverflow.com/questions/1305532/convert-python-dict-to-object
class OptionObject:
    def __init__(self, **entries): 
        self.__dict__.update(entries)

#FASTA Title Handler
def fasta_title_handler(line):
    split_line = line.split(None, 1)
    desc = ""
    if len(split_line) > 1:
        id = split_line[0]
        desc = split_line[1]
    elif (len(split_line) > 0):
        id = split_line[0]
    else:
        fasta_title_handler.noname += 1
        id = "unnamed_sequence_" + str(fasta_title_handler.noname)
    name = id
    return (id, name, desc)
    
fasta_title_handler.noname = 0

#Edited version of Bio.SeqIO.FastaIO.FastaIterator
def FastaIterator(handle, alphabet = single_letter_alphabet, title2ids = fasta_title_handler):
    #Skip any text before the first record (e.g. blank lines, comments)
    line = handle.readline()
    while line and line[0] != ">":
        line = handle.readline()

    while line: 
        if line[0]!=">":
            raise TaskError("Records in Fasta files should start with '>' character")
        if title2ids:
            id, name, descr = title2ids(line[1:].rstrip())
        else:
            descr = line[1:].rstrip()
            id   = descr.split()[0]
            name = id
        
        lines = []
        line = handle.readline()
        while line and line[0] != ">":
        #Remove trailing whitespace, and any internal spaces
        #(and any embedded \r which are possible in mangled files
        #when not opened in universal read lines mode)
            lines.append(line.rstrip().replace(" ","").replace("\r",""))
            line = handle.readline()
        
        #Return the record and then continue...
        yield SeqRecord(Seq(string.join(lines, ''), alphabet), id = id, name = name, description = descr)

def create_logger(logFilepath):
    
    if logFilepath != 'NA':
        def logger(message):
            with open(logFilepath, "a") as logFile:
                logFile.write("[%s] %s\n" % (datetime.now().ctime(), message))
        
    else:
        def logger(message):
            print "[%s] %s" % (datetime.now().ctime(), message)
    
    return logger

def check_fasta_pasta(seq_file):
    #This actually consumes lines from the start of the file
    #The position of the previous line is stored so that when a valid header is found
    #the pointer position can be set back so the FASTA iterator finds it too
    seq_file_last = seq_file.tell()
    seq_line = seq_file.readline()
    while seq_line and seq_line[0] != ">":
        seq_file_last = seq_file.tell()
        seq_line = seq_file.readline()
    
    if not seq_line:
        # Tell Drupal that this node errored out
        raise TaskError("The provided FASTA file did not contain any valid sequences. Sequences must be preceded by a header line with the format '> sequence_name'")
        
    seq_file.seek(seq_file_last)
    
    #if(seq_file.tell() != 0):
    #   logger("Warning: The first valid FASTA header in your input was not at the start of the sequence / file and any characters preceeding it were ignored. Please check that you didn't omit a header by mistake.")

def check_ncbi_sequence(seq_id):
    
    with CachedEntrezFile(seq_id) as ncbi_file:
        check_fasta_pasta(ncbi_file.file)

def check_ncbi_sequence_offtargets():
    
    with CachedEntrezFile(seq_id) as ncbi_file:
        
        check_fasta_pasta(ncbi_file.file)
        
        for record in FastaIterator(ncbi_file.file, alphabet=generic_dna):
            if len(record.seq) > 300000000:
                raise TaskError("Off-Target counting is only supported for NCBI records where all individual sequences are under 300 megabases in size")

class OptParser(OptionParser):
    def error(self, msg):
        """error(msg : string)
        
        Print a usage message incorporating 'msg' to stderr and exit.
        If you override this in a subclass, it should not return -- it
        should either exit or raise an exception.
        """
        #self.print_usage(sys.stderr)
        raise TaskError("%s\n" % (msg))
        #self.exit(2, )

def reverseComplement(sequence):
    
    new_sequence = ""
    
    for i in range(len(sequence)):
        base = sequence[len(sequence) - i - 1]
        if base == 'A':
            new_sequence += 'T'
        elif base == 'C':
            new_sequence += 'G'
        elif base == 'G':
            new_sequence += 'C'
        elif base == 'T':
            new_sequence += 'A'
        else:
            new_sequence += 'N'

    return new_sequence

#https://pypi.python.org/pypi/conditional
class Conditional(object):
    
    """Wrap another context manager and enter it only if condition is true."""

    def __init__(self, condition, contextmanager):
        
        self.condition = condition
        self.contextmanager = contextmanager

    def __enter__(self):
        
        if self.condition:
            return self.contextmanager.__enter__()
        else:
            return False

    def __exit__(self, *args):
        
        if self.condition:
            return self.contextmanager.__exit__(*args)
        else:
            return False
    