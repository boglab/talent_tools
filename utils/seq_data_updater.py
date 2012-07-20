#!/usr/bin/env python

from __future__ import with_statement

import sys
import gzip
import os
import urllib2

seq_dir = "seq_data"
gzip_dir = "seq_data/gzip"

gzip_files = {
    "homo_sapiens": "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/upstream1000.fa.gz",
    "drosophila_melanogaster": "http://hgdownload.cse.ucsc.edu/goldenPath/dm3/bigZips/upstream1000.fa.gz",
    "caenorhabditis_elegans": "http://hgdownload.cse.ucsc.edu/goldenPath/ce6/bigZips/upstream1000.fa.gz",
    "danio_rerio": "http://hgdownload.cse.ucsc.edu/goldenPath/danRer7/bigZips/upstream1000.fa.gz",
    "mus_musculus": "http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/upstream1000.fa.gz",
    "oryza_sativa": "ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_6.1/all.dir/all.1kUpstream.gz",
}

fasta_files = {
    #date is in the filename, can't check for updates
    #"arabidopsis_thaliana": "ftp://ftp.arabidopsis.org/home/tair/Sequences/blast_datasets/TAIR10_blastsets/upstream_sequences/TAIR10_upstream_1000_translation_start_20101028"
}

#process gzipped files
for sequence_name, sequence_url in gzip_files.iteritems():
    
    remote_file = urllib2.urlopen(sequence_url)
    gzipped_filepath = os.path.join(gzip_dir, sequence_name + '.fasta.gz')
    
    #gzip requires the entire file to be present in order to decompress
    with open(gzipped_filepath, 'wb') as gzip_file:
        gzip_file.writelines(remote_file)
    
    #the gzipped file uses 'try' / 'finally' because in python 2.6 gzip files don't have 'with' support
    ungzipped_file = gzip.GzipFile(gzipped_filepath, 'rb')
    try:
        with open(os.path.join(seq_dir, sequence_name + '.fasta'), 'wb') as fasta_file:
            fasta_file.writelines(ungzipped_file)
    finally:
        ungzipped_file.close()

#process normal files
for sequence_name, sequence_url in fasta_files.iteritems():
    
    remote_file = urllib2.urlopen(sequence_url)
    
    with open(os.path.join(seq_dir, sequence_name + '.fasta'), 'wb') as fasta_file:
       fasta_file.writelines(remote_file)
