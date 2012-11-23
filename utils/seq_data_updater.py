#!/usr/bin/env python

from __future__ import with_statement

import sys
import gzip
import os
import urllib2

sys.path.append('/opt/boglab')

from talent.talconfig import BASE_DIR, GENOME_DIR, GENOME_FILE, PROMOTEROME_DIR, PROMOTEROME_FILE

version_dump_filepath = BASE_DIR + "/utils/sequence_versions_dump"

with open(version_dump_filepath, "rb") as versions_file:
    sequence_versions = pickle.load(versions_file)

#genome sequences

genome_urls = {
    "homo_sapiens": "ftp://ftp.ensembl.org/pub/release-{0}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.{0}.dna.toplevel.fa.gz",
    "drosophila_melanogaster": "ftp://ftp.ensembl.org/pub/release-{0}/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP5.{0}.dna.toplevel.fa.gz",
    "caenorhabditis_elegans": "ftp://ftp.ensembl.org/pub/release-{0}/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel215.{0}.dna.toplevel.fa.gz",
    "danio_rerio": "ftp://ftp.ensembl.org/pub/release-{0}/fasta/danio_rerio/dna/Danio_rerio.Zv9.{0}.dna.toplevel.fa.gz",
    "mus_musculus": "ftp://ftp.ensembl.org/pub/release-{0}/fasta/mus_musculus/dna/Mus_musculus.GRCm38.{0}.dna.toplevel.fa.gz",
    "oryza_sativa": "ftp://ftp.ensemblgenomes.org/pub/plants/release-{0}/fasta/oryza_sativa/dna/Oryza_sativa.MSU6.{0}.dna.toplevel.fa.gz",
    "arabidopsis_thaliana": "ftp://ftp.ensemblgenomes.org/pub/plants/release-{0}/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.{0}.dna.toplevel.fa.gz",
}

for sequence_name, sequence_url in genome_urls.iteritems():
    
    try:
        
        remote_file = urllib2.urlopen(sequence_url.format(sequence_versions["genomes"][sequence_name] + 1))
        
        gzipped_filepath = GENOME_DIR + "/gzip/" + sequence_name + '.fasta.gz'
        
        #gzip requires the entire file to be present in order to decompress
        with open(gzipped_filepath, 'wb') as gzip_file:
            gzip_file.writelines(remote_file)
        
        #the gzipped file uses 'try' / 'finally' because in python 2.6 gzip files don't have 'with' support
        ungzipped_file = gzip.GzipFile(gzipped_filepath, 'rb')
        
        try:
            with open(GENOME_FILE % sequence_name, 'wb') as fasta_file:
                fasta_file.writelines(ungzipped_file)
                sequence_versions["genomes"][sequence_name] += 1
        finally:
            ungzipped_file.close()
        
    except urllib2.URLError:
        pass

with open(version_dump_filepath, "wb") as versions_file:
    pickle.dump(sequence_versions, versions_file)

#promoterome sequences

promoterome_urls = {
    "homo_sapiens": "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/bigZips/upstream1000.fa.gz",
    "drosophila_melanogaster": "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Drosophila_melanogaster/bigZips/upstream1000.fa.gz",
    "caenorhabditis_elegans": "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Caenorhabditis_elegans/bigZips/upstream1000.fa.gz",
    "danio_rerio": "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Danio_rerio/bigZips/upstream1000.fa.gz",
    "mus_musculus": "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Mus_musculus/bigZips/upstream1000.fa.gz",
}

for sequence_name, sequence_url in promoterome_urls.iteritems():
    
    try:
        
        remote_file = urllib2.urlopen(sequence_url)
        
        gzipped_filepath = PROMOTEROME_FILE + "/gzip/" + sequence_name + '.fasta.gz'
        
        #gzip requires the entire file to be present in order to decompress
        with open(gzipped_filepath, 'wb') as gzip_file:
            gzip_file.writelines(remote_file)
        
        #the gzipped file uses 'try' / 'finally' because in python 2.6 gzip files don't have 'with' support
        ungzipped_file = gzip.GzipFile(gzipped_filepath, 'rb')
        
        try:
            with open(PROMOTEROME_FILE % sequence_name, 'wb') as fasta_file:
                fasta_file.writelines(ungzipped_file)
        finally:
            ungzipped_file.close()
        
    except urllib2.URLError:
        pass
