#!/usr/bin/env python

import sys
import gzip
import os
import urllib2
import pickle
import subprocess
import warnings

sys.path.append('/opt/boglab')

from talent.talconfig import BASE_DIR, GENOME_DIR, GENOME_FILE, PROMOTEROME_DIR, PROMOTEROME_FILE

version_dump_filepath = BASE_DIR + "/talent/utils/sequence_versions_dump"

with open(version_dump_filepath, "rb") as versions_file:
    sequence_versions = pickle.load(versions_file)

#genome sequences

genome_urls = {
    "homo_sapiens": "ftp://ftp.ensembl.org/pub/release-{0}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.{0}.dna.toplevel.fa.gz",
    "drosophila_melanogaster": "ftp://ftp.ensembl.org/pub/release-{0}/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP5.{0}.dna.toplevel.fa.gz",
    "caenorhabditis_elegans": "ftp://ftp.ensembl.org/pub/release-{0}/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel215.{0}.dna.toplevel.fa.gz",
    "danio_rerio": "ftp://ftp.ensembl.org/pub/release-{0}/fasta/danio_rerio/dna/Danio_rerio.Zv9.{0}.dna.toplevel.fa.gz",
    "mus_musculus": "ftp://ftp.ensembl.org/pub/release-{0}/fasta/mus_musculus/dna/Mus_musculus.GRCm38.{0}.dna.toplevel.fa.gz",
    "rattus_norvegicus": "ftp://ftp.ensembl.org/pub/release-{0}/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_5.0.{0}.dna.toplevel.fa.gz",
    "oryza_sativa": "ftp://ftp.ensemblgenomes.org/pub/plants/release-{0}/fasta/oryza_sativa/dna/Oryza_sativa.MSU6.{0}.dna.toplevel.fa.gz",
    "arabidopsis_thaliana": "ftp://ftp.ensemblgenomes.org/pub/plants/release-{0}/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.{0}.dna.toplevel.fa.gz",
    "brachypodium_distachyon": "ftp://ftp.ensemblgenomes.org/pub/plants/release-{0}/fasta/brachypodium_distachyon/dna/Brachypodium_distachyon.v1.0.{0}.dna.toplevel.fa.gz",
    "solanum_lycopersicum": "ftp://ftp.ensemblgenomes.org/pub/plants/release-{0}/fasta/solanum_lycopersicum/dna/Solanum_lycopersicum.SL2.40.{0}.dna.toplevel.fa.gz"
}

print("Downloading genomes")

subprocess.check_call("mkdir -p %s" % (GENOME_DIR + "/gzip"))

for sequence_name, sequence_url in genome_urls.iteritems():

    print("GENOME: %s" % sequence_name)
    
    destination_file_path = GENOME_FILE % sequence_name
    
    if os.path.isfile(destination_file_path):
        
        # Check for new version of existing
        
        try:
            
            updated_remote_file_path = sequence_url.format(sequence_versions["genomes"][sequence_name] + 1)
            
            #Check if the file exists
            remote_file = urllib2.urlopen(updated_remote_file_path)
            remote_file.close()
            
            gzipped_filepath = GENOME_DIR + "/gzip/" + sequence_name + '.fasta.gz'
            
            subprocess.check_call("wget -O %s %s" % (gzipped_filepath, updated_remote_file_path), shell=True)
            subprocess.check_call("gunzip -c %s > %s" % (gzipped_filepath, (destination_file_path)), shell=True)
            
            sequence_versions["genomes"][sequence_name] += 1
            
        except urllib2.URLError:
            
            pass
        
    else:
        
        # Download missing genome
        
        try:
            
            remote_file_path = sequence_url.format(sequence_versions["genomes"][sequence_name])
            
            #Check if the file exists
            remote_file = urllib2.urlopen(remote_file_path)
            remote_file.close()
            
            gzipped_filepath = GENOME_DIR + "/gzip/" + sequence_name + '.fasta.gz'
            
            subprocess.check_call("wget -O %s %s" % (gzipped_filepath, updated_remote_file_path), shell=True)
            subprocess.check_call("gunzip -c %s > %s" % (gzipped_filepath, (destination_file_path)), shell=True)
            
        except urllib2.URLError:
            
            warnings.warn("Unable to download %s genome, remote file does not exist" % sequence_name)

with open(version_dump_filepath, "wb") as versions_file:
    pickle.dump(sequence_versions, versions_file)

#promoterome sequences

print("Downloading promoteromes")

subprocess.check_call("mkdir -p %s" % (PROMOTEROME_DIR + "/gzip"))

promoterome_urls = {
    "homo_sapiens": "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/bigZips/upstream1000.fa.gz",
    "drosophila_melanogaster": "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Drosophila_melanogaster/bigZips/upstream1000.fa.gz",
    "caenorhabditis_elegans": "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Caenorhabditis_elegans/bigZips/upstream1000.fa.gz",
    "danio_rerio": "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Danio_rerio/bigZips/upstream1000.fa.gz",
    "mus_musculus": "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Mus_musculus/bigZips/upstream1000.fa.gz",
    "rattus_norvegicus": "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Rattus_norvegicus/bigZips/upstream1000.fa.gz",
}

for sequence_name, sequence_url in promoterome_urls.iteritems():
    
    print("PROMOTEROME: %s" % sequence_name)
    
    destination_file_path = PROMOTEROME_FILE % sequence_name
    
    try:
        
        #Check if the file exists
        remote_file = urllib2.urlopen(sequence_url)
        remote_file.close()
        
        gzipped_filepath = PROMOTEROME_DIR + "/gzip/" + sequence_name + '.fasta.gz'
        
        subprocess.check_call("wget -O %s %s" % (gzipped_filepath, sequence_url), shell=True)
        subprocess.check_call("gunzip -c %s > %s" % (gzipped_filepath, (destination_file_path)), shell=True)
        
    except urllib2.URLError:
        
        warnings.warn("Unable to download/update %s promoterome, remote file does not exist" % sequence_name)
