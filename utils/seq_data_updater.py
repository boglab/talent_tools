#!/usr/bin/env python

import sys
import gzip
import os
import urllib2
import pickle
import subprocess
import warnings

from Bio import SeqIO

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
    "gasterosteus_aculeatus": "ftp://ftp.ensembl.org/pub/release-{0}/fasta/gasterosteus_aculeatus/dna/Gasterosteus_aculeatus.BROADS1.{0}.dna.toplevel.fa.gz",
    "mus_musculus": "ftp://ftp.ensembl.org/pub/release-{0}/fasta/mus_musculus/dna/Mus_musculus.GRCm38.{0}.dna.toplevel.fa.gz",
    "rattus_norvegicus": "ftp://ftp.ensembl.org/pub/release-{0}/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_5.0.{0}.dna.toplevel.fa.gz",
    "oryza_sativa": "ftp://ftp.ensemblgenomes.org/pub/plants/release-{0}/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.{0}.dna.toplevel.fa.gz",
    "arabidopsis_thaliana": "ftp://ftp.ensemblgenomes.org/pub/plants/release-{0}/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.{0}.dna.toplevel.fa.gz",
    "brachypodium_distachyon": "ftp://ftp.ensemblgenomes.org/pub/plants/release-{0}/fasta/brachypodium_distachyon/dna/Brachypodium_distachyon.v1.0.{0}.dna.toplevel.fa.gz",
    "solanum_lycopersicum": "ftp://ftp.ensemblgenomes.org/pub/plants/release-{0}/fasta/solanum_lycopersicum/dna/Solanum_lycopersicum.SL2.40.{0}.dna.toplevel.fa.gz",
}

print("Downloading genomes")

subprocess.check_call("mkdir -p %s" % (GENOME_DIR + "/gzip"), shell=True)

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
            
            subprocess.check_call("wget -O %s %s" % (gzipped_filepath, remote_file_path), shell=True)
            subprocess.check_call("gunzip -c %s > %s" % (gzipped_filepath, (destination_file_path)), shell=True)
            
        except urllib2.URLError:
            
            warnings.warn("Unable to download %s genome, remote file does not exist" % sequence_name)

with open(version_dump_filepath, "wb") as versions_file:
    pickle.dump(sequence_versions, versions_file)

#promoterome sequences

print("Downloading promoteromes")

subprocess.check_call("mkdir -p %s" % (PROMOTEROME_DIR + "/gzip"), shell=True)

promoterome_urls = {
    "homo_sapiens": "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/bigZips/upstream1000.fa.gz",
    "drosophila_melanogaster": "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Drosophila_melanogaster/bigZips/upstream1000.fa.gz",
    "caenorhabditis_elegans": "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Caenorhabditis_elegans/bigZips/upstream1000.fa.gz",
    "danio_rerio": "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Danio_rerio/bigZips/upstream1000.fa.gz",
    "mus_musculus": "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Mus_musculus/bigZips/upstream1000.fa.gz",
    "rattus_norvegicus": "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/Rattus_norvegicus/bigZips/upstream1000.fa.gz",
    "oryza_sativa": "http://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_1kb-upstream_2014-06-25.fasta.gz",
    "arabidopsis_thaliana": "ftp://ftp.arabidopsis.org/home/tair/Sequences/blast_datasets/TAIR10_blastsets/upstream_sequences/TAIR10_upstream_1000_translation_start_20101028",
    "brachypodium_distachyon": "ftp://brachypodium.org/brachypodium.org/Annotation/Bdistachyon.MIPS_1_2.promoter.1000.fa.gz",
    "gasterosteus_aculeatus": "ftp://hgdownload.cse.ucsc.edu/goldenPath/gasAcu1/bigZips/upstream1000.fa.gz",
}

for sequence_name, sequence_url in promoterome_urls.iteritems():
    
    print("PROMOTEROME: %s" % sequence_name)
    
    destination_file_path = PROMOTEROME_FILE % sequence_name
    
    try:
        
        #Check if the file exists
        remote_file = urllib2.urlopen(sequence_url)
        remote_file.close()
        
        if sequence_url[-2:] == "gz":
            
            gzipped_filepath = PROMOTEROME_DIR + "/gzip/" + sequence_name + '.fasta.gz'
            
            subprocess.check_call("wget -O %s %s" % (gzipped_filepath, sequence_url), shell=True)
            subprocess.check_call("gunzip -c %s > %s" % (gzipped_filepath, (destination_file_path)), shell=True)
            
        else:
            
            subprocess.check_call("wget -O %s %s" % (destination_file_path, sequence_url), shell=True)
        
    except urllib2.URLError:
        
        warnings.warn("Unable to download/update %s promoterome, remote file does not exist" % sequence_name)


#Generate RAPDB exclusive promoters

subprocess.check_call("wget -O %s %s" % ("/tmp/RAP-MSU.txt.gz", "http://rapdb.dna.affrc.go.jp/download/archive/RAP-MSU.txt.gz"), shell=True)
subprocess.check_call("gunzip -c %s > %s" % ("/tmp/RAP-MSU.txt.gz", "/tmp/RAP-MSU.txt"), shell=True)

with open("/tmp/RAP-MSU.txt", "r") as input_file:
    rapdb_exclusive_ids = set()
    for line in input_file:
        split_line = line.rstrip().split('\t')
        if len(split_line) == 2 and split_line[0] != "None" and split_line[1] == "None":
            rapdb_exclusive_ids.add(split_line[0])

keep_me = [seq for seq in SeqIO.parse(PROMOTEROME_FILE % "oryza_sativa", "fasta") if seq.id.split('-')[0].replace('t', 'g') in rapdb_exclusive_ids]

SeqIO.write(keep_me, PROMOTEROME_FILE % "oryza_sativa_rapdb_only", "fasta")
