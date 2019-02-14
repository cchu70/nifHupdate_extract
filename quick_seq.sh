#!/bin/bash
# Author: Claudia Chu
# Date: 2/14/2019
#
# This is a simple script to run to retrieve every record from the
# NCBI database based on input query.
#
# Be wary that reaching the API limit on queries to the database may cause
# errors depending on the date range and number of entries that contains
# the query term
#
# Files produced:
#
# - xtract_table.txt  table of the records returned by NCBI
# - genome_seq.fasta  fasta file of whole genome records
# - gene_seq.fasta    fasta file of partial sequence records

# Query for nifH gene and pull out accession numbers and additional data
# You can change mindate, maxdate, but NOT the order "Id Caption"
esearch -db nucleotide -query "nifH" | efilter -mindate 05/18/2012 -maxdate 2019 -datetype PDAT | efetch -format docsum | xtract -pattern DocumentSummary -element Id Caption TaxId Organism Title CreateDate > xtract_table.txt

# Takes in xtract file from above command, pulls out the accession number (Caption) and pipes into the sequence extractor
cat xtract_table.txt | grep 'genome' | awk 'BEGIN { ORS="," }; { print $2 }' | python3 nifH_genome_seq_extractor.py > genome_seq.fasta
cat xtract_table.txt | grep 'nifH' | awk 'BEGIN { ORS="," }; { print $2 }' | python3 nifH_genome_seq_extractor.py > gene_seq.fasta