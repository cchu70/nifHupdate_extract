#!/bin/bash
# 2019-04-01 15:55:15.867424

echo Blasting TAXA_2017.nifH.fasta ...
blastn -query TAXA_2017.nifH.fasta -db TAXA_DB -outfmt '6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovhsp' -evalue 0.001 -out TAXA_2017.nifH.blastn.txt

echo Blasting TAXA_2017.genome.fasta ...
blastn -query TAXA_2017.genome.fasta -db TAXA_DB -outfmt '6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovhsp' -evalue 0.001 -out TAXA_2017.genome.blastn.txt

echo Blasting TAXA_2018.nifH.fasta ...
blastn -query TAXA_2018.nifH.fasta -db TAXA_DB -outfmt '6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovhsp' -evalue 0.001 -out TAXA_2018.nifH.blastn.txt

echo Blasting TAXA_2018.genome.fasta ...
blastn -query TAXA_2018.genome.fasta -db TAXA_DB -outfmt '6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovhsp' -evalue 0.001 -out TAXA_2018.genome.blastn.txt

python3 /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract/nifHupdate_launch.py /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract/nifHupdate_config.txt filter_best_alignments /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract

