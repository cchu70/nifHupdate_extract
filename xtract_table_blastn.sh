#!/bin/bash

blastn -query xtract_table_2018.gene.fa -db seqDatabase.fasta -outfmt '6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovhsp' -evalue 0.001 -out xtract_table_2018.gene.blastn.txt
blastn -query xtract_table_2018.genomes.fa -db seqDatabase.fasta -outfmt '6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovhsp' -evalue 0.001 -out xtract_table_2018.genomes.blastn.txt
python3 nifHupdate_launch.py nifHupdate_config.txt filter_alignments
