#!/bin/bash

blastn -query TEST_2018.nifH.fasta -db seqDatabase.fasta -outfmt '6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovhsp' -out TEST_2018.nifH.blastn.txt
python3 nifHupdate_launch.py nifHupdate_config.txt end
