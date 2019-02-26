#!/bin/bash

cat TEST_2018_xtract.txt | grep 'nifH' | awk 'BEGIN { ORS="," }; { print $2 }' | python3 nifHupdate_fasta.py > TEST_2018.nifH.fasta
cat TEST_2018_xtract.txt | grep 'genome' | awk 'BEGIN { ORS="," }; { print $2 }' | python3 nifHupdate_fasta.py > TEST_2018.genome.fasta
python3 nifHupdate_launch.py /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract/nifHupdate_config.txt lenCluster