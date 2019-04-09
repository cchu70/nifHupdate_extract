#!/bin/bash
# 2019-04-07 11:12:24.719591

echo Retrieving TAXA_2017.nifH.fasta ...
cat TAXA_2017.esearch.txt | grep 'nifH' | awk 'BEGIN { ORS="," }; { print $2 }' | python3 /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract/nifHupdate_fasta.py > TAXA_2017.nifH.fasta

echo Retrieving TAXA_2017.genome.fasta ...
cat TAXA_2017.esearch.txt | grep 'genome' | awk 'BEGIN { ORS="," }; { print $2 }' | python3 /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract/nifHupdate_fasta.py > TAXA_2017.genome.fasta

echo Retrieving TAXA_2018.nifH.fasta ...
cat TAXA_2018.esearch.txt | grep 'nifH' | awk 'BEGIN { ORS="," }; { print $2 }' | python3 /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract/nifHupdate_fasta.py > TAXA_2018.nifH.fasta

echo Retrieving TAXA_2018.genome.fasta ...
cat TAXA_2018.esearch.txt | grep 'genome' | awk 'BEGIN { ORS="," }; { print $2 }' | python3 /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract/nifHupdate_fasta.py > TAXA_2018.genome.fasta

python3 /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract/nifHupdate_launch.py /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract/./nifHupdate_config.txt set_db /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract

