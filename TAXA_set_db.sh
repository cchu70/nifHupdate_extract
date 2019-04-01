#!/bin/bash
# 2019-04-01 15:41:54.681575

makeblastdb -in /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract/seqDatabase.fasta -parse_seqids -dbtype nucl -out /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract/TAXA/TAXA_DB
python3 /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract/nifHupdate_launch.py /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract/nifHupdate_config.txt blastn /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract
