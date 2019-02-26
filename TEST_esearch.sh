#!/bin/bash

esearch -db nucleotide -query "nifH [GENE] NOT UNVERIFIED" | efilter -mindate 2018 -maxdate 2019 -datetype PDAT | efetch -format docsum | xtract -pattern DocumentSummary -element Id Caption TaxId Slen Organism Title CreateDate > TEST_2018_xtract.txt
python3 nifHupdate_launch.py /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract/nifHupdate_config.txt fasta