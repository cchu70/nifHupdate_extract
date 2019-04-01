#!/bin/bash
# 2019-03-29 12:43:57.766726

echo esearch for year 2017...
esearch -db nucleotide -query "nifH [GENE] NOT UNVERIFIED" | efilter -mindate 2017 -maxdate 2017 -datetype PDAT | efetch -format docsum | xtract -pattern DocumentSummary -element Id Caption TaxId Slen Organism Title CreateDate > TAXA_2017.esearch.txt


echo esearch for year 2018...
esearch -db nucleotide -query "nifH [GENE] NOT UNVERIFIED" | efilter -mindate 2018 -maxdate 2018 -datetype PDAT | efetch -format docsum | xtract -pattern DocumentSummary -element Id Caption TaxId Slen Organism Title CreateDate > TAXA_2018.esearch.txt


echo Launching next stage fasta
python3 /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract/nifHupdate_launch.py /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract/nifHupdate_config.txt fasta /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract

