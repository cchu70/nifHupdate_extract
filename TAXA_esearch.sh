#!/bin/bash
# 2019-04-11 11:48:59.595149

echo esearch for year 2017...
esearch -db nucleotide -query "nifH [GENE] NOT UNVERIFIED" | efilter -mindate 2017 -maxdate 2017 -datetype PDAT | efetch -format docsum | xtract -pattern DocumentSummary -element Id Caption TaxId Slen CreateDate Organism Title > TAXA_2017.esearch.txt


echo Launching next stage fasta
python3 /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract/nifHupdate_Lib/nifHupdate_launch.py /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract/nifHupdate_config.txt fasta /Users/ClaudiaChu/Documents/Georgia_Tech/GitHub/nifHupdate/nifH_extract log.txt

