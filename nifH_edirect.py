#-db nucleotide -query "nifH" | efilter -mindate 2019 -maxdate 2019 -datetype PDAT | efetch -format docsum | xtract -pattern DocumentSummary -element Id Caption  CreateDate Organism Title Slen TaxId
#subprocess

# input: -batch mode -start date, -end date, -query, -table elements

#!/usr/bin/python
__author__="Claudia Chu"
__date__ ="1/29/19"


from os.path import abspath, join, isfile

from optparse import OptionParser

from sys import stdout, argv, stderr

from os import curdir, system

import subprocess

import time

import datetime

#==============================================================
# DEFAULTS
db_DEF         = "nucleotide"
gene_DEF       = "nifH"
term_DEF       = "nifH [GENE] NOT UNVERIFIED"
datetype_DEF   = "pdat"
startdate_DEF  = 2012
enddate_DEF    = datetime.date.today().year
# enddate_DEF    = 2012
batch_mode_DEF = "none"
elements_DEF   = "Id Caption CreateDate Organism Title Slen TaxId"
file_prefix    = "xtract_table"


# TIME
SLEEP_TIME = 20
QUERY_BUFFER_TIME = 20


#==============================================================
def set_parser():
    usage  = "usage: %prog [options]"

    parser = OptionParser(usage)

    parser.set_description("Scrape nifH from NCBI database using Entrez Direct, or edirect UNIX commands")

    (options, args) = parser.parse_args()

def query(year):

    file_name = "%s_%d" % (file_prefix, year)

    with open(file_name, "w") as f:
            min_date = year
            max_date = year

            # Just for extracting new nifH entries
            if (year == 2012):
                min_date = "05/18/2012" #AFTER or staring at 05/18/2012?
                max_date = "01/01/2013"
            #####

            edirect_cmd = 'esearch -db nucleotide -query \"%s\"' % term_DEF \
                        + ('| efilter -mindate %s -maxdate %s -datetype PDAT' % (min_date, max_date)) \
                        + '| efetch -format docsum | xtract -pattern DocumentSummary -element %s' % elements_DEF
            n = subprocess.Popen( edirect_cmd, shell=True, stdin = subprocess.PIPE, stdout=f)

            print("Querying %s for year %s " % (term_DEF, year))

            #limit = limit + 1
            time.sleep( SLEEP_TIME )
            n.kill()
        #####


#==============================================================
def sort_nifH(year, nf):
    print("Sorting year %s for nifH" % year)
    nifH_cmd = "cat %s_%s | grep 'nifH'" % (file_prefix, year)
    n = subprocess.Popen( nifH_cmd, shell=True, stdin = subprocess.PIPE, stdout=nf)
    time.sleep(2)
    n.kill()


#==============================================================
def sort_genomes(year, gn):
    print("Sorting year %s for genomes" % year)
    genome_cmd = "cat %s_%s | grep 'genome'" % (file_prefix, year)
    p = subprocess.Popen( genome_cmd, shell=True, stdin = subprocess.PIPE, stdout=gn)
    time.sleep(2)

    p.kill()


#==============================================================
def sort_others(year, ot):

    print("Sorting year %s for others" % year)
    others_cmd = "cat %s_%s | grep -v 'nifH' | grep -v 'genome'" % (file_prefix, year)
    k = subprocess.Popen( others_cmd, shell=True, stdin = subprocess.PIPE, stdout=ot)

    time.sleep(2)

    k.kill()


#==============================================================
def wait(subprocess):
    while subprocess.poll() == None:
        time.sleep(2)

#==============================================================
def fasta_genome(year):

# 1) take in the table extracted earlier
# 2) read the accession number
# 3) pull out the fasta sequence: efetch -db nuccore -id [Accession] -format gene_fasta | awk 'BEGIN {RS=">"}/nifH/{print ">"$0}'

    file_name = "%s_%s" % (file_prefix, year)
    fasta_name = "%s_%s.genomes.fa" % (file_prefix, year)

    cmd = """cat %s | grep 'genome' | awk 'BEGIN { ORS="," }; { print $2 }' | python3 nifH_genome_seq_extractor.py > %s""" % (file_name, fasta_name)

    print("Retrieving %s fasta" % year)

    n = subprocess.Popen(cmd, shell=True, stdin = subprocess.PIPE)

    wait(n)
    n.kill()

    print("Finished %s fasta" % year)

#==============================================================
def fasta(year):

# 1) take in the table extracted earlier
# 2) read the accession number
# 3) pull out the fasta sequence: efetch -db nuccore -id [Accession] -format gene_fasta | awk 'BEGIN {RS=">"}/nifH/{print ">"$0}'

    file_name = "%s_%s" % (file_prefix, year)
    fasta_name = "%s_%s.gene.fa" % (file_prefix, year)

    cmd = """cat %s | grep '%s' | awk 'BEGIN { ORS="," }; { print $2 }' | python3 nifH_genome_seq_extractor.py > %s""" % (file_name, gene_DEF, fasta_name)

    print("Retrieving %s %s fasta" % (gene_DEF, year))

    n = subprocess.Popen(cmd, shell=True, stdin = subprocess.PIPE)

    wait(n)
    n.kill()

    print("Finished %s %s fasta" % (gene_DEF, year))

#==============================================================
def accession(year):
    file_name = "%s_%s" % (file_prefix, year)
    acc_list = []
    for line in open(file_name, "r"):
        acc_list.append(line.strip())


#==============================================================
def real_main():
    set_parser()

    # for year in range(startdate_DEF, enddate_DEF + 1):
    #     query(year)
    # #####


    # nifH_file = "%s.nifH_file" % file_prefix
    # genomes_file ="%s.genomes_file" % file_prefix
    # others_file = "%s.others_file" % file_prefix

    # with open(nifH_file, 'a') as nf, open(genomes_file, 'a') as gn, open(others_file, 'a') as ot:
    #     for year in range(startdate_DEF, enddate_DEF + 1):
    #         sort_nifH(year, nf)
    #         sort_genomes(year, gn)
    #         sort_others(year, ot)
    #     #####
    # #####



    for year in range(startdate_DEF, enddate_DEF + 1):
        fasta(year)
    #####

#==============================================================
if ( __name__ == '__main__' ):
    real_main()








