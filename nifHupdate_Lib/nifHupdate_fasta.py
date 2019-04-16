#!/usr/bin/python
__author__="Claudia Chu"
__date__ ="1/29/19"


from os.path import abspath, join, isfile

from optparse import OptionParser

from nifHupdate_lib import wait

from sys import stdout, argv, stderr, stdin

from os import curdir, system

import subprocess

import time

import datetime


#==============================================================
def fasta(accession):

# 1) take in the table extracted earlier
# 2) read the accession number
# 3) pull out the fasta sequence: efetch -db nuccore -id [Accession] -format gene_fasta | awk 'BEGIN {RS=">"}/nifH/{print ">"$0}'


    # extract_cmd = """efetch -db nuccore -id %s -format gene_fasta | awk 'BEGIN {RS=">"}/nifH/{print ">"$0}'""" % (accession)
    extract_cmd = """efetch -db nuccore -id %s -format gene_fasta | awk 'BEGIN {RS=">"}/nifH/{print ">"$0}""" % (accession)
    k = subprocess.Popen( extract_cmd, shell=True, stdin = subprocess.PIPE)

    wait(k)
    k.kill()

#==============================================================
def real_main():

    # take in comma delimeted accession numbers
    fasta(stdin.read()[:-1]) # remove last comma


#==============================================================
if ( __name__ == '__main__' ):
    real_main()
