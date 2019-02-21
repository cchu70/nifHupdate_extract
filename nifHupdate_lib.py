#!/usr/bin/python
__author__="Claudia Chu"
__date__ ="2/18/19"


from os.path import abspath, join, isfile
from Bio import SeqIO
import subprocess
import time

#========================
# Defaults
def_prefix   = "TEST"
def_query    = "nifH"
def_tag      = "GENE"
def_start    = "2012"
def_end      = "2019"
def_datetype = "PDAT"
def_elements = "Id Caption TaxId Slen Organism Title CreateDate"


#========================
def splitbyLength(fastaFile):



#========================
def createShFile(configDict, startDate, endDate):
    query    = configDict["QUERY"]
    prefix   = configDict["PREFIX"]
    query    = configDict["QUERY"]
    datetype = configDict["DATETYPE"]

    esearch  = """esearch -db nucleotide -query \"%s\""""    % (query)
    efilter  = """efilter -mindate %s -maxdate %s -datetype %s""" % (startDate, endDate, datetype)
    efetch_1 = """efetch -format docsum"""
    xtract   = """xtract -pattern DocumentSummary -element %s"""  % (def_elements)
    filters  = """grep -e '%s|genome'""" % (query)

    edirectCmd = """%s | %s | %s | %s > %s_xtract.txt\n""" % (esearch, efilter, efetch_1, xtract, prefix)

    cat = "cat %s_xtract.txt" % (prefix)
    grepNifH = "grep 'nifH'"
    grepGenome = "grep 'genome'"
    acc = """awk 'BEGIN { ORS="," }; { print $2 }'"""

    fastaGeneCmd = """%s | %s | %s | python3 nifH_genome_seq_extractor.py > %s.gene.fasta\n""" % (cat, grepNifH, acc, prefix)
    fastaGenomeCmd = """%s | %s | %s | python3 nifH_genome_seq_extractor.py > %s.genome.fasta\n""" % (cat, grepGenome, acc, prefix)

    shFileName = "%s.%s.sh" % (prefix, startDate)
    sh = open(shFileName, "w")
    sh.write(edirectCmd)
    sh.write(fastaNifHCmd)
    sh.write(fastaGenomeCmd)
    sh.close()
    return shFileName

#========================
def parseConfig(configFile):

    if ( not isfile(configFile)):
        throwError("parse_config did not find %s" % ( configFile ) )

    configDict = {}
    for line in open(configFile, "r"):
        key, val = line.strip().split(None, 1)
        configDict[key] = val

    return configDict


#========================
def launch(shFileName):
    n = subprocess.Popen(shFileName)
    n.poll()


#========================
def parseDate(date):
    nums = date.split('/')

    if len(nums) == 3:
        # MM/DD/YYYY
        year = int(nums[2])
    elif len(nums) == 2:
        year = int(nums[1])
    else:
        year = int(nums[0])
    #####

    return year



