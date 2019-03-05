#!/usr/bin/python
__author__="Claudia Chu"
__date__ ="2/18/19"


from os.path import abspath, join, isfile
from Bio import SeqIO
import subprocess
import time
from shutil import copyfile

#========================
# Defaults
def_prefix   = "TEST"
def_query    = "nifH"
def_tag      = "GENE"
def_start    = "2012"
def_end      = "2019"
def_datetype = "PDAT"
def_elements = "Id Caption TaxId Slen Organism Title CreateDate"

def_blastnOutfmt = '6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovhsp'
def_dbfiles = ["nhr", "nsd", "nin", "nsi", "nsq"]
def_evalue = 0.001

# #========================
# def splitbyLength(fastaFile):



#========================
def createShFile(cmdList, basePath, prefix, stage):

    shFileName = "%s/%s_%s.sh" % (basePath, prefix, stage)

    with open(shFileName, "w") as fh:

        # fh.write("# %s%s%s" % time.localtime(time.time()))
        # print(time.localtime(time.time()))
        fh.write("#!/bin/bash\n")
        fh.write("\n")
        for cmd in cmdList:
            fh.write(cmd)
        #####
    #####

    time.sleep(5)

    return shFileName

#========================
def esearchCmds(configDict, year):

    esearch  = """esearch -db nucleotide -query \"%s\""""    % (configDict["QUERY"])
    efilter  = """efilter -mindate %s -maxdate %s -datetype %s""" % (year, year + 1, configDict["DATETYPE"])
    efetch_1 = """efetch -format docsum"""
    xtract   = """xtract -pattern DocumentSummary -element %s"""  % (def_elements)
    edirectCmd = """%s | %s | %s | %s > %s_%s_xtract.txt\n""" % (esearch, efilter, efetch_1, xtract, configDict["PREFIX"], year)

    return edirectCmd

#========================
def fastaCmds(configDict, year, grepFilter):

    cat = "cat %s_%s_xtract.txt" % (configDict["PREFIX"], year)
    acc = """awk 'BEGIN { ORS="," }; { print $2 }'"""


    fastaCmd = """%s | grep '%s' | %s | python3 nifHupdate_fasta.py > %s_%s.%s.fasta\n""" % (cat, grepFilter, acc, configDict["PREFIX"], year, grepFilter)


    return fastaCmd


#========================
def blastnCmds(configDict, year, source, fmtString = def_blastnOutfmt):

    # fmt = str(outfmtVer) + " ".join([col for col in fmtList])

    fastFileName = "%s_%s.%s.fasta" % (configDict["PREFIX"], year, source)
    outputFile = "%s_%s.%s.blastn.txt" % (configDict["PREFIX"], year, source)

    blastCmd = """blastn -query %s -db %s -outfmt '%s' -evalue %s -out %s\n""" % (fastFileName, configDict["DBNAME"], def_blastnOutfmt, def_evalue, outputFile)
    return blastCmd

#========================
def verifyDb(dbLabel):
    checkFiles = [int(isfile("%s.%s" % (dbname, fileType))) for fileType in def_dbfiles]
    #check if all the database files exist
    if (sum(checkFiles) != len(def_dbfiles)):
        return False
    else:
        return True

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
    n = subprocess.Popen(["bash", shFileName])
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



