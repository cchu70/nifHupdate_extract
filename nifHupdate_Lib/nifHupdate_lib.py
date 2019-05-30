#!/usr/bin/python
__author__="Claudia Chu"
__date__ ="2/18/19"


from os.path import abspath, join, isfile
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import subprocess
import time
from shutil import copyfile
import datetime
from sys import argv, stderr, exit

#========================
# Defaults

def_blastnOutfmt = '6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovhsp sseq'
def_species_pident = 91.9
def_genus_pident = 88.1
def_family_pident = 75
def_dbfiles = ["nhr", "nsd", "nin", "nsi", "nsq"]
def_dbname = "blastnDB"
def_evalue = 0.001
def_minimap_align_len_cutoff = 200
def_blastn_align_len_cutoff = 200

# #========================
# Allowed sets
MINIMAP_LABELS = set(["PREFIX" , "DBFILE", "NUCCORE", "MIN_MINIMAP_ALIGNLEN", "MIN_BLASTN_ALIGNLEN", "PIDENT_CUTOFF"])

minimap_stages = set(['minimap', 'minimap_fasta', 'blastn', 'filter_best_alignments', 'trim_seq', 'cluster', 'deduplicate', 'end', 'rehead_fasta'])

def_validFastaEndings = set(["fa", "fasta", "fna"])
def_expectedFastaFofns = ["minimap_fasta.fofn", "fasta.trim.fofn", "cluster_fasta.fofn", "cluster_fasta_dedup.fofn"]


#========================
def wait(subprocess):
    while subprocess.poll() == None:
        time.sleep(2)

#========================
def createShFile(cmdList, basePath, prefix, stage):

    shFileName = "%s/%s_%s.sh" % (basePath, prefix, stage)

    with open(shFileName, "w") as fh:

        now = datetime.datetime.now()
        fh.write("#!/bin/bash\n")
        fh.write("# %s\n" % str(now))
        fh.write("\n")
        for cmd in cmdList:
            fh.write(cmd + "\n")
            # fh.write("echo " + cmd + "\n")
        #####
    #####

    time.sleep(2)

    return shFileName

#========================

def extractFileName(filePath):
    path = filePath.split("/")[: -1]
    fileName = filePath.split("/")[-1]
    fileNameHead = fileName.split(".")[0]
    return path, fileName, fileNameHead

#========================
def throwError(errorMessage, fh):
    fh.write("\n-----------------\n")
    fh.write("Time: %s\n" % datetime.datetime.now())
    fh.write("%s\n" % errorMessage)
    fh.write("Halting execution\n")
    fh.write("-----------------\n")
    stderr.write("%s\n" % errorMessage)
    assert False


#========================
def parseConfig(configFile, basePath, logFileFh):

    if ( not isfile(configFile)):
        throwError("parse_config did not find config file %s" % ( configFile ) , logFileFh)
    #####

    configDict = {}
    for line in open(configFile, "r"):
        if (line[0] != "#" and line != "\n"):

            # Ignore comments and newlines
            key, val = line.strip().split(None, 1)
            if (key not in MINIMAP_LABELS):
                throwError("%s is not a valid label. Remove it from %s" % (key, configFile), logFileFh)
            else:
                configDict[key] = val

    # check config errors and put in defaults
    # -----------------------------
    # Test prefix
    try:
        prefix = configDict["PREFIX"]
    except:
        throwError("No PREFIX provided in the configuration file %s" % configFile.split("/")[-1], logFileFh)
    #####

     # -----------------------------
    try:
        dbFile = configDict["DBFILE"]
        if (not isfile(dbFile)):
            throwError("Could not find file %s" % dbFile, logFileFh)

    except KeyError:
        throwError("No database fasta file provided in configuration file %s" % configFile.split("/")[-1], logFileFh)
    # -----------------------------

    try:
        blastnDBName = configDict["DBNAME"]
    except:
        configDict["DBNAME"] = def_dbname
    #####

    # -----------------------------
    try:
        nuccore = configDict["NUCCORE"]
        if (not isfile(nuccore)):
            throwError("Could not find database fasta file %s. Make sure you provide the \
                full path to the fofn file or a single fasta file or your database you are searching" % dbFile, logFileFh)
       #####


    except KeyError:
        throwError("No database fasta file provided in configuration file %s" % configFile.split("/")[-1], logFileFh)
    #####

    # Parameters
    # -----------------------------
    try:
        x = configDict["MIN_MINIMAP_ALIGNLEN"]
        configDict["MIN_MINIMAP_ALIGNLEN"] = int(x)
    except KeyError:
        configDict["MIN_MINIMAP_ALIGNLEN"] = def_minimap_align_len_cutoff
    # -----------------------------
    try:
        x = configDict["MIN_BLASTN_ALIGNLEN"]
        configDict["MIN_BLASTN_ALIGNLEN"] = int(x)
    except:
        configDict["MIN_BLASTN_ALIGNLEN"] = def_blastn_align_len_cutoff
    # -----------------------------
    try:
        x = configDict["PIDENT_CUTOFF"]
        configDict["PIDENT_CUTOFF"] = float(x)
    except KeyError:
        configDict["PIDENT_CUTOFF"] = def_family_pident
    #####

    return configDict


#========================
def launch(shFileName):
    print("\nRunning %s\n" % shFileName)
    # n = subprocess.Popen(["bash", "%s >> %s" % (shFileName, outFileName)])
    n = subprocess.Popen(["bash", shFileName])
    wait(n) # wait for the shfile to finish running before proceeding
    n.poll()

#========================
def isFastaFile(fastaFileName):

    fileEnding = fastaFileName.split(".")[-1] # Get the file ending
    isFasta = False
    for ending in def_validFastaEndings:
        if ending == fileEnding:
            isFasta = True
            break
    #####

    return isFasta


####################################################
############# STAGE SPECIFIC FUNCTIONS #############
####################################################



#========================
# blastn stage


def verifyDb(dbLabel):
    checkFiles = [int(isfile("%s.%s" % (dbLabel, fileType))) for fileType in def_dbfiles]
    #check if all the database files exist
    if (sum(checkFiles) != len(def_dbfiles)):
        return False
    else:
        return True

def blastnCmds(dbName, fastaFile, blastOutputFile, fmtString = def_blastnOutfmt):
    blastCmd = """blastn -query %s -db %s -outfmt '%s' -evalue %s -out %s\n""" % (fastaFile, dbName, def_blastnOutfmt, def_evalue, blastOutputFile)
    return blastCmd


def bestAlignment(blastnFile, alignCutOff, pidentCutOff, fh):
    # def_blastnOutfmt = '6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovhsp sstrand'
    seqAlignDict = {}
    count = 0
    for line in open(blastnFile, "r"):
        alignmentData = line.split()
        fastaLabel = alignmentData[0]
        pident     = float(alignmentData[2])
        length = int(alignmentData[3])
        # qcovhsp    = float(alignmentData[14])

        try:
            # Test pident if > 91%
            if (pident > pidentCutOff and length > alignCutOff):
                # seqAlignDict[fastaLabel].append(alignmentData)
                alignmentData[0] += "[%d]" % count # change fasta label
                newLine = "\t".join(alignmentData)
                seqAlignDict[fastaLabel].append(newLine)
                count += 1
        except KeyError:
            count = 0
            # alignmentData[0] += "[%d]" % count # change fasta label
            # newLine = "\t".join(alignmentData)
            seqAlignDict[fastaLabel] = [newLine]
            count += 1
    #####

    # Book keeping
    for label in seqAlignDict:
        for data in seqAlignDict[label]:
            fh.write("%s\n" % data)



custom_fields = ['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'sstrand', 'qcovhsp']
class BlastAlignmentData:
    def __init__(self, qseqid, sseqid, pident, length, qlen, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, sstrand, qcovhsp):
        self.qseqid = qseqid
        self.sseqid = sseqid
        self.pident = pident
        self.length = length
        self.qlen = qlen
        self.mismatch = mismatch
        self.gapopen = gapopen
        self.qstart = qstart
        self.qend = qend
        self.qend = qend
        self.sstart = sstart
        self.send = send
        self.evalue = evalue
        self.bitscore = bitscore
        self.sstrand = sstrand
        self.qcovhsp = qcovhsp
    def show(self):
        print("%s %s" % (self.pident, self.qcovhsp))

#========================

def getBlastnSeq(blastnFile, outputFile):
    fh = open(outputFile, "w")
    # print(blastnFile)
    # print(outputFile)
    for line in open(blastnFile, "r"):
        sseq = line.split()[-1] # sequence is the last line
        blastnId = line.split()[0] # accession number and description
        acc, description = blastnId.split(";", 1)
        seq_obj = Seq(sseq)

        cluster = line.split()[1].split(";")[1]
        record = SeqRecord(seq_obj, "%s;%s;%s" % (acc, cluster, description), '', '')
        # print(record)
        # records.append(record)
        SeqIO.write(record, fh, "fasta")
    ####
    fh.close()


#========================
def reHead_fasta(fastaFileName, outputFileName):
    fh = open(outputFileName, "w")
    for record in SeqIO.parse(fastaFileName, "fasta"):
        # acc = record.id
        headerData = record.description.split(None, 1)
        description = "_".join(headerData[1].split()) # so blastn will keep everything
        headerData[1] = description
        # record.id = acc + ";" + description
        record.id = ";".join(headerData)
        record.description = ""

        SeqIO.write(record, fh, "fasta")
    #####
    fh.close()
    # SeqIO.write(records, "tmp.fasta", "fasta")


#========================
def minimap_filter_alignments(pafFilePath, alignCutOff):
    alignSet = set([])
    for line in open(pafFilePath, "r"):
        alignData = line.split("\t")

        numMismatches = float(alignData[9]) # Col 10
        alignLen = float(alignData[10]) # Col 11
        # print("Mismatches: %d, AlignLen: %d" % (numMismatches, alignLen))
        if (numMismatches / alignLen < .25 and alignLen > alignCutOff):
            # good enough alignment
            alignSet.add(alignData[5])
        #####
    #####
    return alignSet


#========================
def deduplicate(fastaFile):

    # cluster, source, x = fastaFile.split(".", 2)
    # outputFile = "%s.%s.dup.fasta" % (cluster, source)
    header, x = fastaFile.split("/")[-1].split(".", 1)
    outputFile = "%s.dup.fasta" % (header)
    cdHitDupCmd = "cd-hit-dup -i %s -o %s" % (fastaFile, outputFile)
    n = subprocess.Popen(cdHitDupCmd, shell=True)
    n.poll()
    return outputFile
    # return cdHitDupCmd


# TESTING
#========================
def testPrintFile(fileName):
    print("\n======== %s =========" % fileName)
    for line in open(fileName, "r"):
        print(line.strip())
    ######
    print("\n======== %s =========" % fileName)
    assert False

#========================
def test(*statements):
    print("\n==========TESTING===========")
    for statement in statements:
        print(statement)
    #####
    print("============================")
    assert False




