#!/usr/bin/python
__author__="Claudia Chu"
__date__ ="2/18/19"


from os.path import abspath, join, isfile
import os
from Bio import SeqIO
import subprocess
import time
from shutil import copyfile
import datetime
from sys import argv, stderr, exit

#========================
# Defaults
def_query    = "nifH"
def_tag      = "GENE"
def_start    = "2012"
def_end      = "2019"
def_datetype = "PDAT"
def_elements = "Id Caption TaxId Slen CreateDate Organism Title"
def_sorttype = ["nifH", "genome"]

def_blastnOutfmt = '6 qseqid sseqid pident length qlen mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovhsp'
def_dbfiles = ["nhr", "nsd", "nin", "nsi", "nsq"]
def_dbname = "DB"
def_evalue = 0.001
def_minimap_align_len_cutoff = 200

# #========================
# Allowed sets
DATETYPES = set(["PDAT"])
edirect_stages = set(['esearch', 'fasta', 'fasta_rehead','set_db', 'blastn', 'filter_best_alignments', 'trim_seq', 'cluster', 'deduplicate', 'end'])
minimap_stages = set(['minimap', 'minimap_filter', 'set_db', 'blastn', 'filter_best_alignments', 'trim_seq', 'cluster', 'deduplicate', 'end'])
MAX_REQUESTS = 3 # for entrex direct



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
def esearchCmds(configDict, year, outputFile):

    esearch  = """esearch -db nucleotide -query \"%s\""""    % (configDict["QUERY"])
    efilter  = """efilter -mindate %s -maxdate %s -datetype %s""" % (year, year, configDict["DATETYPE"])
    efetch_1 = """efetch -format docsum"""
    xtract   = """xtract -pattern DocumentSummary -element %s"""  % (def_elements)
    edirectCmd = """%s | %s | %s | %s > %s\n""" % (esearch, efilter, efetch_1, xtract, outputFile)

    return edirectCmd

#========================
def minimapCmds(refSeqFileName, querySeqFileName, outputFileName):
    cmd = "minimap2 %s %s > %s" % (refSeqFileName, querySeqFileName, outputFileName)
    return cmd

#========================
def fastaCmds(esearchFile, grepFilter, basePath, fastaFileName):

    cat = "cat %s" % (esearchFile)
    acc = """awk 'BEGIN { ORS="," }; { print $2 }'"""

    fastaCmd = "%s | grep '%s' | %s | python3 %s/nifHupdate_Lib/nifHupdate_fasta.py > %s\n" % (cat, grepFilter, acc, basePath, fastaFileName)

    return fastaCmd

def fasta(esearchFile, sortterm, outputFasta):
    fastaFileHandle = open(outputFasta, "w")

    requestCont = 0
    for line in open(esearchFile.strip(), "r"):
        if (sortterm in line):
            tmpName = "tmp.%s" % outputFasta
            tmpFileHandle = open(tmpName, "w")
            recId, acc, taxId, slen, date, organism, title = line.strip().split("\t")
            fastaCmd = """efetch -db nuccore -id %s -format gene_fasta """ % (acc)


            requestFinished = False
            while (not requestFinished):
                try:
                    print("Retrieving %s" % acc)
                    n = subprocess.Popen(fastaCmd, shell=True, stdin = subprocess.PIPE, stdout = tmpFileHandle)
                    time.sleep(1)
                    wait(n)
                    n.kill()

                    requestFinished = True
                    print("Done!")

                    requestCont += 1 # increment count

                    if (requestCont == MAX_REQUESTS ):
                        time.sleep(10)
                        requestCont = 0
                    #####


                except:
                    requestFinished = False
                    print("Busy...\n")
                    time.sleep(10) # Sleep for 5 minutes
                    print("%s re-running" % acc)
                    # Then Try again
                #####
            #####
            i = 1

            for record in SeqIO.parse(tmpName, "fasta"):
                if ("nifH" in record.description):
                    record.description = "[Acc: %s] [Ver: %d] [Date: %s] [Organism: %s] [Title: %s] [TaxID: %s]" % (record.id, i, date, organism, title, taxId)
                    record.id = acc
                    SeqIO.write(record, fastaFileHandle, "fasta")
                    i += 1
            #####
        #####
    #####

#========================
def blastnCmds(dbName, fastaFile, blastOutputFile, fmtString = def_blastnOutfmt):


    blastCmd = """blastn -query %s -db %s -outfmt '%s' -evalue %s -out %s\n""" % (fastaFile, dbName, def_blastnOutfmt, def_evalue, blastOutputFile)
    return blastCmd

#========================
def bestAlignment(blastnFile, fh):
    seqAlignDict = {}
    for line in open(blastnFile, "r"):
        alignmentData = line.split()
        fastaLabel = alignmentData[0]
        pident     = alignmentData[2]
        qcovhsp    = alignmentData[14]

        try:
            # compare previous qcovhsp, and change if better score
            if (seqAlignDict[fastaLabel][14] < qcovhsp):
                seqAlignDict[fastaLabel] = alignmentData
                # pass in everything about the alignment
            elif (seqAlignDict[fastaLabel][14] == qcovhsp):
                if (seqAlignDict[fastaLabel][2] < pident):
                    seqAlignDict[fastaLabel] = alignmentData
        except:
            seqAlignDict[fastaLabel] = line
        #####
    #####

    # Book keeping
    for label in seqAlignDict:
        # fh.write("\t".join(seqAlignDict[label]) + "\n")
        fh.write(seqAlignDict[label])


#========================
def verifyDb(dbLabel):
    checkFiles = [int(isfile("%s.%s" % (dbLabel, fileType))) for fileType in def_dbfiles]
    #check if all the database files exist
    if (sum(checkFiles) != len(def_dbfiles)):
        return False
    else:
        return True


#========================
def trimSeq(fastaFileName, outputFileName, blastItems):
    # i = 0
    outputFileHandle = open(outputFileName, "w")
    for record in SeqIO.parse(fastaFileName.strip(), "fasta"):
        if record.id in blastItems:

            blastnData = blastItems[record.id]

            # check sequence length
            length = int(blastnData.length)
            qlen = int(blastnData.qlen)
            distance = qlen - length

            # if (distance >= 340 and distance <= (0.5 * qlen)):
            #     record.seq = record.seq[int(blastnData.qstart):int(blastnData.qend)]
            # #####
            record.seq = record.seq[int(blastnData.qstart):int(blastnData.qend)]

            cluster = blastnData.sseqid.split(';')[1]

            # organism = record.description;
            organism = record.description.split(None, 1)[1]; # for nuccore

            header = "%s;%s;%s" % (record.id, cluster, organism)
            print(header)
            record.id = header
            record.description = ""
            SeqIO.write(record, outputFileHandle, "fasta")

#========================
def mapBlast(blastnFofn):
    blastnMap = {}
    for blastnFile in open(blastnFofn, "r"):
        for line in open(blastnFile.strip(), "r"):
            alignmentData = line.split(None)
            qseqid, sseqid, pident, length, qlen, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, sstrand, qcovhsp = line.strip().split(None)
            blastnMap[qseqid] = BlastAlignmentData(qseqid, sseqid, pident, length, qlen, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, sstrand, qcovhsp)
        #####
    #####
    return blastnMap

#========================
def mapEsearch(esearchFofn):
    # "Id Caption TaxId Slen CreateDate Organism Title"
    esearchMap = {}
    for file in open(esearchFofn, "r"):
        for line in open(file.strip(), "r"):
            iD, caption, taxId, slen, createDate, organism, title = line.split('\t')
            esearchMap[caption] = organism
        #####
    #####
    return esearchMap

#========================
def reHead(fastaFileName, esearchMap, outputFile):
    fh = open(outputFile, "w")

    currAcc = ""
    for record in SeqIO.parse(fastaFileName, "fasta"):
        acc = record.id.split(".")[0].split('|')[1]

        try:
            # MAJOR ASSUMPTION THAT THE FIRST ONE WILL ACTUALLY WORK
            organism = esearchMap[acc]
            currAcc = acc
            record.id = acc
            record.description = organism
        except:
            organism = esearchMap[currAcc]
            record.id = currAcc;
            record.description = organism
        #####
        SeqIO.write(record, fh, "fasta")
    #####
    fh.close()

#========================
def whichFastaFofn(configDict):
    fastaFofn = ""
    if (configDict['PATH'] == 'edirect'):
        fastaFofn = "%s.fasta.rehead.fofn" % configDict['PREFIX']
    elif (configDict['PATH'] == 'minimap'):
        fastaFofn = "%s.minimap_filter.fofn" % configDict['PREFIX']

    return fastaFofn

# ====================================================================
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
def deduplicate(fastaFile):

    # cluster, source, x = fastaFile.split(".", 2)
    # outputFile = "%s.%s.dup.fasta" % (cluster, source)
    header, x = fastaFile.split("/")[-1].split(".", 1)
    outputFile = "%s.dup.fasta" % (header)
    cdHitDupCmd = "cd-hit-dup -i %s -o %s" % (fastaFile, outputFile)
    print(cdHitDupCmd)
    n = subprocess.Popen(cdHitDupCmd, shell=True)
    n.poll()
    # return cdHitDupCmd

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

    configDict = {}
    for line in open(configFile, "r"):
        key, val = line.strip().split(None, 1)
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
            throwError("Could not find database fasta file %s" % dbFile, logFileFh)
    except KeyError:
        throwError("No database fasta file provided in configuration file %s" % configFile.split("/")[-1], logFileFh)
    # -----------------------------
    try:
        x = configDict["DBNAME"]
    except KeyError:
        configDict["DBNAME"] = def_dbname


    try:
        errFilePath = configDict["ERRFILE"]
        if (not isfile(esearchFile)):
            # make the file
            fh = open(errFilePath, "w")
            fh.close()
    except:
        configDict["ERRFILE"] = "%s/%s_err.txt" % (basePath, configDict["PREFIX"])
        fh = open(configDict["ERRFILE"], "w")
        fh.close()

    # -----------------------------
    # Test specifying what part of the pipeline to use
    try:
        path = configDict["PATH"]
    except:
        throwError("No PATH provided in the configuration file %s" % configFile.split("/")[-1], logFileFh)

    if (path == 'edirect'):
        try:
            queries = configDict["QUERY"]
        except KeyError:
            throwError("No QUERY provided in the configuration file %s" % configFile.split("/")[-1], logFileFh)
        except:
            throwError("Improper formatting of query terms in the configuration file %s. , logFileFh\
                Please delimit with semicolons ';'", configFile.split("/")[-1])
        # -----------------------------
        # Test DATERANGE, START, END, and DATETYPE
        try:
            useDateRange = configDict["DATERANGE"]
            if (useDateRange.lower() == 'true'):
                try:
                    start = configDict["START"]
                    end   = configDict["END"]
                    if (parseDate(start) > parseDate(end)):
                        throwError("Invalid START and END date in Config File %s. ," +
                            "Check that start date is before end date." % configFile.split("/")[-1])
                except KeyError:
                    throwError("Unable to query in date range without specified start and end time. Please provide START and END tags (MM/DD/YYY, MM/YYYY, or YYYY) in the configuration file %s" % configFile.split("/")[-1], logFileFh)
                #####

                try:
                    datetype = configDict["DATETYPE"]
                except KeyError:
                    throwError("Unable to query data in date range without specified datetype. Options are %s" % ",".join(DATETYPES), logFileFh)
                #####

                # Everything there
                configDict["DATERANGE"] = True
            else:
                throwError("%s is not a valid value. Only valid value for DATERANGE is 'true'. Otherwise, do not include the DATERANGE tag.", useDateRange, logFileFh)
            #####
        except KeyError:
            # check if other tags were used
            try:
                start = configDict["START"]
                throwError("Provided START tag, but DATERANGE is not specified in configFile %s" % configFile.split("/")[-1], logFileFh)
            except:
                pass
            #####

            try:
                end   = configDict["END"]
                throwError("Provided END tag, but DATERANGE is not specified in configFile %s" % configFile.split("/")[-1], logFileFh)
            except:
                pass
            #####

            try:
                datetype = configDict["DATETYPE"]
                throwError("Provided DATETYPE tag, but DATERANGE is not specified in configFile %s" % configFile.split("/")[-1], logFileFh)
            except:
                pass
            #####

            # Default to not use date range
            configDict["DATERANGE"] = False

        # -----------------------------
        try:
            x = configDict["SORTTERMS"]
            configDict["SORTTERMS"] = x.split(None) # store the space delmited sort terms into a list
        except:
            # place default
            configDict["SORTTERMS"] = def_sorttype
        #####

    elif (path == 'minimap'):
        # -----------------------------
        try:
            nuccore = "%s/%s" % (basePath, configDict["NUCCORE"])
            if (not isfile(dbFile)):
                throwError("Could not find database fasta file %s" % dbFile, logFileFh)
            # if (os.stat(nuccore).st_size == 0):
            #     throwError("%s is empty.")
        except KeyError:
            throwError("No database fasta file provided in configuration file %s" % configFile.split("/")[-1], logFileFh)
        #####
    else:
        throwError("%s is not a valid PATH. Choose either edirect or minimap" % path)
    #####
    return configDict


#========================
def launch(shFileName):
    print("Launching %s" % shFileName)
    # n = subprocess.Popen(["bash", "%s >> %s" % (shFileName, outFileName)])
    n = subprocess.Popen(["bash", shFileName])
    wait(n) # wait for the shfile to finish running before proceeding
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




