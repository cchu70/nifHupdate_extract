#!/usr/bin/python
__author__="Claudia Chu"
__date__ ="2/24/19"

from nifHupdate_lib import parseConfig, createShFile, parseDate, launch, esearchCmds, fastaCmds, blastnCmds, bestAlignment, \
 throwError, verifyDb, test, testPrintFile, wait, fasta, trimSeq, mapBlast, deduplicate, minimapCmds, mapEsearch, reHead, \
 def_minimap_align_len_cutoff, whichFastaFofn, extractFileName

from os.path import abspath, join, isfile

from optparse import OptionParser

from os import curdir, system, stat

from sys import stdout

from Bio import SeqIO

import time

import datetime

import subprocess

from sys import argv, stderr

import threading

# ======================================================

def real_main():
    # Read the command line information
    configFile  = argv[1]
    stage       = argv[2]
    basePath    = argv[3]
    logFile     = argv[4]


    # tracking time
    logFileHandle = open(logFile, "a")
    logFileHandle.write("================================\n")
    logFileHandle.write("Stage: %s\n" % stage)
    logFileHandle.write("Time Start: %s\n" % datetime.datetime.now())

    # Read the config file
    configDict = parseConfig(configFile, basePath, logFileHandle)

    if (configDict["PATH"] == 'edirect'):
        # Get numerical year
        startDate = parseDate(configDict["START"])
        endDate = parseDate(configDict["END"])
    #####

    # FREQUENTLY USED CONSTANTS AND VARIABLES
    CMDLIST = []
    PREFIX = configDict["PREFIX"]


# ================= STAGE 1 ===================== #
# esearch: return all accession numbers to records associated with
# the search term.

    if (stage == 'esearch'):
        logFileHandle.write("\n=============== RESTART =================\n")
        # Track esearch text files
        esearchFofn = "%s.esearch.fofn" % PREFIX
        fh = open(esearchFofn, "w")

        for year in range (startDate, endDate):
            outputfile = "%s_%s.esearch.txt\n" % (PREFIX, year)

            # Replace with literal year
            if (year == startDate):
                year = configDict["START"]
            elif (year == endDate):
                year = configDict["END"]
            #####

            CMDLIST.append("echo esearch for year %s..." % year)
            CMDLIST.append(esearchCmds(configDict, year, outputfile))
            fh.write(outputfile)
        #####
        fh.close()

        # Next stage
        nextStage = 'fasta'

# ================= STAGE 2 ===================== #
    elif (stage == 'fasta'):

        fastaFofn = "%s.fasta.fofn" % PREFIX
        fh = open(fastaFofn, "w")
        esearchFofn = "%s/%s/%s.esearch.fofn" % (basePath, PREFIX, PREFIX)
        if (not isfile(esearchFofn)):
            throwError("esearch stage failed: %s not found" % esearchFofn, logFileHandle)
        else:

            threads = []
            for esearchFile in open(esearchFofn, "r"):
                if (not isfile(esearchFile.strip())):
                    throwError("esearch stage failed: %s not found" % file, logFileHandle)
                elif (not stat(esearchFile.strip())):
                    throwError("esearch stage failed: %s empty" % file, logFileHandle)
                #####

                # Make the fasta commands
                prefix = esearchFile.split(".")[0]
                for sortterm in configDict["SORTTERMS"]:
                    fastaFileName = "%s.%s.fasta" % (prefix, sortterm)
                    fh.write("%s\n" % fastaFileName) # record fasta file names

                    CMDLIST.append("echo Retrieving %s ..." % fastaFileName)
                    CMDLIST.append(fastaCmds(esearchFile.strip(), sortterm, basePath, fastaFileName))
                    print("getting %s" % esearchFile)

                #####
            #####
        #####

        fh.close()

        # clean tmp files
        rmCmd = "rm tmp*"
        CMDLIST.append(rmCmd)

# ================ STAGE 2.1 ==================== #
    elif(stage == 'fasta_rehead'):

        esearchFofn = "%s/%s/%s.esearch.fofn" % (basePath, PREFIX, PREFIX)
        esearchMap = mapEsearch(esearchFofn)
        fastaFofn = "%s.fasta.fofn" % PREFIX

        fastaReHeadFofn = "%s.fasta.rehead.fofn" % PREFIX
        fh = open(fastaReHeadFofn, "w")
        if (not isfile(fastaFofn)):
            throwError("fasta stage failed: %s not found" % esearchFofn, logFileHandle)
        else:
            for fastaFile in open(fastaFofn, "r"):
                prefix, source, end = fastaFile.strip().split(".", 2)
                outputFastaFileName = "%s.%s.rehead.fasta" % (prefix, source)
                fh.write("%s\n" % outputFastaFileName)
                print(fastaFile)
                reHead(fastaFile.strip(), esearchMap, outputFastaFileName)
            #####
        #####
        fh.close()

        # Next stage
        nextStage = 'set_db'

# ================= STAGE 3 ===================== #
# Check if all the files exist for a local database for stand alone alignment

    elif (stage == 'set_db'):
        # # TESTING
        # test("In set_db!")
        # ######

        # check database already exists
        if (not verifyDb(configDict["DBFILE"])):
            # move to next stage
            makeDbCmd = "makeblastdb -in %s -parse_seqids -dbtype nucl -out %s/%s/%s" % (configDict["DBFILE"], basePath, PREFIX, configDict["DBNAME"])
            CMDLIST.append(makeDbCmd)
        #####

        nextStage = 'blastn'

# ================= STAGE 4 ===================== #
# Aligning each of the query sequences to the subject sequences from local database.

    elif (stage == 'blastn'):
        fastaFofn = whichFastaFofn(configDict)

        if (not isfile(fastaFofn)):
            throwError("fasta stage failed: %s not found" % esearchFofn, logFileHandle)
        else:

            CMDLIST.append("mkdir blastn")

            fofnFileName = "%s.blastnFiles.fofn" % (PREFIX)
            fh = open(fofnFileName, "w")

            for fastaFile in open(fastaFofn, "r"):

                prefix, source, end = fastaFile.strip().split("/")[-1].split(".", 2)
                outputFileName = "%s.%s.blastn.txt" % (prefix, source)
                fh.write("%s/%s/blastn/%s\n" % (basePath, PREFIX, outputFileName))  # write to blast fofn file
                CMDLIST.append("echo Blasting %s ..." % fastaFile.strip())
                CMDLIST.append(blastnCmds(configDict["DBNAME"], fastaFile.strip(), outputFileName))

            fh.close()
        #####

        CMDLIST.append("mv *.blastn.txt ./blastn")

        nextStage = 'filter_best_alignments'

# ================= STAGE 5 ===================== #
# Parsing through blastn tables to find best alignment for each sequence

    elif (stage == 'filter_best_alignments'):

        fofnFileName = "%s.blastnFiles.fofn" % (PREFIX)
        if (not isfile(fofnFileName)) :
            throwError("%s is not available. Either re-run blastn stage, or cat all your blastn files into this file name." % fofnFileName, logFileHandle)
        #####

        CMDLIST.append("mkdir filter_best_alignments")

        blastnFofn = "%s.blastnFiles.filter.fofn" % (PREFIX)
        ch = open(blastnFofn, "w")

        for blastnFile in open(fofnFileName, "r"):
            prefix, source, end = blastnFile.strip().split("/")[-1].split(".", 2)
            fileName = "%s.%s.blastn.filter.txt" % (prefix, source)
            fh = open(fileName, "w")
            bestAlignment(blastnFile.strip(), fh)
            fh.close()

            # Fofn file
            ch.write("%s/%s/filter_best_alignments/%s\n" % (basePath, PREFIX, fileName))
        #####

        ch.close()

        CMDLIST.append("mv *.blastn.filter.txt ./filter_best_alignments")

        nextStage = 'trim_seq'

# ================= STAGE 6 ===================== #

    elif (stage == 'trim_seq'):

        print('In trim_seq')

        # parse the filtered blastn files into a map
        blastnFofn = "%s.blastnFiles.filter.fofn" % (PREFIX)
        blastnMap = mapBlast(blastnFofn)

        fastaTrimmedFofn = "%s.fasta.trimmed.fofn" % PREFIX
        fh = open(fastaTrimmedFofn, "w")

        # original fasta file
        fastaFofn = whichFastaFofn(configDict)

        CMDLIST.append("mkdir trim_seq")

        for fastaFile in open(fastaFofn, "r"):

            # prefix, source, end = fastaFile.strip().split(".", 2)

            # trimFastaFileName = "%s.%s.trimmed.fasta" % (prefix, source)
            path, fileName, fileNameHead = extractFileName(fastaFile.strip())
            trimFastaFileName = fileNameHead + ".trimmed.fasta"
            fh.write("%s/%s/trim_seq/%s\n" % (basePath, PREFIX, trimFastaFileName))
            trimSeq(fastaFile.strip(), trimFastaFileName, blastnMap)

        #####
        CMDLIST.append("mv *.trimmed.fasta ./trim_seq")

        fh.close()

        nextStage = 'cluster'


# ================= STAGE 7 ===================== #
    elif (stage == 'cluster'):

        clusterFilesFofn = "%s.cluster_fasta.fofn" % PREFIX
        ch = open(clusterFilesFofn, "w")

        clusterFileHandles = {}
        fastaTrimmedFofn = "%s.fasta.trimmed.fofn" % PREFIX

        for fastaTrimmedFileName in open(fastaTrimmedFofn, "r"):
            prefix, source, end = fastaTrimmedFileName.split(".", 2)

            altClusterName = "other.%s.fasta" % source

            for record in SeqIO.parse(fastaTrimmedFileName.strip(), "fasta"):
                cluster = record.id.split(";")[1]
                fileName = cluster + ".%s.fasta" % source
                try:
                    fh = clusterFileHandles[fileName]
                except:
                    print(fileName)
                    try:
                        # new file
                        fh = open(fileName, "a")
                        ch.write("%s\n" % fileName)
                        clusterFileHandles[fileName] = fh
                        # wont work if the cluster name is weirdly formatted
                    except:
                        # rename
                        try:
                            fileName = altClusterName
                            fh = clusterFileHandles[altClusterName]
                        except:
                            # new file
                            fh = open(fileName, "a")
                            ch.write("%s\n" % fileName)
                            clusterFileHandles[fileName] = fh
                    #####
                #####
                SeqIO.write(record, clusterFileHandles[fileName], "fasta")
            #####
        #####
        ch.close()

        for clusterFileHandle in clusterFileHandles:
            clusterFileHandles[clusterFileHandle].close()

        nextStage = 'deduplicate'


# ================= STAGE 8 ===================== #
    elif (stage == 'deduplicate'):
    # cd-hit-dup -i fasta -o output
        print('In deduplication stage!')

        CMDLIST.append("mkdir clusters")

        clusterFastaFofn = "%s.cluster_fasta.fofn" % PREFIX

        dedupFastaFofn = "%s.cluster_fasta_dedup.fofn" % PREFIX
        fh = open(dedupFastaFofn, "w")
        for fastaFile in open(clusterFastaFofn, "r"):
            print(fastaFile)
            CMDLIST.append(deduplicate(fastaFile.strip()))
            fh.write(fastaFile)
        #####
        rmCmd = 'rm -r *.clstr'
        CMDLIST.append(rmCmd)
        CMDLIST.append("mv *.dup.fasta ./clusters")
        nextStage = 'end'

# ================= STAGE 9 ===================== #
    elif (stage == 'end'):
        logFileHandle.write("End Time: %s\n" % datetime.datetime.now())
        logFileHandle.close()

        CMDLIST.append("mkdir ./shFiles")
        CMDLIST.append("mv .sh ./shFiles")
        CMDLIST.append("echo nifHpdate complete!")
        CMDLIST.append("cat %s" % logFile)

        shFileName = createShFile(CMDLIST, basePath, PREFIX, stage)
        launch(shFileName)
        return 0
        # print stats

# ============= ALTERNATIVE APPROACH ============ #
# ================= STAGE A ===================== #
    elif (stage == 'minimap'):
        logFileHandle.write("\n=============== RESTART =================\n")

        print("\nIn minimap stage\n")
        oldSeqDB = configDict["DBFILE"] #full path
        nuccoreDBFofn = "%s/%s" % (basePath, configDict["NUCCORE"]) # path to nuccore files

        miniMapFofn = "%s.minimap.fofn" % PREFIX
        fh = open(miniMapFofn, "w")

        CMDLIST.append("mkdir minimap")

        for nuccoreFile in open(nuccoreDBFofn, "r"):
            outputFileName = "%s.minimap.paf" % (nuccoreFile.split("/")[-1].split(".")[0])
            fh.write("%s/%s/minimap/%s\n" % (basePath, PREFIX, outputFileName))
            CMDLIST.append("echo minimapping %s" % nuccoreFile.split("/")[-1].split(".")[0])
            CMDLIST.append(minimapCmds(nuccoreFile.strip(), oldSeqDB, outputFileName))


        fh.close()

        CMDLIST.append("mv *.paf ./minimap")

        nextStage = 'minimap_filter'

# ================= STAGE B ===================== #
    elif (stage == 'minimap_filter'):
        print("\nIn minimap_filter stage\n")
        # Col 10: Mismatches
        # Col 11: Alignment length

        # if mismatches / alignment length < 25% -> throw out
        # Collect the headers of the fasta files with +75% similarity

        # Filtered fasta file fofn
        fofnFileName = "%s.minimap_filter.fofn" % PREFIX
        fh = open(fofnFileName, "w")

        pafFofnFile = "%s.minimap.fofn" % PREFIX

        # Parse the fofn file
        nuccoreDBFofn = "%s/%s" % (basePath, configDict["NUCCORE"]) # path to nuccore fofn
        nuccoreFilePathDict = {}
        for filePath in open(nuccoreDBFofn, "r"):
            fastaFileID = filePath.split("/")[-1].split(".")[0]
            nuccoreFilePathDict[fastaFileID] = filePath.strip()
        #####

        CMDLIST.append("mkdir minimap_filter")

        # extract the sequences
        for pafFile in open(pafFofnFile, "r"):
            fastaFileID = pafFile.split("/")[-1].split(".")[0] # exclude the file ending
            newFastaFileName = fastaFileID + ".filtered.fasta"
            alignSet = set([])
            for line in open(pafFile.strip(), "r"):
                alignData = line.split("\t")

                numMismatches = float(alignData[9]) # Col 10
                alignLen = float(alignData[10]) # Col 11
                # print("Mismatches: %d, AlignLen: %d" % (numMismatches, alignLen))
                if (numMismatches / alignLen < .25 and alignLen > def_minimap_align_len_cutoff):
                    # good enough alignment
                    alignSet.add(alignData[5])
                #####
            #####
            if alignSet:
                # list is not empty
                extractIDs = "/" + "/,/".join(alignSet) + "/"
                extractCmd = """gunzip -dc %s | awk '%s{p++;print;next} /^>/{p=0} p' > %s""" % (nuccoreFilePathDict[fastaFileID], extractIDs, newFastaFileName)
                CMDLIST.append(extractCmd)
                fh.write("%s/%s/minimap_filter/%s\n" % (basePath, PREFIX, newFastaFileName)) # tracking directory
            #####
            #####
        #####
        CMDLIST.append("mv *.filtered.fasta ./minimap_filter")

        fh.close()

        # Continue with blast
        nextStage = 'set_db'
    #####


    # =============== MOVE TO NEXT STAGE =============== #

    logFileHandle.write("End Time: %s\n" % datetime.datetime.now())
    logFileHandle.close()

    echoCmd = "echo Launching next stage %s" % nextStage
    nextCmd = "python3 %s/nifHupdate_Lib/nifHupdate_launch.py %s %s %s %s\n" % (basePath, configFile, nextStage, basePath, logFile)
    CMDLIST.append(echoCmd)
    CMDLIST.append(nextCmd)
    shFileName = createShFile(CMDLIST, basePath, PREFIX, stage)
    launch(shFileName)


    return 0
#==============================================================
if ( __name__ == '__main__' ):
    real_main()