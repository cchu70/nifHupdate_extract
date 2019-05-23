#!/usr/bin/python
__author__="Claudia Chu"
__date__ ="2/24/19"

<<<<<<< HEAD
from nifHupdate_lib import parseConfig, createShFile, parseDate, launch, blastnCmds, bestAlignment, \
 throwError, verifyDb, test, testPrintFile, wait, deduplicate, \
 def_minimap_align_len_cutoff, extractFileName, reHead_fasta, getBlastnSeq, minimap_filter_alignments, \
=======
from nifHupdate_lib import parseConfig, createShFile, parseDate, launch, esearchCmds, fastaCmds, blastnCmds, bestAlignment, \
 throwError, verifyDb, test, testPrintFile, wait, fasta, trimSeq, mapBlast, deduplicate, minimapCmds, mapEsearch, reHead, \
 def_minimap_align_len_cutoff, whichFastaFofn, extractFileName, reHead_fasta, getBlastnSeq, minimap_filter_alignments, \
>>>>>>> 610dd852ddf62c86371230fb32e4a7cf781260ba
 def_expectedFastaFofns

from os.path import abspath, join, isfile, isdir, getsize

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



    # # tracking time
    logFileHandle = open(logFile, "a")
<<<<<<< HEAD
=======
    # logFileHandle.write("================================\n")
    # logFileHandle.write("Stage: %s\n" % stage)
    # logFileHandle.write("Time Start: %s\n" % datetime.datetime.now())
>>>>>>> 610dd852ddf62c86371230fb32e4a7cf781260ba

    # Read the config file
    configDict = parseConfig(configFile, basePath, logFileHandle)

<<<<<<< HEAD
=======
    # if (configDict["PATH"] == 'edirect'):
    #     # Get numerical year
    #     startDate = parseDate(configDict["START"])
    #     endDate = parseDate(configDict["END"])
    # #####

>>>>>>> 610dd852ddf62c86371230fb32e4a7cf781260ba
    # FREQUENTLY USED CONSTANTS AND VARIABLES
    CMDLIST = []
    PREFIX = configDict["PREFIX"]

    # Is the pipeline at the end
    finished = False


# ================= STAGE 1 ===================== #
<<<<<<< HEAD
    if (stage == 'minimap'):
        # print("\nIn minimap stage\n")
        oldSeqDB = configDict["DBFILE"] #full path
        nuccoreDBFofn = configDict["NUCCORE"] # path to nuccore files

        outputDir = "minimap_output"
        miniMapFofn = "%s.minimap.fofn" % PREFIX
        fh = open(miniMapFofn, "w")

        if (not isdir(outputDir)):
            CMDLIST.append("mkdir %s" % outputDir)
        #####

        for nuccoreFile in open(nuccoreDBFofn, "r"):
            outputFileName = "%s.minimap.paf" % (nuccoreFile.split("/")[-1].split(".")[0])
            fh.write("%s/%s/%s/%s\n" % (basePath, PREFIX, outputDir, outputFileName))
            CMDLIST.append("echo minimapping %s >> %s" % (nuccoreFile.split("/")[-1].split(".")[0], logFile))

            minimapCmd = "minimap2 %s %s > %s" % (nuccoreFile.strip(), oldSeqDB, outputFileName)

            CMDLIST.append(minimapCmd)
            CMDLIST.append("mv %s ./%s" % (outputFileName, outputDir))

        fh.close()



        nextStage = 'minimap_fasta'

# ================= STAGE 2 ===================== #
    elif (stage == 'minimap_fasta'):
        # Col 10: Mismatches
        # Col 11: Alignment length

        # Filtered fasta file fofn
        outputDir = "minimap_filter_output"
        fofnFileName = "%s.minimap_fasta.fofn" % PREFIX
        fh = open(fofnFileName, "w")

        pafFofnFile = "%s.minimap.fofn" % PREFIX

        # Parse the fofn file
        nuccoreDBFofn = configDict["NUCCORE"].strip() # path to nuccore fofn
        nuccoreFilePathDict = {}
        for filePath in open(nuccoreDBFofn, "r"):
            fastaFileID = filePath.split("/")[-1].split(".")[0]
            nuccoreFilePathDict[fastaFileID] = filePath.strip()
        #####

        if (not isdir(outputDir)):
            CMDLIST.append("echo Making new directory %s >> %s" % (outputDir, logFile))
            CMDLIST.append("mkdir %s" % outputDir)

        # extract the sequences
        for pafFilePath in open(pafFofnFile, "r"):
            fastaFileID = pafFilePath.split("/")[-1].split(".")[0] # exclude the file ending
            newFastaFileName = fastaFileID + ".filtered.fasta"
            alignSetAccessions = minimap_filter_alignments(pafFilePath.strip(), configDict["MIN_MINIMAP_ALIGNLEN"])

            if alignSetAccessions:
                # list is not empty
                extractIDs = "|".join(alignSetAccessions) # Regex OR for exact text matches

                CMDLIST.append("echo Retrieving %s >> %s" % (",".join(alignSetAccessions), logFile))
                if (nuccoreFilePathDict[fastaFileID].split(".")[-1] == "gz"):
                    extractCmd = """gunzip -dc %s | awk '/%s/{p++;print;next} /^>/{p=0} p' > %s""" % (nuccoreFilePathDict[fastaFileID], extractIDs, newFastaFileName)
                else:
                    extractCmd = """cat %s | awk '/%s/{p++;print;next} /^>/{p=0} p' > %s""" % (nuccoreFilePathDict[fastaFileID], extractIDs, newFastaFileName)
                #####
                CMDLIST.append(extractCmd)
                CMDLIST.append("mv %s ./%s" % (newFastaFileName, outputDir))
                fh.write("%s/%s/%s/%s\n" % (basePath, PREFIX, outputDir, newFastaFileName)) # tracking directory
            #####
=======
# esearch: return all accession numbers to records associated with
# the search term.

    if (stage == 'esearch'):

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
>>>>>>> 610dd852ddf62c86371230fb32e4a7cf781260ba
            #####
        #####

        fh.close()
<<<<<<< HEAD
        if (getsize(fofnFileName) == 0):
            CMDLIST.append("echo No alignments found. %s is empty. Terminating pipeline." % fofnFileName)
            nextStage = "end"
        else:
            # Continue with next stage
            nextStage = 'rehead_fasta'
        #####

# ================= STAGE 3 ===================== #
    elif (stage == "rehead_fasta"):

        # Filtered fasta files
        fofnFileName = "%s.minimap_fasta.fofn" % PREFIX

        fofnNewFileName = "%s.fasta.rehead.fofn" % PREFIX
        fh = open(fofnNewFileName, "w")



        for fastaFile in open(fofnFileName, "r"):
            # print(fastaFile)
            path, fileName, fileNameHead = extractFileName(fastaFile.strip())
            newFastaFileName = fileNameHead + ".minimap_rehead.fasta"
            fh.write("%s/%s/minimap_rehead/%s\n" % (basePath, PREFIX, newFastaFileName)) # Replace later
            reHead_fasta(fastaFile.strip(), newFastaFileName)

        fh.close()

        if (not isdir("minimap_rehead")):
            CMDLIST.append("mkdir minimap_rehead")
        CMDLIST.append("mv *.minimap_rehead.fasta ./minimap_rehead")

        nextStage = "blastn"
    #####

# ================= STAGE 4 ===================== #

# Input: old database (fasta format), and new sequences to map and trime
# Output: Output format 6 blastn files with the following columns
#     1) qseqid
#     2) sseqid
#     3) pident
#     4) length
#     5) qlen
#     6) mismatch
#     7) gapopen
#     8) qstart
#     9) qend
#     10) sstart
#     11) send
#     12) evalue
#     13) bitscore
#     14) sstrand
#     15) qcovhsp
#     16) sseq
=======

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

# # ================= STAGE 3 ===================== #
# # Check if all the files exist for a local database for stand alone alignment

#     elif (stage == 'set_db'):
#         # # TESTING
#         # test("In set_db!")
#         # ######

#         # check database already exists
#         if (not verifyDb(configDict["DBNAME"])):
#             # move to next stage
#             makeDbCmd = "makeblastdb -in %s -parse_seqids -dbtype nucl -out %s/%s/%s" % (configDict["DBFILE"], basePath, PREFIX, configDict["DBNAME"])
#             CMDLIST.append("echo make blastn db")
#             CMDLIST.append(makeDbCmd)
#         else:
#             CMDLIST.append("DB files already exist")
#         #####

#         nextStage = 'blastn'

# ================= STAGE 4 ===================== #
# Aligning each of the query sequences to the subject sequences from local database.
>>>>>>> 610dd852ddf62c86371230fb32e4a7cf781260ba

    elif (stage == 'blastn'):

        # check database already exists
        if (not verifyDb(configDict["PREFIX"])):
            # move to next stage
            makeDbCmd = "makeblastdb -in %s -parse_seqids -dbtype nucl -out %s/%s/%s" % (configDict["DBFILE"], basePath, PREFIX, configDict["PREFIX"])
            CMDLIST.append("echo make blastn db with name %s >> %s" % (PREFIX, logFile))
            CMDLIST.append(makeDbCmd)
<<<<<<< HEAD
        #####

        outputDir = "blastn_output"
=======
        # else:
        #     CMDLIST.append("echo DB files already exist")
        # #####


        # fastaFofn = whichFastaFofn(configDict)
>>>>>>> 610dd852ddf62c86371230fb32e4a7cf781260ba
        fastaFofn = "%s.fasta.rehead.fofn" % configDict['PREFIX']

        if (not isfile(fastaFofn)):
            throwError("fasta stage failed: %s not found" % esearchFofn, logFileHandle)
        else:

<<<<<<< HEAD
            if (not isdir(outputDir)):
                CMDLIST.append("mkdir " + outputDir)
=======
            if (not isdir("blastn")):
                CMDLIST.append("mkdir blastn")
>>>>>>> 610dd852ddf62c86371230fb32e4a7cf781260ba
            #####

            fofnFileName = "%s.blastnFiles.fofn" % (PREFIX)
            fh = open(fofnFileName, "w")

            for fastaFile in open(fastaFofn, "r"):

                prefix, source, end = fastaFile.strip().split("/")[-1].split(".", 2)
                outputFileName = "%s.%s.blastn.txt" % (prefix, source)
<<<<<<< HEAD
                fh.write("%s/%s/%s/%s\n" % (basePath, PREFIX, outputDir, outputFileName))  # write to blast fofn file
=======
                fh.write("%s/%s/blastn/%s\n" % (basePath, PREFIX, outputFileName))  # write to blast fofn file
>>>>>>> 610dd852ddf62c86371230fb32e4a7cf781260ba
                CMDLIST.append("echo Blasting %s ... >> %s" % (fastaFile.strip(), logFile))
                CMDLIST.append(blastnCmds(configDict["PREFIX"], fastaFile.strip(), outputFileName))

            fh.close()
        #####

<<<<<<< HEAD
        CMDLIST.append("mv *.blastn.txt ./%s" % outputDir)
=======
        CMDLIST.append("mv *.blastn.txt ./blastn")
>>>>>>> 610dd852ddf62c86371230fb32e4a7cf781260ba

        nextStage = 'filter_best_alignments'

# ================= STAGE 5 ===================== #
# Parsing through blastn tables to find best alignment for each sequence

    elif (stage == 'filter_best_alignments'):

        fofnFileName = "%s.blastnFiles.fofn" % (PREFIX)
        if (not isfile(fofnFileName)) :
            throwError("%s is not available. Either re-run blastn stage, or cat all your blastn files into this file name." % fofnFileName, logFileHandle)
        #####

        if (not isdir("filter_best_alignments")):
            CMDLIST.append("mkdir filter_best_alignments")

        blastnFofn = "%s.blastnFiles.filter.fofn" % (PREFIX)
        ch = open(blastnFofn, "w")

        for blastnFile in open(fofnFileName, "r"):
            prefix, source, end = blastnFile.strip().split("/")[-1].split(".", 2)
            fileName = "%s.%s.blastn.filter.txt" % (prefix, source)
            fh = open(fileName, "w")
            bestAlignment(blastnFile.strip(), configDict["MIN_BLASTN_ALIGNLEN"], configDict["PIDENT_CUTOFF"], fh)
            fh.close()

            # Fofn file
            ch.write("%s/%s/filter_best_alignments/%s\n" % (basePath, PREFIX, fileName))
        #####

        ch.close()

        CMDLIST.append("mv *.blastn.filter.txt ./filter_best_alignments")

        nextStage = 'trim_seq'
<<<<<<< HEAD
=======
        # assert False
>>>>>>> 610dd852ddf62c86371230fb32e4a7cf781260ba

# ================= STAGE 6 ===================== #

    elif (stage == 'trim_seq'):

<<<<<<< HEAD

        # parse the filtered blastn files into a map
        blastnFofn = "%s.blastnFiles.filter.fofn" % (PREFIX)
=======
        # print('In trim_seq')

        # parse the filtered blastn files into a map
        blastnFofn = "%s.blastnFiles.filter.fofn" % (PREFIX)
        # blastnMap = mapBlast(blastnFofn)

        # # for b in blastnMap:
        # #     print(blastnMap[b].qseqid, blastnMap[b].sseqid, blastnMap[b].length)
        # # assert False

>>>>>>> 610dd852ddf62c86371230fb32e4a7cf781260ba

        fastaTrimmedFofn = "%s.fasta.trim.fofn" % PREFIX
        fh = open(fastaTrimmedFofn, "w")

        # # original fasta file
        # fastaFofn = whichFastaFofn(configDict)
        if (not isdir("trim_seq")):
            CMDLIST.append("mkdir trim_seq")
<<<<<<< HEAD
        #####

        for blastnFile in open(blastnFofn, "r"):
            path, fileName, fileNameHead = extractFileName(blastnFile.strip())
=======

        # for fastaFile in open(fastaFofn, "r"):

        #     # prefix, source, end = fastaFile.strip().split(".", 2)

        #     # trimFastaFileName = "%s.%s.trimmed.fasta" % (prefix, source)
        #     path, fileName, fileNameHead = extractFileName(fastaFile.strip())

        #     trimFastaFileName = fileNameHead + ".trimmed.fasta"
        #     fh.write("%s/%s/trim_seq/%s\n" % (basePath, PREFIX, trimFastaFileName))
        #     trimSeq(fastaFile.strip(), trimFastaFileName, blastnMap)

        # #####

        for blastnFile in open(blastnFofn, "r"):
            path, fileName, fileNameHead = extractFileName(blastnFile.strip())
            # fileId = blastnFile.split(".")[0]
            # print(fileName)
>>>>>>> 610dd852ddf62c86371230fb32e4a7cf781260ba

            trimFastaFileName = fileNameHead + ".trim.fasta"
            fh.write("%s/%s/trim_seq/%s\n" % (basePath, PREFIX, trimFastaFileName))

            getBlastnSeq(blastnFile.strip(), trimFastaFileName)

        CMDLIST.append("mv *.trim.fasta ./trim_seq")

        fh.close()

        nextStage = 'cluster'


# ================= STAGE 7 ===================== #
    elif (stage == 'cluster'):

        clusterFilesFofn = "%s.cluster_fasta.fofn" % PREFIX
        ch = open(clusterFilesFofn, "w")

        clusterFileHandles = {}
        fastaTrimmedFofn = "%s.fasta.trim.fofn" % PREFIX

        if (not isdir("clusters")):
            CMDLIST.append("mkdir clusters")
        #####

        for fastaTrimmedFileName in open(fastaTrimmedFofn, "r"):
            prefix, source, end = fastaTrimmedFileName.split(".", 2)

            altClusterName = "other.%s.fasta" % source

            for record in SeqIO.parse(fastaTrimmedFileName.strip(), "fasta"):
                cluster = record.id.split(";")[1]
                fileName = cluster + ".%s.fasta" % source
                try:
                    fh = clusterFileHandles[fileName]
                except:
                    # New cluster file found
                    try:
                        # new file
                        fh = open(fileName, "a")
                        # ch.write("%s\n" % fileName)
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
                            clusterFileHandles[fileName] = fh
                    #####
                #####
                SeqIO.write(record, clusterFileHandles[fileName], "fasta")
            #####
        #####

        for clusterFileHandle in clusterFileHandles:
            clusterFileHandles[clusterFileHandle].close()
            ch.write("%s/%s/clusters/%s\n" % (basePath, PREFIX, clusterFileHandle))

        ch.close()
        CMDLIST.append("mv *.trim.fasta clusters/")
        nextStage = 'deduplicate'


# ================= STAGE 8 ===================== #
    elif (stage == 'deduplicate'):
<<<<<<< HEAD
=======
    # cd-hit-dup -i fasta -o output
        # print('In deduplication stage!')
>>>>>>> 610dd852ddf62c86371230fb32e4a7cf781260ba

        if (not isdir("clusters_dedup")):
            CMDLIST.append("mkdir clusters_dedup")

        clusterFastaFofn = "%s.cluster_fasta.fofn" % PREFIX

        dedupFastaFofn = "%s.cluster_fasta_dedup.fofn" % PREFIX
        fh = open(dedupFastaFofn, "w")
        for fastaFile in open(clusterFastaFofn, "r"):
            newFastaFileName = deduplicate(fastaFile.strip())
            fh.write("%s/%s/clusters_dedup/%s\n" % (basePath, PREFIX, newFastaFileName))
        #####
        fh.close()
        rmCmd = 'rm -r *.clstr'
        CMDLIST.append(rmCmd)
        CMDLIST.append("mv *.dup.fasta ./clusters_dedup")
        nextStage = 'end'

# ================= STAGE 9 ===================== #
    elif (stage == 'end'):
<<<<<<< HEAD

=======
        # logFileHandle.write("End Time: %s\n" % datetime.datetime.now())
        # logFileHandle.close()
>>>>>>> 610dd852ddf62c86371230fb32e4a7cf781260ba
        finished = True

        if (not isdir("shFiles")):
            CMDLIST.append("mkdir ./shFiles")
        #####
        CMDLIST.append("mv *.sh ./shFiles")
        CMDLIST.append("echo nifHpdate complete!")

        for fileSuffix in def_expectedFastaFofns:
            fofn = PREFIX + "." + fileSuffix
            if (isfile(fofn)):
                # if the actual fofn file exists
                CMDLIST.append("echo Reading %s >> %s" % (fofn, logFile))
                for file in open(fofn):
                    CMDLIST.append("echo Number of sequences in %s >> %s" % (file.strip(), logFile))
                    CMDLIST.append("""cat %s | grep ">" | wc -l >> %s""" % (file.strip(), logFile))
                #####
            else:
                CMDLIST.append("echo %s is not a file" % fofn)
            #####
        #####

<<<<<<< HEAD
        CMDLIST.append("echo View %s for a brief summary" % logFile)

    # =============== MOVE TO NEXT STAGE =============== #

=======

        # Summary



# ============= ALTERNATIVE APPROACH ============ #
# ================= STAGE A ===================== #
    elif (stage == 'minimap'):
        # logFileHandle.write("\n=============== RESTART =================\n")

        # print("\nIn minimap stage\n")
        oldSeqDB = configDict["DBFILE"] #full path
        nuccoreDBFofn = configDict["NUCCORE"] # path to nuccore files

        miniMapFofn = "%s.minimap.fofn" % PREFIX
        fh = open(miniMapFofn, "w")

        if (not isdir("minimap")):
            CMDLIST.append("mkdir minimap")
        #####

        for nuccoreFile in open(nuccoreDBFofn, "r"):
            outputFileName = "%s.minimap.paf" % (nuccoreFile.split("/")[-1].split(".")[0])
            fh.write("%s/%s/minimap/%s\n" % (basePath, PREFIX, outputFileName))
            CMDLIST.append("echo minimapping %s >> %s" % (nuccoreFile.split("/")[-1].split(".")[0], logFile))
            CMDLIST.append(minimapCmds(nuccoreFile.strip(), oldSeqDB, outputFileName))
            CMDLIST.append("mv %s ./minimap" % outputFileName)

        fh.close()



        nextStage = 'minimap_filter'

# ================= STAGE B ===================== #
    elif (stage == 'minimap_filter'):
        # print("\nIn minimap_filter stage\n")
        # Col 10: Mismatches
        # Col 11: Alignment length

        # if mismatches / alignment length < 25% -> throw out
        # Collect the headers of the fasta files with +75% similarity

        # Filtered fasta file fofn
        fofnFileName = "%s.minimap_filter.fofn" % PREFIX
        fh = open(fofnFileName, "w")

        pafFofnFile = "%s.minimap.fofn" % PREFIX

        # Parse the fofn file
        nuccoreDBFofn = configDict["NUCCORE"].strip() # path to nuccore fofn
        nuccoreFilePathDict = {}
        for filePath in open(nuccoreDBFofn, "r"):
            fastaFileID = filePath.split("/")[-1].split(".")[0]
            nuccoreFilePathDict[fastaFileID] = filePath.strip()
        #####

        if (not isdir("minimap_filter")):
            CMDLIST.append("echo Making new directory minimap_filter >> %s" % logFile)
            CMDLIST.append("mkdir minimap_filter")

        # extract the sequences
        for pafFilePath in open(pafFofnFile, "r"):
            fastaFileID = pafFilePath.split("/")[-1].split(".")[0] # exclude the file ending
            newFastaFileName = fastaFileID + ".filtered.fasta"
            alignSetAccessions = minimap_filter_alignments(pafFilePath.strip(), configDict["MIN_MINIMAP_ALIGNLEN"])
            # for line in open(pafFile.strip(), "r"):
            #     alignData = line.split("\t")

            #     numMismatches = float(alignData[9]) # Col 10
            #     alignLen = float(alignData[10]) # Col 11
            #     # print("Mismatches: %d, AlignLen: %d" % (numMismatches, alignLen))
            #     if (numMismatches / alignLen < .25 and alignLen > def_minimap_align_len_cutoff):
            #         # good enough alignment
            #         alignSetAccessions.add(alignData[5])
            #     #####
            # #####
            if alignSetAccessions:
                # list is not empty
                extractIDs = "|".join(alignSetAccessions) # Regex OR for exact text matches
                if (nuccoreFilePathDict[fastaFileID].split(".")[-1] == "gz"):
                    extractCmd = """gunzip -dc %s | awk '/%s/{p++;print;next} /^>/{p=0} p' > %s""" % (nuccoreFilePathDict[fastaFileID], extractIDs, newFastaFileName)
                else:
                    extractCmd = """cat %s | awk '/%s/{p++;print;next} /^>/{p=0} p' > %s""" % (nuccoreFilePathDict[fastaFileID], extractIDs, newFastaFileName)
                #####
                CMDLIST.append("echo Retrieving %s >> %s" % (",".join(alignSetAccessions), logFile))
                CMDLIST.append(extractCmd)
                CMDLIST.append("mv %s ./minimap_filter" % newFastaFileName)
                fh.write("%s/%s/minimap_filter/%s\n" % (basePath, PREFIX, newFastaFileName)) # tracking directory
            #####
            #####
        #####
        # CMDLIST.append("mv *.filtered.fasta ./minimap_filter")

        fh.close()
        if (getsize(fofnFileName) == 0):
            CMDLIST.append("echo No alignments found. %s is empty. Terminating pipeline." % fofnFileName)
            nextStage = "end"
        else:
            # Continue with next stage
            nextStage = 'rehead_fasta'
        #####

# ================= STAGE C ===================== #
    elif (stage == "rehead_fasta"):

        # Filtered fasta files
        fofnFileName = "%s.minimap_filter.fofn" % PREFIX

        fofnNewFileName = "%s.fasta.rehead.fofn" % PREFIX
        fh = open(fofnNewFileName, "w")



        for fastaFile in open(fofnFileName, "r"):
            # print(fastaFile)
            path, fileName, fileNameHead = extractFileName(fastaFile.strip())
            newFastaFileName = fileNameHead + ".minimap_rehead.fasta"
            fh.write("%s/%s/minimap_rehead/%s\n" % (basePath, PREFIX, newFastaFileName)) # Replace later
            reHead_fasta(fastaFile.strip(), newFastaFileName)

        fh.close()

        if (not isdir("minimap_rehead")):
            CMDLIST.append("mkdir minimap_rehead")
        CMDLIST.append("mv *.minimap_rehead.fasta ./minimap_rehead")
        # CMDLIST.append("rm -r minimap_filter")
        # CMDLIST.append("mv *.minimap_rehead.fofn *.minimap_filter.fofn")
        # CMDLIST.append("mv minimap_rehead ./minimap_filter")

        nextStage = "blastn"
    #####



    # =============== MOVE TO NEXT STAGE =============== #

    # logFileHandle.write("End Time: %s\n" % datetime.datetime.now())
>>>>>>> 610dd852ddf62c86371230fb32e4a7cf781260ba
    logFileHandle.close()

    endState = "echo Ended stage %s >> %s" % (stage, logFile)
    endTime = 'date "+%H:%M:%S   %d/%m/%y" >> ' + ("%s" % logFile)
    CMDLIST.append(endState)
    CMDLIST.append(endTime)

    if not finished:
        nextStageWrite =  "echo Starting stage %s >> %s" % (nextStage, logFile)
        CMDLIST.append(nextStageWrite)

        nextCmd = "python3 %s/nifHupdate_Lib/nifHupdate_launch.py %s %s %s %s\n" % (basePath, configFile, nextStage, basePath, logFile)
        nextStageEcho = "echo Starting stage %s" % nextStage
        CMDLIST.append(nextStageEcho)
        CMDLIST.append(nextCmd)
    #####

<<<<<<< HEAD
=======

>>>>>>> 610dd852ddf62c86371230fb32e4a7cf781260ba
    shFileName = createShFile(CMDLIST, basePath, PREFIX, stage)

    launch(shFileName)

    return 0
#==============================================================
if ( __name__ == '__main__' ):
    real_main()