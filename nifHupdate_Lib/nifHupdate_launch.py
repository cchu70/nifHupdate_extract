#!/usr/bin/python
__author__="Claudia Chu"
__date__ ="2/24/19"

from nifHupdate_lib import parseConfig, createShFile, parseDate, launch, esearchCmds, fastaCmds, blastnCmds, bestAlignment, \
 throwError, verifyDb, test, testPrintFile, wait, fasta, trimSeq, mapBlast

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

    # TESTING
    # test("In Launch File")
    # ######


    # Read the command line information
    configFile  = argv[1]
    stage       = argv[2]
    basePath    = argv[3]
    logFile     = argv[4]

    # TESTING
    # test("ConfigFile %s" % configFile.split("/")[-1], "stage %s" % stage, "basePath %s" % basePath)
    ######


    # for year in range(startDate, endDate):
    #     if (year == startDate):
    #         startDate = configDict["START"]
    #     if (year == endDate):
    #         startDate = configDict["END"]




    # tracking time
    logFileHandle = open(logFile, "a")
    logFileHandle.write("================================\n")
    logFileHandle.write("Stage: %s\n" % stage)
    logFileHandle.write("Time Start: %s\n" % datetime.datetime.now())
    #tic = time.process_time()


    # Read the config file
    configDict = parseConfig(configFile, basePath, logFileHandle)

    # Get numerical year
    startDate = parseDate(configDict["START"])
    endDate = parseDate(configDict["END"])

    # FREQUENTLY USED CONSTANTS AND VARIABLES
    CMDLIST = []
    PREFIX = configDict["PREFIX"]


    # # Move into current running directory
    # cdFolder  = "cd %s" % (configDict["PREFIX"])
    # CMDLIST.append(cdFolder)

    # TESTING
    # test("PREFIX: %s" % PREFIX, "STARTDATE: %s" % startDate, "ENDDATE: %s" % endDate)
    ######


# ================= STAGE 1 ===================== #
# esearch: getting the intial records and accession number

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

        fh.close()
        # Next stage
        nextStage = 'fasta'
        echoCmd = "echo Launching next stage %s" % nextStage
        nextCmd = "python3 %s/nifHupdate_Lib/nifHupdate_launch.py %s %s %s %s\n" % (basePath, configFile, nextStage, basePath, logFile)
        CMDLIST.append(echoCmd)
        CMDLIST.append(nextCmd)

        shFileName = createShFile(CMDLIST, basePath, PREFIX, stage)


        # # TESTING
        # testPrintFile(shFileName)
        # ######

        launch(shFileName)


# ================= STAGE 2 ===================== #
    elif (stage == 'fasta'):

        # # TESTING
        # test("In fasta!")
        # ######

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

                    print("Retrieving %s" % fastaFileName)

                    # # # OLD
                    # if sortterm == "nifH":
                    #     CMDLIST.append("echo Retrieving %s ..." % fastaFileName)
                    #     CMDLIST.append(fastaCmds(esearchFile.strip(), sortterm, basePath, fastaFileName))
                    #     print("getting %s" % esearchFile)
                    # else:
                    #     # just one accession at a time
                    fasta(esearchFile.strip(), sortterm, fastaFileName)
                    # t = threading.Thread(target=fasta, args= (esearchFile.strip(), sortterm, fastaFileName))
                    # threads.append(t)
                    # t.start()
                    # t.join()
                #####
            #####
        #####

        fh.close()

        # testing the time it takes to run
        logFileHandle.write("End Time: %s\n" % time.process_time())


        # clean tmp files
        rmCmd = "rm tmp*"
        # Next stage
        nextStage = 'set_db'
        nextCmd = "python3 %s/nifHupdate_Lib/nifHupdate_launch.py %s %s %s %s\n" % (basePath, configFile, nextStage, basePath, logFile)
        CMDLIST.append(rmCmd)
        CMDLIST.append(nextCmd)
        shFileName = createShFile(CMDLIST, basePath, PREFIX, stage)

        # # TESTING
        # testPrintFile(shFileName)
        # ######
        launch(shFileName)

# ================= STAGE 3 ===================== #
# Check if all the files exist for a local database for stand alone alignment

    elif (stage == 'set_db'):
        # # TESTING
        # test("In set_db!")
        # ######

        # check database already exists
        if (not verifyDb(configDict["DBFILE"])):
            # move to next stage
            makeDbCmd = "makeblastdb -in %s/%s -parse_seqids -dbtype nucl -out %s/%s/%s" % (basePath, configDict["DBFILE"], basePath, PREFIX, configDict["DBNAME"])
            CMDLIST.append(makeDbCmd)
        #####
        logFileHandle.write("End Time: %s\n" % time.process_time())

        nextStage = 'blastn'
        nextCmd = "python3 %s/nifHupdate_Lib/nifHupdate_launch.py %s %s %s %s" % (basePath, configFile, nextStage, basePath, logFile)
        CMDLIST.append(nextCmd)

        shFileName = createShFile(CMDLIST, basePath, PREFIX, stage)
        # # TESTING
        # testPrintFile(shFileName)
        # ######
        launch(shFileName)


# ================= STAGE 4 ===================== #
# Aligning each of the query sequences to the subject sequences from local database.

    elif (stage == 'blastn'):
        # # TESTING
        # test("In blastn!")
        # ######

        fastaFofn = "%s.fasta.fofn" % PREFIX
        if (not isfile(fastaFofn)):
            throwError("fasta stage failed: %s not found" % esearchFofn, logFileHandle)
        else:

            fofnFileName = "%s.blastnFiles.fofn" % (PREFIX)
            fh = open(fofnFileName, "w")

            for fastaFile in open(fastaFofn, "r"):
                print(fastaFile)
                prefix, source, end = fastaFile.strip().split(".", 2)
                outputFile = "%s.%s.blastn.txt" % (prefix, source)
                fh.write(outputFile + "\n") # write to blast fofn file
                CMDLIST.append("echo Blasting %s ..." % fastaFile.strip())
                CMDLIST.append(blastnCmds(configDict["DBNAME"], fastaFile.strip(), outputFile))

            fh.close()
        #####
        logFileHandle.write("End Time: %s\n" % time.process_time())

        nextStage = 'filter_best_alignments'
        nextCmd = "python3 %s/nifHupdate_Lib/nifHupdate_launch.py %s %s %s %s\n" % (basePath, configFile, nextStage, basePath, logFile)

        CMDLIST.append(nextCmd)
        shFileName = createShFile(CMDLIST, basePath, PREFIX, stage)

        launch(shFileName)

# ================= STAGE 5 ===================== #
# Parsing through blastn tables to find best alignment for each sequence

    elif (stage == 'filter_best_alignments'):

        fofnFileName = "%s.blastnFiles.fofn" % (PREFIX)
        if (not isfile(fofnFileName)) :
            throwError("%s is not available. Either re-run blastn stage, or cat all your blastn files into this file name." % fofnFileName, logFileHandle)
        #####


        for blastnFile in open(fofnFileName, "r"):
            prefix, source, end = blastnFile.strip().split(".", 2)
            fileName = "%s.%s.blastn.filter.txt" % (prefix, source)
            fh = open(fileName, "w")
            bestAlignment(blastnFile.strip(), fh)
            fh.close()

        #####

        logFileHandle.write("End Time: %s\n" % time.process_time())

        nextStage = 'trim_seq'
        nextCmd = "python3 %s/nifHupdate_Lib/nifHupdate_launch.py %s %s %s %s\n" % (basePath, configFile, nextStage, basePath, logFile)

        CMDLIST.append(nextCmd)
        shFileName = createShFile(CMDLIST, basePath, configDict["PREFIX"], stage)
        launch(shFileName)

# ================= STAGE 6 ===================== #

    elif (stage == 'trim_seq'):



        # parse the blastn files into a map
        blastnFofn = "%s.blastnFiles.fofn" % (configDict["PREFIX"])
        blastnMap = mapBlast(blastnFofn)

        # for each year and each fasta file (fa_list.txt)
        # 1) get the accession number
        # 2) get associated blastn data: accession number, cluster, type?, cds/genome
        # 3) make new record name with updated data
        # 4) trim query sequence (qstart qend)
        # 5) print to a new file

        # >AF484654;cluster_I;Mesorhizobium_sp._LMG_11892


        fastaTrimmedFofn = "%s.fasta.trimmed.fofn" % PREFIX
        fh = open(fastaTrimmedFofn, "w")

        fastaFofn = "%s.fasta.fofn" % PREFIX
        for fastaFile in open(fastaFofn, "r"):
            prefix, source, end = fastaFile.strip().split(".")

            trimFastaFileName = "%s.%s.trimmed.fasta" % (prefix, source)
            fh.write("%s\n" % trimFastaFileName)
            trimSeq(fastaFile, trimFastaFileName, blastnMap)

        #####
        fh.close()
        logFileHandle.write("End Time: %s\n" % time.process_time())

        nextStage = 'cluster'
        nextCmd = "python3 %s/nifHupdate_Lib/nifHupdate_launch.py %s %s %s %s\n" % (basePath, configFile, nextStage, basePath, logFile)

        CMDLIST.append(nextCmd)
        shFileName = createShFile(CMDLIST, basePath, configDict["PREFIX"], stage)
        launch(shFileName)

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

        for clusterFileHandle in clusterFileHandles:
            clusterFileHandles[clusterFileHandle].close()

        logFileHandle.write("End Time: %s\n" % time.process_time())

        nextStage = 'deduplicate'
        nextCmd = "python3 %s/nifHupdate_Lib/nifHupdate_launch.py %s %s %s %s\n" % (basePath, configFile, nextStage, basePath, logFile)

        CMDLIST.append(nextCmd)
        shFileName = createShFile(CMDLIST, basePath, configDict["PREFIX"], stage)
        launch(shFileName)

# ================= STAGE 8 ===================== #
    elif (stage == 'deduplicate'):
    # cd-hit-dup -i fasta -o output
        clusterFastaFofn = "%s.cluster_fasta.fofn" % PREFIX

        for fastaFile in open(clusterFastaFofn, "r"):
            deduplicate(fastaFile.strip())
        #####
        logFileHandle.write("End Time: %s\n" % time.process_time())

    # Reading stderr


    return 0
#==============================================================
if ( __name__ == '__main__' ):
    real_main()